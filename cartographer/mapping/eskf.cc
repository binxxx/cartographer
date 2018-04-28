#include "cartographer/mapping/eskf.h"

namespace cartographer {
namespace mapping {


Eigen::Matrix3d skew(Eigen::Vector3d w) {
    Eigen::Matrix3d W;
    W <<    0.0, -w.z(),  w.y(),
          w.z(),    0.0, -w.x(),
         -w.y(),  w.x(),    0.0;
    return W;
}

Eigen::Matrix3d angle_axis_to_rotation_matrix(Eigen::Vector3d angle_axis) {

    double scale = 0.5;
    double w = 1.0;
    constexpr double kCutoffAngle = 1e-8;  // We linearize below this angle.
    if (angle_axis.squaredNorm() > kCutoffAngle) {
        const double norm = angle_axis.norm();
        scale = sin(norm / 2.) / norm;
        w = cos(norm / 2.);
    }
    const Eigen::Matrix<double, 3, 1> quaternion_xyz = scale * angle_axis;
    return Eigen::Quaternion<double>(w, quaternion_xyz.x(), quaternion_xyz.y(),
                              quaternion_xyz.z()).toRotationMatrix();
}

Eigen::Matrix3d euler_angle_to_rotation_matrix(Eigen::Vector3d w) {
    Eigen::Matrix3d R;
    R = Eigen::AngleAxisd(w[2], Eigen::Vector3d::UnitZ()) *
        Eigen::AngleAxisd(w[1], Eigen::Vector3d::UnitY()) *
        Eigen::AngleAxisd(w[0], Eigen::Vector3d::UnitX());
    return R;
}

ESKF::ESKF(const mapping::proto::ESKFOptions& options) {

	// TODO(Weikun): initialize params by options
    _imu_freq    = options.imu_frequency();
    _imu_frame   = options.imu_frame();
    _robot_frame = options.robot_frame();
    _imu_enabled = options.imu_enabled();
    _imu_has_quat = options.imu_has_quaternion();
    _smooth_enabled = options.smooth_enabled();
    _smooth_buf_size = options.smooth_buffer_size();
    _smooth_type = options.smooth_type();
    
    _init_roll = options.initial_roll();
    _init_pitch = options.initial_pitch();
    _init_yaw = options.initial_yaw();
    _init_x = options.initial_x();
    _init_y = options.initial_y();
    _init_z = options.initial_z();

    _sigma_acc = options.sigma_accelerometer();
    _sigma_gyr = options.sigma_gyroscope();
    _sigma_bias_acc = options.sigma_bias_accelerometer();
    _sigma_bias_gyr = options.sigma_bias_gyroscope();
    _g = options.gravity();
    _init_bias_acc_x = options.initial_bias_accelerometer_x();
    _init_bias_acc_y = options.initial_bias_accelerometer_y();
    _init_bias_acc_z = options.initial_bias_accelerometer_z();
    _acc_queue_size = options.accelerometer_queue_size();
    _imu_transform = options.imu_transform();

    // initialize nomial states
    _velocity.setZero();
    _quaternion = Eigen::AngleAxisd(_init_yaw,  Eigen::Vector3d::UnitZ())
                * Eigen::AngleAxisd(_init_pitch, Eigen::Vector3d::UnitY())
                * Eigen::AngleAxisd(_init_roll,   Eigen::Vector3d::UnitX());
    _rotation = _quaternion.toRotationMatrix();

    _position << _init_x, _init_y, _init_z;
    _bias_acc << _init_bias_acc_x, _init_bias_acc_y, _init_bias_acc_z;
    _bias_gyr.setZero();

    // initialize error states
    _d_velocity.setZero();
    _d_theta.setZero();
    _d_rotation.setIdentity();
    _d_position.setZero();
    _d_bias_acc.setZero();
    _d_bias_gyr.setZero();

    // initialize imu
    _imu_acceleration.setZero();
    _imu_angular_velocity.setZero();
    _imu_orientation.setIdentity();

    // initialize measurements
    _m_position.setZero();
    _m_theta.setZero();

    // initialize Jacobian matrix;
    _Fx.setZero();
    _Fn.setZero();

    // initialize covariance matrix
     _Sigma.setZero();
     // _Sigma(3,3) = 0.1;
     // _Sigma(4,4) = 0.1;
     // _Sigma(5,5) = 0.1;
    _Q.setZero();

    // gravity
    _gravity << 0.0,0.0,_g;
    // time relatives
    _init_time = true;

    // acc queue
    _acc_queue_count = 0;
    
    //smoother
    _smooth_buf_cnt = 0;
    _x_buf.resize(_smooth_buf_size);
    _y_buf.resize(_smooth_buf_size);
    _z_buf.resize(_smooth_buf_size);
    _vx_buf.resize(_smooth_buf_size);
    _vy_buf.resize(_smooth_buf_size);
    _vz_buf.resize(_smooth_buf_size);
}

ESKF::~ESKF() {}

void ESKF::AddImuData(const sensor::ImuData& imu_data) {

    UpdateTime(imu_data);
    UpdateImu(imu_data);

    PropagateState();
    PropagateCovariance();
}

void ESKF::Update(const Eigen::Matrix<double, 6, 1> &mean_meas,
                  const Eigen::Matrix<double, 6, 6> &cov_meas) {

    for(int i=0; i<mean_meas.size(); i++) {
        if(std::isnan(mean_meas(i))) {
            LOG(WARNING) << "GPF: Recovered measurements contains NaNs.";
            return;
        }
    }
    for(int i=0; i<cov_meas.size(); i++) {
        if(std::isnan(cov_meas(i))) {
            LOG(WARNING) << "GPF: Recovered measurements contains NaNs.";
            return;
        }
    }

    UpdateMeasurementMean(mean_meas);
    UpdateMeasurementCovariance(cov_meas);

    UpdateError();
    UpdateState();
    ResetError();
}

void ESKF::UpdateTime(const sensor::ImuData& imu_data) {

    if(_init_time) {
        _dt = 1.0 / _imu_freq;
        _init_time = false;
    }
    else {
        _dt = common::ToUniversal(imu_data.time) - common::ToUniversal(_imu_time);
        _dt /= 10000000.0;
        // std::cout << _dt << std::endl;
    }
    _imu_time = imu_data.time;

}

void ESKF::UpdateImu(const sensor::ImuData& imu_data) {

    // stacking into a queue
    if(_acc_queue_count < _acc_queue_size) {
        _acc_queue.push_back(imu_data.linear_acceleration);
        _imu_acceleration = imu_data.linear_acceleration;
    }
    else {
        _acc_queue[_acc_queue_count%_acc_queue_size] = imu_data.linear_acceleration;

        Eigen::Vector3d acc_avg;
        acc_avg.setZero();
        for(int i=0; i<_acc_queue_size; i++) {
            acc_avg[0] += _acc_queue[i][0] / (double)_acc_queue_size;
            acc_avg[1] += _acc_queue[i][1] / (double)_acc_queue_size;
            acc_avg[2] += _acc_queue[i][2] / (double)_acc_queue_size;
        }
        _imu_acceleration = acc_avg;
    }
    _acc_queue_count++; 
	_imu_angular_velocity = imu_data.angular_velocity;
    _imu_orientation = Eigen::Quaterniond::Identity();
    _imu_orientation = imu_data.orientation;

    // If imu disabled, set acc, grav to zero
    if(!_imu_enabled) {
        _imu_acceleration.setZero();
        _gravity.setZero();
    } else {
        // If imu has quaternion, remove gravity.
        if(_imu_has_quat) {
            Eigen::Matrix3d imu_rot = _imu_orientation.toRotationMatrix();
            Eigen::Vector3d grav(0.0,0.0,_g);
            _imu_acceleration += imu_rot.transpose() * grav;
            _gravity.setZero();
        }
    }
}
   

void ESKF::PropagateState() {
    Eigen::Vector3d velocity;
    Eigen::Matrix3d rotation;
    Eigen::Vector3d position;
    Eigen::Vector3d bias_acc;
    Eigen::Vector3d bias_gyr;

    // system transition function for nominal state
    velocity = _velocity + (_rotation * (_imu_acceleration - _bias_acc) + _gravity ) * _dt;
    rotation = _rotation * euler_angle_to_rotation_matrix((_imu_angular_velocity - _bias_gyr) * _dt);
    position = _position + _velocity * _dt + 0.5 * (_rotation * (_imu_acceleration - _bias_acc) + _gravity) * _dt * _dt;
    bias_acc = _bias_acc;
    bias_gyr = _bias_gyr;

    // update norminal state to the next step
    _velocity   = velocity;
    _rotation   = rotation;
    _quaternion = Eigen::Quaterniond(_rotation);
    _position   = position;
    _bias_acc   = bias_acc;
    _bias_gyr   = bias_gyr;
}

void ESKF::PropagateCovariance() {
    // compute jacobian
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d Z = Eigen::Matrix3d::Zero();

    Eigen::Matrix<double, 3, 3> R, R_1, R_2;
    R = _rotation;
    R_1 = skew(_imu_acceleration - _bias_acc);
    R_2 = euler_angle_to_rotation_matrix((_imu_angular_velocity - _bias_gyr) * _dt);

    _Fx <<     I,      Z,       -R*R_1*_dt,   -R*_dt,        Z,
           I*_dt,      I,                Z,        Z,        Z,
               Z,      Z,  R_2.transpose(),        Z,   -I*_dt,
               Z,      Z,                Z,        I,        Z,
               Z,      Z,                Z,        Z,        I;

    _Fn << R, Z, Z, Z,
           Z, Z, Z, Z,
           Z, I, Z, Z,
           Z, Z, I, Z,
           Z, Z, Z, I;

    _Q << pow(_sigma_acc * _dt, 2.0) * I, Z, Z, Z,
          Z, pow(_sigma_gyr * _dt, 2.0) * I, Z, Z,
          Z, Z, pow(_sigma_bias_acc * _dt, 2.0) * I, Z,
          Z, Z, Z, pow(_sigma_bias_gyr * _dt, 2.0) * I;

    // update covariance
    _Sigma = _Fx * _Sigma * _Fx.transpose() + _Fn * _Q * _Fn.transpose();
    // std::cout << "dt=" << _Sigma << std::endl;
}

Eigen::Matrix<double, 7, 1> ESKF::GetPoseMean() {
    Eigen::Matrix<double, 7, 1> mean_pose;
    mean_pose[0] = _position.x();
    mean_pose[1] = _position.y();
    mean_pose[2] = _position.z();

    mean_pose[3] = _quaternion.w();
    mean_pose[4] = _quaternion.x();
    mean_pose[5] = _quaternion.y();
    mean_pose[6] = _quaternion.z();
    return mean_pose;
}

transform::Rigid3d ESKF::GetLatestPose() {
    transform::Rigid3d latest_pose(_position, _quaternion);
    return latest_pose; 
}
Eigen::Matrix<double, 6, 6> ESKF::GetPoseCovariance() {
    Eigen::Matrix<double, 6, 6> cov_pose;
    cov_pose = _Sigma.block<6,6>(3,3);

    return cov_pose;
}

Eigen::Quaterniond ESKF::GetGravityOrientation() {
    // return Eigen::Quaterniond::Identity();
    return _quaternion;
}

void ESKF::UpdateMeasurementMean(const Eigen::Matrix<double, 6, 1> &mean_meas) {
    _m_position = mean_meas.block<3,1>(0,0);
    _m_theta = mean_meas.block<3,1>(3,0);
}

void ESKF::UpdateMeasurementCovariance(const Eigen::Matrix<double, 6, 6> &cov_meas) {
    _m_pose_sigma = cov_meas;
}


void ESKF::UpdateError() {
    // assume only pose measurement is used
    Eigen::Matrix<double, 6, 15> H;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d Z = Eigen::Matrix3d::Zero();
    H << Z, I, Z, Z, Z,
         Z, Z, I, Z, Z;

    // measurements
    Eigen::Matrix<double, 6, 1> y;
    y[0] = _m_position.x(); y[1] = _m_position.y(); y[2] = _m_position.z();
    y[3] = _m_theta.x();    y[4] = _m_theta.y();    y[5] = _m_theta.z();

    // kalman gain matrix
    Eigen::Matrix<double, 15, 6> K;
    K = _Sigma * H.transpose() * (H * _Sigma * H.transpose() + 1.0 * _m_pose_sigma).inverse();

    // update error
    Eigen::Matrix<double, 15, 1> x;
    x = K * y;

    _d_velocity << x[0],  x[1],  x[2];
    _d_position << x[3],  x[4],  x[5];
    _d_theta    << x[6],  x[7],  x[8];
    _d_bias_acc << x[9],  x[10], x[11];
    _d_bias_gyr << x[12], x[13], x[14];
    _d_rotation << angle_axis_to_rotation_matrix(_d_theta);
    // std::cout << x.block<6,1>(3,0).transpose() << std::endl;
    // update covariance
    Eigen::Matrix<double, 15, 15> M;
    M = Eigen::MatrixXd::Identity(15,15) - K*H;
    _Sigma = M * _Sigma;// * M.transpose() + K * _m_pose_sigma * K.transpose();

}

void ESKF::UpdateState() {
    _velocity += _d_velocity;
    _position += _d_position;
    _rotation *= _d_rotation;
    _bias_acc += _d_bias_acc;
    _bias_gyr += _d_bias_gyr;
    _quaternion = Eigen::Quaterniond(_rotation);
}

void ESKF::ResetError() {
    _d_velocity.setZero();
    _d_position.setZero();
    _d_theta.setZero();
    _d_rotation.setIdentity();
    _d_bias_acc.setZero();
    _d_bias_gyr.setZero();
}


}  // namespace mapping
}  // namespace cartographer