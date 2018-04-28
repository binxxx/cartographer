
#ifndef CARTOGRAPHER_MAPPING_ESKF_H_
#define CARTOGRAPHER_MAPPING_ESKF_H_

#include <deque>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "cartographer/common/time.h"
#include "cartographer/sensor/imu_data.h"
#include "cartographer/transform/rigid_transform.h"
#include "cartographer/mapping/proto/eskf_options.pb.h"

namespace cartographer {
namespace mapping {

class ESKF {
  public:
    ESKF(const mapping::proto::ESKFOptions& options);
    ~ESKF();

    void AddImuData(const sensor::ImuData& imu_data);
    transform::Rigid3d GetLatestPose();
    
    Eigen::Matrix<double, 7, 1> GetPoseMean();
    Eigen::Matrix<double, 6, 6> GetPoseCovariance();
    Eigen::Quaterniond GetGravityOrientation();
    void UpdateMeasurementMean(const Eigen::Matrix<double, 6, 1> &mean_meas);
    void UpdateMeasurementCovariance(const Eigen::Matrix<double, 6, 6> &cov_meas);
    void Update(const Eigen::Matrix<double, 6, 1> &mean_meas,
                const Eigen::Matrix<double, 6, 6> &cov_meas);

  private:
    void UpdateTime(const sensor::ImuData& imu_data);
    void UpdateImu(const sensor::ImuData& imu_data);
    void PropagateCovariance();
    void PropagateState();
    void UpdateError();
    void UpdateState();
    void ResetError();

    // nominal states
    Eigen::Vector3d    _velocity;
    Eigen::Matrix3d    _rotation;
    Eigen::Quaterniond _quaternion;
    Eigen::Vector3d    _position;
    Eigen::Vector3d    _bias_acc;
    Eigen::Vector3d    _bias_gyr;


    // error states
    Eigen::Vector3d   _d_velocity;
    Eigen::Vector3d   _d_theta;
    Eigen::Matrix3d   _d_rotation;
    Eigen::Vector3d   _d_position;
    Eigen::Vector3d   _d_bias_acc;
    Eigen::Vector3d   _d_bias_gyr;

    // imu measurements
    Eigen::Vector3d    _imu_acceleration;
    Eigen::Vector3d    _imu_angular_velocity;
    Eigen::Quaterniond _imu_orientation;

    // Jacobian matrices
    Eigen::Matrix<double, 15, 15> _Fx;
    Eigen::Matrix<double, 15, 12> _Fn;

    // covariance matrix
    Eigen::Matrix<double, 15, 15> _Sigma;
    Eigen::Matrix<double, 12, 12> _Q;

    // gravity
    double _g;
    Eigen::Vector3d _gravity;

    // time relatives
    common::Time _imu_time; 
    double _dt, _imu_freq;
    bool _init_time;

    // noise params
    double _sigma_acc, _sigma_gyr;
    double _sigma_bias_acc, _sigma_bias_gyr;

    // initialization params
    double _init_bias_acc_x, _init_bias_acc_y, _init_bias_acc_z;
    double _init_roll, _init_pitch, _init_yaw;
    double _init_x, _init_y, _init_z;

    // odometry measurements
    Eigen::Vector3d    _m_theta;
    Eigen::Vector3d    _m_position;

    Eigen::Matrix<double, 6, 6> _m_pose_sigma;

    // a queue to smooth imu accleration measurements
    std::vector<Eigen::Vector3d> _acc_queue;
    int _acc_queue_size;
    int _acc_queue_count;
    
    // frames
    std::string _imu_frame, _robot_frame;
   
    // imu related
    bool _imu_enabled, _imu_has_quat, _imu_transform;

    // smoother
    bool _smooth_enabled;
    int  _smooth_buf_size;
    int  _smooth_buf_cnt;
    std::string _smooth_type;
    std::vector<double> _x_buf;
    std::vector<double> _y_buf;
    std::vector<double> _z_buf;
    std::vector<double> _vx_buf;
    std::vector<double> _vy_buf;
    std::vector<double> _vz_buf;

};

Eigen::Matrix3d skew(Eigen::Vector3d w);
Eigen::Matrix3d angle_axis_to_rotation_matrix(Eigen::Vector3d w);
Eigen::Matrix3d euler_angle_to_rotation_matrix(Eigen::Vector3d w);

}  // namespace mapping
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_ESKF_H_
