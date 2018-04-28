#include "cartographer/mapping/gpf.h"

namespace cartographer {
namespace mapping {

template <typename T>
T add(T p1, T p2) {
    T p;
    p.x = p1.x + p2.x;
    p.y = p1.y + p2.y;
    p.y = p1.z + p2.z;
    return p;
}

template <typename T>
T minus(T p1, T p2) {
    T p;
    p.x = p1.x - p2.x;
    p.y = p1.y - p2.y;
    p.y = p1.z - p2.z;
    return p;
}

template <typename T, typename S>
T multiply(T p, S m) {
    T q;
    q.x = p.x * m;
    q.y = p.y * m;
    q.z = p.z * m;
    return q;
}

template <typename T>
double dot(T p1, T p2) {
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

template <typename T>
T cross(T p1, T p2) {
    T p;
    p.x = p1.y*p2.z - p2.y*p1.z;
    p.y = p2.x*p1.z - p1.x*p2.z;
    p.z = p1.x*p2.y - p2.x*p1.y;
    return p;
}

template <typename T>
double norm(T p) {
    return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}
template <typename T>
T normalize(T p) {

    double n;
    n = norm(p);

    T np = multiply(p, 1/n);
    return np;
}
template <typename T>
std::vector<size_t> sort_index(const std::vector<T> v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

GPF::GPF(){


    // initialize particles pointer
    _ray_sigma = 1.0;
    _cloud_resol = 0.3;
    _set_size = 500;
    _cloud_range = 60.0;
    robot_frame_ = std::string("/base_frame");

    mean_prior_.setZero();
    mean_sample_.setZero();
    mean_posterior_.setZero();
    mean_meas_.setZero();

    Eigen::Matrix<double, 6, 1> sigma;
    sigma << 0.01, 0.01, 0.01, 0.005, 0.005, 0.005;
    cov_prior_ = sigma.asDiagonal();
    cov_sample_.setZero();
    cov_posterior_.setZero();
    cov_meas_.setZero();

    // initialize eskf
    // _eskf_ptr = boost::shared_ptr<ESKF> (new ESKF(nh));

    // initialize particle
    _particles_ptr = boost::shared_ptr<Particles> (new Particles(_set_size, _ray_sigma));

}

GPF::~GPF() {}

void GPF::SetPosePrior(const Eigen::Matrix<double, 7, 1> mean_prior,
                      const Eigen::Matrix<double, 6, 6> cov_prior) {
    mean_prior_ = mean_prior;
    cov_prior_  = cov_prior;
}

void GPF::SetMatchingDistanceMap(const std::shared_ptr<DynamicEDTOctomap> distmap_grid) {
    distmap_grid_ = distmap_grid;
}

void GPF::AddRangeData(const common::Time time,
                       const sensor::RangeData range_data) {
    laser_time_ = time;
    if(distmap_grid_ == nullptr) {
        LOG(WARNING) << "distmap_grid not set";
        return;
    }

    // down sample data with a voxel filter
    const sensor::PointCloud down_sampled_cloud = 
            sensor::VoxelFiltered(range_data.returns,
                                  _cloud_resol);
    if(cloud_ptr_ == nullptr) {
        cloud_ptr_ = std::unique_ptr<sensor::PointCloud>(
            new sensor::PointCloud);
    } else {
        cloud_ptr_->clear();
    }
    for(const Eigen::Vector3f& hit : down_sampled_cloud) {
        const Eigen::Vector3f delta = hit - range_data.origin;
        const float range = delta.norm();
        if (range >= 0.5) {
          if (range <= _cloud_range) {
            cloud_ptr_->push_back(hit);
          }
        }
    }
    // std::cout << "GPF::AddRangeData()::VoxelFiltered size = " << cloud_ptr_->size() << std::endl;
    
    _particles_ptr->SetMean(mean_prior_);
    _particles_ptr->SetCovariance(cov_prior_);
    _particles_ptr->SetCloud(cloud_ptr_);
    _particles_ptr->SetMatchingDistanceMap(distmap_grid_);
    _particles_ptr->Propagate(mean_sample_, cov_sample_,
                              mean_posterior_, cov_posterior_);

    // std::cout << "mean prior:    " << mean_prior_.transpose() << std::endl;
    // std::cout << "mean sample:   " << mean_sample_.transpose() <<std::endl;
    // std::cout << "mean posterior:" << mean_posterior_.transpose() << std::endl;

    // update meas in eskf
    RecoverPoseErrorMeasurements();
    
    // std::cout << "mean meas:     " << mean_meas_.transpose() << std::endl;
}

Eigen::Matrix<double, 6, 1> GPF::GetMeasurementPose() {
    return mean_meas_;
}

Eigen::Matrix<double, 6, 6> GPF::GetMeasurementCovariance() {
    return cov_meas_;
}

void GPF::RecoverPoseErrorMeasurements() {
    Eigen::Matrix<double, 6, 6> K;
    Eigen::Matrix<double, 6, 1> mean_meas;
    Eigen::Matrix<double, 6, 6> cov_meas;

    cov_meas_ = (cov_posterior_.inverse() - cov_sample_.inverse()).inverse();
    CheckPositiveDefinite(cov_meas_);
    K = cov_sample_ * (cov_sample_ + cov_meas_).inverse();
    mean_meas_ = K.inverse() * (mean_posterior_ - mean_sample_) + mean_sample_;

    // check if the recovered pseudo mesure is valid
    for(int i=0; i<mean_meas_.size(); i++) {
        if(std::isnan(mean_meas(i))) {
            LOG(WARNING) << "GPF: Recovered measurements contains NaNs.";
            return;
        }
    }
    for(int i=0; i<cov_meas_.size(); i++) {
        if(std::isnan(cov_meas(i))) {
            LOG(WARNING) << "GPF: Recovered measurements contains NaNs.";
            return;
        }
    }

    // mean_meas_ = mean_meas;
    // cov_meas_ = cov_meas;
}

void GPF::CheckPositiveDefinite(Eigen::Matrix<double, 6, 6> &R) {
    Eigen::Matrix<double, 6, 6> rot;
    Eigen::Matrix<double, 6, 1> scl;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6> >
       eigenSolver(R);
    rot = eigenSolver.eigenvectors();
    scl = eigenSolver.eigenvalues();

    Eigen::Matrix<double, 6, 6> E;
    E.setZero();
    for (int ii=0; ii<6; ++ii) {
        if(scl(ii,0)>0.0) {
            E(ii,ii) = scl(ii,0);
        }
        else {
            E(ii,ii) = 100.0;
        }
    }
    R = rot * E * rot.inverse();
}


}  // namespace mapping
}  // namespace cartographer