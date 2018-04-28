
#ifndef CARTOGRAPHER_MAPPING_GPF_H_
#define CARTOGRAPHER_MAPPING_GPF_H_

#include <deque>
#include <memory>
#include <algorithm>
#include <numeric>
#include <functional>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "cartographer/common/time.h"
#include "cartographer/sensor/range_data.h"
#include "cartographer/sensor/voxel_filter.h"
#include "cartographer/transform/rigid_transform.h"
#include "cartographer/mapping_3d/submaps.h"
#include "cartographer/mapping/pose_estimate.h"
#include "cartographer/mapping/eskf.h"
#include "cartographer/mapping/particles.h"

namespace cartographer {
namespace mapping {


class GPF {
public:
    GPF();
    ~GPF();

    void AddRangeData(const common::Time time,
                      const sensor::RangeData range_data);
    void SetPosePrior(const Eigen::Matrix<double, 7, 1> mean_prior,
                      const Eigen::Matrix<double, 6, 6> cov_prior);
    void SetMatchingDistanceMap(const std::shared_ptr<DynamicEDTOctomap> distmap_grid);

    Eigen::Matrix<double, 6, 1> GetMeasurementPose();
    Eigen::Matrix<double, 6, 6> GetMeasurementCovariance();
private:

    void RecoverPoseErrorMeasurements();
    void CheckPositiveDefinite(Eigen::Matrix<double, 6, 6> &R);

    Eigen::Matrix<double, 7, 1> mean_prior_; // nominal
    Eigen::Matrix<double, 6, 6> cov_prior_;  // error
    Eigen::Matrix<double, 6, 1> mean_sample_;
    Eigen::Matrix<double, 6, 6> cov_sample_;
    Eigen::Matrix<double, 6, 1> mean_posterior_;
    Eigen::Matrix<double, 6, 6> cov_posterior_;
    Eigen::Matrix<double, 6, 1> mean_meas_;
    Eigen::Matrix<double, 6, 6> cov_meas_;

    common::Time laser_time_;
    std::string robot_frame_;

    std::shared_ptr<mapping_3d::Submap> matching_submap_;
    std::shared_ptr<DynamicEDTOctomap> distmap_grid_;
    std::shared_ptr<sensor::PointCloud> cloud_ptr_;
    boost::shared_ptr<ESKF>             _eskf_ptr;
    boost::shared_ptr<Particles>        _particles_ptr;

    double _cloud_resol;
    double _ray_sigma;
    int    _set_size;
    double _cloud_range;

    // nav_msgs::Path _path;
    // std::deque<geometry_msgs::PoseStamped> _pose_deque;
    // tf::TransformBroadcaster _tf_br;
};
}  // namespace mapping
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_GPF_H_
