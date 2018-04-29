
#ifndef CARTOGRAPHER_MAPPING_3D_ESKF_LOCAL_TRAJECTORY_BUILDER_H_
#define CARTOGRAPHER_MAPPING_3D_ESKF_LOCAL_TRAJECTORY_BUILDER_H_

#include <memory>

#include "cartographer/common/time.h"
#include "cartographer/mapping/eskf.h"
#include "cartographer/mapping/gpf.h"
#include "cartographer/mapping/pose_estimate.h"
#include "cartographer/mapping_3d/proto/eskf_local_trajectory_builder_options.pb.h"
#include "cartographer/mapping_3d/scan_matching/ceres_scan_matcher.h"
#include "cartographer/mapping_3d/scan_matching/real_time_correlative_scan_matcher.h"
#include "cartographer/mapping_3d/submaps.h"
#include "cartographer/sensor/imu_data.h"
#include "cartographer/sensor/odometry_data.h"
#include "cartographer/sensor/range_data.h"
#include "cartographer/sensor/voxel_filter.h"
#include "cartographer/transform/rigid_transform.h"

namespace cartographer {
namespace mapping_3d {

// Wires up the local SLAM stack (i.e. pose extrapolator, scan matching, etc.)
// without loop closure.
class ESKFLocalTrajectoryBuilder {
 public:
  struct InsertionResult {
    std::shared_ptr<const mapping::TrajectoryNode::Data> constant_data;
    transform::Rigid3d pose_observation;
    std::vector<std::shared_ptr<const Submap>> insertion_submaps;

    // add new data type
    // std::shared_ptr<const mapping::TrajectoryNode::Data> constant_data_mini;
    // int num_scans;
  };
  explicit ESKFLocalTrajectoryBuilder(
      const proto::ESKFLocalTrajectoryBuilderOptions& options, 
      common::ThreadPool* thread_pool);
  ~ESKFLocalTrajectoryBuilder();

  ESKFLocalTrajectoryBuilder(const ESKFLocalTrajectoryBuilder&) = delete;
  ESKFLocalTrajectoryBuilder& operator=(const ESKFLocalTrajectoryBuilder&) = delete;

  void AddImuData(const sensor::ImuData& imu_data);
  std::unique_ptr<InsertionResult> AddRangeData(
      common::Time time, const sensor::RangeData& range_data);
  void AddOdometerData(const sensor::OdometryData& odometry_data);
  const mapping::PoseEstimate& pose_estimate() const;
  std::shared_ptr<octomap::OcTree> GetMatchingOctomap();
  std::shared_ptr<Submap> GetMatchingSubmap();
  transform::Rigid3d GetMatchingLocalPose();

  /*********************************************************/
  std::unique_ptr<InsertionResult> AddRangeDataMini(common::Time time,
    const sensor::RangeData& range_data);


 private:
  // transform::Rigid3d pose_estimate_copy;
  // Eigen::Quaterniond gravity_alignment_copy;
  std::unique_ptr<InsertionResult> AddAccumulatedRangeData(
      common::Time time, const sensor::RangeData& range_data_in_tracking);
  // std::unique_ptr<InsertionResult> AddAccumulatedRangeData(
  //     common::Time time, const sensor::RangeData& range_data_in_tracking,
  //     const sensor::RangeData& single_range_data);


  std::unique_ptr<InsertionResult> InsertIntoSubmap(
      common::Time time, const sensor::RangeData& range_data_in_tracking,
      // const sensor::RangeData& single_range_data,
      const Eigen::Quaterniond& gravity_alignment,
      const sensor::PointCloud& high_resolution_point_cloud,
      const sensor::PointCloud& low_resolution_point_cloud,
      const transform::Rigid3d& pose_observation);

  common::ThreadPool* thread_pool_;
  const proto::ESKFLocalTrajectoryBuilderOptions options_;
  ActiveSubmaps active_submaps_;

  mapping::PoseEstimate last_pose_estimate_;

  std::unique_ptr<mapping::ESKF> eskf_;
  std::unique_ptr<mapping::GPF> gpf_;

  int num_accumulated_ = 0;
  transform::Rigid3d matching_submap_local_pose_ = transform::Rigid3d::Identity();
  transform::Rigid3f first_pose_estimate_ = transform::Rigid3f::Identity();
  std::unique_ptr<scan_matching::RealTimeCorrelativeScanMatcher>
      real_time_correlative_scan_matcher_;
  std::unique_ptr<scan_matching::CeresScanMatcher> ceres_scan_matcher_;
  sensor::RangeData accumulated_range_data_;
};

}  // namespace mapping_3d
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_3D_LOCAL_TRAJECTORY_BUILDER_H_
