
#include "cartographer/mapping_3d/eskf_local_trajectory_builder.h"

#include <memory>
#include "cartographer/common/make_unique.h"
#include "cartographer/common/time.h"
#include "cartographer/mapping_3d/proto/eskf_local_trajectory_builder_options.pb.h"
#include "cartographer/mapping_3d/proto/submaps_options.pb.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping_3d {

ESKFLocalTrajectoryBuilder::ESKFLocalTrajectoryBuilder(
    const proto::ESKFLocalTrajectoryBuilderOptions& options, 
    common::ThreadPool* thread_pool)
    : thread_pool_(thread_pool) ,
      options_(options),
      active_submaps_(options.submaps_options(), thread_pool),
      real_time_correlative_scan_matcher_(
          common::make_unique<scan_matching::RealTimeCorrelativeScanMatcher>(
              options_.real_time_correlative_scan_matcher_options())),
      ceres_scan_matcher_(common::make_unique<scan_matching::CeresScanMatcher>(
          options_.ceres_scan_matcher_options())),
      accumulated_range_data_{Eigen::Vector3f::Zero(), {}, {}} {}

ESKFLocalTrajectoryBuilder::~ESKFLocalTrajectoryBuilder() {}

void ESKFLocalTrajectoryBuilder::AddImuData(const sensor::ImuData& imu_data) {

  if(eskf_ == nullptr) {
    // TODO: Initialization with options
    eskf_ = std::unique_ptr<mapping::ESKF>(new mapping::ESKF(options_.eskf_options()));
    std::cout << "intializing eskf." << std::endl;
  }

  eskf_->AddImuData(imu_data);
}

std::unique_ptr<ESKFLocalTrajectoryBuilder::InsertionResult>
ESKFLocalTrajectoryBuilder::AddRangeData(const common::Time time,
                                         const sensor::RangeData& range_data) {
  // return nullptr;
  // std::cout << "ESKFLocalTrajectoryBuilder::AddRangeData()" << std::endl;

  if (eskf_ == nullptr) {
    LOG(INFO) << "IMU not yet initialized.";
    return nullptr;
  }
  if (gpf_ == nullptr) {
    // TODO: Initialization with gpf options
    gpf_ = std::unique_ptr<mapping::GPF>(new mapping::GPF());
  }

  std::shared_ptr<mapping_3d::Submap> 
  matching_submap = active_submaps_.GetMatchingSubmap();
  
  // std::cout << "ESKFLocalTrajectoryBuilder::AddRangeData()::Got matching submap" << std::endl;
  gpf_->SetMatchingDistanceMap(matching_submap
                                  ->GetDistanceMap());
  // std::cout << "ESKFLocalTrajectoryBuilder::AddRangeData():: Set distance map for GPF" << std::endl;
  // eskf_->TransformToCurrentDistanceMap(matching_submap
  //                                       ->local_pose()
  //                                       .inverse() *
  //                                      matching_submap_local_pose_);
  matching_submap_local_pose_ = matching_submap->local_pose();

  // std::cout << "ESKFLocalTrajectoryBuilder::AddRangeData()" << std::endl;
  // std::cout << "eskf latest pose:" << eskf_->GetLatestPose().translation() << std::endl;
  // std::cout << "ESKFLocalTrajectoryBuilder::AddRangeData():: Got matching submap local pose" << std::endl;
  // std::cout << matching_submap_local_pose_.rotation() << std::endl;
  // std::cout << matching_submap_local_pose_.translation() << std::endl;  

  // std::cout << "eskf pose mean:\n";
  // std::cout << eskf_->GetPoseMean() << std::endl;
  // std::cout << eskf_->GetPoseCovariance() << std::endl;

  // Transform latest pose into local, then send to gpf
  transform::Rigid3d latest_pose_in_matching_submap = 
      matching_submap_local_pose_.inverse() * 
      eskf_->GetLatestPose();
  Eigen::Matrix<double, 7, 1> pose_mean_in_matching_submap;
  pose_mean_in_matching_submap.block<3,1>(0,0) = 
    latest_pose_in_matching_submap.translation();
  pose_mean_in_matching_submap(3,0) = 
      latest_pose_in_matching_submap.rotation().w();
  pose_mean_in_matching_submap(4,0) = 
      latest_pose_in_matching_submap.rotation().x();
  pose_mean_in_matching_submap(5,0) = 
      latest_pose_in_matching_submap.rotation().y();
  pose_mean_in_matching_submap(6,0) = 
      latest_pose_in_matching_submap.rotation().z();

  // gpf_->SetPosePrior(eskf_->GetPoseMean(), 
  //                    eskf_->GetPoseCovariance());

  gpf_->SetPosePrior(pose_mean_in_matching_submap, 
                     eskf_->GetPoseCovariance());
  // std::cout << "ESKFLocalTrajectoryBuilder::AddRangeData():: Set prior from GPF" << std::endl;
  // std::cout << eskf_->GetPoseMean() << std::endl;

  gpf_->AddRangeData(time, range_data);
  // std::cout << "ESKFLocalTrajectoryBuilder::AddRangeData():: added range to gpf" << std::endl;
  // std::cout << "gpf measurement:\n " << gpf_->GetMeasurementPose() << std::endl;
  // std::cout << "gpf covariance \n" << gpf_->GetMeasurementCovariance() << std::endl;
  
  // Transform gpf error pose from local into global, then update eskf
  eskf_->Update(gpf_->GetMeasurementPose(), 
                gpf_->GetMeasurementCovariance());
  // std::cout << "ESKF::LatestPose: \n";
  // std::cout << eskf_->GetLatestPose().translation().transpose() << std::endl;
  // std::cout << eskf_->GetLatestPose().rotation().w() << " "
  //           << eskf_->GetLatestPose().rotation().x() << " "
  //           << eskf_->GetLatestPose().rotation().y() << " "
  //           << eskf_->GetLatestPose().rotation().z() << std::endl;

  // sensor::RangeData range_data_filtered;
  // for(const Eigen::Vector3f& hit : range_data.returns) {
  //       const Eigen::Vector3f delta = hit - range_data.origin;
  //       const float range = delta.norm();
  //       if (range >= options_.min_range()) {
  //         if (range <= options_.max_range()) {
  //           range_data_filtered.returns.push_back(hit);
  //         }
  //         else {
  //           // range_data_filtered.misses.push_back(hit);
  //         }
  //       }
  //   }
  // return nullptr;
  
  // std::cout << "ESKFLocalTrajectoryBuilder::AddRangeData():: update eskf with measurement" << std::endl;


  if (num_accumulated_ == 0) {
    first_pose_estimate_ = eskf_->GetLatestPose().cast<float>();
    // std::cout << "first pose: " << first_pose_estimate_.translation().transpose() << std::endl;
    // std::cout << "ESKFLocalTrajectoryBuilder::AddRangeData()::Got first pose estimate" << std::endl;
    accumulated_range_data_ =
        sensor::RangeData{Eigen::Vector3f::Zero(), {}, {}};
  }


  const transform::Rigid3f tracking_delta = 
    first_pose_estimate_.inverse() * 
    eskf_->GetLatestPose().cast<float>();
  // std::cout << "eskf latest pose:" << eskf_->GetLatestPose().translation() << std::endl;
  const sensor::RangeData range_data_in_first_tracking =
      sensor::TransformRangeData(range_data, tracking_delta);
  // std::cout << "tracking delta: " << tracking_delta.translation() << std::endl;
  // std::cout << "range data: " << range_data_in_first_tracking.returns[0] << std::endl;
  // std::cout << "ESKFLocalTrajectoryBuilder::AddRangeData()::Got range data in first tracking frame" << std::endl;
  for (const Eigen::Vector3f& hit : range_data_in_first_tracking.returns) {
    const Eigen::Vector3f delta = hit - range_data_in_first_tracking.origin;
    const float range = delta.norm();
    // float range = sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]);
    // std::cout << range << ", ";
    if (range >= options_.min_range()) {
      if (range <= options_.max_range()) {
        accumulated_range_data_.returns.push_back(hit);
      } else {
        // We insert a ray cropped to 'max_range' as a miss for hits beyond the
        // maximum range. This way the free space up to the maximum range will
        // be updated.
        accumulated_range_data_.misses.push_back(
            range_data_in_first_tracking.origin +
            options_.max_range() / range * delta);
      }
    }
  }


  // std::cout << std::endl;
  ++num_accumulated_;

  // last_pose_estimate_ = {
  //     time, eskf_->GetLatestPose(),
  //     sensor::TransformPointCloud(range_data.returns,
  //                                 eskf_->GetLatestPose().cast<float>())};
  last_pose_estimate_ = {
      time, eskf_->GetLatestPose(),
      sensor::TransformPointCloud(accumulated_range_data_.returns,
                                  first_pose_estimate_)};
  // std::cout << "min_range=" << options_.min_range() << std::endl;
  // std::cout << "max_range=" << options_.max_range() << std::endl;

  // std::cout << "accumulated_range_data_ size=" << accumulated_range_data_.returns.size() << std::endl;
  // std::cout << "accumulated_range_data_ size=" << accumulated_range_data_.misses.size() << std::endl;
  // std::cout << "Adding accumulated_range_data_" << std::endl;
  // std::cout << "accumulated scans: " << num_accumulated_ << std::endl;
  if (num_accumulated_ >= options_.scans_per_accumulation()) {
    num_accumulated_ = 0;
    // Accumulated range data in the latest tracking frame
    return AddAccumulatedRangeData(
        time, sensor::TransformRangeData(accumulated_range_data_,
                                         tracking_delta.inverse()));
  }

  // return AddAccumulatedRangeData(
  //       time, sensor::TransformRangeData(accumulated_range_data_,
  //                                        tracking_delta.inverse()),
  //       sensor::TransformRangeData(range_data_in_first_tracking, tracking_delta.inverse()));  
  return nullptr;
}

/************************************************************/

std::unique_ptr<ESKFLocalTrajectoryBuilder::InsertionResult>
ESKFLocalTrajectoryBuilder::AddRangeDataMini(const common::Time time,
                                         const sensor::RangeData& range_data) {


  const transform::Rigid3f tracking_delta = 
    first_pose_estimate_.inverse() * 
    eskf_->GetLatestPose().cast<float>();
  // std::cout << "eskf latest pose:" << eskf_->GetLatestPose().translation() << std::endl;
  const sensor::RangeData range_data_in_first_tracking =
      sensor::TransformRangeData(range_data, tracking_delta);

  std::unique_ptr<ESKFLocalTrajectoryBuilder::InsertionResult> 
  tmp = AddAccumulatedRangeData(
    time, sensor::TransformRangeData(accumulated_range_data_,
                                         tracking_delta.inverse()));

  return std::unique_ptr<InsertionResult>(new InsertionResult{

      std::make_shared<const mapping::TrajectoryNode::Data>(
          mapping::TrajectoryNode::Data{
              tmp->constant_data->time,
              sensor::Compress(sensor::TransformRangeData(range_data_in_first_tracking,
                                         tracking_delta.inverse())),
              tmp->constant_data->gravity_alignment,
              {},  // 'filtered_point_cloud' is only used in 2D.
              tmp->constant_data->high_resolution_point_cloud,
              tmp->constant_data->low_resolution_point_cloud}),
      tmp->pose_observation, tmp->insertion_submaps});  
}

/************************************************************/

std::unique_ptr<ESKFLocalTrajectoryBuilder::InsertionResult>
ESKFLocalTrajectoryBuilder::AddAccumulatedRangeData(
    const common::Time time, const sensor::RangeData& range_data_in_tracking
    // , const sensor::RangeData& single_range_data) {
    ) {
  const sensor::RangeData filtered_range_data = {
      range_data_in_tracking.origin,
      sensor::VoxelFiltered(range_data_in_tracking.returns,
                            0.1),
      sensor::VoxelFiltered(range_data_in_tracking.misses,
                            0.1)};

  if (filtered_range_data.returns.empty()) {
    LOG(WARNING) << "Dropped empty range data.";
    return nullptr;
  }

  // Accessing the second-last element pointer
  std::shared_ptr<const Submap> matching_submap;
  if (active_submaps_.submaps().size() > 1) {
    matching_submap = active_submaps_.submaps().end()[-2]; 
  } 
  else {
    matching_submap = active_submaps_.submaps().front();
  }

  // Get first tracking frame pose in map frame
  const transform::Rigid3d pose_prediction =
      eskf_->GetLatestPose();
  transform::Rigid3d initial_ceres_pose =
      matching_submap->local_pose().inverse() * pose_prediction;
  sensor::AdaptiveVoxelFilter adaptive_voxel_filter(
      options_.high_resolution_adaptive_voxel_filter_options());
  const sensor::PointCloud filtered_point_cloud_in_tracking =
      adaptive_voxel_filter.Filter(filtered_range_data.returns);
  if (options_.use_online_correlative_scan_matching()) {
    // We take a copy since we use 'initial_ceres_pose' as an output argument.
    const transform::Rigid3d initial_pose = initial_ceres_pose;
    real_time_correlative_scan_matcher_->Match(
        initial_pose, filtered_point_cloud_in_tracking,
        matching_submap->high_resolution_hybrid_grid(), &initial_ceres_pose);
  }

  transform::Rigid3d pose_observation_in_submap;
  ceres::Solver::Summary summary;

  sensor::AdaptiveVoxelFilter low_resolution_adaptive_voxel_filter(
      options_.low_resolution_adaptive_voxel_filter_options());
  const sensor::PointCloud low_resolution_point_cloud_in_tracking =
      low_resolution_adaptive_voxel_filter.Filter(filtered_range_data.returns);
  ceres_scan_matcher_->Match(
      matching_submap->local_pose().inverse() * pose_prediction,
      initial_ceres_pose,
      {{&filtered_point_cloud_in_tracking,
        &matching_submap->high_resolution_hybrid_grid()},
       {&low_resolution_point_cloud_in_tracking,
        &matching_submap->low_resolution_hybrid_grid()}},
      &pose_observation_in_submap, &summary);
  const transform::Rigid3d pose_estimate =
      matching_submap->local_pose() * pose_observation_in_submap;

  // TODO: Add matched pose back to eskf
  const transform::Rigid3d error_pose = pose_prediction.inverse() * pose_estimate;
  Eigen::Matrix<double, 6, 1> mean_meas;
  Eigen::Matrix<double, 6, 6> cov_meas;
  Eigen::Matrix<double, 6, 1> cov_meas_vector;
  mean_meas.block<3,1>(0,0) = error_pose.translation();
  mean_meas.block<3,1>(3,0) = 
      transform::RotationQuaternionToAngleAxisVector(
          error_pose.rotation());
  cov_meas_vector << 0.01, 0.01, 0.01, 0.001, 0.001, 0.001;
  cov_meas = Eigen::Matrix<double, 6,6>(cov_meas_vector.asDiagonal());
  // std::cout << "scan matched pose: " << std::endl;
  // std::cout << mean_meas.transpose() << std::endl;
  // std::cout << cov_meas << std::endl;
  // eskf_->Update(mean_meas, cov_meas);
  
  // std::cout << "after update:\n";
  // std::cout << eskf_->GetPoseMean().transpose() << std::endl;
  // last_pose_estimate_ = {
  //     time, pose_estimate,
  //     sensor::TransformPointCloud(filtered_range_data.returns,
  //                                 pose_estimate.cast<float>())};

  const Eigen::Quaterniond gravity_alignment =
      pose_estimate.rotation().cast<double>();

  // std::cout << "Inserting into submap !!!!!!!!!!!!!!!!" << std::endl;
  return InsertIntoSubmap(time, filtered_range_data, gravity_alignment,
                          // filtered_point_cloud_in_tracking,
                          filtered_range_data.returns,
                          low_resolution_point_cloud_in_tracking,
                          pose_estimate);      
  // return InsertIntoSubmap(time, filtered_range_data, single_range_data, gravity_alignment,
  //                         // filtered_point_cloud_in_tracking,
  //                         filtered_range_data.returns,
  //                         low_resolution_point_cloud_in_tracking,
  //                         pose_estimate);
}                             

void ESKFLocalTrajectoryBuilder::AddOdometerData(
    const sensor::OdometryData& odometry_data) {
}

const mapping::PoseEstimate& ESKFLocalTrajectoryBuilder::pose_estimate() const {
  return last_pose_estimate_;
}

std::shared_ptr<octomap::OcTree> ESKFLocalTrajectoryBuilder::GetMatchingOctomap()  {
  std::shared_ptr<mapping_3d::Submap> 
  matching_submap = active_submaps_.GetMatchingSubmap(); 
  return matching_submap->GetOctomap();
}

std::shared_ptr<Submap>  ESKFLocalTrajectoryBuilder::GetMatchingSubmap() {
  return active_submaps_.GetMatchingSubmap();
}


transform::Rigid3d ESKFLocalTrajectoryBuilder::GetMatchingLocalPose() {
  std::shared_ptr<mapping_3d::Submap> 
  matching_submap = active_submaps_.GetMatchingSubmap(); 
  return matching_submap->local_pose();
}

std::unique_ptr<ESKFLocalTrajectoryBuilder::InsertionResult>
ESKFLocalTrajectoryBuilder::InsertIntoSubmap(
    const common::Time time, const sensor::RangeData& range_data_in_tracking,
    // const sensor::RangeData& single_range_data,
    const Eigen::Quaterniond& gravity_alignment,
    const sensor::PointCloud& high_resolution_point_cloud,
    const sensor::PointCloud& low_resolution_point_cloud,
    const transform::Rigid3d& pose_observation) {
  // Querying the active submaps must be done here before calling
  // InsertRangeData() since the queried values are valid for next insertion.
  std::vector<std::shared_ptr<const Submap>> insertion_submaps;
  if(active_submaps_.submaps().size() < 3) {
    for (const std::shared_ptr<Submap>& submap : active_submaps_.submaps()) {
       insertion_submaps.push_back(submap);
    }
  } else {
    for(int i=2; i>0; i--) {
      insertion_submaps.push_back(active_submaps_.submaps().end()[-i]);
    }
  }
  
  // Insert range data in global frame
  active_submaps_.InsertRangeData(
      sensor::TransformRangeData(range_data_in_tracking,
                                 pose_observation.cast<float>()),
      gravity_alignment);

  return std::unique_ptr<InsertionResult>(new InsertionResult{

      std::make_shared<const mapping::TrajectoryNode::Data>(
          mapping::TrajectoryNode::Data{
              time,
              sensor::Compress(range_data_in_tracking),
              gravity_alignment,
              {},  // 'filtered_point_cloud' is only used in 2D.
              high_resolution_point_cloud,
              low_resolution_point_cloud}),
      pose_observation, std::move(insertion_submaps)  });
      // std::make_shared<const mapping::TrajectoryNode::Data>(
      //       mapping::TrajectoryNode::Data{
      //           time,
      //           sensor::Compress(single_range_data),
      //           gravity_alignment,
      //           {},  // 'filtered_point_cloud' is only used in 2D.
      //           high_resolution_point_cloud,
      //           low_resolution_point_cloud}), num_accumulated_});
}

}  // namespace mapping_3d
}  // namespace cartographer
