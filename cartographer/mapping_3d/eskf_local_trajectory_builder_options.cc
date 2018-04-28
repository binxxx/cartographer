
#include "cartographer/mapping_3d/eskf_local_trajectory_builder_options.h"

#include "cartographer/mapping_2d/scan_matching/real_time_correlative_scan_matcher.h"
#include "cartographer/mapping_3d/motion_filter.h"
#include "cartographer/mapping_3d/scan_matching/ceres_scan_matcher.h"
#include "cartographer/mapping_3d/submaps.h"
#include "cartographer/sensor/voxel_filter.h"
#include "cartographer/mapping/eskf_options.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping_3d {

proto::ESKFLocalTrajectoryBuilderOptions CreateESKFLocalTrajectoryBuilderOptions(
    common::LuaParameterDictionary* const parameter_dictionary) {
  proto::ESKFLocalTrajectoryBuilderOptions options;
  options.set_min_range(parameter_dictionary->GetDouble("min_range"));
  options.set_max_range(parameter_dictionary->GetDouble("max_range"));
  options.set_scans_per_accumulation(
      parameter_dictionary->GetInt("scans_per_accumulation"));
  options.set_voxel_filter_size(
      parameter_dictionary->GetDouble("voxel_filter_size"));
  *options.mutable_high_resolution_adaptive_voxel_filter_options() =
      sensor::CreateAdaptiveVoxelFilterOptions(
          parameter_dictionary
              ->GetDictionary("high_resolution_adaptive_voxel_filter")
              .get());
  *options.mutable_low_resolution_adaptive_voxel_filter_options() =
      sensor::CreateAdaptiveVoxelFilterOptions(
          parameter_dictionary
              ->GetDictionary("low_resolution_adaptive_voxel_filter")
              .get());
  options.set_use_online_correlative_scan_matching(
      parameter_dictionary->GetBool("use_online_correlative_scan_matching"));
  *options.mutable_real_time_correlative_scan_matcher_options() =
      mapping_2d::scan_matching::CreateRealTimeCorrelativeScanMatcherOptions(
          parameter_dictionary
              ->GetDictionary("real_time_correlative_scan_matcher")
              .get());
  *options.mutable_ceres_scan_matcher_options() =
      scan_matching::CreateCeresScanMatcherOptions(
          parameter_dictionary->GetDictionary("ceres_scan_matcher").get());
  *options.mutable_motion_filter_options() = CreateMotionFilterOptions(
      parameter_dictionary->GetDictionary("motion_filter").get());
  options.set_imu_gravity_time_constant(
      parameter_dictionary->GetDouble("imu_gravity_time_constant"));
  *options.mutable_submaps_options() = mapping_3d::CreateSubmapsOptions(
      parameter_dictionary->GetDictionary("submaps").get());
  *options.mutable_eskf_options() = mapping::CreateESKFOptions(
      parameter_dictionary->GetDictionary("eskf").get());
  return options;
}

}  // namespace mapping_3d
}  // namespace cartographer
