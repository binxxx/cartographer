/*
 * Copyright 2016 The Cartographer Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef CARTOGRAPHER_MAPPING_TRAJECTORY_BUILDER_H_
#define CARTOGRAPHER_MAPPING_TRAJECTORY_BUILDER_H_

#include <functional>
#include <memory>
#include <string>

#include "cartographer/common/lua_parameter_dictionary.h"
#include "cartographer/common/make_unique.h"
#include "cartographer/common/port.h"
#include "cartographer/common/time.h"
#include "cartographer/mapping/pose_estimate.h"
#include "cartographer/mapping/proto/trajectory_builder_options.pb.h"
#include "cartographer/mapping/submaps.h"
#include "cartographer/mapping_3d/submaps.h"
#include "cartographer/sensor/data.h"
#include "cartographer/sensor/point_cloud.h"
#include "cartographer/sensor/range_data.h"

namespace cartographer {
namespace mapping {

proto::TrajectoryBuilderOptions CreateTrajectoryBuilderOptions(
    common::LuaParameterDictionary* const parameter_dictionary);

// This interface is used for both 2D and 3D SLAM.
class TrajectoryBuilder {
 public:
  using PoseEstimate = mapping::PoseEstimate;

  TrajectoryBuilder() {}
  virtual ~TrajectoryBuilder() {}

  TrajectoryBuilder(const TrajectoryBuilder&) = delete;
  TrajectoryBuilder& operator=(const TrajectoryBuilder&) = delete;

  virtual const PoseEstimate& pose_estimate() const = 0;
  virtual std::shared_ptr<octomap::OcTree> GetMatchingOctomap() = 0;
  virtual std::shared_ptr<mapping::Submap> GetMatchingSubmap() = 0;
  virtual transform::Rigid3d GetMatchingLocalPose() = 0;
  virtual void AddSensorData(const string& sensor_id,
                             std::unique_ptr<sensor::Data> data) = 0;

  void AddRangefinderData(const string& sensor_id, common::Time time,
                          const Eigen::Vector3f& origin,
                          const sensor::PointCloud& ranges) {
    AddSensorData(sensor_id,
                  common::make_unique<sensor::DispatchableRangefinderData>(
                      time, origin, ranges));
  }

  void AddImuData(const string& sensor_id, common::Time time,
                  const Eigen::Vector3d& linear_acceleration,
                  const Eigen::Vector3d& angular_velocity,
                  const Eigen::Quaterniond& orientation) {
    AddSensorData(sensor_id, sensor::MakeDispatchable(sensor::ImuData{
                                 time, linear_acceleration, angular_velocity,
                                 orientation}));
  }

  void AddOdometerData(const string& sensor_id, common::Time time,
                       const transform::Rigid3d& odometer_pose) {
    AddSensorData(sensor_id, sensor::MakeDispatchable(
                                 sensor::OdometryData{time, odometer_pose}));
  }

  void AddFixedFramePoseData(const string& sensor_id, common::Time time,
                             const transform::Rigid3d& fixed_frame_pose) {
    AddSensorData(sensor_id,
                  sensor::MakeDispatchable(
                      sensor::FixedFramePoseData{time, fixed_frame_pose}));
  }
};

}  // namespace mapping
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_TRAJECTORY_BUILDER_H_
