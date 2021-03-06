-- Copyright 2016 The Cartographer Authors
--
-- Licensed under the Apache License, Version 2.0 (the "License"),
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
--
--      http://www.apache.org/licenses/LICENSE-2.0
--
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

MAX_3D_RANGE = 20.

TRAJECTORY_BUILDER_3D_ESKF = {
  min_range = 1.,
  max_range = MAX_3D_RANGE,
  scans_per_accumulation = 1,
  voxel_filter_size = 0.15,

  high_resolution_adaptive_voxel_filter = {
    max_length = 2.,
    min_num_points = 150,
    max_range = 15.,
  },

  low_resolution_adaptive_voxel_filter = {
    max_length = 4.,
    min_num_points = 200,
    max_range = MAX_3D_RANGE,
  },

  use_online_correlative_scan_matching = false,
  real_time_correlative_scan_matcher = {
    linear_search_window = 0.15,
    angular_search_window = math.rad(1.),
    translation_delta_cost_weight = 1e-1,
    rotation_delta_cost_weight = 1e-1,
  },

  ceres_scan_matcher = {
    occupied_space_weight_0 = 1.,
    occupied_space_weight_1 = 6.,
    translation_weight = 5.,
    rotation_weight = 4e2,
    only_optimize_yaw = false,
    ceres_solver_options = {
      use_nonmonotonic_steps = false,
      max_num_iterations = 12,
      num_threads = 1,
    },
  },

  motion_filter = {
    max_time_seconds = 0.5,
    max_distance_meters = 0.1,
    max_angle_radians = 0.004,
  },

  imu_gravity_time_constant = 10.,

  submaps = {
    high_resolution = 0.10,
    high_resolution_max_range = 15.,
    low_resolution = 0.45,
    num_range_data = 160,
    range_data_inserter = {
      hit_probability = 0.55,
      miss_probability = 0.49,
      num_free_space_voxels = 2,
    },
    initial_map_file_name = "map.bt",
    initial_map_resolution = 0.2,
    distance_map_range = 20.0,
    using_eskf = true,
  },

  eskf = {
    imu_frequency = 50,
    imu_frame = "microstrain",
    robot_frame = "microstrain",
    imu_enabled = true,
    imu_has_quaternion = true,
    smooth_enabled = false,
    smooth_buffer_size = 1,
    smooth_type = "mean",
    
    initial_x = 0.0,
    initial_y = 0.0,
    initial_z = 0.0,
    initial_roll = 0.0,
    initial_pitch = 0.0,
    initial_yaw = 0.0,

    sigma_accelerometer = 2.0,
    sigma_gyroscope = 0.01,
    sigma_bias_accelerometer = 0.001,
    sigma_bias_gyroscope = 0.001,

    gravity = 9.806,

    initial_bias_accelerometer_x = 0.0,
    initial_bias_accelerometer_y = 0.0,
    initial_bias_accelerometer_z = 0.0,

    accelerometer_queue_size = 3,
    imu_transform = false,
  },
}
