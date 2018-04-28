#include "cartographer/mapping/eskf_options.h"

namespace cartographer {
namespace mapping {

proto::ESKFOptions CreateESKFOptions(
	common::LuaParameterDictionary* const parameter_dictionary) {

	proto::ESKFOptions options;
  	
  	options.set_initial_x(parameter_dictionary->GetDouble("initial_x"));
  	options.set_initial_y(parameter_dictionary->GetDouble("initial_y"));
  	options.set_initial_z(parameter_dictionary->GetDouble("initial_z"));
  	options.set_initial_roll(parameter_dictionary->GetDouble("initial_roll"));
  	options.set_initial_pitch(parameter_dictionary->GetDouble("initial_pitch"));
  	options.set_initial_yaw(parameter_dictionary->GetDouble("initial_yaw"));

  	options.set_imu_frequency(parameter_dictionary->GetDouble("imu_frequency"));
  	options.set_imu_frame(parameter_dictionary->GetString("imu_frame"));
  	options.set_robot_frame(parameter_dictionary->GetString("robot_frame"));
  	options.set_imu_enabled(parameter_dictionary->GetBool("imu_enabled"));
  	options.set_imu_has_quaternion(parameter_dictionary->GetBool("imu_has_quaternion"));
  	options.set_smooth_enabled(parameter_dictionary->GetBool("smooth_enabled"));
  	options.set_smooth_buffer_size(parameter_dictionary->GetInt("smooth_buffer_size"));
  	options.set_smooth_type(parameter_dictionary->GetString("smooth_type"));

  	options.set_sigma_accelerometer(parameter_dictionary->GetDouble("sigma_accelerometer"));
  	options.set_sigma_gyroscope(parameter_dictionary->GetDouble("sigma_gyroscope"));
  	options.set_sigma_bias_accelerometer(parameter_dictionary->GetDouble("sigma_bias_accelerometer"));
  	options.set_sigma_bias_gyroscope(parameter_dictionary->GetDouble("sigma_bias_gyroscope"));

  	options.set_gravity(parameter_dictionary->GetDouble("gravity"));

  	options.set_initial_bias_accelerometer_x(parameter_dictionary->GetDouble("initial_bias_accelerometer_x"));
  	options.set_initial_bias_accelerometer_y(parameter_dictionary->GetDouble("initial_bias_accelerometer_y"));
  	options.set_initial_bias_accelerometer_z(parameter_dictionary->GetDouble("initial_bias_accelerometer_z"));

  	options.set_accelerometer_queue_size(parameter_dictionary->GetInt("accelerometer_queue_size"));
  	options.set_imu_transform(parameter_dictionary->GetBool("imu_transform"));

  	return options;
}
}
}