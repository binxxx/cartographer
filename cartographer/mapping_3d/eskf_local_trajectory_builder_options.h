
#ifndef CARTOGRAPHER_MAPPING_3D_ESKF_LOCAL_TRAJECTORY_BUILDER_OPTIONS_H_
#define CARTOGRAPHER_MAPPING_3D_ESKF_LOCAL_TRAJECTORY_BUILDER_OPTIONS_H_

#include "cartographer/common/lua_parameter_dictionary.h"
#include "cartographer/mapping_3d/proto/eskf_local_trajectory_builder_options.pb.h"

namespace cartographer {
namespace mapping_3d {

proto::ESKFLocalTrajectoryBuilderOptions CreateESKFLocalTrajectoryBuilderOptions(
    common::LuaParameterDictionary* parameter_dictionary);

}  // namespace mapping_3d
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_3D_LOCAL_TRAJECTORY_BUILDER_OPTIONS_H_
