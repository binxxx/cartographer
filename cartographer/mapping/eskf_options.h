#ifndef CARTOGRAPHER_MAPPING_ESKF_OPTIONS_H_
#define CARTOGRAPHER_MAPPING_ESKF_OPTIONS_H_

#include "cartographer/common/lua_parameter_dictionary.h"
#include "cartographer/mapping/proto/eskf_options.pb.h"

namespace cartographer {
namespace mapping {

proto::ESKFOptions CreateESKFOptions(
	common::LuaParameterDictionary* parameter_dictionary);
}
}

#endif