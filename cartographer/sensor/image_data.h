#ifndef CARTOGRAPHER_SENSOR_IMAGE_DATA_H_
#define CARTOGRAPHER_SENSOR_IMAGE_DATA_H_

// #include "Eigen/Core"
#include "cartographer/common/time.h"
#include "cartographer/sensor/proto/sensor.pb.h"

namespace cartographer {
namespace sensor {

struct ImageData {
  common::Time time;
  size_t width;
  size_t height;
  std::string encoding;

  // pixel data
  std::vector<uint8_t> data;
};

// Converts 'image_data' to a proto::ImageData.
proto::ImageData ToProto(const ImageData& image_data);

// Converts 'proto' to an ImageData.
ImageData FromProto(const proto::ImageData& proto);

}  // namespace sensor
}  // namespace cartographer

#endif  // CARTOGRAPHER_SENSOR_IMAGE_DATA_H_
