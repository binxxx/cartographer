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

#include "cartographer/mapping_3d/submaps.h"

#include <cmath>
#include <limits>

#include "cartographer/common/math.h"
#include "cartographer/sensor/range_data.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping_3d {

namespace {

struct PixelData {
  int min_z = INT_MAX;
  int max_z = INT_MIN;
  int count = 0;
  float probability_sum = 0.f;
  float max_probability = 0.5f;
};

// Filters 'range_data', retaining only the returns that have no more than
// 'max_range' distance from the origin. Removes misses and reflectivity
// information.
sensor::RangeData FilterRangeDataByMaxRange(const sensor::RangeData& range_data,
                                            const float max_range) {
  sensor::RangeData result{range_data.origin, {}, {}};
  for (const Eigen::Vector3f& hit : range_data.returns) {
    if ((hit - range_data.origin).norm() <= max_range) {
      result.returns.push_back(hit);
    }
  }
  return result;
}

std::vector<PixelData> AccumulatePixelData(
    const int width, const int height, const Eigen::Array2i& min_index,
    const Eigen::Array2i& max_index,
    const std::vector<Eigen::Array4i>& voxel_indices_and_probabilities) {
  std::vector<PixelData> accumulated_pixel_data(width * height);
  for (const Eigen::Array4i& voxel_index_and_probability :
       voxel_indices_and_probabilities) {
    const Eigen::Array2i pixel_index = voxel_index_and_probability.head<2>();
    if ((pixel_index < min_index).any() || (pixel_index > max_index).any()) {
      // Out of bounds. This could happen because of floating point inaccuracy.
      continue;
    }
    const int x = max_index.x() - pixel_index[0];
    const int y = max_index.y() - pixel_index[1];
    PixelData& pixel = accumulated_pixel_data[x * width + y];
    ++pixel.count;
    pixel.min_z = std::min(pixel.min_z, voxel_index_and_probability[2]);
    pixel.max_z = std::max(pixel.max_z, voxel_index_and_probability[2]);
    const float probability =
        mapping::ValueToProbability(voxel_index_and_probability[3]);
    pixel.probability_sum += probability;
    pixel.max_probability = std::max(pixel.max_probability, probability);
  }
  return accumulated_pixel_data;
}

// The first three entries of each returned value are a cell_index and the
// last is the corresponding probability value. We batch them together like
// this to only have one vector and have better cache locality.
std::vector<Eigen::Array4i> ExtractVoxelData(
    const HybridGrid& hybrid_grid, const transform::Rigid3f& transform,
    Eigen::Array2i* min_index, Eigen::Array2i* max_index) {
  std::vector<Eigen::Array4i> voxel_indices_and_probabilities;
  const float resolution_inverse = 1.f / hybrid_grid.resolution();

  constexpr float kXrayObstructedCellProbabilityLimit = 0.501f;
  for (auto it = HybridGrid::Iterator(hybrid_grid); !it.Done(); it.Next()) {
    const uint16 probability_value = it.GetValue();
    const float probability = mapping::ValueToProbability(probability_value);
    if (probability < kXrayObstructedCellProbabilityLimit) {
      // We ignore non-obstructed cells.
      continue;
    }

    const Eigen::Vector3f cell_center_submap =
        hybrid_grid.GetCenterOfCell(it.GetCellIndex());
    const Eigen::Vector3f cell_center_global = transform * cell_center_submap;
    const Eigen::Array4i voxel_index_and_probability(
        common::RoundToInt(cell_center_global.x() * resolution_inverse),
        common::RoundToInt(cell_center_global.y() * resolution_inverse),
        common::RoundToInt(cell_center_global.z() * resolution_inverse),
        probability_value);

    voxel_indices_and_probabilities.push_back(voxel_index_and_probability);
    const Eigen::Array2i pixel_index = voxel_index_and_probability.head<2>();
    *min_index = min_index->cwiseMin(pixel_index);
    *max_index = max_index->cwiseMax(pixel_index);
  }
  return voxel_indices_and_probabilities;
}

// Builds texture data containing interleaved value and alpha for the
// visualization from 'accumulated_pixel_data'.
string ComputePixelValues(
    const std::vector<PixelData>& accumulated_pixel_data) {
  string cell_data;
  cell_data.reserve(2 * accumulated_pixel_data.size());
  constexpr float kMinZDifference = 3.f;
  constexpr float kFreeSpaceWeight = 0.15f;
  for (const PixelData& pixel : accumulated_pixel_data) {
    // TODO(whess): Take into account submap rotation.
    // TODO(whess): Document the approach and make it more independent from the
    // chosen resolution.
    const float z_difference = pixel.count > 0 ? pixel.max_z - pixel.min_z : 0;
    if (z_difference < kMinZDifference) {
      cell_data.push_back(0);  // value
      cell_data.push_back(0);  // alpha
      continue;
    }
    const float free_space = std::max(z_difference - pixel.count, 0.f);
    const float free_space_weight = kFreeSpaceWeight * free_space;
    const float total_weight = pixel.count + free_space_weight;
    const float free_space_probability = 1.f - pixel.max_probability;
    const float average_probability = mapping::ClampProbability(
        (pixel.probability_sum + free_space_probability * free_space_weight) /
        total_weight);
    const int delta =
        128 - mapping::ProbabilityToLogOddsInteger(average_probability);
    const uint8 alpha = delta > 0 ? 0 : -delta;
    const uint8 value = delta > 0 ? delta : 0;
    cell_data.push_back(value);                         // value
    cell_data.push_back((value || alpha) ? alpha : 1);  // alpha
  }
  return cell_data;
}

}  // namespace

proto::SubmapsOptions CreateSubmapsOptions(
    common::LuaParameterDictionary* parameter_dictionary) {
  proto::SubmapsOptions options;
  options.set_high_resolution(
      parameter_dictionary->GetDouble("high_resolution"));
  options.set_high_resolution_max_range(
      parameter_dictionary->GetDouble("high_resolution_max_range"));
  options.set_low_resolution(parameter_dictionary->GetDouble("low_resolution"));
  options.set_num_range_data(
      parameter_dictionary->GetNonNegativeInt("num_range_data"));
  *options.mutable_range_data_inserter_options() =
      CreateRangeDataInserterOptions(
          parameter_dictionary->GetDictionary("range_data_inserter").get());
  options.set_initial_map_file_name(
      parameter_dictionary->GetString("initial_map_file_name"));
  options.set_initial_map_resolution(
      parameter_dictionary->GetDouble("initial_map_resolution"));
  options.set_distance_map_range(
      parameter_dictionary->GetDouble("distance_map_range"));
  options.set_using_eskf(
      parameter_dictionary->GetBool("using_eskf"));
  CHECK_GT(options.num_range_data(), 0);
  return options;
}

Submap::Submap(const int max_num_range_data,
               const float high_resolution, const float low_resolution,
               const transform::Rigid3d& local_pose, 
               common::ThreadPool* thread_pool,
               const double initial_map_resolution,
               const double distance_map_range,
               const bool using_octomap)
    : mapping::Submap(local_pose),
      high_resolution_hybrid_grid_(high_resolution),
      low_resolution_hybrid_grid_(low_resolution),
      thread_pool_(thread_pool),
      max_num_range_data_(max_num_range_data),
      distance_map_range_(distance_map_range),
      using_octomap_(using_octomap) {

        // std::cout << "using octomap: " << using_octomap_ << std::endl;
        if(using_octomap_) {
          // octomap grid shares the high resolution
          octomap_grid_ = std::shared_ptr<octomap::OcTree>(
            new octomap::OcTree(initial_map_resolution));
        }
      }

Submap::Submap(const mapping::proto::Submap3D& proto)
    : mapping::Submap(transform::ToRigid3(proto.local_pose())),
      high_resolution_hybrid_grid_(proto.high_resolution_hybrid_grid()),
      low_resolution_hybrid_grid_(proto.low_resolution_hybrid_grid()) {
  SetNumRangeData(proto.num_range_data());
  finished_ = proto.finished();
}

void Submap::ToProto(mapping::proto::Submap* const proto) const {
  auto* const submap_3d = proto->mutable_submap_3d();
  *submap_3d->mutable_local_pose() = transform::ToProto(local_pose());
  submap_3d->set_num_range_data(num_range_data());
  submap_3d->set_finished(finished_);
  *submap_3d->mutable_high_resolution_hybrid_grid() =
      high_resolution_hybrid_grid().ToProto();
  *submap_3d->mutable_low_resolution_hybrid_grid() =
      low_resolution_hybrid_grid().ToProto();
}

void Submap::ToResponseProto(
    const transform::Rigid3d& global_submap_pose,
    mapping::proto::SubmapQuery::Response* const response) const {
  response->set_submap_version(num_range_data());
  // Generate an X-ray view through the 'hybrid_grid', aligned to the xy-plane
  // in the global map frame.
  const float resolution = high_resolution_hybrid_grid_.resolution();
  response->set_resolution(resolution);

  // Compute a bounding box for the texture.
  Eigen::Array2i min_index(INT_MAX, INT_MAX);
  Eigen::Array2i max_index(INT_MIN, INT_MIN);
  const std::vector<Eigen::Array4i> voxel_indices_and_probabilities =
      ExtractVoxelData(high_resolution_hybrid_grid_,
                       global_submap_pose.cast<float>(), &min_index,
                       &max_index);

  const int width = max_index.y() - min_index.y() + 1;
  const int height = max_index.x() - min_index.x() + 1;
  response->set_width(width);
  response->set_height(height);

  const std::vector<PixelData> accumulated_pixel_data = AccumulatePixelData(
      width, height, min_index, max_index, voxel_indices_and_probabilities);
  const string cell_data = ComputePixelValues(accumulated_pixel_data);

  common::FastGzipString(cell_data, response->mutable_cells());
  *response->mutable_slice_pose() = transform::ToProto(
      global_submap_pose.inverse() *
      transform::Rigid3d::Translation(Eigen::Vector3d(
          max_index.x() * resolution, max_index.y() * resolution,
          global_submap_pose.translation().z())));
}

std::shared_ptr<DynamicEDTOctomap> Submap::GetDistanceMap() {
  CHECK(distmap_grid_ != nullptr);
  return distmap_grid_;
}

std::shared_ptr<octomap::OcTree> Submap::GetOctomap() {
  return octomap_grid_;
}

void Submap::InsertRangeData(const sensor::RangeData& range_data,
                             const RangeDataInserter& range_data_inserter,
                             const int high_resolution_max_range) {
  if (finished_) return;
  if (distance_transformed_) return;
  // std::cout << "Submap: InsertRangeData" << std::endl;
  // transform point cloud from global (map) frame to local (submap) frame
  const sensor::RangeData transformed_range_data = sensor::TransformRangeData(
      range_data, local_pose().inverse().cast<float>());
  
  if (using_octomap_) {
    // Schedule a thread to handle octomap and distance map udpate
    // std::cout << "Scheduling distance map update... \n";
    thread_pool_->Schedule([=]() {
      common::MutexLocker locker(&mutex_);
      octomap::Pointcloud octomap_cloud;
      ToOctomapPointCloud(transformed_range_data, octomap_cloud);
      octomap::point3d sensor_origin(transformed_range_data.origin[0],
                                     transformed_range_data.origin[1],
                                     transformed_range_data.origin[2]);
      // octomap_clouds_.push_back(octomap_cloud);
      // sensor_origins_.push_back(sensor_origin);

      octomap_grid_->insertPointCloud(octomap_cloud,
                                      sensor_origin,
                                      octomap::pose6d(0.0,0.0,0.0, 
                                                    0.0,0.0,0.0));
      octomap_grid_->updateInnerOccupancy();
      ComputeDistanceMap();
      std::cout << "distance map update done... \n";

    });
  }
  range_data_inserter.Insert(
      FilterRangeDataByMaxRange(transformed_range_data,
                                high_resolution_max_range),
                                &high_resolution_hybrid_grid_);
  range_data_inserter.Insert(transformed_range_data,
                             &low_resolution_hybrid_grid_);
  SetNumRangeData(num_range_data() + 1);
  
  // high_resolution_max_range_ = (high_resolution_max_range_ > 
  //                               high_resolution_max_range) ? 
  //                              high_resolution_max_range_: 
  //                              high_resolution_max_range;


}
void Submap::ComputeOctomap() {
  for(size_t cloud_index = 0; cloud_index < octomap_clouds_.size();
      ++cloud_index) {
    // Insert range data into octomap_grid_ 
    octomap_grid_->insertPointCloud(octomap_clouds_[cloud_index], 
                                    sensor_origins_[cloud_index],
                                    octomap::pose6d(0.0,0.0,0.0, 
                                                    0.0,0.0,0.0));
  }
  // std::cout << "number of clouds: " << octomap_clouds_.size() << std::endl;
  octomap_grid_->updateInnerOccupancy();
  // std::cout << "update inner occupancy" << std::endl;    
  // octomap_grid_->writeBinaryConst(std::string("/home/aeroscout/carto.bt"));
  // std::cout << "local pose:\n";
  // std::cout << local_pose().translation() << std::endl;
  // std::cout << local_pose().rotation().w() << "\n"
  //           << local_pose().rotation().x() << "\n"
  //           << local_pose().rotation().y() << "\n"
  //           << local_pose().rotation().z() << std::endl;

}

void Submap::LoadOctomap(const std::string map_file_name) {
  // std::string map_file_name = 
  //     std::string("/home/aeroscout/Documents/tunnel_sim/tunnel_sim.bt");
  std::fstream map_file(map_file_name.c_str(), 
                        std::ios_base::binary | 
                        std::ios_base::in);
  if (map_file.is_open()) {
      octomap_grid_->readBinary(map_file_name);
      if(octomap_grid_.get()) {
          if(!octomap_grid_ || octomap_grid_->size() <= 1) {
              exit(-1);
          }
      }
      map_file.close();
  } else {
    LOG(ERROR) << "Octomap file " << map_file_name << " not open.";
    exit(-1);
  }
}

void Submap::ComputeDistanceMap() {
  
  if (distance_transformed_) return;

  if (distmap_grid_ == nullptr) {
    double max_distance = 2.0;
    double min_x, min_y, min_z;
    double max_x, max_y, max_z;
    octomap_grid_->getMetricMin(min_x, min_y, min_z);
    octomap_grid_->getMetricMax(max_x, max_y, max_z);
    min_x = min_y = -distance_map_range_;
    max_x = max_y = distance_map_range_;
    bounding_box_ = {min_x, min_y, min_z, max_x, max_y, max_z};

    octomap::point3d bound_min(min_x - max_distance,
                               min_y - max_distance,
                               min_z - max_distance);
    octomap::point3d bound_max(max_x + max_distance,
                               max_y + max_distance,
                               max_z + max_distance);
    std::cout << "distmap min size: [" << bound_min(0) << ", "
                                       << bound_min(1) << ", "
                                       << bound_min(2) << "]" << std::endl;
    std::cout << "distmap max size: [" << bound_max(0) << ", "
                                       << bound_max(1) << ", "
                                       << bound_max(2) << "]" << std::endl;

    distmap_grid_ = std::shared_ptr<DynamicEDTOctomap>(
                        new DynamicEDTOctomap(max_distance,
                                              octomap_grid_.get(),
                                              bound_min,
                                              bound_max,
                                              false));

  }
  // std::cout<< "distmap defined" << std::endl;
  distmap_grid_->update();
  if (finished_) {
    distance_transformed_ = true;
    // std::cout << "labeling submap transformed.\n";
  }
  // std::cout << "update distance map" << std::endl;
  // CHECK(!distance_transformed_);
  // distance_transformed_ = true;
  // std::cout << "Distance map updated!";
  // std::cout << "Distance map local pose:" << std::endl;  
  // std::cout << local_pose().translation() << std::endl;
  // std::cout << local_pose().rotation().w() << " "
  //           << local_pose().rotation().x() << " "
  //           << local_pose().rotation().y() << " "
  //           << local_pose().rotation().z() << std::endl;
}

std::vector<double> Submap::GetDistanceMapBoundingBox() {
  return bounding_box_;
}
void Submap::ToOctomapPointCloud(const sensor::RangeData& range_data,
                                 octomap::Pointcloud& octomap_cloud) {
  for(auto point : range_data.returns) {
    // TODO: Check range
    octomap_cloud.push_back(point[0], point[1], point[2]);
  }
}

void Submap::Finish() {
  CHECK(!finished_);
  finished_ = true;
}

ActiveSubmaps::ActiveSubmaps(const proto::SubmapsOptions& options,
                             common::ThreadPool* thread_pool)
    : options_(options),
      thread_pool_(thread_pool),
      range_data_inserter_(options.range_data_inserter_options()) {
  // We always want to have at least one submap which we can return and will
  // create it at the origin in absence of a better choice.
  //
  // TODO(whess): Start with no submaps, so that all of them can be
  // approximately gravity aligned.
  // std::cout << "ActiveSubmaps::ActiveSubmaps()::Initializing" << std::endl;
  transform::Rigid3d initial_matching_transform(
        Eigen::Vector3d(0.0,0.0,0.0), 
        Eigen::Quaterniond(1.0,0.0,0.0,0.0));
  if (options_.using_eskf()) {
    matching_submap_ = std::shared_ptr<Submap>(
                        new Submap(options_.num_range_data(),
                                   options_.high_resolution(),
                                   options_.low_resolution(),
                                   initial_matching_transform,
                                   thread_pool_,
                                   options_.initial_map_resolution(),
                                   options_.distance_map_range(),
                                   options_.using_eskf()));
    matching_submap_->LoadOctomap(options_.initial_map_file_name());
    matching_submap_->ComputeDistanceMap();
  }
  AddSubmap(initial_matching_transform);
}

std::vector<std::shared_ptr<Submap>> ActiveSubmaps::submaps() const {
  return submaps_;
}

int ActiveSubmaps::matching_index() const { return matching_submap_index_; }

std::shared_ptr<Submap> ActiveSubmaps::GetMatchingSubmap() {
  // CHECK(matching_submap_ != nullptr);
  
  return matching_submap_;
}

transform::Rigid3d ActiveSubmaps::GetMatchingSubmapLocalPose() {
  if (matching_submap_ != nullptr) {
    return matching_submap_->local_pose();
  } else {
    return transform::Rigid3d::Identity();
  }
}

std::vector<double> ActiveSubmaps::GetMatchingSubmapBoundingBox() {
  if (matching_submap_ != nullptr) {
    return matching_submap_->GetDistanceMapBoundingBox();
  } else {
    std::vector<double> empty_vector;
    return empty_vector;
  }

}
void ActiveSubmaps::InsertRangeData(
    const sensor::RangeData& range_data,
    const Eigen::Quaterniond& gravity_alignment) {

  int cnt = 0;
  for (auto& submap : submaps_) {
    if (!submap->distance_transformed()) {
      // std::cout << cnt << "/" << submaps_.size() << " is not transformed.\n";
      submap->InsertRangeData(range_data, range_data_inserter_,
                              options_.high_resolution_max_range());
    }
    cnt++;
  }

  // std::cout << "active submap size: " << submaps_.size() << std::endl;
  // std::cout << "active submaps end: num_range_data " 
            // << submaps_.back()->num_range_data() << std::endl;
  if (submaps_.back()->num_range_data() == options_.num_range_data()) {
    AddSubmap(transform::Rigid3d(range_data.origin.cast<double>(),
                                 // gravity_alignment));
                                  Eigen::Quaterniond::Identity()));
  }
}

void ActiveSubmaps::AddSubmap(const transform::Rigid3d& local_pose) {
  common::MutexLocker locker(&mutex_);
  if (submaps_.size() > 1) {
    CHECK_EQ(submaps_.size(), 2);
    submaps_.front()->Finish();
    
    if (options_.using_eskf()) {
      // Schedule a thread to handle octomap and distance map udpate
      std::cout << "Scheduling final distance map update..." << std::endl;
      thread_pool_->Schedule([=]() {

        // Wait for all distance updates done then close the submap
        while (true) {
          common::MutexLocker locker(&mutex_);

          if (submaps_.front()->distance_transformed()) {
            submaps_.erase(submaps_.begin());
            std::cout << "old active submap is deleted" << std::endl;

            // Set matching submap to new front
            matching_submap_ = submaps_.front();

            CHECK(matching_submap_ != nullptr);
            CHECK(matching_submap_->GetDistanceMap() != nullptr);
            break;
          } 
          // std::cout << "while loop: waiting ...\n";  
        }
      });
    } else {
      submaps_.erase(submaps_.begin());
    }
    ++matching_submap_index_;
  }
  // std::cout << "using eskf" << options_.using_eskf() << std::endl;
  submaps_.emplace_back(new Submap(options_.num_range_data(),
                                   options_.high_resolution(),
                                   options_.low_resolution(), local_pose, 
                                   thread_pool_,
                                   options_.initial_map_resolution(),
                                   options_.distance_map_range(),
                                   options_.using_eskf()));
  LOG(INFO) << "Added submap " << matching_submap_index_ + submaps_.size();
}

}  // namespace mapping_3d
}  // namespace cartographer
