
#ifndef CARTOGRAPHER_MAPPING_PARTICLES_H_
#define CARTOGRAPHER_MAPPING_PARTICLES_H_

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "cartographer/mapping/multivariate_normal.h"



#include "cartographer/common/time.h"
#include "cartographer/transform/rigid_transform.h"
#include "cartographer/mapping_3d/submaps.h"
#include "cartographer/sensor/range_data.h"

namespace cartographer {
namespace mapping {


struct Particle {
    Eigen::Vector3d translation;
    Eigen::Quaterniond rotation;
    Eigen::Vector3d angle_axis;
    double weight;
    Particle() {
        translation.setZero();
        rotation.setIdentity();
        weight = 0.0;
    }
};

class Particles {
public:

    Particles(int set_size, double ray_sigma);
    ~Particles() {}

    void SetMean(const Eigen::Matrix<double, 7, 1> &mean);
    void SetCovariance(const Eigen::Matrix<double, 6, 6> &cov);
    void SetCloud(const std::shared_ptr<sensor::PointCloud> cloud_ptr);
    void SetMatchingDistanceMap(const std::shared_ptr<DynamicEDTOctomap> distmap_grid);
    void Propagate(Eigen::Matrix<double, 6, 1> &mean_prior,
                   Eigen::Matrix<double, 6, 6> &cov_prior,
                   Eigen::Matrix<double, 6, 1> &mean_posterior,
                   Eigen::Matrix<double, 6, 6> &cov_posterior);

private:
    void DrawParticles();
    void AssignWeightToParticle(Particle &p, sensor::PointCloud &cloud);
    void AssignWeightsToParticles();
    void ComputePosterior();

    std::vector<Particle> pset_;
    std::vector<Particle> d_pset_;

    Eigen::Matrix<double, 7, 1> mean_prior_;
    Eigen::Matrix<double, 6, 1> mean_posterior_;

    Eigen::Matrix<double, 6, 1> d_mean_prior_;
    Eigen::Matrix<double, 6, 6> d_cov_prior_;
    Eigen::Matrix<double, 6, 1> d_mean_sample_;
    Eigen::Matrix<double, 6, 6> d_cov_sample_;
    Eigen::Matrix<double, 6, 1> d_mean_posterior_;
    Eigen::Matrix<double, 6, 6> d_cov_posterior_;

    std::shared_ptr<sensor::PointCloud> cloud_ptr_;
    std::shared_ptr<DynamicEDTOctomap> distmap_grid_;

    int set_size_;
    double ray_sigma_;

};
}  // namespace mapping
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_PARTICLES_H_
