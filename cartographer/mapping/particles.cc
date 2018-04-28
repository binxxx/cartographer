#include "cartographer/mapping/particles.h"
#include "cartographer/mapping/eskf.h"

namespace cartographer {
namespace mapping {

inline double log_likelihood(double x, double sigma) {
    return -0.91893853320467274178 - log(sigma) - 0.5 * x * x / (sigma * sigma);
}
static MultivariateNormal<double, 6> mvn(Eigen::MatrixXd::Zero(6,1),
                                              Eigen::MatrixXd::Identity(6,6));

Particles::Particles(int set_size, double ray_sigma) : 
		set_size_(set_size),
		ray_sigma_(ray_sigma) {

    mean_prior_.setZero();
    mean_posterior_.setZero();
    d_mean_prior_.setZero();
    d_mean_sample_.setZero();
    d_mean_posterior_.setZero();
    d_cov_prior_.setZero();
    d_cov_sample_.setZero();
    d_cov_posterior_.setZero();
    pset_.resize(set_size_);
    d_pset_.resize(set_size_);
    

}

void Particles::SetMean(const Eigen::Matrix<double, 7, 1> &mean) {
    mean_prior_ = mean;
}

void Particles::SetCovariance(const Eigen::Matrix<double, 6, 6> &cov) {
    d_cov_prior_ = cov;
}

void Particles::SetCloud(const std::shared_ptr<sensor::PointCloud> cloud_ptr) {
    cloud_ptr_ = cloud_ptr;
}

void Particles::SetMatchingDistanceMap(const std::shared_ptr<DynamicEDTOctomap> distmap_grid) {
	distmap_grid_ = distmap_grid;
}


void Particles::DrawParticles() {

    mvn.setMean(d_mean_prior_);
    mvn.setCovar(d_cov_prior_);

    for(int i=0; i<set_size_; i++) {
        // random sample error states
        Eigen::Matrix<double, 6, 1> twist;
        mvn.nextSample(twist);
        d_pset_[i].translation = twist.block<3,1>(0,0);
        d_pset_[i].angle_axis = twist.block<3,1>(3,0);

        d_pset_[i].weight = log(1.0/set_size_);
        pset_[i].weight = d_pset_[i].weight;

        // recover nominal states
        pset_[i].translation = mean_prior_.block<3,1>(0,0) + d_pset_[i].translation;

        // recover nominal states: rotation
        Eigen::Matrix3d d_rotation = angle_axis_to_rotation_matrix(d_pset_[i].angle_axis);
        pset_[i].rotation = Eigen::Quaterniond(mean_prior_[3],mean_prior_[4],
                mean_prior_[5],mean_prior_[6]) * Eigen::Quaterniond(d_rotation);
    }

    d_mean_sample_.setZero();
    d_cov_sample_.setZero();

    for(int i=0; i<set_size_; i++) {
        d_mean_sample_.block<3,1>(0,0) += d_pset_[i].translation / set_size_;
        d_mean_sample_.block<3,1>(3,0) += d_pset_[i].angle_axis / set_size_;
    }
    for(int i=0; i<set_size_; i++) {
        Eigen::Matrix<double, 6, 1> twist;
        twist << d_pset_[i].translation - d_mean_sample_.block<3,1>(0,0),
                 d_pset_[i].angle_axis - d_mean_sample_.block<3,1>(3,0);
        d_cov_sample_ += twist*twist.transpose() / set_size_;
    }
}

void Particles::AssignWeightsToParticles() {
//#pragma omp parallel for
    CHECK(cloud_ptr_ != nullptr);
    for(int i=0; i<set_size_; i++) {
        // reproject cloud on to each particle
        sensor::PointCloud cloud_transformed = 
        	sensor::TransformPointCloud(
        		*cloud_ptr_,
        		transform::Rigid3d(pset_[i].translation, 
        						   pset_[i].rotation).cast<float>());
        // weight particle
        AssignWeightToParticle(pset_[i], cloud_transformed);
    }
    // stabilize weights, offset weight values to [-200.0, 0.0] range
    double max_weight(-INFINITY);
    for(int i=0; i<set_size_; i++) {
        if(pset_[i].weight > max_weight) max_weight = pset_[i].weight;
    }
    for(int i=0; i<set_size_; i++) {
        pset_[i].weight -= max_weight;
        if(pset_[i].weight < -200.0) pset_[i].weight = -200.0;
    }

    // normalize weight
    double weight_sum = 0.0;
    for(int i=0; i<set_size_; i++) {
        weight_sum += exp(pset_[i].weight);
    }
    double log_weight_sum = log(weight_sum);
    for(int i=0; i<set_size_; i++) {
        pset_[i].weight -= log_weight_sum;
        d_pset_[i].weight = pset_[i].weight;
    }

}

void Particles::AssignWeightToParticle(Particle &p, sensor::PointCloud &cloud) {
    std::vector<double> weight;
    weight.resize(cloud.size());

    // std::cout << "weighting single particles" << std::endl;
    CHECK(distmap_grid_ != nullptr);
//#pragma omp parallel for
    for(size_t i=0; i<cloud.size(); i++) {
        // the end point of one ray
        octomap::point3d end_pnt(cloud[i][0], cloud[i][1], cloud[i][2]);

        // look up the distance to nearest obstacle
        double dist = distmap_grid_->getDistance(end_pnt);

        // find weight through normal distribution
        // TODO(Weikun): do we really need gridmask?
        // char grid_flag = distmap_grid_->get_gridmask(end_pnt);
        // if(grid_flag != 2) {
        //     if(dist >= 0.0 && dist <= 2.0*ray_sigma_) {
        //         weight[i] = log_likelihood(dist, ray_sigma_);
        //     } else {
        //         weight[i] = log_likelihood(2.0*ray_sigma_, ray_sigma_);
        //     }
        // } else {
        //     if(dist >= 0.0 && dist <= 0.5*ray_sigma_) {
        //         weight[i] = log_likelihood(dist, ray_sigma_);
        //     } else {
        //         weight[i] = log_likelihood(0.5*ray_sigma_, ray_sigma_);
        //     }
        // }
        if(dist >= 0.0 && dist <= 0.5*ray_sigma_) {
            weight[i] = log_likelihood(dist, ray_sigma_);
        } else {
            weight[i] = log_likelihood(0.5*ray_sigma_, ray_sigma_);
        }
    }

    for(size_t i=0; i<cloud.size(); i++) {
        p.weight += weight[i];
    }
}

void Particles::ComputePosterior() {

    d_mean_posterior_.setZero();
    d_cov_posterior_.setZero();

    for(int i=0; i<set_size_; i++) {
        d_mean_posterior_.block<3,1>(0,0) += exp(d_pset_[i].weight) * d_pset_[i].translation;
        d_mean_posterior_.block<3,1>(3,0) += exp(d_pset_[i].weight) * d_pset_[i].angle_axis;

    }
    for(int i=0; i<set_size_; i++) {
        Eigen::Matrix<double, 6, 1> twist;
        twist.block<3,1>(0,0) = d_pset_[i].translation - d_mean_posterior_.block<3,1>(0,0);
        twist.block<3,1>(3,0) = d_pset_[i].angle_axis - d_mean_posterior_.block<3,1>(3,0);
        d_cov_posterior_ += exp(d_pset_[i].weight) * (twist) * (twist).transpose();
    }
}

void Particles::Propagate(Eigen::Matrix<double, 6, 1> &mean_prior,
                          Eigen::Matrix<double, 6, 6> &cov_prior,
                          Eigen::Matrix<double, 6, 1> &mean_posterior,
                          Eigen::Matrix<double, 6, 6> &cov_posterior) {

    // generate particles
    DrawParticles();
    
    // weight each particles
    AssignWeightsToParticles();

    // compute weighted mean and cov
    ComputePosterior();

    mean_prior = d_mean_sample_;
    cov_prior = d_cov_sample_;
    mean_posterior = d_mean_posterior_;
    cov_posterior  = d_cov_posterior_;
}

// std::vector<Particle> Particles::get_pset() {
//     return pset_;
// }

// std::vector<Particle> Particles::get_d_pset() {
//     return d_pset_;
// }


}  // namespace mapping
}  // namespace cartographer