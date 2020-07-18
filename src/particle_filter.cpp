/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * Set the number of particles. Initialize all particles to 
   * first position (based on estimates of x, y, theta and their uncertainties
   * from GPS) and all weights to 1. 
   * Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 10;  // Sets the number of particles
  
  // create gaussians around the initialization parameters with given deviations
  // std::default_random_engine gen;  // now in header file
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  
  // loop over num_particles and fill each particle 
  // with values from the gaussians above
  for (int i = 0; i < num_particles; ++i) {
    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;
    particles.push_back(particle);
  }
  
  // set is_initialized to true, so this routine will not be called again
  is_initialized = true;
  std::cout << "INFO: particle filter initialized" << std::endl;
  
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  // prediction step for each particle
  // std::cout << "INFO: prediction step" << std::endl;
  // std::cout << "velocity: " << velocity << std::endl;
  // std::cout << "yaw_rate: " << yaw_rate << std::endl;
  for (unsigned int i = 0; i < particles.size(); ++i) {
    particles[i].x += (velocity / yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
    particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
    particles[i].theta += yaw_rate * delta_t;
  }
  
  // add gaussian noise to these new values
  for (unsigned int i = 0; i < particles.size(); ++i) {
    // std::default_random_engine gen;  // now in header file
    std::normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
    std::normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
    std::normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
    
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  //struct LandmarkObs {
  //int id;     // Id of matching landmark in the map.
  //double x;   // Local (vehicle coords) x position of landmark observation [m]
  //double y;   // Local (vehicle coords) y position of landmark observation [m]
  //};
  
  // loop over all observed landmarks
  for (unsigned int i = 0; i < observations.size(); ++i) {
    // initialize max_distance with a large value
    double max_distance = 10000;
    
    // now for each observed landmark check the predictions for each particle
    for (unsigned int j = 0; j < predicted.size(); ++j) {
      // and calculate the euclidean distance
      double current_distance = sqrt(pow(observations[i].x - predicted[j].x, 2.0) + pow(observations[i].y - predicted[j].y, 2.0));
      
      // if the euclidead distance is smaller than the current smallest
      // save the best particle id to the observation id and 
      // set the max_distance to this distance
      if (current_distance < max_distance) {
        observations[i].id = predicted[j].id;
        max_distance = current_distance;
      }
    }
    // std::cout << "max_distance: " << max_distance << std::endl;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  // predict measurements to all the map landmarks within sensor range for each particle
  // loop over all particles
  for (unsigned int i = 0; i < particles.size(); ++i) {
    vector<LandmarkObs> predicted;  // create "landmark" vector with all the predicted measurements for this particle
    // look for landmarks in the map within sensor range
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); ++j) {
      double current_distance = sqrt(pow(particles[i].x - map_landmarks.landmark_list[j].x_f, 2.0) + pow(particles[i].y - map_landmarks.landmark_list[j].y_f, 2.0));
      if (current_distance <= sensor_range) {
        LandmarkObs candidate;
        candidate.id = map_landmarks.landmark_list[j].id_i;
        candidate.x = map_landmarks.landmark_list[j].x_f;
        candidate.y = map_landmarks.landmark_list[j].y_f;
        predicted.push_back(candidate);
      }
    }
    // std::cout << "INFO: number of landmarks in range for this particle: " << predicted.size() << std::endl;
    
    // our observations are still in the vehicle frame
    // we need them in the map frame using the current pose of the particle
    vector<LandmarkObs> observations_map_frame;
    for (unsigned int j = 0; j < observations.size(); ++j) {
      LandmarkObs observation_map_frame;
      observation_map_frame.id = observations[j].id;
      observation_map_frame.x = particles[i].x + (cos(particles[i].theta) * observations[j].x) - (sin(particles[i].theta) * observations[j].y);
      observation_map_frame.y = particles[i].y + (sin(particles[i].theta) * observations[j].x) + (cos(particles[i].theta) * observations[j].y);
      observations_map_frame.push_back(observation_map_frame);
    }
    
    // now we have a predicted landmark vector and a observed landmark vector
    // we call dataAssociation here with the transformed observations
    dataAssociation(predicted, observations_map_frame);
    
    // each observation now has the ID of the nearest predicted landmark
    // now, we will need to calculate the multivariate gaussian probabilities
    // for each observed landmark and multiply? them all together (weight of particle)
    // weight init with 1
    particles[i].weight = 1.0;
    for (unsigned int j = 0; j < observations_map_frame.size(); ++j) {
      // we need to find the corresponding predicted values
      unsigned int k = 0;
      bool found_match = false;
      // run the loop as long as we didn't find an id match or the vector is done
      while (not found_match && k < predicted.size()) {
        if (observations_map_frame[j].id == predicted[k].id) {
          //found match
          // std::cout << "INFO: found match" << std::endl;
          found_match = true;
          double scaling_factor = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
          double a = (pow(observations_map_frame[j].x - predicted[k].x, 2.0)) / (2 * (std_landmark[0] * std_landmark[0]));
          double b = (pow(observations_map_frame[j].y - predicted[k].y, 2.0)) / (2 * (std_landmark[1] * std_landmark[1]));
          double proba = scaling_factor * exp(-(a + b));
          // std::cout << "distance x: " << observations_map_frame[j].x - predicted[k].x << std::endl;
          // std::cout << "distance y: " << observations_map_frame[j].y - predicted[k].y << std::endl;
          // std::cout << "proba: " << proba << std::endl;
          particles[i].weight *= proba;
        }
        k++;
      }

    }
    
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  // particles have unnormalized weights
  // we will normalize them
  double max_weight = 0;
  for (unsigned int j = 0; j < particles.size(); ++j) {
    if (particles[j].weight > max_weight) {
      max_weight = particles[j].weight;
    }
  }
  for (unsigned int j = 0; j < particles.size(); ++j) {
    particles[j].weight /= max_weight;
  }
  
  // now we can use the resampling wheel from the lectures
  //p3 = []
  //index = int(random.random() * N)
  // std::default_random_engine gen;  // now in header file
  std::uniform_int_distribution<int> random_distribution(2, particles.size());
  std::cout << "INFO: random index: " << random_distribution(gen) << ", number of particles: " << particles.size() << std::endl;
  //beta = 0.0
  //mw = max(w)

  //for i in range(N):
  //    beta += random.random() * 2.0 * mw
  //    while beta > w[index]:
  //        beta -= w[index]
  //        index = (index + 1) % N
  //    p3.append(p[index])

  //p = p3
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}