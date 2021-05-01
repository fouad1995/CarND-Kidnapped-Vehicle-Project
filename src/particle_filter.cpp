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
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles


  // initalize particles to first position from GPS
  std::default_random_engine gen; // random number generator

    // to make code more readable
  double std_x          = std[0];
  double std_y          = std[1];
  double std_theta      = std[2];

  // create gaussian distribution for GPS position ( ie : adding gaussian noise to GPS data )
  std::normal_distribution<double> gaussian_x(x, std_x);
  std::normal_distribution<double> gaussian_y(y, std_y);
  std::normal_distribution<double> gaussian_theta(theta, std_theta);

  // initializing particles 
  Particle dummy_particle;
  for (int i = 0; i < num_particles; i++) {

      dummy_particle.id     = i;

      // initializing particle position and adding gaussian noise to it
      dummy_particle.x      = gaussian_x(gen);
      dummy_particle.y      = gaussian_y(gen);
      dummy_particle.theta  = gaussian_theta(gen);

      // initalizing weight to 1 
      dummy_particle.weight = 1;


      // save the particle 
      particles.push_back(dummy_particle);

      weights.push_back(1); // initialize all weights to 1 
  }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
    /*
     * prediction Predicts the state for the next time step using the process model.
     * Applying these twoo functions 
     * x = x0 + (v/theta_dot)*[sin(theta+theta_dot*delta_t) - sin(theta)]
     * y = y0 + (v/theta_dot)*[cos(theta) - cos(theta+theta_dot*delta_t)]
     * theta = theta0 + theta_dot*delta_t
     */

     // to make code more readable
    double std_x = std[0];
    double std_y = std[1];
    double std_theta = std[2];
    // adding gaussian noise to velocity and yaw_rate
    std::default_random_engine gen; // random number generator
    // ask about this in code review 
    std::normal_distribution<double> gaussian_x(0, std_x);
    std::normal_distribution<double> gaussian_y(0, std_y);
    std::normal_distribution<double> gaussian_theta(0, std_theta);

    // update the position of all particles using prediction model 

    for (int i = 0; i < num_particles;i++) {
    
        // p here stands for particle 
        double x_0     = particles[i].x     ;
        double y_0     = particles[i].y     ;
        double theta_0 = particles[i].theta ;

        double A    = velocity / yaw_rate; // to prevent redundeny 
        double B    = theta_0 + (yaw_rate * delta_t); // to prevent redundeny 
        // update position 
        particles[i].x      = x_0 + A * (sin( B) - sin(theta_0)) + gaussian_x(gen);
        particles[i].y      = y_0 + A * (cos(theta_0) - cos(B))  + gaussian_y(gen);
        particles[i].theta  = B + gaussian_theta(gen);
        
    }   
        

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

    // predicted   : prediction measurment between one particular particle and all of the map landmarks (what actually the landmark lies on map)

    // observation : actual measurment gatherd from lidar,  transformed sensor measurment to the map coordinate (by particle)
    double minimum_dist = 0.0;
    double minimum_dist_index = 0;
    for (int i = 0; i < observations.size(); i++) {

        for (int j = 0; j < predicted.size(); j++) {
            double distance = dist(observations[i].x, observations[i].y,
                                    predicted[j].x, predicted[j].y);
            if (distance < minimum_dist) {
                minimum_dist = distance;
                minimum_dist_index = j;
            }
                
        }

        // associate the observation with landmark 
        // observation[i] belongs to predicted[minimum_dist_index]
        observations[i].id = predicted[minimum_dist_index].id;
        
        
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

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

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