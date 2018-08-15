/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles.
	num_particles = 100;
	
	// Initialize all particles to first position
	// (based on estimates of x, y, theta and their uncertainties from GPS)
	// Add random Gaussian noise to each particle.
	// create a normal (Gaussian) distribution for x, y and yaw.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
			
	for (int i=0; i<num_particles; i++) {
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1; // and all weights to 1.
		particles.push_back(p);
	}
	
	// Initialization done
	is_initialized = True;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle
	for (int i=0; i<num_particles; i++) {
		particles[i].x = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
		particles[i].y = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
		particles[i].theta = particles[i].theta + yaw_rate*delta_t;
				
		// add random Gaussian noise.
		// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
		//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
		//  http://www.cplusplus.com/reference/random/default_random_engine/
		normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
		particles[i].x = dist_x(gen);
		
		normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
		particles[i].y = dist_y(gen);
		
		normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
		particles[i].theta = dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	for (int i=0; i<observations.size(); i++) {
		double distance = 9999;
		for (int j=0; j<predicted.size(); j++) {
			double d = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			if (d<distance) {
				distance = d;
				observations[i].id = predicted[j].id;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	for (i=0; i<num_particles; i++) {
		
		// find landmarks in range
		vector<LandmarkObs> landmarks_in_range;
		for (int j=0; j<map_landmarks.landmark_list.size(); j++) {
			double dist_landmark = dist(map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f, particles[i].x, particles[i].y);
			if ( dist_landmark <= sensor_range) {
				LandmarkObs obs;
				obs.x = map_landmarks.landmark_list[j].x_f;
				obs.y = map_landmarks.landmark_list[j].y_f;
				obs.id = map_landmarks.landmark_list[j].id_i;
				landmarks_in_range.push_back(obs);
			}
		}
		
		// transformation of observations from VEHICLE'S coordinate system into MAP'S coordinate system
		vector<LandmarkObs> observations_MCS;
		for (int j=0; j<observations.size(); j++) {
			LandmarkObs obs;
        		obs.x = particles[i].x + cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y;
			obs.y = particles[i].y + sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y;
			observations_MCS.push_back(obs);
		}
		
		// Find closest predicted measurements to observed measurements and assign the observed measurements to this particular landmarks
		dataAssociation(landmarks_in_range, observations_MCS);
		
		// calculate weights
		double gauss_norm = ( 1 / (2*M_PI*std_landmark[0]*std_landmark[1]) ); // calculate normalization term
		
		for (int j=0; j<observations_MCS.size(); j++) {
			for (int k=0; k<landmarks_in_range.size(); k++) {
				if (landmarks_in_range[k].id == observations_MCS[j].id) {
					mu_x = landmarks_in_range[k].x;
					mu_y = landmarks_in_range[k].y;
				}
				// ToDo: catch else error!
			}
			// calculate exponent
			double exponent = ((observations_MCS[j].x - mu_x)**2)/(2 * std_landmark[0]**2) + ((observations_MCS[j].y - mu_y)**2)/(2 * std_landmark[1]**2);
			
			// calculate weight using normalization terms and exponent
			particles[i].weight = gauss_norm * exp(-exponent);
		}
	}
	
	// Normalization
	double weights_sum = 0;
	for (int i=0; i<num_particles; i++) {
		weights_sum += particles[i].weight;
	}
	
	for (int i=0; i<num_particles; i++) {
		particles[i].weight /= weights_sum;
	}	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	vector<Particle> particles_resampled(num_particles);
	
	// get all weights
	vector<double> all_weights;
	for (int i=0; i<num_particles; i++) {
		all_weights.push_back(particles[i].weight);
	}
	
	discrete_distribution<> d(all_weights.begin(), all_weights.end());
	
	std::map<int, int> m;
	for (int i=0; i<num_particles; ++i) {
		int j = d(gen);
		particles_resampled[i] = particles[j];
	}
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
