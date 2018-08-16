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
	default_random_engine gen;
    		
	// Set the number of particles.
	num_particles = 100;
	
	// Initialize all particles to first position
	// (based on estimates of x, y, theta and their uncertainties from GPS)
	// Add random Gaussian noise to each particle.
	// create a normal (Gaussian) distribution for x, y and yaw.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	weights = vector<double>(num_particles);
	
	particles = vector<Particle>(num_particles);	
	for (int i=0; i<num_particles; i++) {
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0; // and all weights to 1.
		particles.push_back(p);
		weights[i] = 1.0;
	}
	
	
	// Initialization done
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	default_random_engine gen;
	
	// Add measurements to each particle
	for (int i=0; i<num_particles; i++) {
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;
		
		if ( fabs(yaw_rate) < 0.00001 ) {
			x += velocity*delta_t*cos(theta);
			y += velocity*delta_t*sin(theta);
		}
		else {
			x += (velocity/yaw_rate) * ( sin(theta+yaw_rate*delta_t) - sin(theta));
			y += (velocity/yaw_rate) * ( cos(theta) - cos(theta+yaw_rate*delta_t));
			theta += yaw_rate*delta_t;
		}
				
		// add random Gaussian noise.
		// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
		//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
		//  http://www.cplusplus.com/reference/random/default_random_engine/
		normal_distribution<double> dist_x(x, std_pos[0]);
		particles[i].x = dist_x(gen);
		
		normal_distribution<double> dist_y(y, std_pos[1]);
		particles[i].y = dist_y(gen);
		
		normal_distribution<double> dist_theta(theta, std_pos[2]);
		particles[i].theta = dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
		
	for (int i=0; i<(int)observations.size(); i++) {
		double distance = 9999;
		for (int j=0; j<(int)predicted.size(); j++) {
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
	
	for (int i=0; i<num_particles; i++) {
		
		// find landmarks in range
		vector<LandmarkObs> landmarks_in_range;

		for (int j=0; j<(int)map_landmarks.landmark_list.size(); j++) {
			double mx = map_landmarks.landmark_list[j].x_f;
			double my = map_landmarks.landmark_list[j].y_f;
			double dist_landmark = dist(mx, my, particles[i].x, particles[i].y);
			//double dist_landmark = dist(particles[i].x, particles[i].y, mx, my);
			
			if ( dist_landmark <= sensor_range) {
				LandmarkObs obs;
				obs.id = map_landmarks.landmark_list[j].id_i;
				obs.x = map_landmarks.landmark_list[j].x_f;
				obs.y = map_landmarks.landmark_list[j].y_f;
				landmarks_in_range.push_back(obs);
			}
		}
		
		// transformation of observations from VEHICLE'S coordinate system into MAP'S coordinate system
		vector<LandmarkObs> observations_MCS;
		for (int j=0; j<(int)observations.size(); j++) {
			LandmarkObs obs;
			double p_theta = particles[i].theta;
			double x_obs = observations[j].x;
			double y_obs = observations[j].y;
			obs.id = observations[j].id;
        		obs.x = particles[i].x + cos(p_theta)*x_obs - sin(p_theta)*y_obs;
			obs.y = particles[i].y + sin(p_theta)*x_obs + cos(p_theta)*y_obs;
			observations_MCS.push_back(obs);
		}
		
		
		// Find closest predicted measurements to observed measurements and assign the observed measurements to this particular landmarks
		dataAssociation(landmarks_in_range, observations_MCS);
		
		// calculate weights
		double gauss_norm = ( 1 / (2*M_PI*std_landmark[0]*std_landmark[1]) ); // calculate normalization term
		
		for (int j=0; j<(int)observations_MCS.size(); j++) {
			double mu_x, mu_y;
			bool id_notfound = true;
			for (int k=0; k<(int)landmarks_in_range.size(); k++) {
				
				if (landmarks_in_range[k].id == observations_MCS[j].id) {
					mu_x = landmarks_in_range[k].x;
					mu_y = landmarks_in_range[k].y;
					id_notfound = false;
					break;
				}
				
			}
			if (id_notfound) {
				cout << "no ID found !!!!!!!!!!!!!!!" << endl;
			}
			// calculate exponent
			double ox = observations_MCS[j].x;
			double oy = observations_MCS[j].y;
			double exponent = ( pow((ox-mu_x),2) / (2*pow(std_landmark[0],2)) ) + ( pow((oy-mu_y),2)/ (2*pow(std_landmark[1],2)));
			
			// calculate weight using normalization terms and exponent
			double weight = gauss_norm * exp(-exponent);
			particles[i].weight = weight;
			weights[i] = weight;
		}
	}
	
	// Normalization of the weights
	double weights_sum = 0.0;
	for (int i=0; i<(int)weights.size(); i++) {
		weights_sum += weights[i];
	} 
		
	for (int i=0; i<num_particles; i++) {
		weights[i] /= weights_sum;
	}
			
}

void ParticleFilter::resample() {
	// Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	
	vector<Particle> particles_resampled(num_particles);
	
	discrete_distribution<> d(weights.begin(), weights.end());
	
	for (int i=0; i<num_particles; ++i) {
		particles_resampled[i] = particles[d(gen)];
	}
	
	particles = particles_resampled;

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

	return particle; 
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
