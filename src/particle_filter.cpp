/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 *      Author: Pravin
 *
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
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    num_particles = 50;  // TODO: Set the number of particles

    //Add random Gaussian noise to each particle. Create normal_distribution for x, y, theta.
    default_random_engine rnd;
    normal_distribution<double> norm_dist_x(x, std[0]); //x=mean , std[0]=stddev x
    normal_distribution<double> norm_dist_y(y, std[1]);//y=mean  , std[01]=stddev y
    normal_distribution<double> norm_dist_theta(theta, std[2]);//theta= mean 0, std[2]=stddev theta

    //create particles
    vector<Particle> ps;
    for(int i =0; i < num_particles; i++){
        //create Particle (id=i, weight=1)
        Particle p ;
        p.id=i; //i = particle ID
        p.weight=0; //default weight
        // sample with random Gaussian noise,
        // every particle will get different x,y,theta,
        p.x         =  norm_dist_x(rnd);
        p.y         =  norm_dist_y(rnd);
        p.theta     =  norm_dist_theta(rnd);
        //collect this particle
        ps.push_back(p);
        //weights.push_back(0);
    }

    particles=ps;
    is_initialized=true;

    cout << "INIT DONE" << endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/

    //Add random Gaussian noise to each particle. Create normal_distribution for x, y, theta.
    default_random_engine rnd;

    for( Particle& p : particles){
        //based on given arguments => estimate x,y, theta
        double new_x=0;
        double new_y=0;
        double new_theta = 0;
        //x, y, theta will change velocity and delta_t
        if (fabs(yaw_rate) < 1e-5) {
            //moving straight
            double v_t = velocity * delta_t  ;
            new_theta = p.theta;
            new_x = p.x +  v_t  * cos(p.theta);
            new_y = p.y +  v_t  * sin(p.theta);
        } else {
            //turning : x, y, theta will change based on yaw_rate
            double v_y = velocity / yaw_rate ;
            new_theta = p.theta + yaw_rate * delta_t;
            new_x = p.x +  v_y * (sin(new_theta) - sin(p.theta));
            new_y = p.y +  v_y * (cos(p.theta) - cos(new_theta));
        }

        //add noise
        normal_distribution<double> norm_dist_x(new_x, std_pos[0]); //x=mean  , std_pos[0]=stddev  x
        normal_distribution<double> norm_dist_y(new_y, std_pos[1]);//y=mean, std_pos[1]=stddev  y
        normal_distribution<double> norm_dist_theta(new_theta, std_pos[2]);//theta=mean, std_pos[2]=stddev theta

        p.x         =     norm_dist_x(rnd);
        p.y         =     norm_dist_y(rnd);
        p.theta     =     norm_dist_theta(rnd);
    }

    cout << "prediction -- DONE --" << endl;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.

    // associate each observation to closest predicated landmark
    for(LandmarkObs& o : observations){
        double min_dist = numeric_limits<double>::max();
        for(LandmarkObs& l : predicted) {
            double ob_ecl_dist = dist(o.x, o.y, l.x, l.y);
            if(ob_ecl_dist < min_dist  ){
                //choose the closest (nearest neighbor), with minimum ecludian distance between the observation and landmark  
                min_dist = ob_ecl_dist;
                o.id = l.id;
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

    //total sum of weights -- normalize probability to sum 1
//    double total_weights = 0.0;
    weights.clear();

    for( Particle& p : particles){ //for(unsigned int i=0; i< num_particles; i++) {
//        Particle &p = particles[i];
        double w;
        w = 0;
        //create list of particles close to landmarks
        vector<LandmarkObs> predicted;
        for (Map::single_landmark_s l : map_landmarks.landmark_list) {
            double l_p_d = dist(l.x_f, l.y_f, p.x, p.y);
            //find landmarks if they are close to particle
            if (l_p_d <= sensor_range) {
                LandmarkObs near_landmark;
                near_landmark.id = l.id_i;
                near_landmark.x = l.x_f;
                near_landmark.y = l.y_f;
                predicted.push_back(near_landmark);
            }
        }

        //if no landmarks nearby .... in sensor_range
        if (predicted.size() > 0) {

            // change coordinates: observation coordinate into particle map coordinate
            vector<LandmarkObs> tr_observations;
            const double cos_theta = cos(p.theta);
            const double sin_theta = sin(p.theta);

            vector<double> transformed_x;
            vector<double> transformed_y;
            vector<int> transformed_ids;

            for (const LandmarkObs &o : observations) {
                LandmarkObs tr_observation;
                tr_observation.id = -1;
                tr_observation.x = p.x + (o.x * cos_theta) - (o.y * sin_theta);
                tr_observation.y = p.y + (o.x * sin_theta) + (o.y * cos_theta);
                tr_observations.push_back(tr_observation);
            }

            // Find the predicted measurement that is closest to each observed measurement and assign the
            //   observed measurement to this particular landmark.
            dataAssociation(predicted, tr_observations);

            for (const LandmarkObs &tr_observation : tr_observations) {
                transformed_x.push_back(tr_observation.x);
                transformed_y.push_back(tr_observation.y);
                transformed_ids.push_back(tr_observation.id);
            }
            SetAssociations( p, transformed_ids, transformed_x, transformed_y);
            // update weight
            w = 1;
            double sigma_x = std_landmark[0];
            double sigma_y = std_landmark[1];
            double normal = 1.0 / (2 * M_PI * sigma_x * sigma_y);
            double tx_denom = 2 * pow(sigma_x, 2) ;
            double ty_denom= 2 * pow(sigma_y, 2);
//            for(unsigned int i=0; i <tr_observations.size();i++){
            for (const auto& o : tr_observations) {
//              LandmarkObs& l = predicted[o.id];
//                LandmarkObs  l = predicted[i];
//                LandmarkObs  o =  tr_observations[i];
                double mu_x, mu_y;
                for(const auto& l : predicted){
                    if(l.id==o.id){
                        mu_x = l.x;
                        mu_y = l.y;
                        double expon = pow(o.x - mu_x, 2) / tx_denom + pow(o.y - mu_y, 2) / ty_denom;
                        double calc_w = normal * exp( -  expon  );
                        w *= calc_w;
                        break;
                    }
                }
            }
        }
        p.weight=w;
        weights.push_back(w);
        cout << w << endl;
//        total_weights +=w;
//        weights[i]=w;
//        weights.push_back(w);
    }
    //normalize weight and save cumulative weight
//    weights.clear();
//    double cum=0;
//    for( Particle& p : particles){
//        p.weight /= total_weights;
//        cum +=p.weight;
//        weights.push_back(cum);
//    }

}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<> d(weights.begin(), weights.end()) ;//({40, 10, 10, 40}); //https://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    vector<Particle> samples_ps;
    for (int i = 0; i < num_particles; ++i) { //for (int i = 0; i < particles.size(); ++i) {
        int index= d(rd);
        Particle p = particles[index];
        samples_ps.push_back(p);
//        cout<< " resample " << index << "=" << p.x << "," << p.y << "," << p.theta << "," << p.weight << endl;
    }

    particles.clear();
    weights.clear();

    particles = samples_ps;
}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();


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
