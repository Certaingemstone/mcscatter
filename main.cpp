#include <omp.h>
#include <iostream>
#include <random>
#include <vector>
#include <memory>
#include <Eigen/Dense>
 
using Eigen::Vector3d;

#include "geometry.h"
#include "photon.h"

#define PI 3.14159265

const int N_BATCH = 8;
const int N_PHOTONS_BATCH = 512;
const int PHOTON_LIFETIME = 1000;
const double MIN_ENERGY = 10.0;

Photon spawn_photon_cone(double spread, double energy, double rv) {
    // return a Photon with a direction in a cone centered on z axis,
    // starting at the origin, spread by angle (radian) spread
    // requires a random variable with some distribution on range (0,1)
    double zmin = cos(spread);
    double phi = 2*PI*rv;
    double z = rv * (1 - zmin) + zmin;
    double param = sqrt(1-pow(z, 2));
    Vector3d direction {{param*cos(phi), param*sin(phi), z}};
    Vector3d position{{0, 0, 0}};
    return Photon(position, direction, energy);
}

int main()
{
    // make our cylinders
    std::vector<Object> objects;
    const double P_absorb = 0.1;
    const double P_scatter = 0.1;
    const double theta = PI / 6;
    double sintheta = sin(theta);
    double costheta = cos(theta);
    Eigen::Matrix3d rotX90; // recoil detector
    Vector3d recoil_center = Vector3d{{0, 0, 350}};
    rotX90 << 1, 0, 0,
              0, 0, -1,
              0, 1, 0;
    Cylinder recoil(recoil_center, rotX90, 50.8, 25.4, P_scatter, P_absorb);
    objects.push_back(recoil);

    Eigen::Matrix3d rotY; // scatter detector
    rotY << costheta, 0, -sintheta,
            0,        1,         0,
            sintheta, 0, costheta;
    Vector3d scatter_center = recoil_center + rotY * Vector3d{{0, 0, 250}};
    Cylinder scatter(scatter_center, rotY, 50.8, 25.4, P_scatter, P_absorb);
    objects.push_back(scatter);

    // make our beam blocks...
    //

    const int num_objects = objects.size();

    // have a random number generator 
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0, 1);
    auto rv = std::bind(distribution, generator);

    // do a bunch of photon batches, and write the results of each to file
    for (int batch = 0; batch < N_BATCH; batch++) {

        // spawn photons
        std::vector<Photon> photon_batch;
        for (int i = 0; i < N_PHOTONS_BATCH; i++) {
            photon_batch.push_back(spawn_photon_cone(0.08, 661.6, rv()));
        }
        photon_batch.shrink_to_fit();

        // evolve the photons
#       pragma omp parallel for default(shared)
        for (int i = 0; i < N_PHOTONS_BATCH; i++) {
            // evolve it
            Photon photon = photon_batch[i];
            for (int t = 0; t < PHOTON_LIFETIME; t++) {
                // check if inside
                bool inside = false;
                int obj_id = -1;
                for (int obj = 0; obj < num_objects; obj++) {
                    if (objects[obj].is_inside(photon.get_position())) {
                        inside = true;
                        obj_id = obj;
                    }
                }
                // scatter and absorb if inside something
                if (inside) {
                    if (rv() < objects[obj_id].P_absorb) {
                        photon.absorb();
                        break;
                    }
                    if (rv() < objects[obj_id].P_scatter) {
                        photon.scatter_KN();
                    }
                }
                // move always
                photon.move();
            }
        }   

        // Look at photons, append their interactions to file
        for (Photon photon : photon_batch) {
            
        }
    }
}

