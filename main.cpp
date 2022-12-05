#include <omp.h>
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <memory>
#include <Eigen/Dense>
 
using Eigen::Vector3d;

#include "geometry.h"
#include "photon.h"
#include "xsec.h"
#include "material.h"

#define PI 3.14159265

const int N_BATCH = 4;
const int N_PHOTONS_BATCH = 32768;
const int PHOTON_LIFETIME = 1000;
//const double P_SCATTER = 0.01;
//const double P_ABSORB = 0.005;
const double CONE_ANGLE = 0.08;
const std::string OUT_FILE = "result.csv";

Photon spawn_photon_cone(double spread, double energy, double rv_theta, double rv_phi) {
    // return a Photon with a direction in a cone centered on z axis,
    // starting at the origin, spread by angle (radian) spread
    // requires a random variable with some distribution on range (0,1)
    double zmin = cos(spread);
    double phi = 2*PI*rv_phi;
    double z = rv_theta * (1 - zmin) + zmin;
    double param = sqrt(1-pow(z, 2));
    Vector3d direction {{param*cos(phi), param*sin(phi), z}};
    Vector3d position{{0, 0, 0}};
    return Photon(position, direction.normalized(), energy);
}

void gen_klein_nishina_table(double (*array)[70][181]) {
    // generates 10 keV increments of KN cross section distribution from 10 to 710 keV. 
    // with 1 degree angular resolution from 0 to 180 degrees.
    // i.e. a table 70 by 181 of doubles, representing cumulative distribution of Klein-Nishina
    // probability of scattering, normalized so that at 180 degrees, the value is 1.
    double energy;
    for (int i = 0; i < 70; i++) {
        energy = i * 10 + 10;
        // KN cross section CDF for this energy at each angle
        double cumulate = 0.0;
        for (int angle = 0; angle < 181; angle++) {
            double theta = angle * PI / 180;
            double epsilon = energy / 511;
            double lamb_ratio = 1 / (1 + epsilon * (1 - cos(theta)));
            double dsigdomega = pow(lamb_ratio, 2) * (lamb_ratio + (1/lamb_ratio) - pow(sin(theta), 2));
            cumulate += dsigdomega;
            (*array)[i][angle] = cumulate;
        }
        // normalize the CDF
        double factor = 1.0 / (*array)[i][180];
        for (int angle = 0; angle < 181; angle ++) {
           (*array)[i][angle] = (*array)[i][angle] * factor;
        }
    }
}

int theta_klein_nishina(double (*LUT)[70][181], double energy, double rv) {
    // read the lookup table, given the energy in keV and a random (0,1), get angle, an int (0, 180)
    int idx = std::round(energy / 10);
    if (idx > 0) {
        idx--;
    }
    double (*cdf) = (*LUT)[idx];
    // binary search the cdf for the first value lower than rv, return the corresponding index
    double* p1 = &cdf[0];
    double* p2 = &cdf[180];
    while (p1 < p2) {
        double* mid = p1 + ((p2 - p1) / 2);
        if (rv < *mid) {
            p2 = mid - 1;
        }
        else {
            p1 = mid + 1;
        }
    }
    return p1 - &cdf[0];
}

int main(void) {
    double scatter_angle;
    std::cin >> scatter_angle;

    // prepare file for writing
    // (overwrites existing)
    std::ofstream outfile;
    outfile.open(OUT_FILE);
    outfile << "id,x,y,z,E,obj_id,\n";
    outfile.close();

    // generate lookup table for angle probabilities by energy
    double knLUT[70][181];
    gen_klein_nishina_table(&knLUT);

    // material
    double NaI_e_dens = 9.3883e20; // mm^-3
    int NaI_Z = 38;
    Material NaI = Material(NaI_e_dens, NaI_Z);

    // make our cylinders
    std::vector<Cylinder> objects;
    double sintheta = sin(scatter_angle);
    double costheta = cos(scatter_angle);

    Eigen::Matrix3d rotX90; // recoil detector
    Vector3d recoil_center = Vector3d{{0, 0, 350}};
    rotX90 << 1, 0, 0,
              0, 0, -1,
              0, 1, 0;
    Cylinder recoil(recoil_center, rotX90, 50.8, 25.4, &NaI);
    objects.push_back(recoil);

    Eigen::Matrix3d rotY; // scatter detector
    rotY << costheta, 0, -sintheta,
            0,        1,         0,
            sintheta, 0, costheta;
    Vector3d scatter_center = recoil_center + rotY * Vector3d{{0, 0, 250}};
    Cylinder scatter(scatter_center, rotY, 50.8, 25.4, &NaI);
    objects.push_back(scatter);

    // make our beam blocks...
    //

    const int num_objects = objects.size();

    // have a random number generator 
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0, 1);
    auto rv = std::bind(distribution, generator);

    std::cout << "Setup OK" << std::endl;

    // do a bunch of photon batches, and write the results of each to file
    for (int batch = 0; batch < N_BATCH; batch++) {
        std::cout << "Batch " << batch << " of " << N_BATCH;
        std::cout.flush();
        // spawn photons
        std::vector<Photon> photon_batch;
        for (int i = 0; i < N_PHOTONS_BATCH; i++) {
            photon_batch.push_back(spawn_photon_cone(CONE_ANGLE, 661.6, rv(), rv()));
        }
        photon_batch.shrink_to_fit();

        // for each photon in the batch
#       pragma omp parallel for default(shared)
        for (int i = 0; i < N_PHOTONS_BATCH; i++) {
            // evolve it
            Photon* photon_ptr = &photon_batch[i];
            bool alive = true;
            for (int t = 0; t < PHOTON_LIFETIME; t++) {
                // check if inside
                bool inside = false;
                int obj_id = -1;
                Vector3d * pos_ptr = (*photon_ptr).get_position();
                for (int obj = 0; obj < num_objects; obj++) {
                    if (objects[obj].is_inside(pos_ptr)) {
                        inside = true;
                        obj_id = obj;
                        break;
                    }
                }
                // scatter and absorb if inside something
                if (inside) {
                    Material * mat_ptr = objects[obj_id].material;
                    double P_absorb = (*mat_ptr).absorption_probability((*photon_ptr).get_energy());
                    double P_scatter = (*mat_ptr).scattering_probability((*photon_ptr).get_energy());
                    if (rv() < P_absorb) {
                        (*photon_ptr).absorb(obj_id);
                        alive = false;
                        break;
                    }
                    if (rv() < P_scatter) {
                        double theta = theta_klein_nishina(&knLUT, (*photon_ptr).get_energy(), rv());
                        double phi = 2 * PI * rv();
                        (*photon_ptr).scatter(theta, phi, obj_id);
                    }
                }
                // move always
                (*photon_ptr).move();
            }
            // at end of life, absorb wherever it is
            if (alive) {
                (*photon_ptr).absorb(-1);
            }
        }   

        // Look at photons, append their interactions to file
        outfile.open(OUT_FILE, std::ios::out | std::ios::app);
        
        for (int k = 0; k < N_PHOTONS_BATCH; k++) {
            int id = k + (N_PHOTONS_BATCH * batch);
            Photon * photon_ptr = &photon_batch[k];
            for (std::tuple<Vector3d, double, int> event : *(*photon_ptr).get_events()) {
                Vector3d pos = std::get<0>(event);
                double energy = std::get<1>(event);
                int obj_id = std::get<2>(event);
                outfile << id << "," << pos[0] << "," << pos[1] << "," << pos[2] << "," << energy << "," << obj_id << std::endl; 
            }
        }
        outfile.close();
        std::cout << "\r";
    }
    return 0;
}
