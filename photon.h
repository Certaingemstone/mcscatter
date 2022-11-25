#ifndef PHOTON_H
#define PHOTON_H

#include <cmath>
#include <Eigen/Dense>

using Eigen::Vector3d;

class Photon {
    private:
        Vector3d position;
        Vector3d direction; // is normalized by constructor
        double energy; // in keV
        std::vector<std::tuple<Vector3d, double>> events; // past positions and energies dissipated
    public:
        Photon(Vector3d position, Vector3d direction, double energy);
        void move(); // move one unit in direction
        // void set_position(Vector3d newposition);
        Vector3d* get_position();
        Vector3d* get_direction();
        void scatter(double theta, double phi); // change direction, and change energy according to Compton
        void absorb(); // dump all energy at current position
};

inline Photon::Photon(Vector3d position, Vector3d direction, double energy) {
    position = position;
    direction = direction;
    energy = energy;
    events.push_back(std::tuple<Vector3d, double> {position, 0.0});
}

inline Vector3d* Photon::get_position() {
    return &position;
}
inline Vector3d* Photon::get_direction() {
    return &direction;
}
inline void Photon::move() {
    position = position + direction;
}

inline void Photon::absorb() {
    events.push_back(std::tuple<Vector3d, double> {position, energy});
    energy = 0;
}

inline void Photon::scatter(double theta, double phi) {
    // given an angle (radians), performs Compton scattering by angle theta, in a plane of angle phi
    double Efrac = 1 - (energy * (1-cos(theta) / 511)); // assume electron, 511keV
    events.push_back(std::tuple<Vector3d, double> {position, energy - Efrac * energy});
    energy = Efrac * energy;
    // get orthonormal basis vector including direciton
    Vector3d v1 = Vector3d {{0,0,1}};
    v1 = (v1 - (v1.dot(direction) * direction)).normalized(); // Gram Schmidt
    Eigen::AngleAxisd m1 = Eigen::AngleAxisd(theta, v1); // rotate by theta off current direction
    Eigen::AngleAxisd m2 = Eigen::AngleAxisd(phi, direction); // rotate around old direction by phi
    direction = m2 * m1 * direction;
}

#endif
