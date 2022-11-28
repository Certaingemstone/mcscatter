#ifndef PHOTON_H
#define PHOTON_H

#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using Eigen::Vector3d;

class Photon {
    private:
        Vector3d position;
        Vector3d direction; // is normalized by constructor
        double energy; // in keV
        std::vector<std::tuple<Vector3d, double, int>> events; // past positions energies dissipated and object IDs
    public:
        Photon(Vector3d pos, Vector3d dir, double e);
        void move(); // move one unit in direction
        // void set_position(Vector3d newposition);
        Vector3d* get_position();
        Vector3d* get_direction();
        double get_energy();
        std::vector<std::tuple<Vector3d, double, int>>* get_events();
        void scatter(double theta, double phi, int obj_id); // change direction, and change energy according to Compton
        void absorb(int obj_id); // dump all energy at current position
};

inline Photon::Photon(Vector3d pos, Vector3d dir, double e) {
    position = pos;
    direction = dir.normalized();
    energy = e;
    events.push_back(std::tuple<Vector3d, double, int> {position, 0.0, -1});
}

inline Vector3d* Photon::get_position() {
    return &position;
}
inline Vector3d* Photon::get_direction() {
    return &direction;
}
inline double Photon::get_energy() {
    return energy;
}
inline std::vector<std::tuple<Vector3d, double, int>>* Photon::get_events() {
    return &events;
}
inline void Photon::move() {
    position = position + direction;
}

inline void Photon::absorb(int obj_id) {
    events.push_back(std::tuple<Vector3d, double, int> {position, energy, obj_id});
    energy = 0;
}

inline void Photon::scatter(double theta, double phi, int obj_id) {
    // given an angle (radians), performs Compton scattering by angle theta, in a plane of angle phi
    double Efrac = 1 / (1 + (energy / 511) * (1 - cos(theta))); // assume electron, 511keV
    events.push_back(std::tuple<Vector3d, double, int> {position, energy - Efrac * energy, obj_id});
    energy = Efrac * energy;
    // get orthonormal basis vector including direciton
    Vector3d v1 = Vector3d {{0,0,1}};
    v1 = (v1 - (v1.dot(direction) * direction)).normalized(); // Gram Schmidt
    Eigen::AngleAxisd m1 = Eigen::AngleAxisd(theta, v1); // rotate by theta off current direction
    Eigen::AngleAxisd m2 = Eigen::AngleAxisd(phi, direction); // rotate around old direction by phi
    direction = m2 * m1 * direction;
}

#endif
