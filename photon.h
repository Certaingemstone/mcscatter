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
        void scatter_KN(); // change direction, sampling from appropriate Klein-Nishina prediction
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

inline void Photon::scatter_KN() {

}

#endif
