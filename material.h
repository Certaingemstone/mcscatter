#ifndef MATERIAL_H
#define MATERIAL_H

#include <iostream>
#include "xsec.h"

class Material {
    private:
        static const int table_size = 70;
        double rho_e;
        int Z;
        double scatter_table[table_size];
        double absorb_table[table_size];
    public:
        Material(double electron_density = 1000000., int atomic_number = 1);
        double scattering_probability(double energy);
        double absorption_probability(double energy);
        double read_sca(int idx);
        double get_electron_density();
        int get_Z();
};

inline Material::Material(double electron_density, int atomic_number) {
    rho_e = electron_density;
    Z = atomic_number;
    // Populates array with scattering probabilities per millimeter for energies in intervals of 10 from 10-710 keV.
    for (int idx = 0; idx < table_size; idx++) {
        double energy = idx * 10 + 10;
        //std::cout << energy << std::endl;
        double probability = KN_probability(rho_e, 1, energy);
        scatter_table[idx] = probability;
        //std::cout << s[idx] << std::endl;
        probability = absorb_probability(rho_e, 1, energy, Z);
        absorb_table[idx] = probability;
        //std::cout << a[idx] << std::endl;
    }
}

inline double Material::read_sca(int idx) {
    return scatter_table[idx];
}

inline double Material::scattering_probability(double energy) {
    int key = std::round(energy / 10);
    return scatter_table[key];
}

inline double Material::absorption_probability(double energy) {
    int key = std::round(energy / 10);
    return absorb_table[key];
}

inline double Material::get_electron_density() {
    return rho_e;
}

inline int Material::get_Z() {
    return Z;
}
#endif
