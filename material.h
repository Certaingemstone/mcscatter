#ifndef MATERIAL_H
#define MATERIAL_H

#include "xsec.h"

class Material {
    private:
        double rho_e;
        int Z;
        double (*scatter_table)[70];
        double (*absorb_table)[70];
    public:
        Material();
        Material(double electron_density, int atomic_number);
        double scattering_probability(double energy);
        double absorption_probability(double energy);
};

inline Material::Material(double electron_density, int atomic_number) {
    rho_e = electron_density;
    Z = atomic_number;
    // Populates array with scattering probabilities per millimeter for energies in intervals of 10 from 10-710 keV.
    double s[70];
    double a[70];
    for (int idx = 0; idx < 70; idx++) {
        double energy = idx * 10 + 10;
        double probability = KN_probability(rho_e, 1, energy);
        s[idx] = probability;
        probability = absorb_probability(rho_e, 1, energy, Z);
        a[idx] = probability;
    }
    scatter_table = &s;
    absorb_table = &a;
}

inline double Material::scattering_probability(double energy) {
    int key = std::round(energy / 10);
    return *scatter_table[key];
}

inline double Material::absorption_probability(double energy) {
    int key = std::round(energy / 10);
    return *absorb_table[key];
}

#endif
