#ifndef XSEC_H
#define XSEC_H

#include <cmath>

const double re = 2.81794092e-12; // classical electron radius, in mm

// Implements probabilities based on en.wikipedia.org/wiki/Gamma_ray_cross_section

double KN(double angle, double energy); 
    // Returns Klein-Nishina dsigma/dOmega for a given angle (degrees) and energy in keV, for
    // scattering off of electrons.

double KN_probability(double electron_density, double distance, double energy);
    // Given a distance in mm, and an electron density in number/mm3, and a photon energy in keV
    // Returns the total probability of a scattering event occurring in the interval
    // According to Beer-Lambert attenuation and Klein-Nishina cross section.

double absorb_probability(double electron_density, double distance, double energy, int atomic_number);
    // Same as KN_probability, but uses approximate absorption cross section for photoelectric absorption
    
inline double KN(double angle, double energy) {
    // Klein Nishina differential cross section, angle in degrees, energy in keV
    double theta = angle * 3.14159265 / 180;
    double epsilon = energy / 511;
    double lamb_ratio = 1 / (1 + epsilon * (1 - cos(theta)));
    return 0.5 * pow(re, 2) * pow(lamb_ratio, 2) * (lamb_ratio + (1/lamb_ratio) - pow(sin(theta), 2));
}

inline double KN_probability(double electron_density, double distance, double energy) {
    double total_cross_section = 0;
    // integrate over KN for the given energy
    for (int angle = 0; angle < 180; angle++) {
        total_cross_section += (KN(angle, energy) * 3.14159265 / 180);
    }
    // convert to solid angle integral
    total_cross_section *= (2 * 3.14159265);
    // use cross section and distance to get attenuation in 1mm
    double fraction_scattered = 1 - exp(-electron_density * total_cross_section);
    return fraction_scattered;
}

inline double absorb_probability(double electron_density, double distance, double energy, int atomic_number) {
    double epsilon = energy / 511;
    double total_cross_section = 35.543 * pow(re, 2) * 2.838688e-9 * pow(atomic_number, 5) / pow(epsilon, 3.5);
    double fraction_absorbed = 1 - exp(-electron_density * total_cross_section / atomic_number);
    return fraction_absorbed;
}

#endif
