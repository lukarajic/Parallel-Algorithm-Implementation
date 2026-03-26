#ifndef UTILS_H
#define UTILS_H

#include "nbody.h"

double calculate_kinetic_energy(const System& system, const Config& config);
double calculate_potential_energy(const System& system, const Config& config);
double calculate_total_energy(const System& system, const Config& config);
void calculate_potential_energy_per_particle(System& system, const Config& config);
double calculate_total_mass(const System& system);
Vector3 calculate_center_of_mass(const System& system);
Vector3 calculate_total_momentum(const System& system);
Vector3 calculate_total_angular_momentum(const System& system);
BoundingBox calculate_bounding_box(const System& system);

#endif // UTILS_H
