#ifndef UTILS_H
#define UTILS_H

#include "nbody.h"

double calculate_kinetic_energy(const System& system, const Config& config);
double calculate_potential_energy(const System& system, const Config& config);

#endif // UTILS_H
