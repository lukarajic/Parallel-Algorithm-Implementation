#include <cmath>
#include <omp.h>
#include "utils.h"

double calculate_kinetic_energy(const System& system, const Config& config) {
    double total_ke = 0.0;

    #pragma omp parallel for reduction(+:total_ke)
    for (int i = 0; i < config.num_bodies; ++i) {
        double vx = system.vx[i];
        double vy = system.vy[i];
        double vz = system.vz[i];
        double v_sq = vx * vx + vy * vy + vz * vz;
        total_ke += 0.5 * system.mass[i] * v_sq;
    }

    return total_ke;
}

double calculate_potential_energy(const System& system, const Config& config) {
    double total_pe = 0.0;

    #pragma omp parallel for reduction(+:total_pe)
    for (int i = 0; i < config.num_bodies; ++i) {
        double xi = system.x[i], yi = system.y[i], zi = system.z[i];
        double mi = system.mass[i];

        for (int j = i + 1; j < config.num_bodies; ++j) {
            double dx = system.x[j] - xi;
            double dy = system.y[j] - yi;
            double dz = system.z[j] - zi;
            double dist_sq = dx * dx + dy * dy + dz * dz + config.softening;
            double dist = std::sqrt(dist_sq);
            total_pe -= (config.G * mi * system.mass[j]) / dist;
        }
    }

    return total_pe;
}
