#include <cmath>
#include <omp.h>
#include "utils.h"

double calculate_kinetic_energy(const System& system, const Config& config) {
    double total_ke = 0.0;
    const size_t n = system.size();

    #pragma omp parallel for simd reduction(+:total_ke)
    for (int i = 0; i < (int)n; ++i) {
        const double vx = system.vx[i];
        const double vy = system.vy[i];
        const double vz = system.vz[i];
        const double v_sq = vx * vx + vy * vy + vz * vz;
        total_ke += 0.5 * system.mass[i] * v_sq;
    }

    return total_ke;
}

double calculate_potential_energy(const System& system, const Config& config) {
    double total_pe = 0.0;
    const size_t n = system.size();

    #pragma omp parallel for reduction(+:total_pe)
    for (int i = 0; i < (int)n; ++i) {
        const double xi = system.x[i], yi = system.y[i], zi = system.z[i];
        const double mi = system.mass[i];
        double sub_pe = 0.0;

        #pragma omp simd reduction(+:sub_pe)
        for (int j = i + 1; j < (int)n; ++j) {
            const float r2 = Vector3::dist_sq(xi, yi, zi, system.x[j], system.y[j], system.z[j]) + config.softening;
            const float dist = std::sqrt(r2);
            sub_pe -= (config.G * mi * system.mass[j]) / dist;
        }
        total_pe += sub_pe;
    }

    return total_pe;
}

double calculate_total_mass(const System& system) {
    double total_mass = 0.0;
    const size_t n = system.size();

    #pragma omp parallel for reduction(+:total_mass)
    for (int i = 0; i < (int)n; ++i) {
        total_mass += system.mass[i];
    }

    return total_mass;
}
