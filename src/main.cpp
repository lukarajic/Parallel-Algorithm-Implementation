#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "nbody.h"
#include "init.h"
#include "utils.h"
#include "timer.h"

void compute_forces(System& system, const Config& config) {
    #pragma omp parallel for
    for (int i = 0; i < config.num_bodies; ++i) {
        float fx = 0.0f, fy = 0.0f, fz = 0.0f;
        float xi = system.x[i], yi = system.y[i], zi = system.z[i], mi = system.mass[i];

        for (int j = 0; j < config.num_bodies; ++j) {
            float dx = system.x[j] - xi;
            float dy = system.y[j] - yi;
            float dz = system.z[j] - zi;
            float dist_sq = dx * dx + dy * dy + dz * dz + config.softening;
            float inv_dist = 1.0f / std::sqrt(dist_sq);
            float inv_dist3 = inv_dist * inv_dist * inv_dist;
            float f = config.G * mi * system.mass[j] * inv_dist3;
            fx += dx * f;
            fy += dy * f;
            fz += dz * f;
        }
        system.vx[i] += fx * (config.dt / mi);
        system.vy[i] += fy * (config.dt / mi);
        system.vz[i] += fz * (config.dt / mi);
    }

    #pragma omp parallel for
    for (int i = 0; i < config.num_bodies; ++i) {
        system.x[i] += system.vx[i] * config.dt;
        system.y[i] += system.vy[i] * config.dt;
        system.z[i] += system.vz[i] * config.dt;
    }
}

int main() {
    Config config = {
        5000,    // num_bodies
        10,      // num_steps
        0.01f,   // dt
        6.674e-11f, // G
        1e-9f    // softening
    };

    System system(config.num_bodies);
    init_system(system, config);

    std::cout << "Starting simulation with " << config.num_bodies << " bodies..." << std::endl;
    std::cout << std::fixed << std::setprecision(10);

    double initial_ke = calculate_kinetic_energy(system, config);
    double initial_pe = calculate_potential_energy(system, config);
    std::cout << "Initial Energy: Total=" << initial_ke + initial_pe << " KE=" << initial_ke << " PE=" << initial_pe << std::endl;

    {
        Timer timer("Total simulation");
        for (int step = 0; step < config.num_steps; ++step) {
            compute_forces(system, config);
        }
    }

    double final_ke = calculate_kinetic_energy(system, config);
    double final_pe = calculate_potential_energy(system, config);
    std::cout << "Final Energy:   Total=" << final_ke + final_pe << " KE=" << final_ke << " PE=" << final_pe << std::endl;
    std::cout << "Energy Drift: " << (final_ke + final_pe) - (initial_ke + initial_pe) << std::endl;

    return 0;
}
