#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include "nbody.h"
#include "init.h"
#include "utils.h"
#include "timer.h"

void compute_forces(System& system, const Config& config) {
    const size_t n = system.size();
    #pragma omp parallel for
    for (int i = 0; i < (int)n; ++i) {
        float fx = 0.0f, fy = 0.0f, fz = 0.0f;
        const Vector3 pos_i = system.get_pos(i);
        const float mi = system.mass[i];

        #pragma omp simd reduction(+:fx, fy, fz)
        for (int j = 0; j < (int)n; ++j) {
            const float r2 = Vector3::dist_sq(pos_i.x, pos_i.y, pos_i.z, system.x[j], system.y[j], system.z[j]) + config.softening;
            const float inv_dist = 1.0f / std::sqrt(r2);
            const float inv_dist3 = inv_dist * inv_dist * inv_dist;
            const float f = config.G * mi * system.mass[j] * inv_dist3;
            fx += (system.x[j] - pos_i.x) * f;
            fy += (system.y[j] - pos_i.y) * f;
            fz += (system.z[j] - pos_i.z) * f;
        }
        
        Vector3 vel = system.get_vel(i);
        vel.x += fx * (config.dt / mi);
        vel.y += fy * (config.dt / mi);
        vel.z += fz * (config.dt / mi);
        system.set_vel(i, vel);
    }

    #pragma omp parallel for simd
    for (int i = 0; i < (int)n; ++i) {
        system.x[i] += system.vx[i] * config.dt;
        system.y[i] += system.vy[i] * config.dt;
        system.z[i] += system.vz[i] * config.dt;
    }
}

int main(int argc, char* argv[]) {
    Config config = {
        10000,   // num_bodies
        100,     // num_steps
        0.01f,   // dt
        6.674e-11f, // G
        1e-9f    // softening
    };

    if (argc > 1) {
        config.num_bodies = std::stoi(argv[1]);
    }
    if (argc > 2) {
        config.num_steps = std::stoi(argv[2]);
    }

    if (argc > 3 || (argc > 1 && std::string(argv[1]) == "--help")) {
        std::cout << "Usage: " << argv[0] << " [num_bodies] [num_steps]" << std::endl;
        return 0;
    }

    System system(config.num_bodies);
    init_system(system, config);

    std::cout << "Starting simulation with " << config.num_bodies << " bodies for " << config.num_steps << " steps..." << std::endl;
    std::cout << std::fixed << std::setprecision(10);

    double initial_ke = calculate_kinetic_energy(system, config);
    double initial_pe = calculate_potential_energy(system, config);
    std::cout << "Initial Energy: Total=" << initial_ke + initial_pe << " KE=" << initial_ke << " PE=" << initial_pe << std::endl;

    {
        Timer timer("Total simulation");
        for (int step = 0; step < config.num_steps; ++step) {
            compute_forces(system, config);
            if ((step + 1) % 10 == 0 || step == 0 || step == config.num_steps - 1) {
                std::cout << "Step " << std::setw(4) << step + 1 << " / " << config.num_steps << " completed." << std::endl;
            }
        }
    }

    double final_ke = calculate_kinetic_energy(system, config);
    double final_pe = calculate_potential_energy(system, config);
    std::cout << "Final Energy:   Total=" << final_ke + final_pe << " KE=" << final_ke << " PE=" << final_pe << std::endl;
    std::cout << "Energy Drift: " << (final_ke + final_pe) - (initial_ke + initial_pe) << std::endl;

    return 0;
}
