#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <chrono>
#include "nbody.h"
#include "init.h"
#include "utils.h"
#include "timer.h"
#include "logger.h"
#include "constants.h"

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
        10000,                      // num_bodies
        100,                        // num_steps
        Constants::DEFAULT_DT,      // dt
        Constants::G,               // G
        Constants::DEFAULT_SOFTENING // softening
    };

    bool verbose = false;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--verbose" || arg == "-v") {
            verbose = true;
            Logger::set_level(LogLevel::DEBUG);
        }
    }

    try {
        if (argc > 1 && argv[1][0] != '-') {
            config.num_bodies = std::stoi(argv[1]);
        }
        if (argc > 2 && argv[2][0] != '-') {
            config.num_steps = std::stoi(argv[2]);
        }
        if (argc > 3 && argv[3][0] != '-') {
            config.dt = std::stof(argv[3]);
        }
        if (argc > 4 && argv[4][0] != '-') {
            config.softening = std::stof(argv[4]);
        }
    } catch (const std::exception& e) {
        Logger::error("Parsing command-line arguments: " + std::string(e.what()));
        std::cerr << "Usage: " << argv[0] << " [num_bodies] [num_steps] [dt] [softening] [--verbose]" << std::endl;
        return 1;
    }

    if (argc > 1 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
        std::cout << "Usage: " << argv[0] << " [num_bodies] [num_steps] [dt] [softening] [--verbose]" << std::endl;
        return 0;
    }

    if (!config.validate()) {
        return 1;
    }

    System system(config.num_bodies);
    init_system(system, config);

    config.print();
    std::cout << std::fixed << std::setprecision(10);

    double initial_ke = calculate_kinetic_energy(system, config);
    double initial_pe = calculate_potential_energy(system, config);
    Logger::info("Initial Energy: Total=" + std::to_string(initial_ke + initial_pe));
    Logger::info("Total System Mass: " + std::to_string(calculate_total_mass(system)));
    
    Vector3 initial_cm = calculate_center_of_mass(system);
    Logger::info("Initial CM: (" + std::to_string(initial_cm.x) + ", " + std::to_string(initial_cm.y) + ", " + std::to_string(initial_cm.z) + ")");
    
    Vector3 initial_p = calculate_total_momentum(system);
    Logger::info("Initial Momentum: (" + std::to_string(initial_p.x) + ", " + std::to_string(initial_p.y) + ", " + std::to_string(initial_p.z) + ")");

    Vector3 initial_l = calculate_total_angular_momentum(system);
    Logger::info("Initial Angular Momentum: (" + std::to_string(initial_l.x) + ", " + std::to_string(initial_l.y) + ", " + std::to_string(initial_l.z) + ")");

    BoundingBox initial_bb = calculate_bounding_box(system);
    Logger::info("Initial Extent: Min(" + std::to_string(initial_bb.min.x) + ", " + std::to_string(initial_bb.min.y) + ", " + std::to_string(initial_bb.min.z) + ") Max(" + std::to_string(initial_bb.max.x) + ", " + std::to_string(initial_bb.max.y) + ", " + std::to_string(initial_bb.max.z) + ")");

    double elapsed_seconds = 0.0;
    {
        auto start = std::chrono::high_resolution_clock::now();
        for (int step = 0; step < config.num_steps; ++step) {
            compute_forces(system, config);
            if ((step + 1) % 10 == 0 || step == 0 || step == config.num_steps - 1) {
                Logger::debug("Step " + std::to_string(step + 1) + " / " + std::to_string(config.num_steps) + " completed.");
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        elapsed_seconds = elapsed.count();
    }

    double final_ke = calculate_kinetic_energy(system, config);
    double final_pe = calculate_potential_energy(system, config);
    Logger::info("Final Energy:   Total=" + std::to_string(final_ke + final_pe));
    Logger::info("Energy Drift: " + std::to_string((final_ke + final_pe) - (initial_ke + initial_pe)));

    Vector3 final_cm = calculate_center_of_mass(system);
    Logger::info("Final CM:   (" + std::to_string(final_cm.x) + ", " + std::to_string(final_cm.y) + ", " + std::to_string(final_cm.z) + ")");

    Vector3 final_p = calculate_total_momentum(system);
    Logger::info("Final Momentum:   (" + std::to_string(final_p.x) + ", " + std::to_string(final_p.y) + ", " + std::to_string(final_p.z) + ")");

    Vector3 final_l = calculate_total_angular_momentum(system);
    Logger::info("Final Angular Momentum:   (" + std::to_string(final_l.x) + ", " + std::to_string(final_l.y) + ", " + std::to_string(final_l.z) + ")");

    BoundingBox final_bb = calculate_bounding_box(system);
    Logger::info("Final Extent:   Min(" + std::to_string(final_bb.min.x) + ", " + std::to_string(final_bb.min.y) + ", " + std::to_string(final_bb.min.z) + ") Max(" + std::to_string(final_bb.max.x) + ", " + std::to_string(final_bb.max.y) + ", " + std::to_string(final_bb.max.z) + ")");

    config.summary(elapsed_seconds);

    return 0;
}
