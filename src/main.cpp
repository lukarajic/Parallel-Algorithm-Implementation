#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <chrono>
#include <algorithm>
#include "nbody.h"
#include "init.h"
#include "utils.h"
#include "timer.h"
#include "logger.h"
#include "constants.h"

void compute_forces_direct(System& system, const Config& config) {
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
}

void compute_forces_barnes_hut(System& system, const Config& config) {
    static OctreePool pool(config.num_bodies * 20);
    pool.reset();

    BoundingBox boundary = calculate_bounding_box(system);
    boundary.min = boundary.min - Vector3(0.1f, 0.1f, 0.1f);
    boundary.max = boundary.max + Vector3(0.1f, 0.1f, 0.1f);

    OctreeNode* root = pool.allocate(boundary);
    const size_t n = system.size();
    for (int i = 0; i < (int)n; ++i) {
        root->insert(i, system, pool);
    }
    
    #pragma omp parallel
    {
        #pragma omp single
        root->update_properties(system);
    }

    #pragma omp parallel for
    for (int i = 0; i < (int)n; ++i) {
        Vector3 force = root->compute_force(i, system, config.G, config.softening, config.theta);
        const float mi = system.mass[i];
        Vector3 vel = system.get_vel(i);
        vel.x += force.x * (config.dt / mi);
        vel.y += force.y * (config.dt / mi);
        vel.z += force.z * (config.dt / mi);
        system.set_vel(i, vel);
    }
}

void compute_forces(System& system, const Config& config) {
    if (config.use_barnes_hut) {
        compute_forces_barnes_hut(system, config);
    } else {
        compute_forces_direct(system, config);
    }

    const size_t n = system.size();
    #pragma omp parallel for simd
    for (int i = 0; i < (int)n; ++i) {
        system.x[i] += system.vx[i] * config.dt;
        system.y[i] += system.vy[i] * config.dt;
        system.z[i] += system.vz[i] * config.dt;
    }
}

void report_min_pe(const System& system) {
    auto it = std::min_element(system.potential_energy.begin(), system.potential_energy.end());
    size_t idx = std::distance(system.potential_energy.begin(), it);
    Logger::info("Min Potential Energy: " + std::to_string(*it) + " at index " + std::to_string(idx));
}

int main(int argc, char* argv[]) {
    Config config = {
        10000,                      // num_bodies
        100,                        // num_steps
        Constants::DEFAULT_DT,      // dt
        Constants::G,               // G
        Constants::DEFAULT_SOFTENING, // softening
        0.5f,                       // theta
        false                       // use_barnes_hut
    };

    bool verbose = false;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--verbose" || arg == "-v") {
            verbose = true;
            Logger::set_level(LogLevel::DEBUG);
        } else if (arg == "--barnes-hut" || arg == "-bh") {
            config.use_barnes_hut = true;
        } else if (arg == "--theta" && i + 1 < argc) {
            config.theta = std::stof(argv[++i]);
        }
    }

    try {
        if (argc > 1 && argv[1][0] != '-') {
            config.num_bodies = std::stoi(argv[1]);
        }
        if (argc > 2 && argv[2][0] != '-') {
            config.num_steps = std::stoi(argv[2]);
        }
    } catch (const std::exception& e) {
        Logger::error("Parsing command-line arguments: " + std::string(e.what()));
        return 1;
    }

    if (argc > 1 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
        std::cout << "Usage: " << argv[0] << " [num_bodies] [num_steps] [--barnes-hut] [--theta T] [--verbose]" << std::endl;
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
    
    Vector3 initial_cm = calculate_center_of_mass(system);
    Logger::info("Initial CM: (" + std::to_string(initial_cm.x) + ", " + std::to_string(initial_cm.y) + ", " + std::to_string(initial_cm.z) + ")");

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

    config.summary(elapsed_seconds);

    return 0;
}
