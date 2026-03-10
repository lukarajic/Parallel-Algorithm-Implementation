#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include "nbody.h"

void init_system(System& system, int n) {
    std::mt19937 gen(42);
    std::uniform_real_distribution<float> pos_dist(-100.0f, 100.0f);
    std::uniform_real_distribution<float> mass_dist(1.0f, 10.0f);

    for (int i = 0; i < n; ++i) {
        system.x[i] = pos_dist(gen);
        system.y[i] = pos_dist(gen);
        system.z[i] = pos_dist(gen);
        system.vx[i] = 0.0f;
        system.vy[i] = 0.0f;
        system.vz[i] = 0.0f;
        system.mass[i] = mass_dist(gen);
    }
}

void compute_forces(System& system, int n, float dt) {
    const float G = 6.674e-11f;
    const float softening = 1e-9f;

    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        float fx = 0.0f, fy = 0.0f, fz = 0.0f;
        float xi = system.x[i], yi = system.y[i], zi = system.z[i], mi = system.mass[i];

        for (int j = 0; j < n; ++j) {
            float dx = system.x[j] - xi;
            float dy = system.y[j] - yi;
            float dz = system.z[j] - zi;
            float dist_sq = dx * dx + dy * dy + dz * dz + softening;
            float inv_dist = 1.0f / std::sqrt(dist_sq);
            float inv_dist3 = inv_dist * inv_dist * inv_dist;
            float f = G * mi * system.mass[j] * inv_dist3;
            fx += dx * f;
            fy += dy * f;
            fz += dz * f;
        }
        system.vx[i] += fx * (dt / mi);
        system.vy[i] += fy * (dt / mi);
        system.vz[i] += fz * (dt / mi);
    }

    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        system.x[i] += system.vx[i] * dt;
        system.y[i] += system.vy[i] * dt;
        system.z[i] += system.vz[i] * dt;
    }
}

int main() {
    const int num_bodies = 5000;
    const int num_steps = 10;
    const float dt = 0.01f;

    System system(num_bodies);
    init_system(system, num_bodies);

    std::cout << "Starting simulation with " << num_bodies << " bodies (SoA)..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    for (int step = 0; step < num_steps; ++step) {
        compute_forces(system, num_bodies, dt);
    }
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Simulation finished in " << elapsed.count() << " seconds." << std::endl;

    return 0;
}
