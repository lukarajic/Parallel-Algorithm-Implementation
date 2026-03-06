#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include "nbody.h"

void init_bodies(std::vector<Body>& bodies, int n) {
    std::mt19937 gen(42);
    std::uniform_real_distribution<float> pos_dist(-100.0f, 100.0f);
    std::uniform_real_distribution<float> mass_dist(1.0f, 10.0f);

    for (int i = 0; i < n; ++i) {
        bodies[i] = {
            {pos_dist(gen), pos_dist(gen), pos_dist(gen)},
            {0.0f, 0.0f, 0.0f},
            mass_dist(gen)
        };
    }
}

void compute_forces(std::vector<Body>& bodies, float dt) {
    const float G = 6.674e-11f; // Gravitational constant
    const float softening = 1e-9f; // To avoid division by zero

    for (size_t i = 0; i < bodies.size(); ++i) {
        Vector3 force = {0.0f, 0.0f, 0.0f};
        for (size_t j = 0; j < bodies.size(); ++j) {
            if (i == j) continue;

            Vector3 r = bodies[j].position - bodies[i].position;
            float dist_sq = r.x * r.x + r.y * r.y + r.z * r.z + softening;
            float inv_dist = 1.0f / std::sqrt(dist_sq);
            float inv_dist3 = inv_dist * inv_dist * inv_dist;

            float f = G * bodies[i].mass * bodies[j].mass * inv_dist3;
            force += r * f;
        }
        // Update velocity: F = ma => a = F/m
        bodies[i].velocity += force * (dt / bodies[i].mass);
    }

    // Update positions
    for (auto& b : bodies) {
        b.position += b.velocity * dt;
    }
}

int main() {
    const int num_bodies = 1000;
    const int num_steps = 10;
    const float dt = 0.01f;

    std::vector<Body> bodies(num_bodies);
    init_bodies(bodies, num_bodies);

    std::cout << "Starting simulation with " << num_bodies << " bodies..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    for (int step = 0; step < num_steps; ++step) {
        compute_forces(bodies, dt);
    }
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Simulation finished in " << elapsed.count() << " seconds." << std::endl;

    return 0;
}
