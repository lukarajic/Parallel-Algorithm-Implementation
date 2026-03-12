#include <random>
#include "init.h"

void init_system(System& system, const Config& config) {
    std::mt19937 gen(42);
    std::uniform_real_distribution<float> pos_dist(-100.0f, 100.0f);
    std::uniform_real_distribution<float> mass_dist(1.0f, 10.0f);

    for (int i = 0; i < config.num_bodies; ++i) {
        system.x[i] = pos_dist(gen);
        system.y[i] = pos_dist(gen);
        system.z[i] = pos_dist(gen);
        system.vx[i] = 0.0f;
        system.vy[i] = 0.0f;
        system.vz[i] = 0.0f;
        system.mass[i] = mass_dist(gen);
    }
}
