#include <cmath>
#include <omp.h>
#include <limits>
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

void calculate_potential_energy_per_particle(System& system, const Config& config) {
    const size_t n = system.size();

    #pragma omp parallel for
    for (int i = 0; i < (int)n; ++i) {
        const double xi = system.x[i], yi = system.y[i], zi = system.z[i];
        const double mi = system.mass[i];
        double pe = 0.0;

        for (int j = 0; j < (int)n; ++j) {
            if (i == j) continue;
            const float r2 = Vector3::dist_sq(xi, yi, zi, system.x[j], system.y[j], system.z[j]) + config.softening;
            const float dist = std::sqrt(r2);
            pe -= (config.G * mi * system.mass[j]) / dist;
        }
        system.potential_energy[i] = (float)pe;
    }
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

Vector3 calculate_center_of_mass(const System& system) {
    double cm_x = 0.0, cm_y = 0.0, cm_z = 0.0;
    double total_mass = 0.0;
    const size_t n = system.size();

    #pragma omp parallel for reduction(+:cm_x, cm_y, cm_z, total_mass)
    for (int i = 0; i < (int)n; ++i) {
        const double m = system.mass[i];
        cm_x += system.x[i] * m;
        cm_y += system.y[i] * m;
        cm_z += system.z[i] * m;
        total_mass += m;
    }

    if (total_mass > 0) {
        return Vector3((float)(cm_x / total_mass), (float)(cm_y / total_mass), (float)(cm_z / total_mass));
    }
    return Vector3(0.0f, 0.0f, 0.0f);
}

Vector3 calculate_total_momentum(const System& system) {
    double px = 0.0, py = 0.0, pz = 0.0;
    const size_t n = system.size();

    #pragma omp parallel for reduction(+:px, py, pz)
    for (int i = 0; i < (int)n; ++i) {
        const double m = system.mass[i];
        px += system.vx[i] * m;
        py += system.vy[i] * m;
        pz += system.vz[i] * m;
    }

    return Vector3((float)px, (float)py, (float)pz);
}

Vector3 calculate_total_angular_momentum(const System& system) {
    double lx = 0.0, ly = 0.0, lz = 0.0;
    const size_t n = system.size();

    #pragma omp parallel for reduction(+:lx, ly, lz)
    for (int i = 0; i < (int)n; ++i) {
        const double m = system.mass[i];
        const double rx = system.x[i];
        const double ry = system.y[i];
        const double rz = system.z[i];
        const double vx = system.vx[i];
        const double vy = system.vy[i];
        const double vz = system.vz[i];

        lx += m * (ry * vz - rz * vy);
        ly += m * (rz * vx - rx * vz);
        lz += m * (rx * vy - ry * vx);
    }

    return Vector3((float)lx, (float)ly, (float)lz);
}

BoundingBox calculate_bounding_box(const System& system) {
    float min_x = std::numeric_limits<float>::max();
    float min_y = std::numeric_limits<float>::max();
    float min_z = std::numeric_limits<float>::max();
    float max_x = std::numeric_limits<float>::lowest();
    float max_y = std::numeric_limits<float>::lowest();
    float max_z = std::numeric_limits<float>::lowest();
    const size_t n = system.size();

    #pragma omp parallel for reduction(min:min_x, min_y, min_z) reduction(max:max_x, max_y, max_z)
    for (int i = 0; i < (int)n; ++i) {
        if (system.x[i] < min_x) min_x = system.x[i];
        if (system.y[i] < min_y) min_y = system.y[i];
        if (system.z[i] < min_z) min_z = system.z[i];
        if (system.x[i] > max_x) max_x = system.x[i];
        if (system.y[i] > max_y) max_y = system.y[i];
        if (system.z[i] > max_z) max_z = system.z[i];
    }

    return {{min_x, min_y, min_z}, {max_x, max_y, max_z}};
}
