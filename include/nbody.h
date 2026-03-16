#ifndef NBODY_H
#define NBODY_H

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

struct Vector3 {
    float x, y, z;

    Vector3() : x(0.0f), y(0.0f), z(0.0f) {}
    Vector3(float x, float y, float z) : x(x), y(y), z(z) {}

    Vector3 operator+(const Vector3& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    Vector3& operator+=(const Vector3& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    Vector3 operator-(const Vector3& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    Vector3 operator*(float scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    float length_sq() const {
        return x * x + y * y + z * z;
    }

    float length() const {
        return std::sqrt(length_sq());
    }

    static inline float dist_sq(float x1, float y1, float z1, float x2, float y2, float z2) {
        float dx = x2 - x1;
        float dy = y2 - y1;
        float dz = z2 - z1;
        return dx * dx + dy * dy + dz * dz;
    }
};

struct Body {
    Vector3 position;
    Vector3 velocity;
    float mass;
};

struct System {
    std::vector<float> x, y, z;
    std::vector<float> vx, vy, vz;
    std::vector<float> mass;

    System(int n) : x(n), y(n), z(n), vx(n), vy(n), vz(n), mass(n) {}
    size_t size() const { return x.size(); }

    Vector3 get_pos(size_t i) const { return {x[i], y[i], z[i]}; }
    Vector3 get_vel(size_t i) const { return {vx[i], vy[i], vz[i]}; }
    
    void set_pos(size_t i, const Vector3& v) { x[i] = v.x; y[i] = v.y; z[i] = v.z; }
    void set_vel(size_t i, const Vector3& v) { vx[i] = v.x; vy[i] = v.y; vz[i] = v.z; }
};

struct Config {
    int num_bodies;
    int num_steps;
    float dt;
    float G;
    float softening;

    void print() const {
        std::cout << "--- Simulation Configuration ---" << std::endl;
        std::cout << "Bodies:     " << num_bodies << std::endl;
        std::cout << "Steps:      " << num_steps << std::endl;
        std::cout << "Time Step:  " << dt << std::endl;
        std::cout << "G:          " << std::scientific << (double)G << std::endl;
        std::cout << "Softening:  " << std::scientific << (double)softening << std::defaultfloat << std::endl;
        std::cout << "---------------------------------" << std::endl;
    }
};

#endif // NBODY_H
