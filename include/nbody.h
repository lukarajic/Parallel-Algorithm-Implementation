#ifndef NBODY_H
#define NBODY_H

#include <vector>

struct Vector3 {
    float x, y, z;

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
};

#endif // NBODY_H
