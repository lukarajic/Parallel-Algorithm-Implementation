#ifndef NBODY_H
#define NBODY_H

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

#endif // NBODY_H
