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

    Vector3 cross(const Vector3& other) const {
        return {
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        };
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

struct BoundingBox {
    Vector3 min;
    Vector3 max;

    Vector3 center() const {
        return (min + max) * 0.5f;
    }
};

struct System;
struct OctreePool;

struct OctreeNode {
    BoundingBox boundary;
    float total_mass;
    Vector3 center_of_mass;
    int particle_idx;
    OctreeNode* children[8];

    OctreeNode();
    OctreeNode(const BoundingBox& box);
    ~OctreeNode();

    void reset(const BoundingBox& box);
    bool is_leaf() const;
    int get_octant(const Vector3& pos) const;
    BoundingBox create_child_boundary(int octant) const;
    void insert(int new_particle_idx, const System& system, OctreePool& pool);
    void update_properties(const System& system);
    Vector3 compute_force(int target_idx, const System& system, float G, float softening, float theta) const;
};

struct OctreePool {
    std::vector<OctreeNode> nodes;
    size_t next_free;

    OctreePool(size_t capacity) : nodes(capacity), next_free(0) {}

    OctreeNode* allocate(const BoundingBox& box) {
        if (next_free >= nodes.size()) return nullptr;
        OctreeNode* node = &nodes[next_free++];
        node->reset(box);
        return node;
    }

    void reset() {
        next_free = 0;
    }
};

struct System {
    std::vector<float> x, y, z;
    std::vector<float> vx, vy, vz;
    std::vector<float> mass;
    std::vector<float> potential_energy;

    System(int n) : x(n), y(n), z(n), vx(n), vy(n), vz(n), mass(n), potential_energy(n) {}
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
    float theta;
    bool use_barnes_hut;

    void print() const {
        std::cout << "--- Simulation Configuration ---" << std::endl;
        std::cout << "Bodies:     " << num_bodies << std::endl;
        std::cout << "Steps:      " << num_steps << std::endl;
        std::cout << "Time Step:  " << dt << std::endl;
        std::cout << "G:          " << std::scientific << (double)G << std::endl;
        std::cout << "Softening:  " << std::scientific << (double)softening << std::endl;
        std::cout << "Algorithm:  " << (use_barnes_hut ? "Barnes-Hut (O(N log N))" : "Direct Sum (O(N^2))") << std::defaultfloat << std::endl;
        if (use_barnes_hut) {
            std::cout << "Theta:      " << theta << std::endl;
        }
        std::cout << "---------------------------------" << std::endl;
    }

    bool validate() const {
        if (num_bodies <= 0) {
            std::cerr << "Error: Number of bodies must be positive." << std::endl;
            return false;
        }
        if (num_steps <= 0) {
            std::cerr << "Error: Number of steps must be positive." << std::endl;
            return false;
        }
        if (dt <= 0.0f) {
            std::cerr << "Error: Time step must be positive." << std::endl;
            return false;
        }
        return true;
    }

    void summary(double elapsed_seconds) const {
        std::cout << "--- Simulation Summary ---" << std::endl;
        std::cout << "Bodies:      " << num_bodies << std::endl;
        std::cout << "Steps:       " << num_steps << std::endl;
        std::cout << "Total Time:  " << elapsed_seconds << " s" << std::endl;
        if (elapsed_seconds > 0) {
            double interactions = (double)num_bodies * num_bodies * num_steps;
            std::cout << "Performance: " << (interactions / elapsed_seconds) / 1e6 << " million interactions/s" << std::endl;
        }
        std::cout << "--------------------------" << std::endl;
    }
};

#endif // NBODY_H
