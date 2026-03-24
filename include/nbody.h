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

struct OctreeNode {
    BoundingBox boundary;
    float total_mass;
    Vector3 center_of_mass;
    
    // Index of the particle if this is a leaf node, -1 otherwise
    int particle_idx;
    
    // Pointers to the 8 children (octants)
    OctreeNode* children[8];

    OctreeNode(const BoundingBox& box) : boundary(box), total_mass(0.0f), center_of_mass(0.0f, 0.0f, 0.0f), particle_idx(-1) {
        for (int i = 0; i < 8; ++i) children[i] = nullptr;
    }

    ~OctreeNode() {
        for (int i = 0; i < 8; ++i) {
            if (children[i]) delete children[i];
        }
    }

    bool is_leaf() const {
        for (int i = 0; i < 8; ++i) {
            if (children[i]) return false;
        }
        return true;
    }

    int get_octant(const Vector3& pos) const {
        Vector3 center = boundary.center();
        int octant = 0;
        if (pos.x >= center.x) octant |= 1;
        if (pos.y >= center.y) octant |= 2;
        if (pos.z >= center.z) octant |= 4;
        return octant;
    }

    BoundingBox create_child_boundary(int octant) const {
        Vector3 center = boundary.center();
        BoundingBox child_box;
        child_box.min.x = (octant & 1) ? center.x : boundary.min.x;
        child_box.max.x = (octant & 1) ? boundary.max.x : center.x;
        child_box.min.y = (octant & 2) ? center.y : boundary.min.y;
        child_box.max.y = (octant & 2) ? boundary.max.y : center.y;
        child_box.min.z = (octant & 4) ? center.z : boundary.min.z;
        child_box.max.z = (octant & 4) ? boundary.max.z : center.z;
        return child_box;
    }

    void insert(int new_particle_idx, const System& system) {
        Vector3 pos = system.get_pos(new_particle_idx);

        // Case 1: Internal node
        if (!is_leaf()) {
            int octant = get_octant(pos);
            if (!children[octant]) {
                children[octant] = new OctreeNode(create_child_boundary(octant));
            }
            children[octant]->insert(new_particle_idx, system);
            return;
        }

        // Case 2: Empty leaf
        if (particle_idx == -1) {
            particle_idx = new_particle_idx;
            return;
        }

        // Case 3: Leaf already contains a particle - must split
        int existing_idx = particle_idx;
        particle_idx = -1; // Mark as internal

        // Move existing particle
        int oct1 = get_octant(system.get_pos(existing_idx));
        if (!children[oct1]) children[oct1] = new OctreeNode(create_child_boundary(oct1));
        children[oct1]->insert(existing_idx, system);

        // Insert new particle
        int oct2 = get_octant(pos);
        if (!children[oct2]) children[oct2] = new OctreeNode(create_child_boundary(oct2));
        children[oct2]->insert(new_particle_idx, system);
    }

    void update_properties(const System& system) {
        if (is_leaf()) {
            if (particle_idx != -1) {
                total_mass = system.mass[particle_idx];
                center_of_mass = system.get_pos(particle_idx);
            }
            return;
        }

        total_mass = 0.0f;
        center_of_mass = {0.0f, 0.0f, 0.0f};

        for (int i = 0; i < 8; ++i) {
            if (children[i]) {
                children[i]->update_properties(system);
                total_mass += children[i]->total_mass;
                center_of_mass += children[i]->center_of_mass * children[i]->total_mass;
            }
        }

        if (total_mass > 0) {
            center_of_mass = center_of_mass * (1.0f / total_mass);
        }
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

    void print() const {
        std::cout << "--- Simulation Configuration ---" << std::endl;
        std::cout << "Bodies:     " << num_bodies << std::endl;
        std::cout << "Steps:      " << num_steps << std::endl;
        std::cout << "Time Step:  " << dt << std::endl;
        std::cout << "G:          " << std::scientific << (double)G << std::endl;
        std::cout << "Softening:  " << std::scientific << (double)softening << std::defaultfloat << std::endl;
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
