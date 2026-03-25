#include <cmath>
#include "nbody.h"

OctreeNode::OctreeNode() 
    : boundary({{0,0,0}, {0,0,0}}), total_mass(0.0f), center_of_mass(0.0f, 0.0f, 0.0f), particle_idx(-1) {
    for (int i = 0; i < 8; ++i) children[i] = nullptr;
}

OctreeNode::OctreeNode(const BoundingBox& box) 
    : boundary(box), total_mass(0.0f), center_of_mass(0.0f, 0.0f, 0.0f), particle_idx(-1) {
    for (int i = 0; i < 8; ++i) children[i] = nullptr;
}

OctreeNode::~OctreeNode() {
    // No recursive delete, pool manages memory
}

void OctreeNode::reset(const BoundingBox& box) {
    boundary = box;
    total_mass = 0.0f;
    center_of_mass = {0.0f, 0.0f, 0.0f};
    particle_idx = -1;
    for (int i = 0; i < 8; ++i) children[i] = nullptr;
}

bool OctreeNode::is_leaf() const {
    for (int i = 0; i < 8; ++i) {
        if (children[i]) return false;
    }
    return true;
}

int OctreeNode::get_octant(const Vector3& pos) const {
    Vector3 center = boundary.center();
    int octant = 0;
    if (pos.x >= center.x) octant |= 1;
    if (pos.y >= center.y) octant |= 2;
    if (pos.z >= center.z) octant |= 4;
    return octant;
}

BoundingBox OctreeNode::create_child_boundary(int octant) const {
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

void OctreeNode::insert(int new_particle_idx, const System& system, OctreePool& pool) {
    Vector3 pos = system.get_pos(new_particle_idx);

    if (!is_leaf()) {
        int octant = get_octant(pos);
        if (!children[octant]) {
            children[octant] = pool.allocate(create_child_boundary(octant));
        }
        children[octant]->insert(new_particle_idx, system, pool);
        return;
    }

    if (particle_idx == -1) {
        particle_idx = new_particle_idx;
        return;
    }

    int existing_idx = particle_idx;
    particle_idx = -1;

    int oct1 = get_octant(system.get_pos(existing_idx));
    if (!children[oct1]) children[oct1] = pool.allocate(create_child_boundary(oct1));
    children[oct1]->insert(existing_idx, system, pool);

    int oct2 = get_octant(pos);
    if (!children[oct2]) children[oct2] = pool.allocate(create_child_boundary(oct2));
    children[oct2]->insert(new_particle_idx, system, pool);
}

void OctreeNode::update_properties(const System& system) {
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
            #pragma omp task shared(system)
            children[i]->update_properties(system);
        }
    }
    #pragma omp taskwait

    for (int i = 0; i < 8; ++i) {
        if (children[i]) {
            total_mass += children[i]->total_mass;
            center_of_mass += children[i]->center_of_mass * children[i]->total_mass;
        }
    }

    if (total_mass > 0) {
        center_of_mass = center_of_mass * (1.0f / total_mass);
    }
}

Vector3 OctreeNode::compute_force(int target_idx, const System& system, float G, float softening, float theta) const {
    if (particle_idx == target_idx) return {0.0f, 0.0f, 0.0f};

    Vector3 pos_target = system.get_pos(target_idx);
    Vector3 r = center_of_mass - pos_target;
    float dist_sq = r.length_sq() + softening;
    float dist = std::sqrt(dist_sq);

    float size = boundary.max.x - boundary.min.x;
    if (is_leaf() || (size / dist < theta)) {
        if (total_mass == 0) return {0.0f, 0.0f, 0.0f};
        float inv_dist = 1.0f / dist;
        float inv_dist3 = inv_dist * inv_dist * inv_dist;
        float f = G * system.mass[target_idx] * total_mass * inv_dist3;
        return r * f;
    }

    Vector3 total_f = {0.0f, 0.0f, 0.0f};
    for (int i = 0; i < 8; ++i) {
        if (children[i]) {
            total_f += children[i]->compute_force(target_idx, system, G, softening, theta);
        }
    }
    return total_f;
}
