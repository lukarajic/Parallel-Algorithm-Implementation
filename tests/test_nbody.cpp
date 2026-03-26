#include <iostream>
#include <cassert>
#include <cmath>
#include "nbody.h"
#include "utils.h"

void test_vector3_ops() {
    Vector3 v1(1.0f, 2.0f, 3.0f);
    Vector3 v2(4.0f, 5.0f, 6.0f);

    assert(v1.x == 1.0f);
    Vector3 v_sum = v1 + v2;
    assert(v_sum.x == 5.0f);
    assert(v_sum.y == 7.0f);
    assert(v_sum.z == 9.0f);

    Vector3 v_diff = v2 - v1;
    assert(v_diff.x == 3.0f);
    assert(v_diff.y == 3.0f);
    assert(v_diff.z == 3.0f);

    Vector3 v_mul = v1 * 2.0f;
    assert(v_mul.x == 2.0f);
    assert(v_mul.y == 4.0f);
    assert(v_mul.z == 6.0f);

    Vector3 v_len(3.0f, 4.0f, 0.0f);
    assert(v_len.length_sq() == 25.0f);
    assert(std::abs(v_len.length() - 5.0f) < 1e-6f);

    std::cout << "test_vector3_ops passed!" << std::endl;
}

void test_dist_sq() {
    float x1 = 0.0f, y1 = 0.0f, z1 = 0.0f;
    float x2 = 1.0f, y2 = 2.0f, z2 = 2.0f;
    float d2 = Vector3::dist_sq(x1, y1, z1, x2, y2, z2);
    assert(std::abs(d2 - 9.0f) < 1e-6f);

    std::cout << "test_dist_sq passed!" << std::endl;
}

void test_system_utils() {
    System system(2);
    system.x[0] = 0.0f; system.y[0] = 0.0f; system.z[0] = 0.0f;
    system.vx[0] = 0.0f; system.vy[0] = 0.0f; system.vz[0] = 0.0f;
    system.mass[0] = 10.0f;

    system.x[1] = 10.0f; system.y[1] = 0.0f; system.z[1] = 0.0f;
    system.vx[1] = 1.0f; system.vy[1] = 0.0f; system.vz[1] = 0.0f;
    system.mass[1] = 10.0f;

    // Total Mass
    assert(std::abs(calculate_total_mass(system) - 20.0) < 1e-6);

    // Center of Mass
    Vector3 cm = calculate_center_of_mass(system);
    assert(std::abs(cm.x - 5.0f) < 1e-6f);
    assert(std::abs(cm.y - 0.0f) < 1e-6f);
    assert(std::abs(cm.z - 0.0f) < 1e-6f);

    // Total Momentum
    Vector3 p = calculate_total_momentum(system);
    assert(std::abs(p.x - 10.0f) < 1e-6f);
    assert(std::abs(p.y - 0.0f) < 1e-6f);
    assert(std::abs(p.z - 0.0f) < 1e-6f);

    // Bounding Box
    BoundingBox bb = calculate_bounding_box(system);
    assert(bb.min.x == 0.0f);
    assert(bb.max.x == 10.0f);
    assert(bb.min.y == 0.0f);
    assert(bb.max.y == 0.0f);

    std::cout << "test_system_utils passed!" << std::endl;
}

void test_energy() {
    System system(2);
    Config config = {2, 1, 0.01f, 1.0f, 0.0f, 0.5f, false};
    
    // Body 1: At origin, moving at 1m/s in X, mass 1kg
    system.x[0] = 0.0f; system.y[0] = 0.0f; system.z[0] = 0.0f;
    system.vx[0] = 1.0f; system.vy[0] = 0.0f; system.vz[0] = 0.0f;
    system.mass[0] = 1.0f;

    // Body 2: At (1,0,0), stationary, mass 1kg
    system.x[1] = 1.0f; system.y[1] = 0.0f; system.z[1] = 0.0f;
    system.vx[1] = 0.0f; system.vy[1] = 0.0f; system.vz[1] = 0.0f;
    system.mass[1] = 1.0f;

    // KE = 0.5 * 1 * 1^2 + 0.5 * 1 * 0^2 = 0.5
    // PE = -G * m1 * m2 / dist = -1.0 * 1 * 1 / 1 = -1.0
    // Total = 0.5 - 1.0 = -0.5
    
    double ke = calculate_kinetic_energy(system, config);
    double pe = calculate_potential_energy(system, config);
    double total = calculate_total_energy(system, config);

    assert(std::abs(ke - 0.5) < 1e-6);
    assert(std::abs(pe - (-1.0)) < 1e-6);
    assert(std::abs(total - (-0.5)) < 1e-6);

    std::cout << "test_energy passed!" << std::endl;
}

void test_momentum() {
    System system(2);
    
    // Body 1: Mass 2kg, velocity (1, 2, 3)
    system.x[0] = 0.0f; system.y[0] = 0.0f; system.z[0] = 0.0f;
    system.vx[0] = 1.0f; system.vy[0] = 2.0f; system.vz[0] = 3.0f;
    system.mass[0] = 2.0f;

    // Body 2: Mass 3kg, velocity (-1, 0, 1)
    system.x[1] = 1.0f; system.y[1] = 0.0f; system.z[1] = 0.0f;
    system.vx[1] = -1.0f; system.vy[1] = 0.0f; system.vz[1] = 1.0f;
    system.mass[1] = 3.0f;

    // P = m1*v1 + m2*v2
    // Px = 2*1 + 3*(-1) = 2 - 3 = -1
    // Py = 2*2 + 3*0 = 4
    // Pz = 2*3 + 3*1 = 6 + 3 = 9
    
    Vector3 p = calculate_total_momentum(system);

    assert(std::abs(p.x - (-1.0f)) < 1e-6f);
    assert(std::abs(p.y - 4.0f) < 1e-6f);
    assert(std::abs(p.z - 9.0f) < 1e-6f);

    std::cout << "test_momentum passed!" << std::endl;
}

int main() {
    test_vector3_ops();
    test_dist_sq();
    test_system_utils();
    test_energy();
    test_momentum();
    std::cout << "All tests passed!" << std::endl;
    return 0;
}
