#include <iostream>
#include <cassert>
#include <cmath>
#include "nbody.h"

void test_vector3_ops() {
    Vector3 v1(1.0f, 2.0f, 3.0f);
    Vector3 v2(4.0f, 5.0f, 6.0f);

    // Addition
    Vector3 v_sum = v1 + v2;
    assert(v_sum.x == 5.0f);
    assert(v_sum.y == 7.0f);
    assert(v_sum.z == 9.0f);

    // Subtraction
    Vector3 v_diff = v2 - v1;
    assert(v_diff.x == 3.0f);
    assert(v_diff.y == 3.0f);
    assert(v_diff.z == 3.0f);

    // Scalar multiplication
    Vector3 v_mul = v1 * 2.0f;
    assert(v_mul.x == 2.0f);
    assert(v_mul.y == 4.0f);
    assert(v_mul.z == 6.0f);

    // Length and length_sq
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

int main() {
    test_vector3_ops();
    test_dist_sq();
    std::cout << "All tests passed!" << std::endl;
    return 0;
}
