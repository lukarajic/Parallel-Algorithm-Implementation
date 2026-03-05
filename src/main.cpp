#include <iostream>
#include <omp.h>

int main() {
    std::cout << "Parallel N-Body Simulation" << std::endl;

    #ifdef _OPENMP
        std::cout << "OpenMP is supported!" << std::endl;
        #pragma omp parallel
        {
            #pragma omp single
            std::cout << "Running with " << omp_get_num_threads() << " threads." << std::endl;
        }
    #else
        std::cout << "OpenMP is NOT supported." << std::endl;
    #endif

    return 0;
}
