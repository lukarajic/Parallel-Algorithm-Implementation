# Parallel N-Body Simulation: Architectural Design

This project implements a high-performance N-Body simulation in C++ using OpenMP. The architecture is designed to demonstrate key skills in High-Performance Computing (HPC) and parallel algorithm optimization.

## 1. Memory Layout: Structure of Arrays (SoA)

Unlike a traditional "Array of Structures" (AoS) where each body is an object `Body { x, y, z, vx, vy, vz, mass }`, this project uses a **Structure of Arrays (SoA)** layout.

- **Design**: All X-coordinates are in one contiguous vector, all Y-coordinates in another, and so on.
- **Benefit**: This is highly cache-friendly. When calculating forces, the CPU pre-fetches contiguous blocks of data into the L1/L2 caches. It also allows for efficient **SIMD Vectorization**, as data can be loaded directly into wide registers (AVX-2/AVX-512) without shuffling.

## 2. Algorithms

### Direct Sum ($O(N^2)$)
The baseline implementation calculates forces between every pair of particles.
- **Optimization**: Parallelized using `omp parallel for` and explicit `omp simd` hints for the inner loops.

### Barnes-Hut ($O(N \log N)$)
For larger simulations, we implement the Barnes-Hut optimization.
- **Spatial Partitioning**: A recursive **Octree** divides 3D space into octants.
- **Multipole Expansion**: Distant clusters of particles are approximated as a single point mass at their center of mass.
- **Complexity**: Reduces the number of interactions from $N^2$ to $N \log N$, allowing simulations with 100,000+ particles.

## 3. Memory Management: Octree Pool

A common bottleneck in tree-based algorithms is frequent heap allocation (`new`/`delete`). 

- **Optimization**: We implement an `OctreePool` that pre-allocates a large contiguous block of `OctreeNode` objects.
- **Result**: Tree construction is significantly faster as it only involves incrementing a pointer and resetting node state, completely eliminating heap fragmentation and allocation overhead during the simulation loop.

## 4. Parallelization Strategy

- **OpenMP Tasks**: Used for the recursive tree property updates (mass and center-of-mass calculation).
- **Data Parallelism**: Used for the force calculation and position integration phases.
- **Thread Safety**: RAII-based `Timer` and a thread-safe `Logger` ensure correct behavior in multi-threaded environments.
