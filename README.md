# High-Performance Parallel N-Body Simulation

A high-performance C++ implementation of the N-Body problem, optimized for multi-core CPUs. This project demonstrates advanced skills in **Parallel Computing**, **Algorithmic Optimization**, and **Memory Architecture**.

## 🚀 Key Technical Highlights

-   **Algorithmic Optimization**: Implements both the $O(N^2)$ Direct Sum and the $O(N \log N)$ **Barnes-Hut** algorithm using a recursive Octree.
-   **Multi-threading (OpenMP)**: Leverages OpenMP for data-parallel force calculations and task-based parallel tree property updates.
-   **Memory Layout (SoA)**: Utilizes a **Structure of Arrays (SoA)** memory layout instead of AoS to maximize L1/L2 cache locality and enable seamless SIMD vectorization.
-   **SIMD Vectorization**: Explicitly utilizes `#pragma omp simd` hints to leverage wide CPU registers (AVX-2 / AVX-512).
-   **Custom Memory Management**: Implements an `OctreePool` (Memory Pool) to eliminate the overhead of heap allocations during the simulation loop.
-   **Physical Validation**: Includes robust monitoring for conservation of **Energy**, **Linear Momentum**, and **Angular Momentum**.

## 🛠 Build and Run

### Prerequisites
-   CMake (3.10+)
-   C++17 Compiler (GCC, Clang, or AppleClang)
-   OpenMP Library

### Building with Makefile
The project includes a professional Makefile wrapper:
```bash
make build
```

### Running Benchmarks
Compare the performance of Direct Sum vs. Barnes-Hut:
```bash
make benchmark
```

### Manual Execution
```bash
./build/nbody [num_bodies] [num_steps] [--barnes-hut] [--theta T] [--verbose]
```

## 📊 Performance Comparison (20,000 Bodies)

| Algorithm | Complexity | Total Time (50 steps) |
| :--- | :--- | :--- |
| **Direct Sum** | $O(N^2)$ | ~2.08s |
| **Barnes-Hut** | $O(N \log N)$ | **~1.40s** |

*Benchmarks conducted on an 8-core CPU. See `docs/ARCHITECTURE.md` for more details on optimizations.*

## 🧪 Testing
The project includes a suite of unit tests for core vector math and physical validation:
```bash
make test
```

## 📂 Project Structure
-   `include/`: Header files (Data structures, Utilities, Constants)
-   `src/`: Implementation files (Simulation logic, Tree logic, Init)
-   `tests/`: Unit testing suite
-   `docs/`: Detailed architectural documentation
-   `scripts/`: Benchmarking and analysis tools
