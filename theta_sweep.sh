#!/bin/bash

# Build the project first
echo "Building project..."
mkdir -p build && cd build && cmake .. > /dev/null && make -j > /dev/null
cd ..

# Define parameters
N=20000
STEPS=20
THETAS=(0.3 0.5 0.7 0.9)

echo "--------------------------------------------------------"
echo "| Theta  | Total Time (s) | Performance (M/s) |"
echo "--------------------------------------------------------"

for T in "${THETAS[@]}"; do
    OUTPUT=$(./build/nbody $N $STEPS --barnes-hut --theta $T)
    TIME=$(echo "$OUTPUT" | grep "Total Time" | awk '{print $3}')
    PERF=$(echo "$OUTPUT" | grep "Performance" | awk '{print $2}')
    printf "| %6.1f | %14s | %17s |\n" "$T" "$TIME" "$PERF"
done
echo "--------------------------------------------------------"
