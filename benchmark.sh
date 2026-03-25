#!/bin/bash

if [[ "$1" == "--clean" ]]; then
    echo "Cleaning build directory..."
    rm -rf build
    exit 0
fi

# Build the project first
echo "Building project..."
mkdir -p build && cd build && cmake .. > /dev/null && make -j > /dev/null
cd ..

# Define parameters
STEPS=20
BODIES=(5000 10000 20000)

echo "--------------------------------------------------------"
echo "| Bodies | Algorithm | Total Time (s) | Performance (M/s) |"
echo "--------------------------------------------------------"

for N in "${BODIES[@]}"; do
    # Run Direct Sum
    OUTPUT=$(./build/nbody $N $STEPS)
    TIME=$(echo "$OUTPUT" | grep "Total Time" | awk '{print $3}')
    PERF=$(echo "$OUTPUT" | grep "Performance" | awk '{print $2}')
    printf "| %6d | Direct    | %14s | %17s |\n" "$N" "$TIME" "$PERF"

    # Run Barnes-Hut
    OUTPUT=$(./build/nbody $N $STEPS --barnes-hut)
    TIME=$(echo "$OUTPUT" | grep "Total Time" | awk '{print $3}')
    PERF=$(echo "$OUTPUT" | grep "Performance" | awk '{print $2}')
    printf "| %6d | Barnes-Hut| %14s | %17s |\n" "$N" "$TIME" "$PERF"
    echo "--------------------------------------------------------"
done
