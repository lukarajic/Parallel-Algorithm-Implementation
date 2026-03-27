.PHONY: all build test clean benchmark sweep

all: build

build:
	mkdir -p build && cd build && cmake .. && make -j

test: build
	cd build && ctest --output-on-failure

benchmark: build
	./benchmark.sh

sweep: build
	./theta_sweep.sh

clean:
	rm -rf build
	./benchmark.sh --clean
