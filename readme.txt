Project Overview: Molecular Dynamics Simulation Framework

This project implements a molecular dynamics simulation framework, including components for managing the simulation context, applying physical forces, integrating equations of motion, and outputting results.

Compiling the Project with Make

The project includes a Makefile to simplify the compilation process. You can use the make command to build the executables and manage dependencies.

make

Performance Benchmarking

Available Benchmark Scripts

The project provides several benchmark scripts for performance analysis:

- benchmark_bonds.sh
- benchmark_dispersion.sh
- benchmark_electrostatic.sh

Running Benchmarks

# General benchmark execution format
./benchmark_<type> -name

Benchmark Output

- Results are saved to results/data
- Performance visualizations available in results/plot_speedups.ipynb

Executable Options

Electrostatic Interactions Executable

Parallel Execution
./bin/benchmark_electrostatic -n <system_size> -s <steps> -v <verbose> -w <write> -e <threads>

Sequential Execution
./bin/benchmark_electrostatic -n <system_size> -s <steps> -v <verbose> -w <write>
./bin/benchmark_electrostatic -n <system_size> -s <steps> -v <verbose> -w <write> -e 0

Dispersion Interactions Executable

Parallel Execution
./bin/benchmark_dispersion -n <system_size> -s <steps> -v <verbose> -w <write> -d <threads>

Sequential Execution
./bin/benchmark_dispersion -n <system_size> -s <steps> -v <verbose> -w <write>
./bin/benchmark_dispersion -n <system_size> -s <steps> -v <verbose> -w <write> -d 0

Bonded Interactions Executable

Parallel Execution
./bin/benchmark_bonds -n <system_size> -s <steps> -v <verbose> -w <write> -h <threads>

Sequential Execution
./bin/benchmark_bonds -n <system_size> -s <steps> -v <verbose> -w <write>
./bin/benchmark_bonds -n <system_size> -s <steps> -v <verbose> -w <write> -h 0

Command-Line Option Details

Option    Description
-n        System size (number of particles)
-s        Number of simulation steps
-v        Verbose output mode
-w        Write output to file
-e        Threads for electrostatic calculations
-d        Threads for dispersion calculations
-h        Threads for bonded interactions

Example Usage

Parallel Electrostatic Simulation
# Run with 1000 particles, 500 steps, 4 threads
./bin/benchmark_electrostatic -n 1000 -s 500 -e 4

Sequential Electrostatic Simulation
# Run with 1000 particles, 500 steps
./bin/benchmark_electrostatic -n 1000 -s 500
