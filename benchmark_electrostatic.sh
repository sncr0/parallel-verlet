#!/bin/bash


if [ -z "$1" ]; then
  echo "Usage: $0 <name>"
  exit 1
fi
name="$1"

perf='/usr/lib/linux-tools/5.15.0-126-generic/perf'

# Configuration
EXECUTABLE="./bin/molecular_simulation"    # Path to your executable
OUTPUT_FILE="results/data/${name}_benchmark_results.csv"   # CSV output file
E_VALUES=(0 1 2 4 8 12 16 32 64 128 256)                      # Values for the -e flag (threads)
N_VALUES=(1000 5000 10000 25000 50000 100000 1000000)                 # Values for the -n flag (size)
REPEATS=3                                 # Number of repetitions per configuration

# Initialize CSV file with headers
echo -n "Size/Threads" > $OUTPUT_FILE
for e in "${E_VALUES[@]}"; do
  echo -n ",$e" >> $OUTPUT_FILE
done
echo >> $OUTPUT_FILE

# Run benchmarks and populate the grid
for n in "${N_VALUES[@]}"; do
  echo -n "$n" >> $OUTPUT_FILE
  for e in "${E_VALUES[@]}"; do
    echo "Benchmarking for -n $n and -e $e..."


    # Run perf with repeats and capture elapsed times
    # perf_output=$(/usr/lib/linux-tools/5.15.0-126-generic/perf stat -e cache-misses,cache-references --repeat $REPEATS $EXECUTABLE -n $n -s 1 -e $e 2>&1)
    perf_output=$(perf stat -e cache-misses,cache-references --repeat $REPEATS $EXECUTABLE -n $n -s 1 -e $e 2>&1)

    # Extract average elapsed time and append it to the grid
    avg_time=$(echo "$perf_output" | awk '/seconds user/ {sum += $1; count++} END {if (count > 0) print sum / count}')
    echo -n ",$avg_time" >> $OUTPUT_FILE
  done
  echo >> $OUTPUT_FILE
done

echo "Grid benchmarking complete. Results saved to $OUTPUT_FILE."