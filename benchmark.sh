#!/bin/bash

# Script to benchmark molecular simulation with varying threads

# Configuration
EXECUTABLE="./bin/molecular_simulation"  # Path to your executable
OUTPUT_FILE="benchmark_results.csv"      # CSV output file
THREAD_COUNTS=(1 2 4 8 16)               # Number of threads to test
NUM_RUNS=5                               # Number of runs per thread count

# Ensure the output file is empty
echo "Threads,Average_Time(s)" > $OUTPUT_FILE

# Run benchmarks
for threads in "${THREAD_COUNTS[@]}"; do
  total_time=0

  echo "Benchmarking with -e $threads threads..."

  for ((run=1; run<=NUM_RUNS; run++)); do
    # Capture runtime using /usr/bin/time
    runtime=$(/usr/bin/time -f "%e" $EXECUTABLE -n 10000 -n 1 -e $threads 2>&1 | tail -n1)
    echo "Run $run for $threads threads: ${runtime}s"
    total_time=$(echo "$total_time + $runtime" | bc)
  done

  # Calculate average time
  avg_time=$(echo "scale=4; $total_time / $NUM_RUNS" | bc)
  echo "Average runtime for $threads threads: ${avg_time}s"

  # Append results to CSV file
  echo "$threads,$avg_time" >> $OUTPUT_FILE
done

echo "Benchmarking complete. Results saved to $OUTPUT_FILE."
