#ifndef CHRONOMETER_H
#define CHRONOMETER_H

#include <chrono>
#include <unordered_map>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <shared_mutex> // For std::shared_mutex
#include <mutex>

struct Timing {
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::nanoseconds duration_ns;
    std::chrono::microseconds duration_us;   
    std::chrono::milliseconds duration_ms;
    std::chrono::seconds duration_s;

    Timing() : start(std::chrono::high_resolution_clock::now()) {}
};

class Chronometer {
public:
    void start(const std::string& tag, int thread_id = 0);

    void end(const std::string& tag, int thread_id = 0);

    void printTiming(const std::string& tag, const std::string& timescale, int thread_id = 0) const;

private:
    // through timepoints and over threads, should maybe be dynamically allocated for number of threads
    // map < tag, map < thread_id, vector < Timing > > >
    std::unordered_map<std::string, std::unordered_map<int, std::vector<Timing>>> timings_map;
    mutable std::shared_mutex mtx_; // Shared mutex for thread-safe read/write


};

#endif // CHRONOMETER_H
