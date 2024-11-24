#include "Chronometer.h"


// Start timing a section of code
void Chronometer::start(const std::string& tag, int thread_id) {
    std::unique_lock lock(mtx_); // Write lock

    // First, ensure the tag exists and get a reference to its thread map
    auto& thread_map = timings_map[tag];
    // Then, add the timing to the specific thread's vector
    thread_map[thread_id].push_back(Timing());
}


void Chronometer::end(const std::string& tag, int thread_id) {
    std::unique_lock lock(mtx_); // Write lock

    auto it = timings_map.find(tag);
    if (it == timings_map.end()) {
        std::cerr << "Warning: Tag \"" << tag << "\" was not started!\n";
        return;
    }

    auto thread_it = it->second.find(thread_id);
    if (thread_it == it->second.end() || thread_it->second.empty()) {
        std::cerr << "Warning: No timing data for tag \"" << tag << "\" on thread " << thread_id << "!\n";
        return;
    }

    auto& timing = thread_it->second.back();
    timing.end = std::chrono::high_resolution_clock::now();
    
    // Calculate durations
    auto duration = timing.end - timing.start;
    timing.duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
    timing.duration_us = std::chrono::duration_cast<std::chrono::microseconds>(duration);
    timing.duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
    timing.duration_s = std::chrono::duration_cast<std::chrono::seconds>(duration);
}

void Chronometer::printTiming(const std::string& tag, const std::string& timescale, int thread_id) const {
    // Check if the tag exists in the hashmap
    auto tag_it = timings_map.find(tag);
    if (tag_it == timings_map.end()) {
        // Tag not found, print an error and return
        std::cerr << "Error: No timing information found for tag \"" << tag << "\".\n";
        return;
    }

    // Get the vector of Timing objects for the tag
    const auto& thread_timings = tag_it->second;
    if (thread_timings.empty()) {
        std::cerr << "Error: No recorded timings for tag \"" << tag << "\".\n";
        return;
    }

    auto thread_it = thread_timings.find(thread_id);

    if (thread_it == thread_timings.end()) {
        // Tag not found, print an error and return
        std::cerr << "Error: No timing information found for tag \"" << tag << "\" on thread " << thread_id << ".\n";
        return;
    }

    const auto& timings = thread_it->second;
    if (timings.empty()) {
        std::cerr << "Error: No recorded timings for tag \"" << tag << "\" on thread " << thread_id << ".\n";
        return;
    }

    const auto& latest_timing = timings.back();


    // Get the most recent Timing object
    // if (thread_id < thread_timings.size()) {
    //     const auto& latest_timing = thread_timings[thread_id].back();
    //     // Do something with latest_timing
    // }
    // // const auto& latest_timing = timings[thread_id].back();

    // Determine the correct time unit
    long long print_time = 0;
    if (timescale == "ns") {
        print_time = latest_timing.duration_ns.count();
    } else if (timescale == "us") {
        print_time = latest_timing.duration_us.count();
    } else if (timescale == "ms") {
        print_time = latest_timing.duration_ms.count();
    } else if (timescale == "s") {
        print_time = latest_timing.duration_s.count();
    } else {
        // Invalid timescale, default to milliseconds
        std::cerr << "Warning: Invalid timescale \"" << timescale << "\", defaulting to \"ms\".\n";
        print_time = latest_timing.duration_ms.count();
    }

    // Print the timing information
    printf("Timing for section \"%s\" (thread %d): %lld %s\n", tag.c_str(), thread_id, print_time, timescale.c_str());
}

