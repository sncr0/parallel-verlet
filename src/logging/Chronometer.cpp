#include "Chronometer.h"

// Start timing a section of code
void Chronometer::start(const std::string& tag, int thread_id) {
    timings_map[tag][thread_id].push_back(Timing());
}

// End timing a section of code
void Chronometer::end(const std::string& tag, int thread_id) {
    auto endTime = std::chrono::high_resolution_clock::now();

    // Check if the tag exists in startTimes
    if (timings_map.find(tag) != timings_map.end() &&
        timings_map[tag].find(thread_id) != timings_map[tag].end() &&
        !timings_map[tag][thread_id].empty()) {
            
        timings_map[tag][thread_id].back().end = endTime;
        timings_map[tag][thread_id].back().duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - timings_map[tag][thread_id].back().start);
        timings_map[tag][thread_id].back().duration_us = std::chrono::duration_cast<std::chrono::microseconds>(endTime - timings_map[tag][thread_id].back().start);
        timings_map[tag][thread_id].back().duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - timings_map[tag][thread_id].back().start);
        timings_map[tag][thread_id].back().duration_s = std::chrono::duration_cast<std::chrono::seconds>(endTime - timings_map[tag][thread_id].back().start);
    } else {
        std::cerr << "Warning: Tag \"" << tag << "\" was not started!\n";
    }
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

