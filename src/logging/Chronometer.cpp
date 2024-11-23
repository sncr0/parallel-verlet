#include "Chronometer.h"

// Start timing a section of code
void Chronometer::start(const std::string& tag) {
    timings_map[tag] = Timing();
}

// End timing a section of code
void Chronometer::end(const std::string& tag) {
    auto endTime = std::chrono::high_resolution_clock::now();

    // Check if the tag exists in startTimes
    if (timings_map.find(tag) != timings_map.end()) {
        timings_map[tag].end = endTime;
        timings_map[tag].duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - timings_map[tag].start);
        timings_map[tag].duration_us = std::chrono::duration_cast<std::chrono::microseconds>(endTime - timings_map[tag].start);
        timings_map[tag].duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - timings_map[tag].start);
        timings_map[tag].duration_s = std::chrono::duration_cast<std::chrono::seconds>(endTime - timings_map[tag].start);
    } else {
        std::cerr << "Warning: Tag \"" << tag << "\" was not started!\n";
    }
}

void Chronometer::printTiming(const std::string& tag, const std::string& timescale) const {
    long print_time;
    for (auto& entry : timings_map) {
        if (timescale == "ns") {
            print_time = entry.second.duration_ns.count();
        }
        else if (timescale == "us") {
            print_time = entry.second.duration_us.count();
        }
        else if (timescale == "ms") {
            print_time = entry.second.duration_ms.count();
        }
        else if (timescale == "s") {
            print_time = entry.second.duration_s.count();
        }
        else {
            print_time = entry.second.duration_ms.count();
        }
    }

        printf("\nTiming for section \"%s\": %ld %s\n", tag.c_str(), print_time, timescale.c_str());
}


