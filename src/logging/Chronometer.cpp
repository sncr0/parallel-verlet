#include "Chronometer.h"

// Start timing a section of code
void Chronometer::start(const std::string& tag) {
    timings_map[tag] = std::chrono::high_resolution_clock::now();
}

// End timing a section of code
void Chronometer::end(const std::string& tag) {
    auto endTime = std::chrono::high_resolution_clock::now();

    // Check if the tag exists in startTimes
    if (timings_map.find(tag) != timings_map.end()) {
        timings_map[tag].end = endTime;
        timings_map[tag].duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - timings_map[tag].start).count();
    } else {
        std::cerr << "Warning: Tag \"" << tag << "\" was not started!\n";
    }
}

// Print all recorded timings
void Chronometer::printTimings() const {
    printf("\nRecorded Timings:\n");
    printf("%-20s%-15s\n", "Section", "Time (ms)");
    printf("%s\n", std::string(35, '-').c_str());

    for (const auto& entry : timings_map) {
        printf("%-20s%-15lld\n", entry.first.c_str(), entry.second.duration);
    }
}

void Chronometer::printTiming(const std::string& tag) const {
    if (timings_map.find(tag) != timings_map.end()) {
        printf("\nTiming for section \"%s\":\n", tag.c_str());
        printf("%-15s%-15s\n", "Start", "End");
        printf("%s\n", std::string(30, '-').c_str());
        printf("%-15lld%-15lld\n", timings_map.at(tag).start.time_since_epoch().count(), timings_map.at(tag).end.time_since_epoch().count());
    } else {
        std::cerr << "Warning: Tag \"" << tag << "\" was not started!\n";
    }
}


