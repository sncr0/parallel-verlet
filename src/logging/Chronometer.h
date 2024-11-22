#ifndef CHRONOMETER_H
#define CHRONOMETER_H

#include <chrono>
#include <unordered_map>
#include <string>
#include <iostream>
#include <iomanip>

struct timing {
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    long long duration;

    timing(std::chrono::high_resolution_clock::time_point start) : start(start) {}
};

class Chronometer {
public:
    // Start timing a section of code with a given tag
    void start(const std::string& tag);

    // End timing a section of code with a given tag
    void end(const std::string& tag);

    // Print all recorded timings
    void printTimings() const;

    void printTiming(const std::string& tag) const; 

private:
    // Map to store start times for tags
    std::unordered_map<std::string, timing> timings_map;

};

#endif // CHRONOMETER_H
