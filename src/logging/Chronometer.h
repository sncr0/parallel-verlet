#ifndef CHRONOMETER_H
#define CHRONOMETER_H

#include <chrono>
#include <unordered_map>
#include <string>
#include <iostream>
#include <iomanip>

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
    void start(const std::string& tag);

    void end(const std::string& tag);

    void printTiming(const std::string& tag, const std::string& timescale) const;

private:
    std::unordered_map<std::string, Timing> timings_map;

};

#endif // CHRONOMETER_H
