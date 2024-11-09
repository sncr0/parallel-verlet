// Verbose.h
#ifndef VERBOSE_H
#define VERBOSE_H

#include <iostream>

extern int vflag;  // Declare the verbosity flag (defined in Verbose.cpp)

// #define VERBOSE(fmt...) do { if (vflag) { std::printf(fmt); std::fflush(stdout); } } while(0)
#define VERBOSE(fmt, ...) do { if (vflag) { std::printf(fmt, ##__VA_ARGS__); std::fflush(stdout); } } while(0)

void setVerboseFlag(int flag);  // Function to set verbosity flag

#endif // VERBOSE_H