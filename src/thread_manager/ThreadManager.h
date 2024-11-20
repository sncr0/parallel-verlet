#ifndef THREAD_MANAGER_H
#define THREAD_MANAGER_H

// Struct to manage the number of threads for various tasks
struct ThreadManager {
    int harmonic_bond_threads;   // Number of threads for harmonic bonds
    int lennard_jones_threads;   // Number of threads for Lennard-Jones interactions

    // Default constructor
    ThreadManager()
        : harmonic_bond_threads(1), lennard_jones_threads(1) {}

    // Constructor with initial values
    ThreadManager(int hbThreads, int ljThreads)
        : harmonic_bond_threads(hbThreads), lennard_jones_threads(ljThreads) {}
};

#endif // THREAD_MANAGER_H
