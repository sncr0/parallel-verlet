#ifndef THREAD_MANAGER_H
#define THREAD_MANAGER_H

// Struct to manage the number of threads for various tasks
struct ThreadManager {
    int harmonic_bond_threads;   // Number of threads for harmonic bonds
    int lennard_jones_threads;   // Number of threads for Lennard-Jones interactions
    int electrostatic_threads;   // Number of threads for electrostatic forces

    // Default constructor
    ThreadManager()
        : harmonic_bond_threads(1), lennard_jones_threads(1), electrostatic_threads(1) {}

    // Constructor with initial values
    ThreadManager(int hb_threads, int lj_threads, int es_threads)
        : harmonic_bond_threads(hb_threads), lennard_jones_threads(lj_threads), electrostatic_threads(es_threads) {}
};

#endif // THREAD_MANAGER_H
