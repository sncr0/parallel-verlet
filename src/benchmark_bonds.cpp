// main.cpp

#include "system/System.h"
#include "context/Context.h"
#include "integrator/VerletIntegrator.h"
#include "forces/LennardJonesForce.h"
#include "forces/HarmonicBondForce.h"
#include "forces/ElectrostaticForce.h"
#include <iostream>
#include "io/XYZWriter.h"
#include "getopt.h"
#include "logging/Verbose.h"
#include "thread_manager/ThreadManager.h"
#include <bits/stdc++.h>
#include "logging/Chronometer.h"



int main(int argc, char **argv) {
    int c;
    int harmonic_bond_threads = 0;
    int electrostatic_bond_threads = 0;
    int dispersion_force_threads = 0;
    int num_atoms = 1;
    int num_steps = 1;
    int write_flag = 0;

    while ((c = getopt(argc, argv, "vh:e:n:d:s:w")) != -1) {
        switch (c) {
            case 'v':
                setVerboseFlag(1);
                break;
            case 'h': {
                try {
                    int value = std::stoi(optarg);
                    if (value < 0) {
                        throw std::invalid_argument("non-positive value");
                    }
                    harmonic_bond_threads = value;
                    printf("harmonic_bond_threads: %d\n", harmonic_bond_threads);
                } catch (const std::exception& e) {
                    std::cerr << "Invalid number of threads for -h: " << optarg << "\n";
                    return EXIT_FAILURE;
                }
                break;
            }
            case 'e': {
                try {
                    int value = std::stoi(optarg);
                    if (value < 0) {
                        throw std::invalid_argument("non-positive value");
                    }
                    electrostatic_bond_threads = value;
                    printf("electrostatic_bond_threads: %d\n", electrostatic_bond_threads);
                } catch (const std::exception& e) {
                    std::cerr << "Invalid number of threads for -e: " << optarg << "\n";
                    return EXIT_FAILURE;
                }
                break;
            }
            case 'd': {
                try {
                    int value = std::stoi(optarg);
                    if (value < 0) {
                        throw std::invalid_argument("non-positive value");
                    }
                    dispersion_force_threads = value;
                    printf("dispersion_force_threads: %d\n", dispersion_force_threads);
                } catch (const std::exception& e) {
                    std::cerr << "Invalid number of threads for -d: " << optarg << "\n";
                    return EXIT_FAILURE;
                }
                break;
            }
            case 'n': {
                try {
                    int value = std::stoi(optarg);
                    if (value <= 0) {
                        throw std::invalid_argument("non-positive value");
                    }
                    num_atoms = value;
                    printf("num_atoms: %d\n", num_atoms);
                } catch (const std::exception& e) {
                    std::cerr << "Invalid number of atoms for -n: " << optarg << "\n";
                    return EXIT_FAILURE;
                }
                break;
            }
            case 's': {
                try {
                    int value = std::stoi(optarg);
                    if (value <= 0) {
                        throw std::invalid_argument("non-positive value");
                    }
                    num_steps = value;
                    printf("num_steps: %d\n", num_steps);
                } catch (const std::exception& e) {
                    std::cerr << "Invalid number of steps for -s: " << optarg << "\n";
                    return EXIT_FAILURE;
                }
                break;
            }
            case 'w':
                write_flag = 1;
                break;
            default:
                std::cerr << "Unknown option: " << c << "\n";
                return EXIT_FAILURE;
        }
    }


    printf("Benchmarking harmonic bond system\n");
    Chronometer chronometer;
    ThreadManager thread_manager(harmonic_bond_threads, dispersion_force_threads, electrostatic_bond_threads);

    VERBOSE("Creating the system\n");
    System system;

    VerletIntegrator integrator(0.01, thread_manager, chronometer);

    VERBOSE("Adding forces to the system\n");
    auto hbForce1 = std::make_shared<HarmonicBondForce>(1.0, 1.0);
    hbForce1->num_threads = harmonic_bond_threads;

    VERBOSE("Adding particles to the system\n");
    for (int i = 0; i < num_atoms; ++i) {
        int x = rand()%500;
        int y = rand()%500;
        int z = rand()%500;
        system.addParticle(1.0, -1.0, x, y, z);
        system.addParticle(1.0, 1.0, x+1, y+1, z+1);
        hbForce1->addBond(0, 1);
        VERBOSE("Particle %d added at position (%d, %d, %d)\n", i, x, y, z);
        VERBOSE("Particle %d added at position (%d, %d, %d)\n", i+1, x+1, y+1, z+1);
        VERBOSE("Bond added between particles %d and %d\n", i, i+1);
    }

    integrator.addForce(hbForce1);

    Context context(system, integrator);

    VERBOSE("Starting the simulation with %zu particles\n", system.getNumParticles());

    if (write_flag) {
        VERBOSE("Writing trajectory to trajectory.xyz\n");
        XYZWriter trajectoryWriter("trajectory.xyz");
        const int output_interval = 100; // Output every 100 steps
        for (int step = 0; step < num_steps; ++step) {
            context.runSimulation(1);
            if (step % output_interval == 0) {
                trajectoryWriter.writeFrame(system); // Write the current frame
            }
        }
        trajectoryWriter.close();
    } else {
        for (int step = 0; step < num_steps; ++step) {
            context.runSimulation(1);
        }
    }

    return 0;
}
