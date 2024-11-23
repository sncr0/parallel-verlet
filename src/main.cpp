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
    int harmonic_bond_threads = 1;
    int electrostatic_bond_threads = 1;

    while ((c = getopt(argc, argv, "vh:e:")) != -1) {
        switch (c) {
            case 'v':
                setVerboseFlag(1);
                break;
            case 'h': {
                try {
                    int value = std::stoi(optarg);
                    if (value <= 0) {
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
                    if (value <= 0) {
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
            default:
                std::cerr << "Unknown option: " << c << "\n";
                return EXIT_FAILURE;
        }
    }

    Chronometer chronometer;

    // chronometer.start("system_creation");
    ThreadManager thread_manager(harmonic_bond_threads, 1, electrostatic_bond_threads);

    System system;
    // system.addParticle(1.0, 1.0, 0.0, 0.0, 0.0);
    // system.addParticle(1.0, 1.0, 1.0, 1.0, 1.0);

    // system.addParticle(1.0, 2.0, 2.0, 2.0);
    // system.addParticle(1.0, 3.0, 3.0, 3.0);
    VerletIntegrator integrator(0.01, thread_manager, chronometer);

    auto electrostatic_force = std::make_shared<ElectrostaticForce>(1.0);
    integrator.addForce(electrostatic_force);
    // auto hbForce1 = std::make_shared<HarmonicBondForce>(1.0, 1.0);


    for (int i = 0; i < 10000; ++i) {
        // Placing particles along a 1D line (e.g., x-axis)
        // system.addParticle(1.0, 1.0, i * 3.0, 1.0, 0.0);  // (mass, x, y, z)
        int x = rand()%500;
        int y = rand()%500;
        int z = rand()%500;
        // printf("x: %d, y: %d, z: %d\n", x, y, z);
        system.addParticle(1.0, -1.0, x, y, z);  // (mass, x, y, z)
        // system.addParticle(rand()%100, rand()%100, rand()%100, 1.0, 0.0);  // (mass, x, y, z)
        // show coords:;
    }

    // integrator.addForce(hbForce1);


    // // auto ljForce = std::make_shared<LennardJonesForce>(0.1, 1.0);
    // hbForce1->addBond(0, 1);



    // auto hbForce2 = std::make_shared<HarmonicBondForce>(1.0, 1.0);
    // hbForce2->addBond(2, 3);

    // // integrator.addForce(ljForce);
    // integrator.addForce(hbForce2);

    Context context(system, integrator);



    // chronometer.end("system_creation");
    // chronometer.printTiming("system_creation", "us");
    // chronometer.printTimings();
    VERBOSE("Starting the simulation with %zu particles\n", system.getNumParticles());

    // Run the simulation
    // context.runSimulation(100000);
    // XYZWriter trajectoryWriter("trajectory.xyz");

    // Run the simulation and write trajectory
    const int numSteps = 1;
    const int outputInterval = 100; // Output every 100 steps

    for (int step = 0; step < numSteps; ++step) {
        // context.step();
         context.runSimulation(1);
        if (step % outputInterval == 0) {
            // trajectoryWriter.writeFrame(system); // Write the current frame
        }
    }

    // Close the trajectory writer
    // trajectoryWriter.close();
// 
    // Output the final positions of particles
    // for (size_t i = 0; i < system.getNumParticles(); ++i) {
    //     auto pos = system.getParticle(i).getPosition();
    //     std::cout << "Particle " << i << " final position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
    // }

    return 0;
}
