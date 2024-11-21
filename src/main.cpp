// main.cpp

#include "system/System.h"
#include "context/Context.h"
#include "integrator/VerletIntegrator.h"
#include "forces/LennardJonesForce.h"
#include "forces/HarmonicBondForce.h"
#include <iostream>
#include "io/XYZWriter.h"
#include "getopt.h"
#include "logging/Verbose.h"
#include "thread_manager/ThreadManager.h"
#include <chrono>


int main(int argc, char **argv) {
    int c;
    int harmonic_bond_threads = 1;
    while ((c = getopt(argc,argv,"vh:")) != -1 ){
        switch(c) {
            case 'v':
                setVerboseFlag(1);
                break;
            case 'h':
                char* end;
                long value = std::strtol(optarg, &end, 10);
                if (*end != '\0' || value <= 0) {
                    std::cerr << "Invalid number of threads: " << optarg << "\n";
                    return EXIT_FAILURE;
                }
                harmonic_bond_threads = static_cast<int>(value);
                printf("harmonic_bond_threads: %d\n", harmonic_bond_threads);
                break;
        }
    }

    std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
    ThreadManager thread_manager(harmonic_bond_threads, 1);

    System system;
    // system.addParticle(1.0, 0.0, 0.0, 0.0);
    // system.addParticle(1.0, 1.0, 1.0, 1.0);

    // system.addParticle(1.0, 2.0, 2.0, 2.0);
    // system.addParticle(1.0, 3.0, 3.0, 3.0);
    VerletIntegrator integrator(0.01, thread_manager);
    auto hbForce1 = std::make_shared<HarmonicBondForce>(1.0, 1.0);


    for (int i = 0; i < 10000000; ++i) {
        // Placing particles along a 1D line (e.g., x-axis)
        system.addParticle(1.0, i * 1.1, 0.0, 0.0);  // (mass, x, y, z)
        system.addParticle(1.0, i * 1.0, 1.0, 0.0);  // (mass, x, y, z)

        hbForce1->addBond(2*i, 2*i + 1);

    }

    integrator.addForce(hbForce1);


    // // auto ljForce = std::make_shared<LennardJonesForce>(0.1, 1.0);
    // hbForce1->addBond(0, 1);

    // auto hbForce2 = std::make_shared<HarmonicBondForce>(1.0, 1.0);
    // hbForce2->addBond(2, 3);

    // // integrator.addForce(ljForce);
    // integrator.addForce(hbForce2);

    Context context(system, integrator);


    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);
    printf("Time to create system: %f\n", elapsed.count());

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

    // Output the final positions of particles
    // for (size_t i = 0; i < system.getNumParticles(); ++i) {
    //     auto pos = system.getParticle(i).getPosition();
    //     std::cout << "Particle " << i << " final position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
    // }

    return 0;
}
