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



int main(int argc, char **argv) {
    int c;
    while ((c = getopt(argc,argv,"v")) != -1 ){
        switch(c) {
            case 'v':
                setVerboseFlag(1);
                break;
        }
    }


    System system;
    // system.addParticle(1.0, 0.0, 0.0, 0.0);
    // system.addParticle(1.0, 1.0, 1.0, 1.0);

    // system.addParticle(1.0, 2.0, 2.0, 2.0);
    // system.addParticle(1.0, 3.0, 3.0, 3.0);
    VerletIntegrator integrator(0.01);
    auto hbForce1 = std::make_shared<HarmonicBondForce>(1.0, 1.0);


    for (int i = 0; i < 1000000; ++i) {
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

    VERBOSE("Starting the simulation with %zu particles\n", system.getNumParticles());

    // Run the simulation
    // context.runSimulation(100000);
    XYZWriter trajectoryWriter("trajectory.xyz");

    // Run the simulation and write trajectory
    const int numSteps = 1;
    const int outputInterval = 100; // Output every 100 steps

    for (int step = 0; step < numSteps; ++step) {
        // context.step();
         context.runSimulation(1);
        if (step % outputInterval == 0) {
            trajectoryWriter.writeFrame(system); // Write the current frame
        }
    }

    // Close the trajectory writer
    trajectoryWriter.close();

    // Output the final positions of particles
    // for (size_t i = 0; i < system.getNumParticles(); ++i) {
    //     auto pos = system.getParticle(i).getPosition();
    //     std::cout << "Particle " << i << " final position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
    // }

    return 0;
}
