// main.cpp

#include "system/System.h"
#include "context/Context.h"
#include "integrator/VerletIntegrator.h"
#include "forces/LennardJonesForce.h"
#include "forces/HarmonicBondForce.h"
#include <iostream>
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
    system.addParticle(1.0, 0.0, 0.0, 0.0);
    system.addParticle(1.0, 1.0, 1.0, 1.0);

    VerletIntegrator integrator(0.01);


    auto ljForce = std::make_shared<LennardJonesForce>(0.1, 1.0);
    auto hbForce = std::make_shared<HarmonicBondForce>(1.0, 1.0);
    hbForce->addBond(0, 1);

    integrator.addForce(ljForce);
    integrator.addForce(hbForce);

    Context context(system, integrator);

    VERBOSE("Starting the simulation with %zu particles\n", system.getNumParticles());

    // Run the simulation
    context.runSimulation(100000);

    // Output the final positions of particles
    for (size_t i = 0; i < system.getNumParticles(); ++i) {
        auto pos = system.getParticle(i).getPosition();
        std::cout << "Particle " << i << " final position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
    }

    return 0;
}
