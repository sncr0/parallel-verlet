#include <iostream>
#include <openmm/System.h>
#include <openmm/Context.h>
#include <openmm/Integrator.h>
#include <openmm/VerletIntegrator.h>
#include <openmm/NonbondedForce.h>
#include <openmm/Platform.h>
#include <cstdlib>
#include <ctime>

using namespace OpenMM;

int main(int argc, char** argv) {
    // Set up simulation parameters
    const double timeStep = 0.01; // Time step in ps
    const int numParticles = 1000; // Number of particles
    const double boxSize = 100.0; // Box dimensions in nm
    const int numSteps = 1000; // Number of steps

    // Seed RNG
    std::srand(std::time(0));

    // Create the system
    System system;

    // Add particles with random positions inside the box
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0); // Mass = 1 amu
    }

    // Create the NonbondedForce for Electrostatics and van der Waals interactions
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::PME); // Use PME for long-range electrostatics
    nonbonded->setCutoffDistance(1.0); // Cutoff for short-range interactions in nm

    // Add particles to the NonbondedForce
    for (int i = 0; i < numParticles; i++) {
        double charge = ((i % 2) == 0) ? 1.0 : -1.0; // Alternate charges
        double sigma = 0.3; // Lennard-Jones sigma in nm
        double epsilon = 0.1; // Lennard-Jones epsilon in kJ/mol
        nonbonded->addParticle(charge, sigma, epsilon);

        // Randomize positions within the box
        double x = (std::rand() / (double)RAND_MAX) * boxSize;
        double y = (std::rand() / (double)RAND_MAX) * boxSize;
        double z = (std::rand() / (double)RAND_MAX) * boxSize;
        system.setDefaultPeriodicBoxVectors({boxSize, 0, 0}, {0, boxSize, 0}, {0, 0, boxSize});
    }

    // Add the NonbondedForce to the system
    system.addForce(nonbonded);

    // Create a VerletIntegrator
    VerletIntegrator integrator(timeStep);

    // Create a context
    Context context(system, integrator, Platform::getPlatformByName("CPU"));

    // Initialize particle positions randomly
    std::vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++) {
        positions[i] = Vec3((std::rand() / (double)RAND_MAX) * boxSize,
                            (std::rand() / (double)RAND_MAX) * boxSize,
                            (std::rand() / (double)RAND_MAX) * boxSize);
    }
    context.setPositions(positions);

    // Run the simulation
    for (int step = 0; step < numSteps; step++) {
        integrator.step(1);
        if (step % 100 == 0) {
            std::cout << "Step: " << step << " completed." << std::endl;
        }
    }

    // Clean up
    delete nonbonded;

    return 0;
}
