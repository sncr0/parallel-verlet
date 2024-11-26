// integrator/VerletIntegrator.h

#ifndef VERLET_INTEGRATOR_H
#define VERLET_INTEGRATOR_H

#include "../system/System.h"
#include <array>
#include <vector>
#include <memory>
#include "../forces/Force.h"
#include "../thread_manager/ThreadManager.h"

class VerletIntegrator {
public:
    explicit VerletIntegrator(double timestep);//, ThreadManager& thread_manager, Chronometer& chronometer);
    void step(System& system);
    void addForce(std::shared_ptr<Force> force);


private:
    double timestep;
    
    // Parameters for Lennard-Jones potential
    const double epsilon = 1.0; // Depth of the potential well
    const double sigma = 1.0;   // Distance at which the potential is minimum
    // ThreadManager& thread_manager;
    // Chronometer& chronometer;

    std::vector<std::shared_ptr<Force>> forces; // Collection of forces
    std::vector<std::array<double, 3>> old_forces;
    bool first_step = true;


};

#endif
