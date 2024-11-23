#ifndef ELECTROSTATIC_FORCE_H
#define ELECTROSTATIC_FORCE_H

#include "Force.h"
#include <utility>
#include <vector>

/**
 * Represents an electrostatic force between pairs of particles.
 */
class ElectrostaticForce : public Force {
public:
    ElectrostaticForce(double coulombConstant);

    /**
     * Adds a bond between two particles by their indices.
     * @param particle1Index - Index of the first particle.
     * @param particle2Index - Index of the second particle.
     */
    void addBond(size_t particle1Index, size_t particle2Index);

    /**
     * Compute the electrostatic forces and update the force array.
     * @param system - The system containing particles.
     * @param forces - A reference to the forces array to update.
     */
    void compute(System& system, std::vector<std::array<double, 3>>& forces, 
                    ThreadManager& thread_manager,
                    Chronometer& chronometer) const override;

    void compute_parallel(System& system, std::vector<std::array<double, 3>>& forces, 
                            ThreadManager& thread_manager,
                            Chronometer& chronometer) const override;

private:
    double coulombConstant;
    std::vector<std::pair<size_t, size_t>> bonds; // List of particle pairs for interactions
    std::vector<std::chrono::duration<long long>> thread_times;
};

#endif // ELECTROSTATIC_FORCE_H