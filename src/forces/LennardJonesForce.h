#ifndef LENNARD_JONES_FORCE_H
#define LENNARD_JONES_FORCE_H

#include "Force.h"
#include <cmath>

class LennardJonesForce : public Force {
public:
    LennardJonesForce(double epsilon, double sigma);

    void compute(System& system, std::vector<std::array<double, 3>>& forces,
                    ThreadManager& thread_manager,
                    Chronometer& chronometer) const override;

private:
    double epsilon; // Depth of the potential well
    double sigma;   // Finite distance where inter-particle potential is zero
};

#endif // LENNARD_JONES_FORCE_H
