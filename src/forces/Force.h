#ifndef FORCE_H
#define FORCE_H

#include <vector>
#include <array>
#include "../system/System.h"
#include "../thread_manager/ThreadManager.h"

class Force {
public:
    virtual ~Force() = default;

    /**
     * Compute forces and update the force array.
     * @param system - The system containing particles.
     * @param forces - A reference to the forces array to update.
     */
    virtual void compute(System& system,
                            std::vector<std::array<double, 3>>& forces,
                            ThreadManager& thread_manager) const = 0;
};

#endif // FORCE_H
