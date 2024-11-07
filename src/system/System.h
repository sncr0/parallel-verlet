// system/System.h

#ifndef SYSTEM_H
#define SYSTEM_H

#include "Particle.h"
#include <vector>
#include <cstddef> 

class System {
public:
    System();
    ~System();
    
    void addParticle(double mass, double x, double y, double z);
    size_t getNumParticles() const;
    Particle& getParticle(size_t index);

private:
    std::vector<Particle> particles;
};

#endif