// system/System.h

#ifndef SYSTEM_H
#define SYSTEM_H

#include "Particle.h"
#include "ParticleData.h"
#include <vector>
#include <cstddef> 

class System {
public:
    System();
    ~System();
    
    void addParticle(double mass, double charge, 
                            double x, double y, double z,
                            double vx=0, double vy=0, double vz=0);
    size_t getNumParticles() const;
    std::shared_ptr<Particle> getParticle(size_t index);    
    const std::shared_ptr<const Particle> getParticle(size_t index) const;


private:
    ParticleData particle_data;
};

#endif