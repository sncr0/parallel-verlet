// system/Particle.h

#ifndef PARTICLE_H
#define PARTICLE_H

#include <array>

class Particle {
public:
    Particle(double mass, double x, double y, double z);

    double getMass() const;
    const std::array<double, 3>& getPosition() const;
    void setPosition(double x, double y, double z);

private:
    double mass;
    std::array<double, 3> position;
};

#endif
