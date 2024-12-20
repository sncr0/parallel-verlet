// system/Particle.h

#ifndef PARTICLE_H
#define PARTICLE_H

#include <array>

class Particle {
public:
    Particle(double mass, double charge, double x, double y, double z,
             double vx, double vy, double vz); 

    double getMass() const;
    double getCharge() const;
    const std::array<double, 3>& getPosition() const;
    void setPosition(double x, double y, double z);
    const std::array<double, 3>& getVelocity() const;
    void setVelocity(double x, double y, double z);

private:
    double mass;
    double charge;
    std::array<double, 3> position;
    std::array<double, 3> velocity;
};

#endif
