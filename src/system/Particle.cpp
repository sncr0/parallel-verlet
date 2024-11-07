// system/Particle.cpp

#include "Particle.h"

Particle::Particle(double mass, double x, double y, double z)
    : mass(mass), position({x, y, z}) {}

double Particle::getMass() const {
    return mass;
}

const std::array<double, 3>& Particle::getPosition() const {
    return position;
}

void Particle::setPosition(double x, double y, double z) {
    position = {x, y, z};
}
