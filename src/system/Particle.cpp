// system/Particle.cpp

#include "Particle.h"

Particle::Particle(double mass, 
                    double charge,
                    double x, double y, double z, 
                    double vx, double vy, double vz) : mass(mass), charge(charge), position({x, y, z}), velocity({vx, vy, vz}) {}

double Particle::getMass() const {
    return mass;
}

double Particle::getCharge() const {
    return charge;
}

const std::array<double, 3>& Particle::getPosition() const {
    return position;
}

void Particle::setPosition(double x, double y, double z) {
    position = {x, y, z};
}

const std::array<double, 3>& Particle::getVelocity() const {
    return velocity;
}

void Particle::setVelocity(double x, double y, double z) {
    velocity = {x, y, z};
}
