// system/System.cpp

#include "System.h"

System::System() {}

System::~System() {}

void System::addParticle(double mass, double x, double y, double z) {
    particles.emplace_back(mass, x, y, z, 0, 0, 0);
}

size_t System::getNumParticles() const {
    return particles.size();
}

Particle& System::getParticle(size_t index) {
    return particles.at(index);
}

const Particle& System::getParticle(size_t index) const {
    return particles.at(index); // Const version
}
