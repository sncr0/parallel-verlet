// system/System.cpp

#include "System.h"

System::System() {}

System::~System() {}

void System::addParticle(double mass, double charge, 
                            double x, double y, double z,
                            double vx, double vy, double vz) {

    particle_data.insertData(mass, charge, 
                                x, y, z,
                                vx, vy, vz);

}

size_t System::getNumParticles() const {
    return particle_data.getNumData();
}

std::shared_ptr<Particle> System::getParticle(size_t index) {
    return particle_data.getParticle(index);
}

const std::shared_ptr<const Particle> System::getParticle(size_t index) const {
    return particle_data.getParticle(index); // Const version
}
