// system/ParticleData.cpp

#include "ParticleData.h"

ParticleData::ParticleData() {}

ParticleData::~ParticleData() {}

void ParticleData::insertData(double x, double y, double z, 
                                double vx, double vy, double vz, 
                                double charge, double mass) {

x_coordinates.emplace_back(x);
y_coordinates.emplace_back(x);
z_coordinates.emplace_back(x);

x_velocities.emplace_back(x);
y_velocities.emplace_back(x);
z_velocities.emplace_back(x);

charges.emplace_back(x);
masses.emplace_back(x);
};

std::shared_ptr<Particle> ParticleData::getParticle(size_t index) {
    std::shared_ptr<double> x = std::make_shared<double>(x_coordinates[index]);
    std::shared_ptr<double> y = std::make_shared<double>(x_coordinates[index]);
    std::shared_ptr<double> z = std::make_shared<double>(x_coordinates[index]);

    std::shared_ptr<double> vx = std::make_shared<double>(x_velocities[index]);
    std::shared_ptr<double> vy = std::make_shared<double>(y_velocities[index]);
    std::shared_ptr<double> vz = std::make_shared<double>(z_velocities[index]);

    std::shared_ptr<double> charge = std::make_shared<double>(charges[index]);
    std::shared_ptr<double> mass = std::make_shared<double>(masses[index]);

    return std::make_shared<Particle>(x, y, z, vx, vy, vz, charge, mass);
};

const std::shared_ptr<const Particle> ParticleData::getParticle(size_t index) const {
    std::shared_ptr<double> x = std::make_shared<double>(x_coordinates[index]);
    std::shared_ptr<double> y = std::make_shared<double>(x_coordinates[index]);
    std::shared_ptr<double> z = std::make_shared<double>(x_coordinates[index]);

    std::shared_ptr<double> vx = std::make_shared<double>(x_velocities[index]);
    std::shared_ptr<double> vy = std::make_shared<double>(y_velocities[index]);
    std::shared_ptr<double> vz = std::make_shared<double>(z_velocities[index]);

    std::shared_ptr<double> charge = std::make_shared<double>(charges[index]);
    std::shared_ptr<double> mass = std::make_shared<double>(masses[index]);

    return std::make_shared<const Particle>(x, y, z, vx, vy, vz, charge, mass);
};

size_t ParticleData::getNumData() const { return x_coordinates.size(); }