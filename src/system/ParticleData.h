// system/System.h

#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H

#include <vector>
#include <cstddef>
#include "Particle.h"


class ParticleData {
public:
    ParticleData();
    ~ParticleData();
    void insertData(double x, double y, double z, 
                        double vx, double vy, double vz, 
                        double charge, double mass);
    size_t getNumData() const;

    std::shared_ptr<Particle> getParticle(size_t index);
    const std::shared_ptr<const Particle> getParticle(size_t index) const;

private:
    std::vector<double> x_coordinates;
    std::vector<double> y_coordinates;
    std::vector<double> z_coordinates;

    std::vector<double> x_velocities;
    std::vector<double> y_velocities;
    std::vector<double> z_velocities;

    std::vector<double> charges;
    std::vector<double> masses;
};

#endif