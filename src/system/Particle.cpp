// system/Particle.cpp

#include "Particle.h"
#include <stdexcept>

Particle::Particle(
    std::shared_ptr<double> x, 
    std::shared_ptr<double> y, 
    std::shared_ptr<double> z,
    std::shared_ptr<double> vx, 
    std::shared_ptr<double> vy, 
    std::shared_ptr<double> vz,
    std::shared_ptr<double> c, 
    std::shared_ptr<double> m
) : 
    x_position(x), 
    y_position(y), 
    z_position(z),
    x_velocity(vx), 
    y_velocity(vy), 
    z_velocity(vz),
    charge(c), 
    mass(m)
{
    // Validate that all pointers are non-null
    if (!x || !y || !z || !vx || !vy || !vz || !c || !m) {
        throw std::invalid_argument("All pointer arguments must be non-null");
    }
}

const double& Particle::get_x_position() const { return *x_position; }
const double& Particle::get_y_position() const { return *x_position; }
const double& Particle::get_z_position() const { return *x_position; }

const double& Particle::get_x_velocity() const { return *x_velocity; }
const double& Particle::get_y_velocity() const { return *y_velocity; }
const double& Particle::get_z_velocity() const { return *z_velocity; }

const double& Particle::get_charge() const { return *charge; }
const double& Particle::get_mass() const { return *mass; }


void Particle::set_x_position(double value) { *x_position = value; }
void Particle::set_y_position(double value) { *y_position = value; }
void Particle::set_z_position(double value) { *z_position = value; }

void Particle::set_x_velocity(double value) { *x_velocity = value; }
void Particle::set_y_velocity(double value) { *y_velocity = value; }
void Particle::set_z_velocity(double value) { *z_velocity = value; }

void Particle::set_charge(double value) { *charge = value; }
void Particle::set_mass(double value) { *mass = value; }
