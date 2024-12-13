// system/Particle.h

#ifndef PARTICLE_H
#define PARTICLE_H

#include <array>
#include <memory>

class Particle {
public:
    Particle(
        std::shared_ptr<double> x,
        std::shared_ptr<double> y,
        std::shared_ptr<double> z,
        std::shared_ptr<double> vx,
        std::shared_ptr<double> vy,
        std::shared_ptr<double> vz,
        std::shared_ptr<double> c,
        std::shared_ptr<double> m
    );

    const double& get_x_position() const;
    const double& get_y_position() const;
    const double& get_z_position() const;

    const double& get_x_velocity() const;
    const double& get_y_velocity() const;
    const double& get_z_velocity() const;

    const double& get_charge() const;
    const double& get_mass() const;


    void set_x_position(double value);
    void set_y_position(double value);
    void set_z_position(double value);

    void set_x_velocity(double value);
    void set_y_velocity(double value);
    void set_z_velocity(double value);

    void set_charge(double value);
    void set_mass(double value);


private:
    std::shared_ptr<double> x_position;
    std::shared_ptr<double> y_position;
    std::shared_ptr<double> z_position;

    std::shared_ptr<double> x_velocity;
    std::shared_ptr<double> y_velocity;
    std::shared_ptr<double> z_velocity;

    std::shared_ptr<double> charge;
    std::shared_ptr<double> mass;
};

#endif
