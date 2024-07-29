#ifndef RAT_HPP
#define RAT_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <fstream>

// this struct will be used to store information about the components of position, velocity and acceleration 
struct gp {
    double a;
    double b;
};

void print(gp const& p);

gp operator*(double c, const gp &stru);
gp operator+(const gp &lhs, const gp &rhs);
gp& operator+=(gp &lhs, const gp &rhs);

class particle {
    static constexpr double charge_unit = 1.6E-19; // absolute value of an electron's charge
    double _mass;
    double _charge;
    double _init_y;
    gp _position;
    gp _velocity;
    gp _acceleration;

public:
    particle(double mass, double charge, double init_y, gp position, gp velocity);
    particle() : _mass(3.26628E-25), _charge(79 * charge_unit), _init_y(0), _position{0., 0.}, _velocity{0., 0.}, _acceleration{0., 0.} {}

    // Methods to get vectors position, velocity, and acceleration
    gp getposition() const;
    gp getvelocity() const;
    gp getacceleration() const;
    double getcharge() const;
    double getmass() const;
    double getinity() const;

    // Methods to set vectors position, velocity, and acceleration
    void set_xposition(double a);
    void set_yposition(double b);
    void set_xvelocity(double a);
    void set_yvelocity(double b);

    // Method that returns the distance from the origin
    double O_distance() const;

    // Method that returns the distance between two particles
    double P_distance(const particle &pArticle) const;

    // Method that returns the vector difference between two particles
    gp vector_distance(const particle &pArticle) const;

    // Returns the unit vector of distance
    gp versor_distance(const particle &pArticle) const;

    // Method that returns the acceleration unit vector
    gp acceleration_versor(const particle &pArticle) const;

    // Method to get the magnitude of the acceleration
    double acceleration_module(const particle &pArticle);

    gp acceleration_vector(const particle &pArticle);

    // Method to evolve the particle using the Runge-Kutta method
    void evolve(const particle& pArticle, double rilevator_radius, std::ofstream& angular_file);
};

// Define the arctangent to return values between 0 and 2*pi
double atangent(double y, double x);

#endif
