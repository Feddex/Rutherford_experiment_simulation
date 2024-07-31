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
gp operator/(const gp &stru, double c);
gp operator+(const gp &lhs, const gp &rhs);
gp& operator+=(gp &lhs, const gp &rhs);


//classes definition 
class particle {
    protected:
    static constexpr double charge_unit = 1.6E-19;
    const double _charge;
    gp _position;

public:
    particle(double charge, gp position);
    gp get_position() const;
    double get_charge() const;
};

class goldenNucleus : public particle {
    static constexpr double charge_unit = 1.6E-19;
    
    public:
    goldenNucleus() : particle( 79., {0., 0.}) {} //we can put an integer as charge 
    

 };



class alpha : public particle {
    const double _mass;
    double _init_y;
    gp _velocity;
    gp _acceleration;

public:
    alpha(double mass, double charge, double init_y, gp position, gp velocity);
    
    // Methods to get vectors position, velocity, and acceleration
    double getmass() const;
    double getinity() const;
    gp get_velocity() const;
    gp get_acceleration() const;
    

    // Methods to set vectors position, velocity, and acceleration
    void set_position(const gp a);
    void set_velocity(const gp a);
    void set_acceleration(const gp acceleration);

    // // Method that returns the distance from the origin, therefore the distance from golden nucleus
    double O_distance() const;

    // // Method that returns the distance between two particles
    // double P_distance(const particle &pArticle) const;

    // // Method that returns the vector difference between two particles
    // gp vector_distance(const particle &pArticle) const;

    // // Returns the unit vector of distance
    // gp versor_distance(const particle &pArticle) const;

    // // Method that returns the acceleration unit vector
    // gp acceleration_versor(const particle &pArticle) const;

    // // Method to get the magnitude of the acceleration
    // double acceleration_module(const particle &pArticle);

    // gp acceleration_vector(const particle &pArticle);

    // // Method to evolve the particle using the Runge-Kutta method
    // void evolve(const particle& pArticle, double rilevator_radius, std::ofstream& angular_file);
};

// Define the arctangent to return values between 0 and 2*pi
double atangent(double y, double x);

#endif
