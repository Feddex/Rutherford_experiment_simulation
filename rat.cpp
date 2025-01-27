#include "rat.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <fstream>

void print(gp const &p)
{
    std::cout << "(" << p.a << ", " << p.b << ")" << '\n';
}

gp operator*(double c, const gp &stru) { return {c * stru.a, c * stru.b}; }
gp operator/(const gp &stru, double c) { return {stru.a / c, stru.b / c}; }
gp operator+(const gp &lhs, const gp &rhs) { return {lhs.a + rhs.a, lhs.b + rhs.b}; }
gp &operator+=(gp &lhs, const gp &rhs)
{
    lhs.a += rhs.a;
    lhs.b += rhs.b;
    return lhs;
}
//pourposely made for this specific configuration
double atangent(double y, double x)
{
    double pi = 3.14159265359;
    if (x == 0 && y > 0)
        return pi / 2;
    if (x == 0 && y < 0)
        return 3 * pi / 2;
    if (x > 0 && y == 0)
        return 0; // ok
    if (x < 0 && y == 0)
        return pi; // ok
    if (x > 0 && y > 0)
        return atan(y / x);
    if (x > 0 && y < 0)
        return  2*pi + atan(y / x); // there was 2 * pi + atan(y / x);
    if (x < 0 && y > 0)
        return pi + atan(y / x);
    if (x < 0 && y < 0)
        return pi + atan(y / x);
    return 0;
}

particle::particle(double charge, gp position): _charge(charge*charge_unit), _position(position) {}

alpha::alpha(double mass, double charge, double init_y, gp position, gp velocity)
    : _mass(mass), particle::particle(charge, position), _init_y(init_y), _velocity(velocity)
{
    if (mass < 0)
        throw std::invalid_argument("Mass must be positive");
}


double particle::get_charge() const { return _charge; }
gp particle::get_position() const { return _position; }


gp alpha::get_velocity() const { return _velocity; }
gp alpha::get_acceleration() const { return _acceleration; }
double alpha::getmass() const { return _mass; }
double alpha::getinity() const { return _init_y; }

void alpha::set_position(const gp a) { _position = a; }
void alpha::set_velocity(const gp a) { _velocity = a; }
void alpha::set_acceleration(const gp acceleration) { _acceleration = acceleration; }

double alpha::O_distance() const
{
    return std::sqrt(std::pow(_position.a, 2) + std::pow(_position.b, 2));
}


// // Yoshida coefficients
// const double w1 = 1.0 / (2.0 - std::pow(2.0, 1.0 / 3.0));
// const double w0 =  1.0 - 2.0 * w1;

// const double c1 = w1 / 2.0;
// const double c2 = (w0 + w1) / 2.0;
// const double c3 = c2;
// const double c4 = c1;

// const double d1 = w1;
// const double d2 = w0;
// const double d3 = w1;

// void alpha::evolve(const particle &pArticle, double rilevator_radius, std::ofstream &angular_file)
// {
//     double t = 0;
//     double dt = 1e-22;
    

//     std::vector<gp> transit;
//     transit.push_back(_position);

//     // //FOR LP
//     // // Initial half-step for velocity
//     // gp initial_acc = acceleration_vector(pArticle);
//     // _velocity += 0.5 * dt * initial_acc;

//     while (O_distance() < rilevator_radius)
//     {
//         // //YOSHIDA
//         // Yoshida Step 1
//         _position += c1 * dt * _velocity;
//         _velocity += d1 * dt * acceleration_vector(pArticle);
       
//         // Yoshida Step 2
//         _position += c2 * dt * _velocity;
//         _velocity += d2 * dt * acceleration_vector(pArticle);
        
//         // Yoshida Step 3
//         _position += c3 * dt * _velocity;
//          _velocity += d3 * dt * acceleration_vector(pArticle);

//         // Yoshida Step 4
//         _position += c4 * dt * _velocity;

//         // // Record the position for analysis
//         // transit.push_back(_position);

        
        //RK
        // gp k1_v = acceleration_vector(pArticle);
        // gp k1_p = _velocity;

        // gp k2_v = acceleration_vector(pArticle) + (dt / 2) * k1_v;
        // gp k2_p = _velocity + (dt / 2) * k1_v;

        // gp k3_v = acceleration_vector(pArticle) + (dt / 2) * k2_v;
        // gp k3_p = _velocity + (dt / 2) * k2_v;

        // gp k4_v = acceleration_vector(pArticle) + dt * k3_v;
        // gp k4_p = _velocity + dt * k3_v;

        // _position += (dt / 6) * (k1_p + 2 * k2_p + 2 * k3_p + k4_p);
        // _velocity += (dt / 6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);


        //LP
        //  // Update position
        // _position += dt * _velocity;

        // // Compute new acceleration
        // gp new_acc = acceleration_vector(pArticle);

        // // Update velocity with the new acceleration
        // _velocity += 0.5 * dt * new_acc;

                

    //     // Record the position for analysis
    //     transit.push_back(_position);
        
    //     // std::cout << "Position: ";
    //     //     print(_position);
    //         // std::cout << "Velocity: ";
    //         // print(_velocity);
    //         // std::cout << "acc: ";
    //         // print(acceleration_vector(pArticle));
    // }




///ANGOLO NON TOCCARE 
//     // Determine the angle of the last position relative to the rilevator_radius
//     auto it = (transit.end() - 1);
//     auto is = it - 1;
//     double diff_lastpos_radius = std::abs(O_distance() - rilevator_radius);
//     double diff_penlastpos_radius = std::abs(std::sqrt(std::pow((*is).a, 2) + std::pow((*is).b, 2)) - rilevator_radius);

//     if (diff_lastpos_radius < diff_penlastpos_radius)
//     {
//         double angular = atangent((*it).b-getinity(), (*it).a);
//         angular_file << angular << '\n';
//     }
//     else
//     {
//         double angular = atangent((*is).b-getinity(), (*is).a);
//         angular_file << angular << '\n';
//     }
// }