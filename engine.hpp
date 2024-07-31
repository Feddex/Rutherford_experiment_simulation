#ifndef ENGINE_HPP
#define ENGINE_HPP

#include "rat.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <fstream>

class engine
{
    //protected:
public:
    const double _k = 8.9875517873681764e9;
    double particle_distance(alpha const &alpha, goldenNucleus const &gold);
    gp vector_distance(alpha const &alpha, goldenNucleus const &gold);
    gp versor_distance(alpha const &alpha, goldenNucleus const &gold);

    // this method sets the acceleration of the alpha particle
    void applyForce(alpha& alpha, goldenNucleus const &gold);
    void angle_calc(double const &init_y, double const &rilevator_radius, std::vector<gp> &transit, std::ofstream &angle_file);
    virtual void evolve(alpha& alpha, goldenNucleus const &gold, double rilevator_radius, std::ofstream &angle_file, double dt);

    // // each method integrates the trajecotry and returns the record of positions of the single particle
    // std::vector<gp> leapfrog(alpha& alpha, goldenNucleus const &gold, double dt, double rilevator_radius);
    // std::vector<gp> RK(alpha& alpha, goldenNucleus const &gold, double dt, double rilevator_radius);
   
    
};

class yoshidaEngine : public engine {
    private: 
    // Yoshida coefficients
    const double w1 = 1.0 / (2.0 - std::pow(2.0, 1.0 / 3.0));
    const double w0 = 1.0 - 2.0 * w1;

    const double c1 = w1 / 2.0;
    const double c2 = (w0 + w1) / 2.0;
    const double c3 = c2;
    const double c4 = c1;

    const double d1 = w1;
    const double d2 = w0;
    const double d3 = w1;

    public:
    void evolve(alpha& alpha, goldenNucleus const &gold, double rilevator_radius, std::ofstream &angle_file, double dt) override;
    //YOSHIDA INTEGRATION BLOCK
    std::vector<gp> Yoshida(alpha& alpha, goldenNucleus const &gold, double dt, double rilevator_radius);
    void yoshida_steps_package (alpha& alpha, const goldenNucleus gold, const double dt);
    void yoshida_step1(alpha& alpha, const goldenNucleus gold, const double dt);
    void yoshida_step2(alpha& alpha, const goldenNucleus gold, const double dt);
    void yoshida_step3(alpha& alpha, const goldenNucleus gold, const double dt);
    void yoshida_step4(alpha& alpha, const goldenNucleus gold, const double dt);
};

class leapFrogEngine : public engine {

    public:
    void evolve(alpha& alpha, goldenNucleus const &gold, double rilevator_radius, std::ofstream &angle_file, double dt) override;
    //LEAPFROG INTEGRATION BLOCK
    std::vector<gp> Leapfrog(alpha& alpha, goldenNucleus const &gold, double dt, double rilevator_radius);
    void leapfrog_steps_package (alpha& alpha, const goldenNucleus gold, const double dt);
    void leapfrog_step1(alpha& alpha, const goldenNucleus gold, const double dt);
    void leapfrog_step2(alpha& alpha, const goldenNucleus gold, const double dt);
};

class RKEngine : public engine {

    public:
    void evolve(alpha& alpha, goldenNucleus const &gold, double rilevator_radius, std::ofstream &angle_file, double dt) override;
    //RK INTEGRATION BLOCK
    std::vector<gp> RK(alpha& alpha, goldenNucleus const &gold, double dt, double rilevator_radius);
    void RK_steps_package (alpha& alpha, const goldenNucleus gold, const double dt);
    
};



    
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

    // RK
    //  gp k1_v = acceleration_vector(pArticle);
    //  gp k1_p = _velocity;

    // gp k2_v = acceleration_vector(pArticle) + (dt / 2) * k1_v;
    // gp k2_p = _velocity + (dt / 2) * k1_v;

    // gp k3_v = acceleration_vector(pArticle) + (dt / 2) * k2_v;
    // gp k3_p = _velocity + (dt / 2) * k2_v;

    // gp k4_v = acceleration_vector(pArticle) + dt * k3_v;
    // gp k4_p = _velocity + dt * k3_v;

    // _position += (dt / 6) * (k1_p + 2 * k2_p + 2 * k3_p + k4_p);
    // _velocity += (dt / 6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);

    // LP
    //   // Update position
    //  _position += dt * _velocity;

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




#endif
