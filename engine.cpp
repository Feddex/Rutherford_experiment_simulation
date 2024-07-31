
#include "engine.hpp"
#include "rat.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <fstream>

// BASE ENGINE CLASS FUNCTION
double engine::particle_distance(alpha const &alpha, goldenNucleus const &gold)
{
    double delta_x = alpha.get_position().a - gold.get_position().a;
    double delta_y = alpha.get_position().b - gold.get_position().b;
    return std::sqrt(std::pow(delta_x, 2) + std::pow(delta_y, 2));
}

gp engine::vector_distance(alpha const &alpha, goldenNucleus const &gold)
{
    double delta_x = alpha.get_position().a - gold.get_position().a;
    double delta_y = alpha.get_position().b - gold.get_position().b;
    return {delta_x, delta_y};
};

gp engine::versor_distance(alpha const &alpha, goldenNucleus const &gold)
{
    gp vector_distance = engine::vector_distance(alpha, gold);
    double distance = engine::particle_distance(alpha, gold);
    return vector_distance / distance;
};

// this method sets the acceleration of the alpha particle
void engine::applyForce(alpha& alpha, goldenNucleus const &gold)
{
    gp versor_distance = engine::versor_distance(alpha, gold);
    double distance = particle_distance(alpha, gold);
    double acceleration_module = (_k * alpha.get_charge() * gold.get_charge()) / (std::pow(distance, 2) * alpha.getmass()); ///CHECK EXPRESSION

    gp acceleration_vector = acceleration_module * versor_distance;
    alpha.set_acceleration(acceleration_vector);
};

// this function return the corrected angle of deflection.
void engine::angle_calc(double const &init_y, double const &rilevator_radius, std::vector<gp> &transit, std::ofstream &angle_file)
{
    auto it = (transit.end() - 1);
    auto is = it - 1;

    // this values are required to determine which position between the last one or the previous one in trasnit vector was the closest to the rilevator screen.
    // the closes is used to determine the angle of deflection  //POSSBILE BUG
    double diff_lastpos_radius = std::abs(std::sqrt(std::pow((*it).a, 2) + std::pow((*it).b, 2)) - rilevator_radius);    // c'era O_distance() al posto di (*it)
    double diff_penlastpos_radius = std::abs(std::sqrt(std::pow((*is).a, 2) + std::pow((*is).b, 2)) - rilevator_radius); // questi misurano la distanza dal centro, se il protone non è al centro non va bene

    if (diff_lastpos_radius < diff_penlastpos_radius)
    {
        double angular = atangent((*it).b - init_y, (*it).a);
        angle_file << angular << '\n';
    }
    else
    {
        double angular = atangent((*is).b - init_y, (*is).a);
        angle_file << angular << '\n';
    }
};

void engine::evolve(alpha& alpha, goldenNucleus const &gold, double rilevator_radius, std::ofstream &angle_file, double dt)
{
    double init_y = alpha.getinity();
    // FILL THE INTEGRATOR
}


// DERIVED INTEGRATOR CLASSES:

// 1) YOSHIDA ENGINE
//  this function is the core of the engine class
void yoshidaEngine::evolve(alpha& alpha, goldenNucleus const &gold, double rilevator_radius, std::ofstream &angle_file, double dt)
{
    double init_y = alpha.getinity();
    std::vector<gp> transit = Yoshida(alpha, gold, rilevator_radius, dt);
    angle_calc(init_y, rilevator_radius, transit, angle_file);
}

// YOSHIDA INTEGRATION BLOCK
std::vector<gp> yoshidaEngine::Yoshida(alpha& alpha, goldenNucleus const &gold, double rilevator_radius, double dt)
{
    // definition of transit vector: this will be returned
    std::vector<gp> transit;
    transit.push_back(alpha.get_position());
    //
    std::cout << (alpha.O_distance() < rilevator_radius) <<'\n';
    // while the distance of the alpha particle from the origin (center of the cirular rilevator screen) is less then the radius of the rilevator, the integration continues.
    while (alpha.O_distance() < rilevator_radius)
    {
        // apply the yoshida step function
        yoshida_steps_package(alpha, gold, dt);
        // std::cout << "position: " ;
        // print(alpha.get_position());
        // Record the position for analysis
        transit.push_back(alpha.get_position());
    }

    return transit;
};
void yoshidaEngine::yoshida_steps_package(alpha& alpha, const goldenNucleus gold, const double dt)
{
    // Yoshida Step 1
    yoshida_step1(alpha, gold, dt);
    // Yoshida Step 2
    yoshida_step2(alpha, gold, dt);
    // Yoshida Step 3
    yoshida_step3(alpha, gold, dt);
    // Yoshida Step 4
    yoshida_step4(alpha, gold, dt);
}
void yoshidaEngine::yoshida_step1(alpha& alpha, const goldenNucleus gold, const double dt)
{
    gp newAlphaPos1 = alpha.get_position() + c1 * dt * alpha.get_velocity();
    alpha.set_position(newAlphaPos1);
    applyForce(alpha, gold); // RECALCULATING ALPHA'S ACCELERATION
    gp newAlphaVel1 = alpha.get_velocity() + d1 * dt * alpha.get_acceleration();
    alpha.set_velocity(newAlphaVel1);
};
void yoshidaEngine::yoshida_step2(alpha& alpha, const goldenNucleus gold, const double dt)
{
    gp newAlphaPos2 = alpha.get_position() + c2 * dt * alpha.get_velocity();
    alpha.set_position(newAlphaPos2);
    applyForce(alpha, gold); // RECALCULATING ALPHA'S ACCELERATION
    gp newAlphaVel2 = alpha.get_velocity() + d2 * dt * alpha.get_acceleration();
    alpha.set_velocity(newAlphaVel2);
};
void yoshidaEngine::yoshida_step3(alpha& alpha, const goldenNucleus gold, const double dt)
{
    gp newAlphaPos3 = alpha.get_position() + c3 * dt * alpha.get_velocity();
    alpha.set_position(newAlphaPos3);
    applyForce(alpha, gold); // RECALCULATING ALPHA'S ACCELERATION
    gp newAlphaVel3 = alpha.get_velocity() + d3 * dt * alpha.get_acceleration();
    alpha.set_velocity(newAlphaVel3);
};
void yoshidaEngine::yoshida_step4(alpha& alpha, const goldenNucleus gold, const double dt)
{
    gp newAlphaPos4 = alpha.get_position() + c4 * dt * alpha.get_velocity();
    alpha.set_position(newAlphaPos4);
};

// 2) LEAPFROG ENGINE
void leapFrogEngine::evolve(alpha& alpha, goldenNucleus const &gold, double rilevator_radius, std::ofstream &angle_file, double dt)
{
    double init_y = alpha.getinity();
    std::vector<gp> transit = Leapfrog(alpha, gold, rilevator_radius, dt);
    angle_calc(init_y, rilevator_radius, transit, angle_file);
}

std::vector<gp> leapFrogEngine::Leapfrog(alpha& alpha, goldenNucleus const &gold, double dt, double rilevator_radius)
{
    // definition of transit vector: this will be returned
    std::vector<gp> transit;
    transit.push_back(alpha.get_position());

    // fisrt half step
    leapfrog_step1(alpha, gold, dt);

    // while the distance of the alpha particle from the origin (center of the cirular rilevator screen) is less then the radius of the rilevator, the integration continues.
    while (alpha.O_distance() < rilevator_radius)
    {

        leapfrog_steps_package(alpha, gold, dt);

        // Record the position for analysis
        transit.push_back(alpha.get_position());
    }

    return transit;
}
// LEAPFROG INTEGRATION BLOCK
void leapFrogEngine::leapfrog_step1(alpha& alpha, const goldenNucleus gold, const double dt)
{
    applyForce(alpha, gold); // RECALCULATING ACCELERATION
    gp newAlphaVel1 = alpha.get_velocity() + 0.5 * dt * alpha.get_acceleration();
    alpha.set_velocity(newAlphaVel1);
}

void leapFrogEngine::leapfrog_step2(alpha& alpha, const goldenNucleus gold, const double dt)
{
    gp newAlphaPos1 = alpha.get_position() + dt * alpha.get_velocity();
    alpha.set_position(newAlphaPos1);
}

void leapFrogEngine::leapfrog_steps_package(alpha& alpha, const goldenNucleus gold, const double dt)
{
    leapfrog_step2(alpha, gold, dt); //step 2 before the 1, beacuse the sequense is 1 before cycle and then 2, 1, 2,1 ,2 ,1 ...
    leapfrog_step1(alpha, gold, dt);
}


// 3) RK ENGINE
void RKEngine::evolve(alpha& alpha, goldenNucleus const &gold, double rilevator_radius, std::ofstream &angle_file, double dt)
{
    double init_y = alpha.getinity();
    std::vector<gp> transit = RK(alpha, gold, rilevator_radius, dt);
    angle_calc(init_y, rilevator_radius, transit, angle_file);
}

std::vector<gp> RKEngine::RK(alpha& alpha, goldenNucleus const &gold, double dt, double rilevator_radius)
{
    // definition of transit vector: this will be returned
    std::vector<gp> transit;
    transit.push_back(alpha.get_position());

    // while the distance of the alpha particle from the origin (center of the cirular rilevator screen) is less then the radius of the rilevator, the integration continues.
    while (alpha.O_distance() < rilevator_radius)
    {

        //RK STEPS
        RK_steps_package(alpha, gold, dt);

        // Record the position for analysis
        transit.push_back(alpha.get_position());
    }

    return transit;
}
//  RK INTEGRATION BLOCK
void RKEngine::RK_steps_package (alpha& alpha, const goldenNucleus gold, const double dt){
    // RK coefficients 
        applyForce(alpha, gold);
        gp k1_v = alpha.get_acceleration();
        gp k1_p = alpha.get_velocity();

       gp k2_v = alpha.get_acceleration() + (dt / 2) * k1_v;
       gp k2_p = alpha.get_velocity() + (dt / 2) * k1_v;
       
       gp k3_v = alpha.get_acceleration() + (dt / 2) * k2_v;
       gp k3_p = alpha.get_velocity() + (dt / 2) * k2_v;

       gp k4_v = alpha.get_acceleration() + dt * k3_v;
       gp k4_p = alpha.get_velocity() + dt * k3_v;

       gp newAlphaPos = alpha.get_position() + (dt / 6) * (k1_p + 2 * k2_p + 2 * k3_p + k4_p);
       gp newAlphaVel = alpha.get_velocity() + (dt / 6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
       alpha.set_position(newAlphaPos);
       alpha.set_velocity(newAlphaVel);
};
