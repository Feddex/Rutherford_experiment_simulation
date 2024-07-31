#include "rat.hpp"
#include <iostream>
#include <random>
#include <fstream>

int main() {
    double charge_unit = 1.6E-19;
    double alpha_mass = 6.64424E-27;
    double rilevator_radius = 2.11E-7;
    
     //this is the gold atom nucleus 
    const particle gold_nucleus;
    //std::cout<<"mass: "<<gold_nucleus.getmass()<<" charge: "<<gold_nucleus.getcharge()<<'\n';
    std::vector<double> angular_vector;

    //let's define the distribution that will generate the initial position of the alpha particle 
    std::random_device r;
    std::default_random_engine eng(r());
    std::normal_distribution<> d(0, 3*1.44E-10);

    std::ofstream angular_file("angulars.txt");
    
    //this is the numer of alpha particle fired
    int num_shot = 10000;

    for (int i = 0; i < num_shot; ++i) {
       double init_y = d(eng);
        particle alpha{alpha_mass, 2 * charge_unit, init_y ,{-4.5e-14, init_y}, {1778866, 0}}; //ther was 1.53E7
        alpha.evolve(gold_nucleus, rilevator_radius, angular_file);
    }

    // particle alpha{alpha_mass, 2 * charge_unit, {-2.108E-7, 0}, {282006, 0}}; //ther was 1.53E7
    // alpha.evolve(gold_nucleus, rilevator_radius, angular_file);

    return 0;
}
