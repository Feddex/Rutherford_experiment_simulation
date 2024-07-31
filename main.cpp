#include "rat.hpp"
#include "engine.hpp"
#include <iostream>
#include <random>
#include <fstream>
#include <cmath>

int main()
{
    double rilevator_radius = 1.46E-10;
    goldenNucleus gold;
    std::vector<double> angle_vector;

    std::random_device r;
    std::default_random_engine eng(r());
    double a = -1.44E-10;
    double b = 1.44E-10;

    // CHOSE THE INITIAL POSITION DISTRIBUTION
     std::uniform_real_distribution<> d(a, b);
    //std::normal_distribution<> d(0, 1.44E-10);

    std::ofstream angle_file("angulars.txt");
    std::ofstream init_pos_dis_file("dist_initial_pos.txt");

    yoshidaEngine engine;

    // //this is the numer of alpha particle fired
    int num_shot = 1000000;
    double dt = 1E-22;

    double unitc = 1.6E-19;
    double eps= 8.8541878176E-12;
    double amass= 6.64E-27;
    double vel= 1.57E7;
    double en = 0.5*amass*vel*vel;
    double n = (2*79*unitc*unitc)/(8*M_PI*eps*en);
    std::cout<<"n: "<< n<<'\n';

    for (int i = 0; i < num_shot; ++i)
    {
        //std::cout << i << '\n';
        double init_y = d(eng);
        init_pos_dis_file << init_y << '\n';

        double theorical_outcam = 2 * atangent(n,init_y);
        angle_file<<theorical_outcam<<'\n';
        // alpha alpha{6.64E-27,2, init_y, {-4.5e-14, init_y}, {1778866, 0}}; // ther was 1.53E7
        // engine.evolve(alpha, gold, rilevator_radius, angle_file, dt);
    }

    // alpha alpha{6.64E-27,2, 0, {-4.5e-14, 0}, {1778866, 0}}; // ther was 1.53E7
    // engine.evolve(alpha, gold, rilevator_radius, angle_file, dt);

    //alpha alpha{6.64E-27, 2, 0, {-4.5e-14, 0}, {1778866, 0}}; // ther was 1.53E7
    // engine.applyForce(alpha, gold);
    // std::cout<<"position alpha: ";
    // print(alpha.get_position());
    // std::cout<<"position gold: ";
    // print(gold.get_position());
    // std::cout<<"charges: "<<alpha.get_charge()<<", "<<gold.get_charge();
    // std::cout<<"acceleration: ";
    // print(alpha.get_acceleration());


    // engine.evolve(alpha, gold, rilevator_radius, angle_file, dt);

    return 0;
}
