#include "rat.hpp"

void particle::evolve(const particle &pArticle, double rilevator_radius, std::ofstream &angular_file)
{
    double dt = 1e-17; // Time step
    

    std::vector<gp> transit;
    transit.push_back(_position);

    // Initial half-step for velocity
    gp initial_acc = acceleration_vector(pArticle);
    _velocity += 0.5 * dt * initial_acc;

    double t = 0;

    while (O_distance() < rilevator_radius)
    {
        // Update position
        _position += dt * _velocity;

        // Compute new acceleration
        gp new_acc = acceleration_vector(pArticle);

        // Update velocity with the new acceleration
        _velocity += 0.5 * dt * new_acc;

        // Store position for later analysis
        transit.push_back(_position);

        // Update time and iteration count
        t += dt;
       
       
            // std::cout << "Time: " << t << " dt: " << dt << '\n';
            // std::cout << "Position: ";
            // print(_position);
            // std::cout << "Velocity: ";
            // print(_velocity);
            // std::cout << "acc: ";
            // print(acceleration_vector(pArticle));
        
    }

    // Determine the angular position where the particle exits the rilevator_radius
    auto it = (transit.end() - 1);
    auto is = it - 1;
    double diff_lastpos_radius = std::abs(O_distance() - rilevator_radius);
    double diff_penlastpos_radius = std::abs(std::sqrt(std::pow((*is).a, 2) + std::pow((*is).b, 2)) - rilevator_radius);

    if (diff_lastpos_radius < diff_penlastpos_radius)
    {
        double angular = atangent((*it).b, (*it).a);
        angular_file << angular << '\n';
    }
    else
    {
        double angular = atangent((*is).b, (*is).a);
        angular_file << angular << '\n';
    }
}
