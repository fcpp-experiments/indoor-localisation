// Copyright © 2026 Giorgio Audrito and Leonardo Bertolino. All Rights Reserved.

/**
 * @file graphic.cpp
 * @brief Runs a single execution of the cooperative indoor localisation case study with a graphical user interface.
 */

#include "lib/localisation.hpp"

using namespace fcpp;

//! @brief The main function.
int main(int argc, char *argv[]) {
    using namespace fcpp;

    option::plot_t p;
    std::cout << "/*\n";
    {
        // Simulation parameters: communication radius and general variance
        real_t comm_radius = 100;
        real_t variance = 0.2;
        if (argc == 3) {
            comm_radius = std::atof(argv[1]);
            variance = std::atof(argv[2]);
        }
        std::weibull_distribution<real_t> distr = distribution::make<std::weibull_distribution>(real_t(1.0), variance);
        real_t half_radius = 100 - 100*variance;
        // The network object type (interactive simulator with given options).
        using net_t = component::interactive_simulator<option::list>::net;
        // The initialisation values (simulation name).
        auto init_v = common::make_tagged_tuple<option::name, option::plotter, option::radius, option::half_radius, option::variance, option::random>("Cooperative Indoor Localisation", &p, comm_radius, half_radius, variance, distr);
        // Construct the network object.
        net_t network{init_v};
        // Run the simulation until exit.
        network.run();
    }
    std::cout << "*/\n";
    std::cout << plot::file("localisation", p.build());
    return 0;
}
