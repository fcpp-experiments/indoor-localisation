// Copyright © 2026 Giorgio Audrito and Leonardo Bertolino. All Rights Reserved.

/**
 * @file graphic.cpp
 * @brief Runs a single execution of the aggregate indoor localisation case study with a graphical user interface.
 */

#include "lib/localisation.hpp"

using namespace fcpp;

//! @brief The main function.
int main(int argc, char *argv[]) {
    using namespace fcpp;

    // The plotter object.
    option::gui_plot p;
    std::cout << "/*\n";
    {
        // Calculate simulation parameters.
        real_t comm_radius = option::def_rad;
        real_t variance = option::def_var;
        real_t speed = option::def_v;
        std::string algo = "mlcoop_real";
        if (argc == 5) {
            comm_radius = std::atof(argv[1]);
            variance = std::atof(argv[2]);
            speed = std::atof(argv[3]);
            algo = argv[4];
        }
        real_t half_radius = 100 - variance;
        variance /= 100;
        std::weibull_distribution<real_t> distr = distribution::make<std::weibull_distribution>(real_t(1.0), variance);
        // The network object type (interactive simulator with given options).
        using net_t = component::interactive_simulator<option::list<false>>::net;
        // The initialisation values (simulation name).
        auto init_v = common::make_tagged_tuple_t(
            option::name{},         "Cooperative Indoor Localisation",
            option::plotter{},      &p,
            option::radius{},       comm_radius,
            option::half_radius{},  half_radius,
            option::variance{},     variance,
            option::random{},       distr,
            option::speed{},        speed,
            option::display{},      algo
        );
        // Construct the network object.
        net_t network{init_v};
        // Run the simulation until exit.
        network.run();
    }
    std::cout << "*/\n";
    // Builds the resulting plots.
    std::cout << plot::file("graphic", p.build(), {{"MAX_CROP", "1.05"}, {"LOG_LIN", "10"}, {"SIGMA", "0.1"}});
    return 0;
}
