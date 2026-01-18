// Copyright © 2026 Giorgio Audrito and Leonardo Bertolino. All Rights Reserved.

/**
 * @file graphic.cpp
 * @brief Runs a single execution of the cooperative indoor localisation case study with a graphical user interface.
 */

#include "lib/localisation.hpp"

using namespace fcpp;

//! @brief The main function.
int main() {
    using namespace fcpp;

    option::plot_t p;
    std::cout << "/*\n";
    {
        //! @brief The network object type (interactive simulator with given options).
        using net_t = component::interactive_simulator<option::list>::net;
        //! @brief The initialisation values (simulation name).
        auto init_v = common::make_tagged_tuple<option::name, option::plotter>("Cooperative Indoor Localisation", &p);
        //! @brief Construct the network object.
        net_t network{init_v};
        //! @brief Run the simulation until exit.
        network.run();
    }
    std::cout << "*/\n";
    std::cout << plot::file("localisation", p.build());
    return 0;
}
