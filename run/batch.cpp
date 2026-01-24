// Copyright © 2026 Giorgio Audrito and Leonardo Bertolino. All Rights Reserved.

/**
 * @file batch.cpp
 * @brief Runs a batch of executions of the aggregate indoor localisation case study.
 */

#include "lib/localisation.hpp"

using namespace fcpp;

//! @brief The main function.
int main() {
    using namespace fcpp;

    // The plotter object.
    option::batch_plot p;
    // The component type (batch simulator with given options).
    using comp_t = component::batch_simulator<option::list<true>>;
    // The list of initialisation values to be used for simulations.
    auto init_list = batch::make_tagged_tuple_sequence(
        batch::arithmetic<option::seed       >( 0,  99,  1),                        // 100 different random seeds
        batch::arithmetic<option::half_radius>(50,  98,  2, (int)option::def_hr),   //  26 different variances
        batch::arithmetic<option::radius     >(50, 300, 10, (int)option::def_rad),  //  26 different communication radiuses
        // generate output file name for the run
        batch::stringify<option::output>("output/batch", "txt"),
        // computes half radius from variance
        batch::formula<option::variance, real_t>([](auto const& x) {
            return real_t(100 - common::get<option::half_radius>(x)) / 100;
        }),
        // computes random distribution from variance
        batch::formula<option::random, std::weibull_distribution<real_t>>([](auto const& x) {
            return distribution::make<std::weibull_distribution>(real_t(1.0), (real_t)common::get<option::variance>(x));
        }),
        batch::constant<option::plotter>(&p) // reference to the plotter object
    );
    // Runs the given simulations.
    batch::run(comp_t{}, init_list);
    // Builds the resulting plots.
    std::cout << plot::file("batch", p.build());
    return 0;
}
