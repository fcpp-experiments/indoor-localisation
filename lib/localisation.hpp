// Copyright © 2026 Giorgio Audrito and Leonardo Bertolino. All Rights Reserved.

/**
 * @file localisation.cpp
 * @brief Cooperative indoor localisation case study.
 */

#ifndef LOCALISATION_H_
#define LOCALISATION_H_

#include "lib/fcpp.hpp"
#include "lib/dv.hpp"
#include "lib/ksource.hpp"
#include "lib/coop.hpp"

/**
 * @brief Namespace containing all the objects in the FCPP library.
 */
namespace fcpp {

//! @brief Dummy ordering between positions (allows positions to be used as secondary keys in ordered tuples).
template <size_t n>
bool operator<(vec<n> const&, vec<n> const&) {
    return false;
}

//! @brief Namespace containing the libraries of coordination routines.
namespace coordination {

//! @brief Tags used in the node storage.
namespace tags {
    //! @brief General debugging information.
    struct debug {};
    //! @brief Whether the node is a anchor or not.
    struct is_anchor {};
    //! @brief The level of variance to consider in the simulation.
    struct variance {};
    //! @brief A weibull distribution with a mean of 1 and given variance.
    struct random {};

    //! @brief Color of the current node.
    struct node_color {};
    //! @brief Size of the current node.
    struct node_size {};
    //! @brief Shape of the current node.
    struct node_shape {};

    //! @brief dv real algorithm
    struct dv_real {};
    //! @brief dv hop algorithm
    struct dv_hop {};
    //! @brief ksource real algorithm
    struct ksource12_real {};
    //! @brief ksource hop algorithm
    struct ksource12_hop {};
    //! @brief ksource real algorithm
    struct ksource6_real {};
    //! @brief ksource hop algorithm
    struct ksource6_hop {};
    //! @brief coop algorithm
    struct coop_real {};

    //! @brief estimated position for an algorithm
    template <typename T>
    struct pos {};

    //! @brief distance error for an algorithm
    template <typename T>
    struct error {};
    
    //! @brief message size for an algorithm
    template <typename T>
    struct msg_size {};
}


//! @brief Runs an algorithm and saves monitoring data.
GEN(A, F) void monitor_algorithm(ARGS, A, F&& fun) { CODE
    using namespace tags;
    size_t msiz_pre = node.cur_msg_size();
    node.storage(pos<A>{}) = std::forward<F>(fun)();
    node.storage(error<A>{}) = distance(node.position(), node.storage(pos<A>{}));
    node.storage(msg_size<A>{}) = node.cur_msg_size() - msiz_pre;
}
//! @brief Storage list for function monitor_algorithm.
GEN_EXPORT(A) monitor_algorithm_s = storage_list<
    tags::pos<A>,       vec<2>,
    tags::error<A>,     double,
    tags::msg_size<A>,  size_t
>;
//! @brief Aggregator list for function monitor_algorithm.
GEN_EXPORT(A) monitor_algorithm_a = storage_list<
    tags::error<A>,     aggregator::mean<double>,
    tags::msg_size<A>,  aggregator::mean<double>
>;


// @brief Main function.
MAIN() {
    // import tag names in the local scope.
    using namespace tags;
    field<real_t> nbr_dist = map_hood([&](real_t d){
        return d * node.storage(random{})(node.generator());
    }, node.nbr_dist());

    monitor_algorithm(CALL, dv_real{}, [&](){
        return dv(CALL, node.storage(is_anchor{}), nbr_dist, 80);
    });
    monitor_algorithm(CALL, dv_hop{}, [&](){
        return dv(CALL, node.storage(is_anchor{}), 1, 0);
    });
    monitor_algorithm(CALL, ksource12_real{}, [&](){
        return ksource(CALL, 12, node.storage(is_anchor{}), nbr_dist, 80);
    });
    monitor_algorithm(CALL, ksource12_hop{}, [&](){
        return ksource(CALL, 12, node.storage(is_anchor{}), 1, 0);
    });
    monitor_algorithm(CALL, ksource6_real{}, [&](){
        return ksource(CALL, 6, node.storage(is_anchor{}), nbr_dist, 80);
    });
    monitor_algorithm(CALL, ksource6_hop{}, [&](){
        return ksource(CALL, 6, node.storage(is_anchor{}), 1, 0);
    });
    monitor_algorithm(CALL, coop_real{}, [&](){
        return nb_coop(CALL, node.storage(is_anchor{}), nbr_dist);
    });

    // usage of node storage
    node.storage(node_size{})  = node.storage(is_anchor{}) ? 12 : 8;
    node.storage(node_shape{}) = node.storage(is_anchor{}) ? shape::cube : shape::sphere;
    node.storage(node_color{}) = node.storage(is_anchor{}) ? color(RED) : color(GREEN);
}
//! @brief Export list for the main function.
FUN_EXPORT main_t = export_list<dv_t, ksource_t, nb_coop_t>;
//! @brief Storage list for the main function.
FUN_EXPORT main_s = storage_list<
    tags::debug,        std::string,
    tags::random,       std::weibull_distribution<real_t>,
    tags::is_anchor,    bool,
    tags::node_color,   color,
    tags::node_size,    double,
    tags::node_shape,   shape,
    monitor_algorithm_s<tags::dv_real>,
    monitor_algorithm_s<tags::dv_hop>,
    monitor_algorithm_s<tags::ksource12_real>,
    monitor_algorithm_s<tags::ksource12_hop>,
    monitor_algorithm_s<tags::ksource6_real>,
    monitor_algorithm_s<tags::ksource6_hop>,
    monitor_algorithm_s<tags::coop_real>
>;
//! @brief Aggregator list for the main function.
FUN_EXPORT main_a = storage_list<
    monitor_algorithm_a<tags::dv_real>,
    monitor_algorithm_a<tags::dv_hop>,
    monitor_algorithm_a<tags::ksource12_real>,
    monitor_algorithm_a<tags::ksource12_hop>,
    monitor_algorithm_a<tags::ksource6_real>,
    monitor_algorithm_a<tags::ksource6_hop>,
    monitor_algorithm_a<tags::coop_real>
>;

} // namespace coordination


// SYSTEM SETUP

//! @brief Namespace for component options.
namespace option {

//! @brief Import tags to be used for component options.
using namespace component::tags;
//! @brief Import tags used by aggregate functions.
using namespace coordination::tags;

//! @brief Plot of error over time.
using error_plot_t = plot::plotter<coordination::main_a, plot::time, error>;
//! @brief Plot of message size over time.
using msize_plot_t = plot::plotter<coordination::main_a, plot::time, msg_size>;
//! @brief Plotter class for all plots.
using plot_t = plot::join<error_plot_t, msize_plot_t>;

//! @brief End of simulated time.
constexpr size_t end_time = 100;
//! @brief Description of the round schedule.
using round_s = sequence::periodic<
    distribution::interval_n<times_t, 0, 1>,        // uniform time in the [0,1] interval for start
    distribution::weibull<                          // weibull-distributed time for interval (1±variance)
        distribution::constant_n<times_t, 1>,
        distribution::constant_n<times_t, 1, 10, variance>
    >,
    distribution::constant_n<times_t, end_time+2>   // the constant end_time+2 number for end
>;
//! @brief The sequence of network snapshots (one every simulated second).
using log_s = sequence::periodic_n<1, 0, 1, end_time>;

//! @brief The connection predicate (100% at 0m, 50% at 80m, 0% at 100m).
using connect_t = connect::radial<80, connect::fixed<100>>;

//! @brief The sequence of anchor generation events (20 devices all generated at time 0).
using anchor_spawn_s = sequence::multiple_n<20, 0>;
//! @brief The distribution of initial anchor positions (random in a 500x500 square).
using anchor_pos_d = sequence::rectangle_n<1, 0, 0, 500, 500, 20>;
//! @brief The sequence of device generation events (100 devices all generated at time 0).
using device_spawn_s = sequence::multiple_n<100, 0>;
//! @brief The distribution of initial device positions (random in a 500x500 square).
using device_pos_d = distribution::rect_n<1, 0, 0, 500, 500>;

//! @brief The general simulation options.
DECLARE_OPTIONS(list,
    parallel<true>,      // multithreading enabled on node rounds
    synchronised<false>, // optimise for asynchronous networks
    program<coordination::main>,            // program to be run (refers to MAIN above)
    exports<coordination::main_t>,          // export type list (types used in messages)
    node_store<coordination::main_s>,       // the contents of the node storage
    aggregators<coordination::main_a>,      // the tags and corresponding aggregators to be logged
    plot_type<plot_t>,                      // the plotter object
    connector<connect_t>,                   // connection predicate
    retain<metric::retain<5,1>>,            // messages are kept for 5 seconds before expiring
    round_schedule<round_s>,                // the sequence generator for round events on nodes
    log_schedule<log_s>,                    // the sequence generator for log events on the network
    spawn_schedule<anchor_spawn_s>,         // the sequence generator of anchor creation events on the network
    init<
        random,     distribution::constant_i<std::weibull_distribution<real_t>, random>,
        variance,   distribution::constant_i<real_t, variance>,
        is_anchor,  distribution::constant_n<bool, true>,
        x,          anchor_pos_d
    >,
    spawn_schedule<device_spawn_s>,         // the sequence generator of device creation events on the network
    init<
        random,     distribution::constant_i<std::weibull_distribution<real_t>, random>,
        variance,   distribution::constant_i<real_t, variance>,
        is_anchor,  distribution::constant_n<bool, false>,
        x,          device_pos_d
    >,
    dimension<2>,           // dimensionality of the space
    shape_tag<node_shape>,  // the shape of a node is read from this tag in the store
    size_tag<node_size>,    // the size  of a node is read from this tag in the store
    color_tag<node_color>   // the color of a node is read from this tag in the store
);

} // namespace option

} // namespace fcpp

#endif // LOCALISATION_H_
