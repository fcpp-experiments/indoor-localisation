// Copyright © 2026 Giorgio Audrito and Leonardo Bertolino. All Rights Reserved.

/**
 * @file localisation.cpp
 * @brief Aggregate indoor localisation case study.
 */

#ifndef LOCALISATION_H_
#define LOCALISATION_H_

#include "lib/fcpp.hpp"
#include "lib/dv.hpp"
#include "lib/coop.hpp"

/**
 * @brief Namespace containing all the objects in the FCPP library.
 */
namespace fcpp {

//! @brief The time at which part of the devices are dead.
constexpr size_t bad_time = 50;

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
    struct dv_all_real {};
    //! @brief dv hop algorithm
    struct dv_all_hop {};
    //! @brief ksource real algorithm
    struct dv_6close_real {};
    //! @brief ksource hop algorithm
    struct dv_6close_hop {};
    //! @brief nbcoop real algorithm
    struct nbcoop_real {};
    //! @brief mlcoop real algorithm
    struct mlcoop_real {};
    //! @brief wmlcoop real algorithm
    struct wmlcoop_real {};

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
    tags::error<A>,     real_t,
    tags::msg_size<A>,  size_t
>;
//! @brief Aggregator list for function monitor_algorithm.
GEN_EXPORT(A) monitor_algorithm_a = storage_list<
    tags::error<A>,     aggregator::mean<real_t>,
    tags::msg_size<A>,  aggregator::mean<real_t>
>;


// @brief Main function.
MAIN() {
    // import tag names in the local scope.
    using namespace tags;
    // usage of node storage
    node.storage(node_size{})  = node.storage(is_anchor{}) ? 12 : 8;
    node.storage(node_shape{}) = node.storage(is_anchor{}) ? shape::cube : shape::sphere;
    node.storage(node_color{}) = color(GREEN);
    // 1/4 of the nodes are down between times bad_time and 2*bad_time
    if (node.uid % 4 == 0 and node.current_time() > bad_time and node.next_time() < 2*bad_time) {
        node.storage(node_color{}) = color(RED);
        return;
    }
    // distances with error
    field<real_t> nbr_dist = map_hood([&](real_t d){
        return d * node.storage(random{})(node.generator());
    }, node.nbr_dist());
    // initial random position
    vec<2> init = make_vec(node.next_real(0,500), node.next_real(0,500));

    monitor_algorithm(CALL, dv_all_real{}, [&](){
        return dv_all(CALL, init, node.storage(is_anchor{}), nbr_dist, 80, 1000);
    });
    monitor_algorithm(CALL, dv_all_hop{}, [&](){
        int max_dist = 150000 / (node.net.storage(component::tags::half_radius{})*node.net.storage(component::tags::radius{}));
        return dv_all(CALL, init, node.storage(is_anchor{}), 1, 1, max_dist);
    });
    monitor_algorithm(CALL, dv_6close_real{}, [&](){
        return dv_kclose(CALL, 6, init, node.storage(is_anchor{}), nbr_dist, 80);
    });
    monitor_algorithm(CALL, dv_6close_hop{}, [&](){
        return dv_kclose(CALL, 6, init, node.storage(is_anchor{}), 1, 1);
    });
    monitor_algorithm(CALL, nbcoop_real{}, [&](){
        return nb_coop(CALL, init, node.storage(is_anchor{}), nbr_dist);
    });
    monitor_algorithm(CALL, mlcoop_real{}, [&](){
        return ml_coop(CALL, init, node.storage(is_anchor{}), nbr_dist);
    });
    /*
    monitor_algorithm(CALL, wmlcoop_real{}, [&](){
        real_t aw = 15000 / (node.net.storage(tags::variance{})*node.net.storage(component::tags::half_radius{})*node.net.storage(component::tags::radius{}));
        return wml_coop(CALL, init, node.storage(is_anchor{}), nbr_dist, aw, 0.005);
    });
    */
}
//! @brief Export list for the main function.
FUN_EXPORT main_t = export_list<dv_all_t, dv_kclose_t, nb_coop_t, ml_coop_t, wml_coop_t>;
//! @brief Storage list for the main function.
FUN_EXPORT main_s = storage_list<
    tags::debug,        std::string,
    tags::random,       std::weibull_distribution<real_t>,
    tags::is_anchor,    bool,
    tags::node_color,   color,
    tags::node_size,    real_t,
    tags::node_shape,   shape,
    monitor_algorithm_s<tags::dv_all_real>,
    monitor_algorithm_s<tags::dv_all_hop>,
    monitor_algorithm_s<tags::dv_6close_real>,
    monitor_algorithm_s<tags::dv_6close_hop>,
    monitor_algorithm_s<tags::nbcoop_real>,
    monitor_algorithm_s<tags::mlcoop_real>
>;
//! @brief Aggregator list for the main function.
FUN_EXPORT main_a = storage_list<
    monitor_algorithm_a<tags::dv_all_real>,
    monitor_algorithm_a<tags::dv_all_hop>,
    monitor_algorithm_a<tags::dv_6close_real>,
    monitor_algorithm_a<tags::dv_6close_hop>,
    monitor_algorithm_a<tags::nbcoop_real>,
    monitor_algorithm_a<tags::mlcoop_real>
>;

} // namespace coordination


// SYSTEM SETUP

//! @brief Namespace for component options.
namespace option {

//! @brief Import tags to be used for component options.
using namespace component::tags;
//! @brief Import tags used by aggregate functions.
using namespace coordination::tags;

//! @brief End of simulated time.
constexpr size_t end_time = 3*bad_time;

//! @brief Generic plot given X axis, Y axis and filter description Fs
template<typename X, template<class> class Y, typename... Fs>
using general_plot = plot::filter<Fs..., plot::plotter<coordination::main_a, X, Y, common::type_sequence<aggregator::stats<real_t>>>>;

//! @brief The simulation time after which we measure performance.
constexpr size_t mean_time = 0;
//! @brief The default variance simulation parameter.
constexpr size_t def_var = 20;
//! @brief The default half-radius simulation parameter.
constexpr size_t def_hr = 100 - def_var;
//! @brief The default radius simulation parameter.
constexpr size_t def_rad = 150;
//! @brief Plot of error over time.
using error_time_plot = general_plot<plot::time, error,    half_radius, filter::equal<100-def_var>, radius, filter::equal<def_rad>>;
//! @brief Plot of message size over time.
using msize_time_plot = general_plot<plot::time, msg_size, half_radius, filter::equal<100-def_var>, radius, filter::equal<def_rad>>;
//! @brief Plot of error over variance.
using error_var_plot = general_plot<variance, error,    plot::time, filter::above<mean_time>, radius, filter::equal<def_rad>>;
//! @brief Plot of message size over variance.
using msize_var_plot = general_plot<variance, msg_size, plot::time, filter::above<mean_time>, radius, filter::equal<def_rad>>;
//! @brief Plot of error over radius.
using error_rad_plot = general_plot<radius, error,    plot::time, filter::above<mean_time>, half_radius, filter::equal<100-def_var>>;
//! @brief Plot of message size over radius.
using msize_rad_plot = general_plot<radius, msg_size, plot::time, filter::above<mean_time>, half_radius, filter::equal<100-def_var>>;
//! @brief Plotter class for all batch plots.
using batch_plot = plot::join<error_time_plot, msize_time_plot, error_var_plot, msize_var_plot, error_rad_plot, msize_rad_plot>;

//! @brief Plot of error over time.
using error_plot = general_plot<plot::time, error>;
//! @brief Plot of message size over time.
using msize_plot = general_plot<plot::time, msg_size>;
//! @brief Plotter class for all GUI plots.
using gui_plot = plot::join<error_plot, msize_plot>;

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
template <bool batch>
DECLARE_OPTIONS(list,
    parallel<false>,     // multithreading disabled on node rounds
    synchronised<false>, // optimise for asynchronous networks
    program<coordination::main>,            // program to be run (refers to MAIN above)
    exports<coordination::main_t>,          // export type list (types used in messages)
    node_store<coordination::main_s>,       // the contents of the node storage
    net_store<                              // the contents of the net storage
        variance,       real_t,
        half_radius,    size_t,
        radius,         size_t
    >,
    aggregators<coordination::main_a>,      // the tags and corresponding aggregators to be logged
    extra_info<                             // general parameters to use for plotting
        variance,       real_t,
        half_radius,    size_t,
        radius,         size_t
    >,
    plot_type<                              // the plotter object
        std::conditional_t<batch, batch_plot, gui_plot>
    >,
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
