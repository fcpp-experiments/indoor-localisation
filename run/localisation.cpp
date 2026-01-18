// Copyright © 2021 Giorgio Audrito. All Rights Reserved.

/**
 * @file exercises.cpp
 * @brief Quick-start aggregate computing exercises.
 */

// [INTRODUCTION]
//! Importing the FCPP library.
#include "lib/fcpp.hpp"
#include "algorithms/dvHop.hpp"
#include "algorithms/bis_ksource_broadcast.hpp"
#include "algorithms/non_bayesian_cooperative_loocalization.hpp"
#include "anchor/anchorLayout.hpp"

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

    //! @brief Color of the current node.
    struct node_color {};
    //! @brief Size of the current node.
    struct node_size {};
    //! @brief Shape of the current node.
    struct node_shape {};
    //! @brief is anchor or not
    struct is_anchor {};

    //! @brief dv real algorithm
    struct dv_real {};
    //! @brief dv hop algorithm
    struct dv_hop {};
    //! @brief ksource real algorithm
    struct ksource_real {};
    //! @brief ksource hop algorithm
    struct ksource_hop {};
    //! @brief coop algorithm
    struct coop {};
    //! @brief ideal reference
    struct ideal {};
    
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

//! @brief The maximum communication range between nodes.
constexpr size_t communication_range = 100;


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
    
    // GRIGLIA or PERIMETRO
    AnchorLayout anchor_layout = PERIMETRO;

    double side = 500.0;
    double step = 100.0;
    int rows = 4;
    int cols = 5;
    int total_anchors = rows * cols;

    if (node.uid < total_anchors) {
        node.position() = positionAnchor(
            CALL,
            node.uid,
            anchor_layout,
            side, 
            step, 
            rows,
            cols
        );
        node.storage(is_anchor{}) = true;
        node.storage(node_color{}) = color(RED);

    } else {
        node.storage(is_anchor{}) = false;
        node.storage(node_color{}) = color(GREEN);
    }

    monitor_algorithm(CALL, dv_hop{}, [&](){
        return dv(CALL, node.storage(is_anchor{}), 1, 0);
    });
    monitor_algorithm(CALL, dv_real{}, [&](){
        return dv(CALL, node.storage(is_anchor{}), node.nbr_dist(), 80);
    });
    monitor_algorithm(CALL, ksource_hop{}, [&](){
        return ksource(CALL, node.storage(is_anchor{}), 1, 0);
    });
    monitor_algorithm(CALL, ksource_real{}, [&](){
        return ksource(CALL, node.storage(is_anchor{}), node.nbr_dist(), 80);
    });
    monitor_algorithm(CALL, coop{}, [&](){
        return nBayesianCoop(CALL, node.storage(is_anchor{}));
    });
    monitor_algorithm(CALL, ideal{}, [&](){
        return node.position();
    });

    // usage of node storage
    node.storage(node_size{})  = 10;
    node.storage(node_shape{}) = shape::sphere;   

}
//! @brief Export list for the main function.
FUN_EXPORT main_t = export_list<dv_t, ksource_t, nBayesianCoop_t>;
//! @brief Storage list for the main function.
FUN_EXPORT main_s = storage_list<
    tags::debug,        std::string,
    tags::is_anchor,    bool,
    tags::node_color,   color,
    tags::node_size,    double,
    tags::node_shape,   shape,
    monitor_algorithm_s<tags::dv_hop>,
    monitor_algorithm_s<tags::dv_real>,
    monitor_algorithm_s<tags::ksource_hop>,
    monitor_algorithm_s<tags::ksource_real>,
    monitor_algorithm_s<tags::coop>,
    monitor_algorithm_s<tags::ideal>
>;
//! @brief Aggregator list for the main function.
FUN_EXPORT main_a = storage_list<
    monitor_algorithm_a<tags::dv_hop>,
    monitor_algorithm_a<tags::dv_real>,
    monitor_algorithm_a<tags::ksource_hop>,
    monitor_algorithm_a<tags::ksource_real>,
    monitor_algorithm_a<tags::coop>,
    monitor_algorithm_a<tags::ideal>
>;

} // namespace coordination

// [SYSTEM SETUP]

//! @brief Namespace for component options.
namespace option {

//! @brief Import tags to be used for component options.
using namespace component::tags;
//! @brief Import tags used by aggregate functions.
using namespace coordination::tags;

//! @brief Number of people in the area.
constexpr int node_num = 100;
//! @brief Dimensionality of the space.
constexpr size_t dim = 2;

//! @brief Plot of error over time.
using error_plot_t = plot::plotter<coordination::main_a, plot::time, error>;
//! @brief Plot of message size over time.
using msize_plot_t = plot::plotter<coordination::main_a, plot::time, msg_size>;
//! @brief Plotter class for all plots.
using plot_t = plot::join<error_plot_t, msize_plot_t>;

//! @brief Description of the round schedule.
using round_s = sequence::periodic<
    distribution::interval_n<times_t, 0, 1>,    // uniform time in the [0,1] interval for start
    distribution::weibull_n<times_t, 10, 1, 10> // weibull-distributed time for interval (10/10=1 mean, 1/10=0.1 deviation)
>;
//! @brief The sequence of network snapshots (one every simulated second).
using log_s = sequence::periodic_n<1, 0, 1>;
//! @brief The sequence of node generation events (node_num devices all generated at time 0).
using spawn_s = sequence::multiple_n<node_num, 0>;
//! @brief The distribution of initial node positions (random in a 500x500 square).
using rectangle_d = distribution::rect_n<1, 0, 0, 500, 500>;

//! @brief The general simulation options.
DECLARE_OPTIONS(list,
    parallel<true>,      // multithreading enabled on node rounds
    synchronised<false>, // optimise for asynchronous networks
    program<coordination::main>,        // program to be run (refers to MAIN above)
    exports<coordination::main_t>,      // export type list (types used in messages)
    retain<metric::retain<2,1>>,        // messages are kept for 2 seconds before expiring
    round_schedule<round_s>,            // the sequence generator for round events on nodes
    log_schedule<log_s>,                // the sequence generator for log events on the network
    spawn_schedule<spawn_s>,            // the sequence generator of node creation events on the network
    tuple_store<coordination::main_s>,  // the contents of the node storage
    aggregators<coordination::main_a>,  // the tags and corresponding aggregators to be logged
    plot_type<plot_t>, // the plotter object
    init<
        x,      rectangle_d // initialise position randomly in a rectangle for new nodes
    >,
    dimension<dim>, // dimensionality of the space
    connector<connect::fixed<100, 1, dim>>, // connection allowed within a fixed comm range
    shape_tag<node_shape>, // the shape of a node is read from this tag in the store
    size_tag<node_size>,   // the size  of a node is read from this tag in the store
    color_tag<node_color>  // the color of a node is read from this tag in the store
);

} // namespace option

} // namespace fcpp


//! @brief The main function.
int main() {
    using namespace fcpp;

    option::plot_t p;
    std::cout << "/*\n";
    {
        //! @brief The network object type (interactive simulator with given options).
        using net_t = component::interactive_simulator<option::list>::net;
        //! @brief The initialisation values (simulation name).
        auto init_v = common::make_tagged_tuple<option::name, option::plotter>("Exercises", &p);
        //! @brief Construct the network object.
        net_t network{init_v};
        //! @brief Run the simulation until exit.
        network.run();
    }
    std::cout << "*/\n";
    std::cout << plot::file("localisation", p.build());
    return 0;
}
