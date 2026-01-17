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

    // usage of node storage
    node.storage(node_size{})  = 10;
    node.storage(node_shape{}) = shape::sphere;   

}
//! @brief Export types used by the main function (update it when expanding the program).
FUN_EXPORT main_t = export_list<dv_t, ksource_t, nBayesianCoop_t>;

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
//! @brief The contents of the node storage as tags and associated types.
using store_t = tuple_store<
    debug,                      std::string,
    node_color,                 color,
    node_size,                  double,
    node_shape,                 shape,
    is_anchor,                  bool,
    pos<dv_hop>,                vec<2>,
    pos<dv_real>,               vec<2>,
    pos<ksource_hop>,           vec<2>,
    pos<ksource_real>,          vec<2>,
    pos<coop>,                  vec<2>,
    error<dv_hop>,              double,
    error<dv_real>,             double,
    error<ksource_hop>,         double,
    error<ksource_real>,        double,
    error<coop>,                double,
    msg_size<dv_hop>,           size_t,
    msg_size<dv_real>,          size_t,
    msg_size<ksource_hop>,      size_t,
    msg_size<ksource_real>,     size_t,
    msg_size<coop>,             size_t
>;
//! @brief The tags and corresponding aggregators to be logged (change as needed).
using aggregator_t = aggregators<
    error<dv_hop>,          aggregator::mean<double>,
    error<dv_real>,         aggregator::mean<double>,
    error<ksource_hop>,     aggregator::mean<double>,
    error<ksource_real>,    aggregator::mean<double>,
    error<coop>,            aggregator::mean<double>,
    msg_size<dv_hop>,       aggregator::mean<double>,
    msg_size<dv_real>,      aggregator::mean<double>,
    msg_size<ksource_hop>,  aggregator::mean<double>,
    msg_size<ksource_real>, aggregator::mean<double>,
    msg_size<coop>,         aggregator::mean<double>
>;
//! @brief Plot of error over time.
using error_plot_t = plot::plotter<aggregator_t, plot::time, error>;
//! @brief Plot of message size over time.
using msize_plot_t = plot::plotter<aggregator_t, plot::time, msg_size>;
//! @brief Plotter class for all plots.
using plot_t = plot::join<error_plot_t, msize_plot_t>;

//! @brief The general simulation options.
DECLARE_OPTIONS(list,
    parallel<true>,      // multithreading enabled on node rounds
    synchronised<false>, // optimise for asynchronous networks
    program<coordination::main>,   // program to be run (refers to MAIN above)
    exports<coordination::main_t, std::unordered_set<int, fcpp::common::hash<int>>, std::unordered_map<int, int, fcpp::common::hash<int>>, tuple<int, int>>, // export type list (types used in messages)
    retain<metric::retain<2,1>>,   // messages are kept for 2 seconds before expiring
    round_schedule<round_s>, // the sequence generator for round events on nodes
    log_schedule<log_s>,     // the sequence generator for log events on the network
    spawn_schedule<spawn_s>, // the sequence generator of node creation events on the network
    store_t,       // the contents of the node storage
    aggregator_t,  // the tags and corresponding aggregators to be logged
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
