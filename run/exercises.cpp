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
    //! @brief Color of the current node.
    struct node_color {};
    //! @brief Size of the current node.
    struct node_size {};
    //! @brief Shape of the current node.
    struct node_shape {};
    //! @brief is anchor or not
    struct is_anchor {};

    //! @brief pos stimata
    struct pos_dv_real {};
    //! @brief pos stimata
    struct pos_dv_hop {};
    //! @brief pos stimata
    struct pos_ksource_real {};
    //! @brief pos stimata
    struct pos_ksource_hop {};
    //! @brief pos stimata
    struct pos_coop {};

    //! @brief error misurazione
    struct dist_dv_real {};
    //! @brief error misurazione
    struct dist_dv_hop {};
    //! @brief error misurazione
    struct dist_ksource_real {};
    //! @brief error misurazione
    struct dist_ksource_hop {};
    //! @brief error misurazione
    struct dist_coop {};
    
    template <typename T>
    struct error {};
    
}

//! @brief The maximum communication range between nodes.
constexpr size_t communication_range = 100;

// @brief Main function.
MAIN() {
    // import tag names in the local scope.
    using namespace tags;

    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_int_distribution<> dis(0, 500); 
    
    // GRIGLIA or PERIMETRO
    AnchorLayout anchor_layout = PERIMETRO;

    int id = node.uid;
    double side = 500.0;
    double step = 100.0;
    int rows = 4;
    int cols = 5;
    int total_anchors = rows * cols;

    if (id < total_anchors) {
        node.position() = positionAnchor(
            CALL,
            id,
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
         if (node.current_time() == 0){
            node.storage(pos_coop{})[0] = dis(gen);
            node.storage(pos_coop{})[1] = dis(gen);
        }  
    }

    std::vector<int> my_anchor_keys;
    if (node.storage(is_anchor{}))
        my_anchor_keys = { id };

    node.storage(pos_dv_hop{}) = dvHop(CALL, id, node.storage(is_anchor{}), my_anchor_keys, 1, 0);
    node.storage(pos_dv_real{}) = dvHop(CALL, id, node.storage(is_anchor{}), my_anchor_keys, node.nbr_dist(), 80);
    node.storage(pos_ksource_hop{}) = bis_ksource(CALL, node.storage(is_anchor{}), 1, 0); 
    node.storage(pos_ksource_real{}) = bis_ksource(CALL, node.storage(is_anchor{}), node.nbr_dist(), 80);
    node.storage(pos_coop{}) = nBayesianCoop(CALL, node.storage(pos_coop{})[0], node.storage(pos_coop{})[1], node.storage(is_anchor{}));

    node.storage(error<dist_dv_hop>{}) = distance(node.position(), node.storage(pos_dv_hop{}));
    node.storage(error<dist_dv_real>{}) = distance(node.position(), node.storage(pos_dv_real{}));
    node.storage(error<dist_ksource_real>{}) = distance(node.position(), node.storage(pos_ksource_real{}));
    node.storage(error<dist_ksource_hop>{}) = distance(node.position(), node.storage(pos_ksource_hop{}));
    node.storage(error<dist_coop>{}) = distance(node.position(), node.storage(pos_coop{}));
    
    // usage of node storage
    node.storage(node_size{})  = 10;
    node.storage(node_shape{}) = shape::sphere;   

}
//! @brief Export types used by the main function (update it when expanding the program).
FUN_EXPORT main_t = export_list<dvHop_t, bis_ksource_t, nBayesianCoop_t>;

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
    node_color,                 color,
    node_size,                  double,
    node_shape,                 shape,
    is_anchor,                  bool,
    pos_dv_hop,                 vec<2>,
    pos_dv_real,                vec<2>,
    pos_ksource_real,           vec<2>,
    pos_ksource_hop,            vec<2>,
    pos_coop,                   vec<2>,
    error<dist_ksource_real>,   double,
    error<dist_dv_hop>,         double,
    error<dist_coop>,           double,
    error<dist_ksource_hop>,    double,
    error<dist_dv_real>,        double
    
>;
//! @brief The tags and corresponding aggregators to be logged (change as needed).
using aggregator_t = aggregators<
    error<dist_ksource_real>,   aggregator::mean<double>,
    error<dist_coop>,           aggregator::mean<double>,
    error<dist_dv_hop>,         aggregator::mean<double>,
    error<dist_ksource_hop>,    aggregator::mean<double>,
    error<dist_dv_real>,        aggregator::mean<double>
>;

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

    //! @brief The network object type (interactive simulator with given options).
    using net_t = component::interactive_simulator<option::list>::net;
    //! @brief The initialisation values (simulation name).
    auto init_v = common::make_tagged_tuple<option::name>("Exercises");
    //! @brief Construct the network object.
    net_t network{init_v};
    //! @brief Run the simulation until exit.
    network.run();
    return 0;
}
