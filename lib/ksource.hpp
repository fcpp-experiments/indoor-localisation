#ifndef KSOURCE_H_
#define KSOURCE_H_

#include "lib/coordination/spreading.hpp"
#include "lib/data/vec.hpp"

#include "multilateration.hpp"

/**
 * @brief Namespace containing all the objects in the FCPP library.
 */
namespace fcpp {

//! @brief Namespace containing the libraries of coordination routines.
namespace coordination {

//! @brief Estimates the node position by multilateration with the 20 closest anchors.
FUN vec<2> ksource(ARGS, bool is_anchor, field<real_t> nbr_dist, real_t info_speed) { CODE
    std::unordered_map<device_t,real_t> anchor_distance_map;
    std::unordered_map<device_t,vec<2>> anchor_pos_map;

    old(CALL, 1.0, [&](real_t correction){
        auto hop_map_all = bis_ksource_broadcast(CALL, is_anchor, make_tuple(node.position(), correction), 6, 1, info_speed, [&](){
            return nbr_dist;
        });
        real_t apx_dist = 0;
        real_t true_dist = 0;
        for (auto const& [id, t] : hop_map_all) {
            real_t dist = get<0>(t);
            vec<2> pos = get<0>(get<2>(t));
            real_t corr = get<1>(get<2>(t));
            if (is_anchor) {
                true_dist += distance(node.position(), pos);
                apx_dist += dist;
            } else {
                anchor_pos_map[id] = pos;
                anchor_distance_map[id] = dist * corr;
            }
        }
        if (is_anchor && true_dist != 0 && apx_dist != 0)
            correction = true_dist/apx_dist;
        return correction;
    });
    if (is_anchor) return node.position();
    return multilateration(anchor_pos_map, anchor_distance_map);
}
//! @brief Export list for ksource.
FUN_EXPORT ksource_t = export_list<bis_ksource_broadcast_t<tuple<vec<2>, real_t>>, real_t>;

} // namespace coordination

} // namespace fcpp

#endif
