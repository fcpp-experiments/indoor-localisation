#ifndef DV_H_
#define DV_H_

#include "lib/common/option.hpp"
#include "lib/coordination/spreading.hpp"
#include "lib/data/vec.hpp"

#include "multilateration.hpp"

/**
 * @brief Namespace containing all the objects in the FCPP library.
 */
namespace fcpp {

//! @brief Namespace containing the libraries of coordination routines.
namespace coordination {

//! @brief Estimates the node position by multilateration with every other anchor.
FUN vec<2> dv(ARGS, bool is_anchor, field<real_t> nbr_dist, real_t info_speed){ CODE
    std::unordered_map<device_t,real_t> anchor_distance_map;
    std::unordered_map<device_t,vec<2>> anchor_pos_map;
    common::option<device_t> my_anchor_keys;
    if (is_anchor) my_anchor_keys = node.uid;

    old(CALL, 1.0, [&](real_t correction){
        auto hop_map_all = spawn(CALL, [&](device_t anchor_id){
            real_t dist = bis_distance(CALL, node.uid == anchor_id, 1, info_speed, [&](){
                return nbr_dist;
            });
            auto t = broadcast(CALL, dist, make_tuple(node.position(), correction));
            return make_tuple(make_tuple(dist, t), true);
        }, my_anchor_keys);
        real_t apx_dist = 0;
        real_t true_dist = 0;
        for (auto const& [id, t] : hop_map_all) {
            real_t dist = get<0>(t);
            vec<2> pos = get<0>(get<1>(t));
            real_t corr = get<1>(get<1>(t));
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
//! @brief Export list for dv.
FUN_EXPORT dv_t = export_list<real_t, spawn_t<device_t, bool>, broadcast_t<real_t, tuple<vec<2>, real_t>>, bis_distance_t>;

} // namespace coordination

} // namespace fcpp

#endif
