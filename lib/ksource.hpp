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
FUN vec<2> ksource(ARGS, int k, vec<2> init, bool is_anchor, field<real_t> nbr_dist, real_t info_speed) { CODE
    std::vector<tuple<vec<2>,real_t>> anchors;

    return old(CALL, init, [&](vec<2> pos){
        old(CALL, 1.0, [&](real_t correction){
            auto anchor_map = bis_ksource_broadcast(CALL, is_anchor, make_tuple(node.position(), correction), k, 1, info_speed, [&](){
                return nbr_dist;
            });
            real_t apx_dist = 0;
            real_t true_dist = 0;
            for (auto const& t : anchor_map) {
                real_t dist = get<0>(t.second);
                vec<2> pos = get<0>(get<2>(t.second));
                real_t corr = get<1>(get<2>(t.second));
                if (is_anchor) {
                    true_dist += distance(node.position(), pos);
                    apx_dist += dist;
                } else {
                    anchors.emplace_back(pos, dist * corr);
                }
            }
            if (is_anchor && true_dist != 0 && apx_dist != 0)
                correction = true_dist/apx_dist;
            return correction;
        });
        if (is_anchor) return node.position();
        return multilateration(pos, anchors);
    });
}
//! @brief Export list for ksource.
FUN_EXPORT ksource_t = export_list<vec<2>, bis_ksource_broadcast_t<tuple<vec<2>, real_t>>, real_t>;

} // namespace coordination

} // namespace fcpp

#endif // KSOURCE_H_
