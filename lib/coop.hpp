#ifndef COOP_H_
#define COOP_H_

#include "lib/coordination.hpp"
#include "lib/data.hpp"

/**
 * @brief Namespace containing all the objects in the FCPP library.
 */
namespace fcpp {

//! @brief Namespace containing the libraries of coordination routines.
namespace coordination {

//! @brief Non-bayesian cooperative localization.
FUN vec<2> nb_coop(ARGS, bool is_anchor, field<real_t> nbr_dist){ CODE
    vec<2> init = make_vec(node.next_real(0,500), node.next_real(0,500));
    return nbr(CALL, init, [&](field<vec<2>> nbr_pos) {
        auto nbr_pos_dist = list_hood(CALL, std::vector<tuple<vec<2>, double>>{}, make_tuple(nbr_pos, nbr_dist));
        if (is_anchor) return node.position();
        real_t alpha = 0.1;
        vec<2> pos = self(CALL, nbr_pos);
        for (auto const& t : nbr_pos_dist) {
            vec<2> pos_stim = get<0>(t);
            real_t sensed_dist = get<1>(t);
            real_t pos_dist = distance(pos, pos_stim);
            real_t delta = sensed_dist - pos_dist;
            vec<2> diff = pos - pos_stim;
            real_t length = norm(diff);
            if (length > 1e-6) pos += (alpha * delta / length) * diff;
        }
        return pos;
    });
}
//! @brief Export list for coop.
FUN_EXPORT nb_coop_t = export_list<vec<2>>;

} // namespace coordination

} // namespace fcpp

#endif // COOP_H_
