#ifndef COOP_H_
#define COOP_H_

#include "lib/coordination.hpp"
#include "lib/data.hpp"

#include "multilateration.hpp"

/**
 * @brief Namespace containing all the objects in the FCPP library.
 */
namespace fcpp {

//! @brief Namespace containing the libraries of coordination routines.
namespace coordination {

//! @brief Non-bayesian cooperative localization.
FUN vec<2> nb_coop(ARGS, vec<2> init, bool is_anchor, field<real_t> nbr_dist){ CODE
    return nbr(CALL, init, [&](field<vec<2>> nbr_pos) {
        auto nbr_pos_dist = list_hood(CALL, std::vector<tuple<vec<2>, real_t>>{}, make_tuple(nbr_pos, nbr_dist), tags::nothing{});
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


//! @brief Cooperative localization based on multilateration.
FUN vec<2> ml_coop(ARGS, vec<2> init, bool is_anchor, field<real_t> nbr_dist){ CODE
    return nbr(CALL, init, [&](field<vec<2>> nbr_pos) {
        auto pos_dist_vec = list_hood(CALL, std::vector<tuple<vec<2>, real_t>>{}, make_tuple(nbr_pos, nbr_dist), tags::nothing{});
        if (is_anchor) return node.position();
        return multilateration(self(CALL, nbr_pos), pos_dist_vec);
    });
}
//! @brief Export list for ml_coop.
FUN_EXPORT ml_coop_t = export_list<vec<2>>;


/**
 * @brief Cooperative localization based on weighted multilateration.
 *
 * The anchor_weight should be the inverse standard deviation of nbr_dist measurements.
 * If it is a uniform distribution with given standard deviation s multiplied by true nbr_dist measurements,
 * its standard deviation can be approximated as s * half_radius * 2/3.
 * The device_weight should be the inverse standard deviation of init positions.
 * If it is a uniform distribution between (0,0) and (S,S) its standard deviation is S/√6.
 */
FUN vec<2> wml_coop(ARGS, vec<2> init, bool is_anchor, field<real_t> nbr_dist, real_t anchor_weight, real_t device_weight){ CODE
    return nbr(CALL, init, [&](field<vec<2>> nbr_pos) {
        vec<2> pos = node.position();
        nbr(CALL, device_weight, [&](field<real_t> nbr_weights){
            auto pos_dist_vec = list_hood(CALL, std::vector<tuple<vec<2>, real_t>>{}, make_tuple(nbr_pos, nbr_dist), tags::nothing{});
            auto weights_vec = list_hood(CALL, std::vector<real_t>{}, nbr_weights, tags::nothing{});
            if (is_anchor) return anchor_weight;
            real_t weight;
            pos = multilateration(self(CALL, nbr_pos), pos_dist_vec, weights_vec, anchor_weight, weight);
            return weight;
        });
        return pos;
    });
}
//! @brief Export list for wml_coop.
FUN_EXPORT wml_coop_t = export_list<vec<2>, real_t>;

} // namespace coordination

} // namespace fcpp

#endif // COOP_H_
