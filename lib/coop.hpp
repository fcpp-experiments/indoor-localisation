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
    return nbr(CALL, init, [&](field<vec<2>> prev) {
        field<bool> nbr_anchor = nbr(CALL, is_anchor);
        tuple<field<double>, field<vec<2>>, field<bool>> map = make_tuple(nbr_dist, prev, nbr_anchor);
        auto hop_map_all = list_hood(CALL, std::vector<tuple<double, vec<2>, bool>>{}, map);

        double real_distance, supposed_distance;
        double alpha = 0.1;
        vec<2> searchPos = self(CALL, prev);

        for (auto const& t : hop_map_all) {
            bool anc = get<2>(t);
            if (!is_anchor){
                vec<2> pos_stim = get<1>(t);

                real_distance = get<0>(t);
                supposed_distance = distance(searchPos, pos_stim);
                double ex = real_distance - supposed_distance;

                //Vettore direzione dal nodo al vicino
                vec<2> diff = searchPos - pos_stim;
                double length = norm(diff); // lunghezza del vettore differenza
                if (length > 0) {
                    vec<2> dir = diff / length; // direzione unitaria
                                              // 4. Aggiorna posizione
                    searchPos += alpha * ex * dir;
                }
            }
        }
        return searchPos;
    });
}
//! @brief Export list for coop.
FUN_EXPORT nb_coop_t = export_list<vec<2>, bool>;

} // namespace coordination

} // namespace fcpp

#endif // COOP_H_
