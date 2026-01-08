#ifndef FCPP_NON_BAYESIAN_COOPERATIVE_LOOCALIZATION_H_
#define FCPP_NON_BAYESIAN_COOPERATIVE_LOOCALIZATION_H_
#include "lib/coordination.hpp"
#include "lib/data.hpp"

namespace fcpp{
namespace coordination{
    FUN vec<2> nBayesianCoop(ARGS, double x_stimato, double y_stimato, bool is_anchor){ CODE
        field<tuple<vec<2>, vec<2>, bool>> map = nbr(CALL, make_tuple(node.position(), make_vec(x_stimato, y_stimato), is_anchor));
        auto hop_map_all = list_hood(CALL,  std::vector< tuple<vec<2>, vec<2>, bool> >{}, map);

        int real_distance, supposed_distance;
        double alpha = 0.1;

        vec<2> searchPos = make_vec(x_stimato, y_stimato);

        for (auto const& t : hop_map_all) {
            bool anc = get<2>(t);
            if (!is_anchor){
                vec<2> pos = get<0>(t);
                vec<2> pos_stim = get<1>(t);
                
                real_distance = distance(node.position(), pos);
                supposed_distance = distance(searchPos, pos_stim);

                int ex = real_distance - supposed_distance;

                //Vettore direzione dal nodo al vicino
                vec<2> diff = searchPos - pos_stim;
                double norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1]); // lunghezza del vettore differenza
                if (norm > 0) {
                    vec<2> dir = diff / norm; // direzione unitaria
                    // 4. Aggiorna posizione
                    searchPos += alpha * ex * dir;
                    x_stimato = searchPos[0];
                    y_stimato = searchPos[1];
                }
            }
        }

        return make_vec(x_stimato, y_stimato);
    }
    FUN_EXPORT nBayesianCoop_t = export_list<tuple<vec<2>, vec<2>, bool>, double, int, broadcast_t<int, tuple<vec<2>, double>>>;
}
}

#endif