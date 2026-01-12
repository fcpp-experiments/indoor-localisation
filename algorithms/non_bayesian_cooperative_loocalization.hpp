#ifndef FCPP_NON_BAYESIAN_COOPERATIVE_LOOCALIZATION_H_
#define FCPP_NON_BAYESIAN_COOPERATIVE_LOOCALIZATION_H_
#include "lib/coordination.hpp"
#include "lib/data.hpp"

namespace fcpp{
namespace coordination{
    FUN vec<2> nBayesianCoop(ARGS, bool is_anchor){ CODE

        auto state = old(CALL, vec<2>{node.next_real(0,501), node.next_real(0,501)}, [&](vec<2> prev){ 
            field<tuple<vec<2>, vec<2>, bool>> map = nbr(CALL, make_tuple(node.position(), prev, is_anchor));
            auto hop_map_all = list_hood(CALL,  std::vector< tuple<vec<2>, vec<2>, bool> >{}, map);

            double real_distance, supposed_distance;
            double alpha = 0.1;

            vec<2> searchPos = prev;

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
                    }
                }
            }
            return searchPos; 
        });

    return state;
    }
    FUN_EXPORT nBayesianCoop_t = export_list<tuple<vec<2>, vec<2>, bool>, vec<2>>;
}
}

#endif