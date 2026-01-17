#ifndef TRILATERATION_H_
#define TRILATERATION_H_
#include "lib/coordination.hpp"
#include "lib/data.hpp"

namespace fcpp{
namespace coordination{
    FUN vec<2> trilaterazione(ARGS, bool is_anchor, std::unordered_map<int,double> anchor_distance_map, std::unordered_map<int,int> anchor_x_map, std::unordered_map<int,int> anchor_y_map){ CODE
        // Se è un'ancora, restituisce (1,1) o la sua posizione reale se disponibile
        if (is_anchor) return make_vec(1, 1);

        auto& dist_map = anchor_distance_map;
        std::vector<int> anchor_ids;
        for (auto const& [id, _] : dist_map) anchor_ids.push_back(id);

        // Variabili di output con valore di default
        double x_est = 0;
        double y_est = 0;
        if (anchor_ids.size() >= 3) {
            int a_ref = anchor_ids[0];
            double x1 = anchor_x_map[a_ref];
            double y1 = anchor_y_map[a_ref];
            double r1 = dist_map[a_ref];

            double ATA11 = 0, ATA12 = 0, ATA22 = 0;
            double ATb1 = 0, ATb2 = 0;

            for (size_t i = 1; i < anchor_ids.size(); ++i) {
                int ai = anchor_ids[i];
                double xi = anchor_x_map[ai];
                double yi = anchor_y_map[ai];
                double ri = dist_map[ai];

                double Ai1 = 2 * (xi - x1);
                double Ai2 = 2 * (yi - y1);
                double bi  = r1*r1 - ri*ri + xi*xi - x1*x1 + yi*yi - y1*y1;

                ATA11 += Ai1 * Ai1;
                ATA12 += Ai1 * Ai2;
                ATA22 += Ai2 * Ai2;
                ATb1  += Ai1 * bi;
                ATb2  += Ai2 * bi;
            }

            double det = ATA11 * ATA22 - ATA12 * ATA12;

            if (std::abs(det) > 1e-9) {
                // Caso normale: sistema risolvibile
                x_est = (ATb1 * ATA22 - ATA12 * ATb2) / det;
                y_est = (ATA11 * ATb2 - ATA12 * ATb1) / det;
            } else {
                // Caso degenerato (ancore allineate): calcoliamo il baricentro delle ancore
                double sum_x = 0, sum_y = 0;
                for (int id : anchor_ids) {
                    sum_x += anchor_x_map[id];
                    sum_y += anchor_y_map[id];
                }
                x_est = sum_x / anchor_ids.size();
                y_est = sum_y / anchor_ids.size();
            }
            return make_vec(x_est, y_est);
        }

        return make_vec(1,1);
    }
    FUN_EXPORT trilaterazione_t = export_list<>;
}
}

#endif
