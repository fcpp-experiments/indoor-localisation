#ifndef FCPP_BIS_KSOURCE_BROADCAST_H_
#define FCPP_BIS_KSOURCE_BROADCAST_H_
#include "lib/coordination.hpp"
#include "lib/data.hpp"

namespace fcpp{
namespace coordination{
    FUN vec<2> bis_ksource(ARGS, bool is_anchor){CODE
        std::unordered_map<int, int, fcpp::common::hash<int>> hop_map;
        std::unordered_map<int,int, fcpp::common::hash<int>> anchor_x_map;
        std::unordered_map<int,int, fcpp::common::hash<int>> anchor_y_map;
        std::unordered_map<int,int, fcpp::common::hash<int>> anchor_correction_map;
        std::unordered_map<int,double, fcpp::common::hash<int>> anchor_distance_map;

        double x_est;
        double y_est;

        double correction = old(CALL, 0, [&](double correction){
            auto hop_map_all = bis_ksource_broadcast(CALL, is_anchor, make_tuple(node.position(), correction), 5, 1, 80);

            for (auto const& [id, t] : hop_map_all){
                double dist = get<0>(t);
                auto tuple = get<2>(t);
                vec<2> pos = get<0>(tuple);
                double corr = get<1>(tuple);
                hop_map[id] = dist;
                anchor_x_map[id] = pos[0];
                anchor_y_map[id] = pos[1];
                anchor_correction_map[id] = corr;
                if (is_anchor)
                    anchor_distance_map[id] = dist;
                else
                    anchor_distance_map[id] = dist;
            }

            int hop = 0;
            double distance = 0;
            if (is_anchor){
                for (auto const& [key, value] : anchor_distance_map){
                    distance += value;
                    hop += hop_map[key];
                }   
                if (distance != 0 && hop != 0)
                    correction = distance/hop;
            }

            return correction;
        });

        if (!is_anchor){
            auto& dist_map = anchor_distance_map;
            auto& x_map = anchor_x_map;
            auto& y_map = anchor_y_map;

            std::vector<int> anchor_ids;
            for (auto const& [id, _] : dist_map)
                anchor_ids.push_back(id);

            int a1 = -1, a2 = -1, a3 = -1;

            if (anchor_ids.size() >= 3) {
                // Prendo la prima ancora come riferimento
                int a_ref = anchor_ids[0];
                double x1 = x_map[a_ref];
                double y1 = y_map[a_ref];
                double r1 = dist_map[a_ref];

                // Inizializzo i coefficienti della matrice normale
                double ATA11 = 0, ATA12 = 0, ATA22 = 0;
                double ATb1 = 0, ATb2 = 0;

                // Costruisci e accumula (AᵀA) e (Aᵀb)
                for (size_t i = 1; i < anchor_ids.size(); ++i) {
                    int ai = anchor_ids[i];
                    double xi = x_map[ai];
                    double yi = y_map[ai];
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

                // Determinante del sistema normale
                double det = ATA11 * ATA22 - ATA12 * ATA12;

                if (fabs(det) > 1e-12) {
                    x_est = (ATb1 * ATA22 - ATA12 * ATb2) / det;
                    y_est = (ATA11 * ATb2 - ATA12 * ATb1) / det;
                }

                return make_vec(x_est, y_est);
            }
        }
        return make_vec(0,0);
    }
    FUN_EXPORT bis_ksource_t = export_list<bis_ksource_broadcast_t<tuple<vec<2>, double>>>;
}
}

#endif