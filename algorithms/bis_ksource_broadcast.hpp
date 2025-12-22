#ifndef FCPP_BIS_KSOURCE_BROADCAST_H_
#define FCPP_BIS_KSOURCE_BROADCAST_H_
#include "lib/coordination.hpp"
#include "lib/data.hpp"
#include "trilateration.hpp"

namespace fcpp{
namespace coordination{
    FUN vec<2> bis_ksource(ARGS, bool is_anchor, field<real_t> nbr_dist, real_t info_speed){CODE
        std::unordered_map<int, real_t, fcpp::common::hash<int>> hop_map;
        std::unordered_map<int,int, fcpp::common::hash<int>> anchor_x_map;
        std::unordered_map<int,int, fcpp::common::hash<int>> anchor_y_map;
        std::unordered_map<int,int, fcpp::common::hash<int>> anchor_correction_map;
        std::unordered_map<int,double, fcpp::common::hash<int>> anchor_distance_map;

        double x_est;
        double y_est;

        double correction = old(CALL, 1, [&](double correction){
            auto hop_map_all = bis_ksource_broadcast(CALL, is_anchor, make_tuple(node.position(), correction), 5, 1, info_speed, [&](){
                return nbr_dist;
            });

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
                    anchor_distance_map[id] = distance(node.position(), pos);
                else
                    anchor_distance_map[id] = dist * corr;
            }

            real_t hop = 0;
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

        vec<2> pos = trilaterazione(CALL, is_anchor, anchor_distance_map, anchor_x_map, anchor_y_map);
        return pos;

    }
    FUN_EXPORT bis_ksource_t = export_list<bis_ksource_broadcast_t<tuple<vec<2>, double>>, trilaterazione_t>;
}
}

#endif