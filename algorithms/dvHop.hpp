#ifndef FCPP_DV_HOP_H_
#define FCPP_DV_HOP_H_
#include "lib/coordination.hpp"
#include "lib/data.hpp"
#include "trilateration.hpp"

namespace fcpp{
namespace coordination{
    
    using dvhop_state_t = tuple<
            double, // correction
            std::unordered_map<int, real_t,    fcpp::common::hash<int>>, // hop_map
            std::unordered_map<int, double,    fcpp::common::hash<int>>, // correction map
            std::unordered_map<int, double, fcpp::common::hash<int>>, // dist_map
            std::unordered_map<int, int,    fcpp::common::hash<int>>, // x_map
            std::unordered_map<int, int,    fcpp::common::hash<int>>  // y_map
        >;

    FUN vec<2> dvHop(ARGS, int nodeid, bool is_anchor, std::vector<int> my_anchor_keys, field<real_t> nbr_dist, real_t info_speed){ CODE

        double x_est;
        double y_est;

        auto state = old(CALL, dvhop_state_t{0.0, {}, {}, {}, {}, {}}, [&](dvhop_state_t prev){
            auto next = prev;

            auto& correction = get<0>(next);
            auto& hop_map = get<1>(next);
            auto& anchor_correction_map = get<2>(next);
            auto& anchor_distance_map = get<3>(next);
            auto& anchor_x_map = get<4>(next);
            auto& anchor_y_map = get<5>(next);

            auto hop_map_all = spawn(CALL, [&](int nodeid){
                using fcpp::coordination::abf_hops;
                bool is_source = (node.uid == nodeid);
                //int d = abf_hops(CALL, is_source);
                real_t d = bis_distance(CALL, is_anchor, 1, info_speed, [&](){
                    return nbr_dist;
                });
                auto r = broadcast(CALL, d, make_tuple(node.position(), correction));
                return make_tuple(make_tuple(d, get<0>(r), get<1>(r)), true);
            }, my_anchor_keys); 

            for (auto const& [id, t] : hop_map_all) {
                real_t hop = get<0>(t);
                vec<2> pos = get<1>(t);
                double corr = get<2>(t);
                hop_map[id] = hop;
                anchor_x_map[id] = pos[0];
                anchor_y_map[id] = pos[1];
                anchor_correction_map[id] = corr;
                if (is_anchor)
                    anchor_distance_map[id] = distance(node.position(), pos);
                else
                    anchor_distance_map[id] = corr * hop;
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

            return next;
        });

        vec<2> pos = trilaterazione(CALL, is_anchor, get<3>(state), get<4>(state), get<5>(state));
        return pos;

    }
    FUN_EXPORT dvHop_t = export_list<dvhop_state_t, broadcast_t<real_t, tuple<vec<2>, double>>, bis_distance_t, trilaterazione_t, double, int, broadcast_t<int, tuple<vec<2>, double>>>;

}
}

#endif