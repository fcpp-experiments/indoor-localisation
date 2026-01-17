#ifndef MULTILATERATION_H_
#define MULTILATERATION_H_

#include "lib/data/vec.hpp"

namespace fcpp {

namespace coordination {

vec<2> multilateration(std::unordered_map<device_t, vec<2>> anchor_pos, std::unordered_map<device_t, real_t> anchor_dist) {
    // Initial approximation: anchors centroid
    vec<2> p{0,0};
    for (const auto& a : anchor_pos) p += a.second;
    p /= anchor_pos.size();
    double lambda = 1e-3; // normal equation parameter
    
    for (int iter = 0; iter < 100; ++iter) {
        // Accumulate JᵀJ and Jᵀr
        double H00 = 0, H01 = 0, H11 = 0;
        vec<2> g{0,0};
        double cost = 0;
        
        for (const auto& a : anchor_pos) {
            vec<2> delta = p - a.second;
            double r = norm(delta);
            
            if (r < 1e-8) continue; // avoid singularity
            
            double ri = r - anchor_dist[a.first];
            vec<2> J = delta / r;
            
            H00 += J[0] * J[0];
            H01 += J[0] * J[1];
            H11 += J[1] * J[1];
            
            g += J * ri;
            cost += ri * ri;
        }
        
        // Damped Hessian
        H00 += lambda;
        H11 += lambda;
        
        // Solve 2x2 system
        double det = H00 * H11 - H01 * H01;
        if (std::abs(det) < 1e-12)
            break;
        
        vec<2> dp{-H11 * g[0] + H01 * g[1], H01 * g[0] - H00 * g[1]};
        dp /= det;
        
        vec<2> pNew = p + dp;
        
        // Evaluate new cost
        double newCost = 0;
        for (const auto& a : anchor_pos) {
            vec<2> delta = pNew - a.second;
            double r = norm(delta);
            double ri = r - anchor_dist[a.first];
            newCost += ri * ri;
        }
        
        // Accept or reject step
        if (newCost < cost) {
            p = pNew;
            lambda *= 0.3;
            if (norm(dp) < 1e-6) break; // tolerance
        } else {
            lambda *= 2.0;
        }
    }
    return p;
}

} // namespace coordination

} // namespace fcpp

#endif
