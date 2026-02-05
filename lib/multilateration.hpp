#ifndef MULTILATERATION_H_
#define MULTILATERATION_H_

#include <vector>

#include "lib/data/vec.hpp"

/**
 * @brief Namespace containing all the objects in the FCPP library.
 */
namespace fcpp {

//! @brief Namespace containing the libraries of coordination routines.
namespace coordination {

//! @brief Weighted nonlinear least-squares 2D multilateration with the Levenberg–Marquardt method, given an approximated position and a base weight for anchors.
vec<2> multilateration(vec<2> pos, std::vector<tuple<vec<2>, real_t>> const& anchors, std::vector<real_t> const& weights, real_t base_weight, real_t& weight) {
    // Handle special cases
    if (anchors.size() == 0) return pos;
    if (anchors.size() == 1) {
        vec<2> diff = pos - get<0>(anchors[0]);
        real_t len = norm(diff);
        if (len < 1e-8) return pos;
        diff *= get<1>(anchors[0]) / len;
        return get<0>(anchors[0]) + diff;
    }
    real_t lambda = 1e-3; // normal equation parameter
    for (int iter = 0; iter < 100; ++iter) {
        // Accumulate JᵀJ and Jᵀr
        real_t H00 = 0, H01 = 0, H11 = 0;
        vec<2> g{0,0};
        real_t cost = 0;
        for (int i=0; i<anchors.size(); ++i) {
            vec<2> delta = pos - get<0>(anchors[i]);
            real_t r = norm(delta);
            if (r < 1e-8) continue; // avoid singularity
            real_t ri = r - get<1>(anchors[i]);
            vec<2> J = delta / r;
            if (weights.size()) {
                ri *= weights[i];
                J *= weights[i];
            }
            H00 += J[0] * J[0];
            H01 += J[0] * J[1];
            H11 += J[1] * J[1];
            g += J * ri;
            cost += ri * ri;
        }
        // Estimated variance
        real_t vNew;
        if (weights.size()) {
            real_t det = H00 * H11 - H01 * H01;
            vNew = (H11+H00) / det;
        }
        // Damped Hessian
        H00 += lambda;
        H11 += lambda;
        // Solve 2x2 system
        real_t det = H00 * H11 - H01 * H01;
        if (std::abs(det) < 1e-12)
            break;
        vec<2> dp{-H11 * g[0] + H01 * g[1], H01 * g[0] - H00 * g[1]};
        dp /= det;
        vec<2> pNew = pos + dp;
        // Evaluate new cost
        real_t newCost = 0;
        for (int i=0; i<anchors.size(); ++i) {
            vec<2> delta = pNew - get<0>(anchors[i]);
            real_t r = norm(delta);
            real_t ri = r - get<1>(anchors[i]);
            if (weights.size()) ri *= weights[i];
            newCost += ri * ri;
        }
        // Accept or reject step
        if (newCost < cost) {
            if (weights.size()) weight = 1 / std::sqrt(1 / (base_weight*base_weight) + vNew);
            pos = pNew;
            lambda *= 0.3;
            if (norm(dp) < 1e-6) break; // tolerance
        } else {
            lambda *= 2.0;
        }
    }
    return pos;
}
//! @brief Nonlinear least-squares 2D multilateration with the Levenberg–Marquardt method, given an approximated position.
inline vec<2> multilateration(vec<2> pos, std::vector<tuple<vec<2>, real_t>> anchors) {
    real_t weight;
    return multilateration(pos, anchors, {}, 0, weight);
}

} // namespace coordination

} // namespace fcpp

#endif // MULTILATERATION_H_
