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

//! @brief Nonlinear least-squares 2D multilateration with the Levenberg–Marquardt method, given an approximated position.
vec<2> multilateration(vec<2> pos, std::vector<tuple<vec<2>, real_t>> anchors) {
    // Handle special cases
    if (anchors.size() == 0) return pos;
    if (anchors.size() == 1) {
        vec<2> diff = pos - get<0>(anchors[0]);
        diff *= get<1>(anchors[0]) / max(norm(diff), 1e-3);
        return get<0>(anchors[0]) + diff;
    }
    double lambda = 1e-3; // normal equation parameter
    for (int iter = 0; iter < 100; ++iter) {
        // Accumulate JᵀJ and Jᵀr
        double H00 = 0, H01 = 0, H11 = 0;
        vec<2> g{0,0};
        double cost = 0;
        for (const auto& a : anchors) {
            vec<2> delta = pos - get<0>(a);
            double r = norm(delta);
            if (r < 1e-8) continue; // avoid singularity
            double ri = r - get<1>(a);
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
        vec<2> pNew = pos + dp;
        // Evaluate new cost
        double newCost = 0;
        for (const auto& a : anchors) {
            vec<2> delta = pNew - get<0>(a);
            double r = norm(delta);
            double ri = r - get<1>(a);
            newCost += ri * ri;
        }
        // Accept or reject step
        if (newCost < cost) {
            pos = pNew;
            lambda *= 0.3;
            if (norm(dp) < 1e-6) break; // tolerance
        } else {
            lambda *= 2.0;
        }
    }
    return pos;
}

} // namespace coordination

} // namespace fcpp

#endif // MULTILATERATION_H_
