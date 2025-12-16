#ifndef FCPP_NON_BAYESIAN_COOPERATIVE_LOOCALIZATION_H_
#define FCPP_NON_BAYESIAN_COOPERATIVE_LOOCALIZATION_H_
#include "lib/coordination.hpp"
#include "lib/data.hpp"

namespace fcpp{
namespace coordination{
    FUN vec<2> nBayesianCoop(ARGS){ CODE
        field<tuple<vec<2>, vec<2>, bool>> map = nbr(CALL, make_tuple(node.position(), make_vec(node.storage(x_stimato{}), node.storage(y_stimato{})), node.storage(is_anchor{})));
        auto hop_map_all = list_hood(CALL,  std::vector< tuple<vec<2>, vec<2>, bool> >{}, map);

    }

}
}

#endif