#ifndef FCPP_ANCHOR_LAYOUT_H_
#define FCPP_ANCHOR_LAYOUT_H_
#include "lib/coordination.hpp"
#include "lib/data.hpp"

namespace fcpp{
namespace coordination{
    
    enum AnchorLayout {
        PERIMETRO = 1,
        GRIGLIA   = 2
    };

    FUN vec<2> positionAnchor(ARGS, int id, AnchorLayout layout, double side, double step, int rows, int cols){CODE
        if (layout == PERIMETRO) {
        int anchors_per_side = static_cast<int>(side / step);
        int pos = id;

        double x = 0.0, y = 0.0;

        if (pos <= anchors_per_side) {
            x = pos * step;
            y = side;
        }
        else if (pos <= anchors_per_side * 2) {
            x = side;
            y = side - (pos - anchors_per_side) * step;
        }
        else if (pos <= anchors_per_side * 3) {
            x = side - (pos - anchors_per_side * 2) * step;
            y = 0;
        }
        else {
            x = 0;
            y = (pos - anchors_per_side * 3) * step;
        }

        return make_vec(x, y);
        }
        else {
            double step_x = side / (cols - 1);
            double step_y = side / (rows - 1);

            int row = id / cols;
            int col = id % cols;

            return make_vec(col * step_x, row * step_y);
        }
    }
    FUN_EXPORT positionAnchor_t = export_list<>;
}
}

#endif