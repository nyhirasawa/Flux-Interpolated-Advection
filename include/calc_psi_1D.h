#pragma once

#include "grid_1d.h"
#include "define_float_type.h"

namespace smoke_simulation {
    void calc_psi_on_cell_face_from_density_on_cell_center_1D(
        const Grid_1D &all_grid,
        std::vector<MY_FLOAT_TYPE> &psi_on_cell_face_y,
        const std::vector<MY_FLOAT_TYPE> &density_on_cell_center,
        const int Grid_num_y,
        const MY_FLOAT_TYPE cell_length,
        const std::string interpolation_method
    );
}
