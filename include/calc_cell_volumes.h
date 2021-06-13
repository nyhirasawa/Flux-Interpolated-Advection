#pragma once

#include "grid.h"

namespace smoke_simulation{
    void calc_all_cell_volumes(
        Grid &all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const bool use_MacCormack_scheme,
        const std::string interpolation_method_for_velocity
    );

    void calc_cell_volume_contribution(
        const cell_face &face,
        MY_FLOAT_TYPE &cell_volume
    );
}
