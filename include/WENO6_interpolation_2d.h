#pragma once

#include <vector>
#include "define_float_type.h"
#include "grid.h"
namespace smoke_simulation{

MY_FLOAT_TYPE WENO6_interpolation_2D_cell_center_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_center_values,
    const Grid &all_grid
);

MY_FLOAT_TYPE WENO6_interpolation_2D_cell_face_x_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
    const Grid &all_grid
);

MY_FLOAT_TYPE WENO6_interpolation_2D_cell_face_y_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_y_values,
    const Grid &all_grid
);

}//namespace smoke_simulation
