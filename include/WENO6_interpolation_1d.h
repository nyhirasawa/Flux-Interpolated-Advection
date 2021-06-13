#pragma once

#include <vector>
#include "define_float_type.h"
#include "grid.h"
#include "grid_3d.h"

namespace smoke_simulation{
    MY_FLOAT_TYPE WENO6_interpolation_1Dy_cell_center_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid &all_grid
    );

    MY_FLOAT_TYPE WENO6_interpolation_1Dy_cell_face_x_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid &all_grid
    );

    MY_FLOAT_TYPE WENO6_interpolation_1Dy_cell_face_y_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_y_values,
        const Grid &all_grid
    );

    MY_FLOAT_TYPE WENO6_interpolation_1Dy_cell_center_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_3D &all_grid
    );

    MY_FLOAT_TYPE WENO6_interpolation_1Dy_cell_face_x_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid_3D &all_grid
    );

    MY_FLOAT_TYPE WENO6_interpolation_1Dy_cell_face_y_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_y_values,
        const Grid_3D &all_grid
    );

    MY_FLOAT_TYPE WENO6_interpolation_1Dy_cell_face_z_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_z_values,
        const Grid_3D &all_grid
    );
}//namespace smoke_simulation
