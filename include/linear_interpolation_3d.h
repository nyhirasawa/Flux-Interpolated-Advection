#pragma once

#include <vector>
#include "define_float_type.h"
#include "grid_3d.h"
namespace smoke_simulation{
    MY_FLOAT_TYPE linear_interpolation_3D(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_3D &all_grid,
        const VEC3_TYPE origin_pos,
        const int grid_num_x,
        const int grid_num_y,
        const int grid_num_z,
        const MY_FLOAT_TYPE cell_length
    );

    MY_FLOAT_TYPE linear_interpolation_3D_cell_center_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_3D &all_grid
    );

    MY_FLOAT_TYPE linear_interpolation_3D_cell_face_values(
        const int dim,
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid_3D &all_grid
    );

    MY_FLOAT_TYPE linear_interpolation_3D_cell_face_x_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid_3D &all_grid
    );

    MY_FLOAT_TYPE linear_interpolation_3D_cell_face_y_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_y_values,
        const Grid_3D &all_grid
    );

    MY_FLOAT_TYPE linear_interpolation_3D_cell_face_z_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_z_values,
        const Grid_3D &all_grid
    );

    MY_FLOAT_TYPE linear_interpolation_3D_psi_velocity_x(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_3D &all_grid
    );

    MY_FLOAT_TYPE linear_interpolation_3D_psi_velocity_y(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_3D &all_grid
    );

    MY_FLOAT_TYPE linear_interpolation_3D_psi_velocity_z(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_3D &all_grid
    );
}//namespace smoke_simulation
