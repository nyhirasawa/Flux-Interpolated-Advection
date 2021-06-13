#pragma once

#include <vector>
#include "define_float_type.h"
#include "grid.h"
#include "grid_3d.h"
#include "grid_1d.h"

namespace smoke_simulation{
    //////////////////////////////////////////////////
    ////////// 1次元の補間
    //////////////////////////////////////////////////
    MY_FLOAT_TYPE linear_interpolation_1D(
        MY_FLOAT_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_1D &all_grid,
        const MY_FLOAT_TYPE origin_pos,
        const int grid_num_y,
        const MY_FLOAT_TYPE cell_length
    );

    MY_FLOAT_TYPE linear_interpolation_1D_cell_center_values(
        MY_FLOAT_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_1D &all_grid
    );

    //////////////////////////////////////////////////
    ////////// 2次元の補間
    //////////////////////////////////////////////////
    MY_FLOAT_TYPE linear_interpolation_1Dy(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid &all_grid,
        const VEC3_TYPE origin_pos,
        const int grid_num_x,
        const int grid_num_y,
        const MY_FLOAT_TYPE cell_length
    );
    MY_FLOAT_TYPE linear_interpolation_1Dy_cell_center_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid &all_grid
    );
    MY_FLOAT_TYPE linear_interpolation_1Dy_cell_face_x_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid &all_grid
    );
    MY_FLOAT_TYPE linear_interpolation_1Dy_cell_face_y_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_y_values,
        const Grid &all_grid
    );
    MY_FLOAT_TYPE linear_interpolation_1Dy_psi_velocity_x(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid &all_grid
    );
    MY_FLOAT_TYPE linear_interpolation_1Dy_psi_velocity_y(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid &all_grid
    );
    //////////////////////////////////////////////////
    ////////// 3次元の補間
    //////////////////////////////////////////////////
    MY_FLOAT_TYPE linear_interpolation_1Dy(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_3D &all_grid,
        const VEC3_TYPE origin_pos,
        const int grid_num_x,
        const int grid_num_y,
        const MY_FLOAT_TYPE cell_length
    );
    MY_FLOAT_TYPE linear_interpolation_1Dy_cell_center_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_3D &all_grid
    );
    MY_FLOAT_TYPE linear_interpolation_1Dy_cell_face_x_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid_3D &all_grid
    );
    MY_FLOAT_TYPE linear_interpolation_1Dy_cell_face_y_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_y_values,
        const Grid_3D &all_grid
    );
    MY_FLOAT_TYPE linear_interpolation_1Dy_cell_face_z_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_z_values,
        const Grid_3D &all_grid
    );
    MY_FLOAT_TYPE linear_interpolation_1Dy_psi_velocity_x(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid_3D &all_grid
    );
    MY_FLOAT_TYPE linear_interpolation_1Dy_psi_velocity_y(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid_3D &all_grid
    );
    MY_FLOAT_TYPE linear_interpolation_1Dy_psi_velocity_z(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid_3D &all_grid
    );
}//namespace smoke_simulation
