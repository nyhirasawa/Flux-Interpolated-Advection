#pragma once

#include <vector>
#include "define_float_type.h"
#include "grid.h"
namespace smoke_simulation{
    MY_FLOAT_TYPE linear_interpolation_2D(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid &all_grid,
        const VEC3_TYPE origin_pos,
        const int grid_num_x,
        const int grid_num_y,
        const MY_FLOAT_TYPE cell_length
    );
    MY_FLOAT_TYPE linear_interpolation_2D_cell_center_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid &all_grid
    );
    MY_FLOAT_TYPE linear_interpolation_2D_cell_face_x_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid &all_grid
    );
    MY_FLOAT_TYPE linear_interpolation_2D_cell_face_y_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_y_values,
        const Grid &all_grid
    );
    MY_FLOAT_TYPE linear_interpolation_2D_psi_velocity_x(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid &all_grid
    );

    MY_FLOAT_TYPE linear_interpolation_2D_psi_velocity_y(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid &all_grid
    );

    // psi_density を積分で計算するとき用
    MY_FLOAT_TYPE linear_interpolation_2D_cell_face_y_values_use_integral(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &psi_cell_face_y,
        const std::vector<MY_FLOAT_TYPE> &density_cell_center,
        const Grid &all_grid,
        const int num_gauss_quadrature_point_for_integrate_density
    );

}//namespace smoke_simulation
