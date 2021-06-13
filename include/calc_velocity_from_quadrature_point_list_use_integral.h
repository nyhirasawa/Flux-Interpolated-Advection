#pragma once

#include "gauss_quadrature_points.h"
#include "grid_3d.h"

namespace smoke_simulation{
    // i_thread 番目のスレッドで計算する求積点のリスト quadrature_point_list_of_i_thread から 移流後の質量密度場 velocity_x_after_advectへの寄与を計算する
    void calc_velocity_from_quadrature_point_list_use_integral(
        const int dim,
        const quadrature_point_vector &quadrature_point_list_for_interpolate_velocity,
        const quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &velocity_after_advect,
        const int i_thread,
        const std::string interpolation_method,
        const std::vector<MY_FLOAT_TYPE> &advected_values
    );

    // i_thread 番目のスレッドで計算する求積点のリスト quadrature_point_list_of_i_thread から 移流後の質量密度場 velocity_x_after_advectへの寄与を計算する
    void calc_velocity_x_from_quadrature_point_list_use_integral(
        const quadrature_point_vector &quadrature_point_list_for_interpolate_velocity_x,
        const quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &velocity_x_after_advect,
        const int i_thread,
        const std::string interpolation_method,
        const std::vector<MY_FLOAT_TYPE> &advected_values
    );

    // i_thread 番目のスレッドで計算する求積点のリスト quadrature_point_list_of_i_thread から 移流後の質量密度場 velocity_x_after_advectへの寄与を計算する
    void calc_velocity_y_from_quadrature_point_list_use_integral(
        const quadrature_point_vector &quadrature_point_list_for_interpolate_velocity_y,
        const quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &velocity_y_after_advect,
        const int i_thread,
        const std::string interpolation_method,
        const std::vector<MY_FLOAT_TYPE> &advected_values
    );

    // i_thread 番目のスレッドで計算する求積点のリスト quadrature_point_list_of_i_thread から 移流後の質量密度場 velocity_x_after_advectへの寄与を計算する
    void calc_velocity_z_from_quadrature_point_list_use_integral(
        const quadrature_point_vector &quadrature_point_list_for_interpolate_velocity_z,
        const quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &velocity_z_after_advect,
        const int i_thread,
        const std::string interpolation_method,
        const std::vector<MY_FLOAT_TYPE> &advected_values
    );
} //namespace smoke_simulation
