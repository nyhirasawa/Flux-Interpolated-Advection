#pragma once

#include "grid_3d.h"

#include "gauss_quadrature_points.h"
#include "parallelize_functions.h"

namespace smoke_simulation{
//--------------------------------------------------------------------------
// 移流関連
//--------------------------------------------------------------------------
    //advect項の計算
    //速度場を時間 -dt だけバックトレースしてadvect項を計算する
    void advect_velocity_flux_advection_3D_std_thread(
        const int dim,
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const bool use_MacCormack_scheme,
        const bool use_clamping_in_MacCormack_scheme,
        const int num_gauss_quad_boundary,
        const int num_gauss_quad_bulk,
    //		const bool enable_eliminate_zero_velocity_false_diffusion,
        const std::string split_method_velocity,
        const std::string integral_method,
        const int num_gauss_quadrature_points_on_one_triangle,
        std::vector<MY_FLOAT_TYPE> &advected_values,
        const std::string interpolation_method,
        const std::string interpolation_method_in_calc_psi,
        const bool use_zero_velocity_for_backtrace,
        const bool use_integral,
        const int num_gauss_quadrature_point_for_integrate_density,
        const bool enable_cell_volume_correction,
        const MY_FLOAT_TYPE minimal_area
    );

    //advect項の計算
    //速度場を時間 -dt だけバックトレースしてadvect項を計算する
    void advect_velocity_x_flux_advection_3D_std_thread(
        const std::vector<MY_FLOAT_TYPE> &before_advect_values,
        std::vector<MY_FLOAT_TYPE> &after_advect_values,
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const bool use_MacCormack_scheme,
        const bool use_clamping_in_MacCormack_scheme,
        const int num_gauss_quad_boundary,
        const int num_gauss_quad_bulk,
        const bool enable_eliminate_zero_velocity_false_diffusion,
    //    const std::string split_method,
        const std::string integral_method,
        const int num_gauss_quadrature_points_on_one_triangle,
        const std::string interpolation_method,
    //    const std::string interpolation_method_in_calc_psi,
        bool use_integral,
        const int num_gauss_quadrature_point_for_integrate_density,
        const bool enable_cell_volume_correction,
        const MY_FLOAT_TYPE minimal_area
    );

    //advect項の計算
    //速度場を時間 -dt だけバックトレースしてadvect項を計算する
    void advect_velocity_y_flux_advection_3D_std_thread(
        const std::vector<MY_FLOAT_TYPE> &before_advect_values,
        std::vector<MY_FLOAT_TYPE> &after_advect_values,
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const bool use_MacCormack_scheme,
        const bool use_clamping_in_MacCormack_scheme,
        const int num_gauss_quad_boundary,
        const int num_gauss_quad_bulk,
        const bool enable_eliminate_zero_velocity_false_diffusion,
    //    const std::string split_method,
        const std::string integral_method,
        const int num_gauss_quadrature_points_on_one_triangle,
        const std::string interpolation_method,
    //    const std::string interpolation_method_in_calc_psi,
        bool use_integral,
        const int num_gauss_quadrature_point_for_integrate_density,
        const bool enable_cell_volume_correction,
        const MY_FLOAT_TYPE minimal_area
    );

    //advect項の計算
    //速度場を時間 -dt だけバックトレースしてadvect項を計算する
    void advect_velocity_z_flux_advection_3D_std_thread(
        const std::vector<MY_FLOAT_TYPE> &before_advect_values,
        std::vector<MY_FLOAT_TYPE> &after_advect_values,
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const bool use_MacCormack_scheme,
        const bool use_clamping_in_MacCormack_scheme,
        const int num_gauss_quad_boundary,
        const int num_gauss_quad_bulk,
        const bool enable_eliminate_zero_velocity_false_diffusion,
    //    const std::string split_method,
        const std::string integral_method,
        const int num_gauss_quadrature_points_on_one_triangle,
        const std::string interpolation_method,
    //    const std::string interpolation_method_in_calc_psi,
        bool use_integral,
        const int num_gauss_quadrature_point_for_integrate_density,
        const bool enable_cell_volume_correction,
        const MY_FLOAT_TYPE minimal_area
    );

    //advect項の計算
    //速度場を時間 -dt だけバックトレースしてadvect項を計算する
    void advect_velocity_xyz_flux_advection_3D_std_thread(
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const bool use_MacCormack_scheme,
        const bool use_clamping_in_MacCormack_scheme,
        const int num_gauss_quad_boundary,
        const int num_gauss_quad_bulk,
        const bool enable_eliminate_zero_velocity_false_diffusion,
    //    const std::string split_method,
        const std::string integral_method,
        const int num_gauss_quadrature_points_on_one_triangle,
        const std::string interpolation_method,
    //    const std::string interpolation_method_in_calc_psi,
        bool use_integral,
        const int num_gauss_quadrature_point_for_integrate_density,
        const bool enable_cell_volume_correction,
        const MY_FLOAT_TYPE minimal_area
    );

//--------------------------------------------------------------------------
// 並列化関連
//--------------------------------------------------------------------------
    //calc_density_in_cell_from_quadrature_point_list()をstd::thread で並列に実行する関数
    void parallelize_calc_velocity_from_quadrature_point_list (
        auto func,
        const quadrature_point_vector &quadrature_point_list,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &velocity_after_advect,
        const std::string interpolation_method
    );

    //calc_density_in_cell_from_quadrature_point_list()をstd::thread で並列に実行する関数
    void parallelize_calc_velocity_from_quadrature_point_list (
        auto func,
        const int dim,
        const quadrature_point_vector &quadrature_point_list,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &velocity_after_advect,
        const std::string interpolation_method
    );

    //backtrace_and_split_cell_faces()をstd::thread で並列に実行する関数
    void parallelize_backtrace_and_calc_all_quadrature_point_list_for_integrate_velocity(
        auto func,
        size_t size,
        const int dim,
        quadrature_point_vector &quadrature_point_list_for_interpolate_density,
        quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
        quadrature_point_vector &quadrature_point_list_for_calc_cell_volume,
        bool use_MacCormack_scheme,
        const Grid_3D &all_grid,
        MY_FLOAT_TYPE time_step_length,
//        bool set_zero_normal_velocity_at_boundary,
        std::string split_method,
//        bool is_zero_velocity_false_diffusion_correction_mode,
        const std::string integral_method,
        const int num_gauss_quadrature_points,
        const std::string interpolation_method,
        const bool use_zero_velocity_for_backtrace,
        const int num_gauss_quadrature_point_for_integrate_density,
        const bool enable_cell_volume_correction,
        const std::string backtrace_usage,
        const MY_FLOAT_TYPE minimal_area
    );

    //calc_velocity_from_quadrature_point_list_use_integral()をstd::thread で並列に実行する関数
    void parallelize_calc_velocity_from_quadrature_point_list_use_integral (
        auto func,
        const int dim,
        const quadrature_point_vector &quadrature_point_list_for_interpolate_density,
        const quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &substance_density_after_advect,
        const std::string interpolation_method,
        const std::vector<MY_FLOAT_TYPE> &advected_values,
        const MY_FLOAT_TYPE minimal_area
    );

} // namespace smoke_simulation
