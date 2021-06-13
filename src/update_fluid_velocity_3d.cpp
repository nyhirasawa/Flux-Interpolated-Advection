#include "update_fluid_velocity.h"

#include <math.h>
#include <iostream>
#include <chrono>//時間計測用
#include <fstream>//ファイル書き出し用
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

#include "grid_3d.h"
#include "linear_solver.h"
#include "physical_const.h"
#include "sparse_matrix.h"
#include "utils.h"
#include "define_float_type.h"
#include "linear_interpolation_3d.h"
#include "calc_psi_3D.h"
#include "gauss_quadrature_points.h"
#include "backtrace_and_calc_all_quadrature_point_list.h"
#include "backtrace_and_calc_all_quadrature_point_list_for_integrate_density.h"
#include "move_substances_3d.h"
#include "calc_density_in_cell_from_quadrature_point_list.h"
#include "calc_density_in_cell_from_quadrature_point_list_for_integrate_density.h"
#include "calc_velocity_from_quadrature_point_list.h"
#include "calc_velocity_from_quadrature_point_list_use_integral.h"
#include "advect_velocity_semi_lagrangian_3d.h"
#include "calc_pressure_3d.h"
#include "advect_velocity_flux_advection_3d.h"

namespace smoke_simulation{
//advect項の計算
//速度場を時間 -dt だけバックトレースしてadvect項を計算する
void advect_fluid_3D(
    Grid_3D& all_grid,
    int i_frame,
    const MY_FLOAT_TYPE time_step_length,
    const bool use_flux_advection,
    const bool use_MacCormack_scheme,
    const bool use_clamping_in_MacCormack_scheme,
    const int num_gauss_quad_boundary,
    const int num_gauss_quad_bulk,
    const bool enable_eliminate_zero_velocity_false_diffusion,
    const std::string split_method,
    const std::string integral_method,
    const std::string interpolation_method,
    const int num_gauss_quadrature_point_for_integrate_density,
    const bool enable_cell_volume_correction,
    const bool use_integral,
    const MY_FLOAT_TYPE minimal_area
){
    //移流項の計算
    if (use_flux_advection) {
        advect_velocity_xyz_flux_advection_3D_std_thread(
            all_grid,
            time_step_length,
            use_MacCormack_scheme,
            use_clamping_in_MacCormack_scheme,
            num_gauss_quad_boundary,
            num_gauss_quad_bulk,
            enable_eliminate_zero_velocity_false_diffusion,
//            split_method,
            integral_method,
            /*num_gauss_quadrature_points_on_one_triangle = */num_gauss_quad_bulk,
            interpolation_method,
//            interpolation_method_in_calc_psi,
            use_integral,
            num_gauss_quadrature_point_for_integrate_density,
            enable_cell_volume_correction,
            minimal_area
        );
    }
    else {
        advect_velocity_semi_lagrangian_3D(
            all_grid,
            time_step_length,
            interpolation_method
        );
    }
}

//流体の 1 time step
void update_fluid_velocity_3D(
    Grid_3D& all_grid,
    int i_frame,
    const MY_FLOAT_TYPE time_step_length,
    const bool use_flux_advection,
    const bool use_MacCormack_scheme,
    const bool use_clamping_in_MacCormack_scheme,
    const int num_gauss_quad_boundary,
    const int num_gauss_quad_bulk,
    const bool enable_eliminate_zero_velocity_false_diffusion,
    const std::string split_method,
    const std::string integral_method,
    const std::string interpolation_method,
    const bool enable_cell_volume_correction,
    const bool fix_velocity,
    const int num_gauss_quadrature_point_for_integrate_density,
    const bool use_integral,
    const MY_FLOAT_TYPE minimal_area
){
    ////時間計測開始
    auto start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
//    add_force_fluid_3D(all_grid, time_step_length);
    ////時間計測終了
    auto end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
    auto duration_time = end_time - start_time;        // 要した時間を計算
    auto duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
    // 要した時間をミリ秒（1/1000秒）に変換して表示
//    std::cout <<"add_force_fluid_3D: "<< duration_time_msec << " milli sec \n";
    if (i_frame == 0) {
        ////時間計測開始
        start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
        calc_pressure_gradient_term_3D(all_grid);
        ////時間計測終了
        end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
        duration_time = end_time - start_time;        // 要した時間を計算
        duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
        // 要した時間をミリ秒（1/1000秒）に変換して表示
        std::cout <<"calc_pressure_gradient_term_3D: "<< duration_time_msec << " milli sec \n";

        ////時間計測開始
        start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
        advect_fluid_3D(
            all_grid,
            i_frame,
            time_step_length,
            use_flux_advection,
            use_MacCormack_scheme,
            use_clamping_in_MacCormack_scheme,
            num_gauss_quad_boundary,
            num_gauss_quad_bulk,
            enable_eliminate_zero_velocity_false_diffusion,
            split_method,
            integral_method,
            interpolation_method,
            num_gauss_quadrature_point_for_integrate_density,
            enable_cell_volume_correction,
            use_integral,
            minimal_area
        );

        ////時間計測終了
        end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
        duration_time = end_time - start_time;        // 要した時間を計算
        duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
        // 要した時間をミリ秒（1/1000秒）に変換して表示
        std::cout <<"advect_fluid_3D: "<< duration_time_msec << " milli sec \n";
    }
    else {
        ////時間計測開始
        start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
        advect_fluid_3D(
            all_grid,
            i_frame,
            time_step_length,
            use_flux_advection,
            use_MacCormack_scheme,
            use_clamping_in_MacCormack_scheme,
            num_gauss_quad_boundary,
            num_gauss_quad_bulk,
            enable_eliminate_zero_velocity_false_diffusion,
            split_method,
            integral_method,
            interpolation_method,
            num_gauss_quadrature_point_for_integrate_density,
            enable_cell_volume_correction,
            use_integral,
            minimal_area
        );
        ////時間計測終了
        end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
        duration_time = end_time - start_time;        // 要した時間を計算
        duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
        // 要した時間をミリ秒（1/1000秒）に変換して表示
        std::cout <<"advect_fluid_3D: "<< duration_time_msec << " milli sec \n";

        ////時間計測開始
        start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
        calc_pressure_gradient_term_3D(all_grid);
        ////時間計測終了
        end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
        duration_time = end_time - start_time;        // 要した時間を計算
        duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
        // 要した時間をミリ秒（1/1000秒）に変換して表示
        std::cout <<"calc_pressure_gradient_term_3D: "<< duration_time_msec << " milli sec \n";
    }
}

}
