#ifndef UPDATE_FLUID_VELOCITY_3D_H
#define UPDATE_FLUID_VELOCITY_3D_H
#include <vector>

#include "grid_3d.h"
#include "physical_const.h"
#include "gauss_quadrature_points.h"

#include "define_float_type.h"

namespace smoke_simulation{

//advect項の計算
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
);

//Poisson eq. を解くことによる圧力の計算
//void calc_pressure_3D(Grid_3D& all_grid);
//pressure gradient termの計算
//void calc_pressure_gradient_term_3D(Grid_3D& all_grid);
//上記4ステップをまとめただけの関数
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
);

}

#endif //UPDATE_FLUID_VELOCITY_H
