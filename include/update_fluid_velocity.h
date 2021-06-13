#ifndef UPDATE_FLUID_VELOCITY_H
#define UPDATE_FLUID_VELOCITY_H
#include <vector>

#include "grid.h"
#include "physical_const.h"

#include "define_float_type.h"

namespace smoke_simulation {
	void advect_velocity(
		Grid& all_grid,
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
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density,
		const bool enable_cell_volume_correction
	);

	//advect項の計算
	//速度場を時間 -dt だけバックトレースしてadvect項を計算する
	void advect_velocity_flux_advection(
		Grid& all_grid,
		const MY_FLOAT_TYPE time_step_length,
		const bool use_MacCormack_scheme,
		const bool use_clamping_in_MacCormack_scheme,
		const int num_gauss_quad_boundary,
		const int num_gauss_quad_bulk,
		const bool enable_eliminate_zero_velocity_false_diffusion,
		const std::string split_method,
		const std::string integral_method,
		const std::string interpolation_method,
		const std::string interpolation_method_in_calc_psi,
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density,
		const bool enable_cell_volume_correction
	);

	//advect項の計算
	//速度場を時間 -dt だけバックトレースしてadvect項を計算する
	void advect_velocity_flux_advection_x(
		Grid& all_grid,
		const MY_FLOAT_TYPE time_step_length,
		const bool use_MacCormack_scheme,
		const bool use_clamping_in_MacCormack_scheme,
		const int num_gauss_quad_boundary,
		const int num_gauss_quad_bulk,
//		const bool enable_eliminate_zero_velocity_false_diffusion,
		const std::string split_method,
		const std::string integral_method,
		std::vector<MY_FLOAT_TYPE> &advected_values,
		const std::string interpolation_method,
		const std::string interpolation_method_in_calc_psi,
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density,
		const bool enable_cell_volume_correction
//		std::vector<MY_FLOAT_TYPE> &velocity_after_advect_x
	);

	//advect項の計算
	//速度場を時間 -dt だけバックトレースしてadvect項を計算する
	void advect_velocity_flux_advection_y(
		Grid& all_grid,
		const MY_FLOAT_TYPE time_step_length,
		const bool use_MacCormack_scheme,
		const bool use_clamping_in_MacCormack_scheme,
		const int num_gauss_quad_boundary,
		const int num_gauss_quad_bulk,
//		const bool enable_eliminate_zero_velocity_false_diffusion,
		const std::string split_method,
		const std::string integral_method,
		std::vector<MY_FLOAT_TYPE> &advected_values,
		const std::string interpolation_method,
		const std::string interpolation_method_in_calc_psi,
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density,
		const bool enable_cell_volume_correction
//		std::vector<MY_FLOAT_TYPE> &velocity_after_advect_x
	);

	void correct_velocity_x_and_volume_by_pressure_solve(
		Grid& all_grid
	);
	void correct_velocity_y_and_volume_by_pressure_solve(
		Grid& all_grid
	);

	//advect項の計算
	void advect_velocity_semi_lagrangian(
		Grid& all_grid,
		const MY_FLOAT_TYPE time_step_length,
		const std::string interpolation_method
	);

	//advect項の計算
	//速度場を時間 -dt だけバックトレースしてadvect項を計算する
	void advect_velocity_semi_lagrangian_x(
		Grid& all_grid,
		const MY_FLOAT_TYPE time_step_length,
		const std::string interpolation_method,
		std::vector<MY_FLOAT_TYPE> &velocity_after_advect_x
	);

	//advect項の計算
	//速度場を時間 -dt だけバックトレースしてadvect項を計算する
	void advect_velocity_semi_lagrangian_y(
		Grid& all_grid,
		const MY_FLOAT_TYPE time_step_length,
		const std::string interpolation_method,
		std::vector<MY_FLOAT_TYPE> &velocity_after_advect_y
	);

	//上記4ステップをまとめただけの関数
	void update_fluid_velocity(
		Grid& all_grid,
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
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density
	);

	// velocity の数値がnanを含むかをチェック
	void check_nan_in_velocity(const Grid& all_grid);

	MY_FLOAT_TYPE integrate_normal_component_of_psi_velocity_x(
		const Grid& all_grid,
		cell_face face,
		const int num_gauss_quad_boundary,
		const int num_gauss_quad_bulk,
		const std::string split_method,
		const std::string interpolation_method,
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density
	);

	MY_FLOAT_TYPE integrate_normal_component_of_psi_velocity_y(
		const Grid& all_grid,
		cell_face face,
		const int num_gauss_quad_boundary,
		const int num_gauss_quad_bulk,
		const std::string split_method,
		const std::string interpolation_method,
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density
	);

	// density の MacCormack 法でのclamping の処理
	MY_FLOAT_TYPE clamp_integral_of_normal_component_of_psi_velocity_x_of_MacCormack(
		const Grid& all_grid,
		MY_FLOAT_TYPE integral_of_normal_component_of_psi_after_advect,
		VEC3_TYPE after_backtrace_position,
		const std::string interpolation_method
	);

	// density の MacCormack 法でのclamping の処理
	MY_FLOAT_TYPE clamp_integral_of_normal_component_of_psi_velocity_y_of_MacCormack(
		const Grid& all_grid,
		MY_FLOAT_TYPE integral_of_normal_component_of_psi_after_advect,
		VEC3_TYPE after_backtrace_position,
		const std::string interpolation_method
	);
}

#endif //UPDATE_FLUID_VELOCITY_H
