#ifndef MOVING_SUBSTANCES_H
#define MOVING_SUBSTANCES_H
#include <vector>
#include "grid.h"
#include "physical_const.h"
#include "cell_vertex.h"
#include "cell_face.h"
#include "calc_backtraced_face.h"

#include "define_float_type.h"

namespace smoke_simulation {
	//密度場の移流項を semi-Lagrangian で計算
	void advect_density_semi_lagrangian(
		Grid& all_grid,
		const MY_FLOAT_TYPE time_step_length,
		const bool use_MacCormack_scheme,
		const bool use_clamping_in_MacCormack_scheme,
		const int num_gauss_quad_boundary,
		const int num_gauss_quad_bulk,
		const std::string interpolation_method
	);
	//密度場の移流項を Lentine の方法で計算
	void advect_density_Lentine(
		Grid& all_grid,
		const MY_FLOAT_TYPE time_step_length,
		const bool enable_cell_volume_correction
	);
	//密度場の移流項を flux advection で計算
	void advect_density_flux_advection(
		Grid& all_grid,
		const MY_FLOAT_TYPE time_step_length,
//		const bool set_zero_normal_velocity_at_boundary,
		const bool use_MacCormack_scheme,
		const bool use_clamping_in_MacCormack_scheme,
		const int num_gauss_quad_boundary,
		const int num_gauss_quad_bulk,
		const bool enable_eliminate_zero_velocity_false_diffusion,
		const std::string split_method,
		const std::string integral_method,
		std::vector<MY_FLOAT_TYPE> &advected_values,
		const std::string interpolation_method,
	    const std::string interpolation_method_in_calc_psi,
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density,
		const bool enable_cell_volume_correction
	);
	//上の4ステップをまとめただけの関数(substance densityの1時間ステップ分の更新に相当)
	void move_substances(
		Grid& all_grid,
		int i_frame,
		const MY_FLOAT_TYPE time_step_length,
		const bool use_flux_advection,
//		const bool set_zero_normal_velocity_at_boundary,
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
		const bool enable_cell_volume_correction,
		const bool use_lentines_advection
	);
	void correct_mass_and_volume_by_volume_gradient_diffusion(
		Grid& all_grid,
		const int iteration_num
	);
	void correct_mass_and_volume_by_pressure_solve(
		std::vector<MY_FLOAT_TYPE> &corrected_values,
		const std::vector<MY_FLOAT_TYPE> &cell_volume,
		const Grid& all_grid
	);
	// psi から計算した密度が正しいかのチェック
	void check_psi_density(Grid& all_grid);
	void check_negative_dinsity(Grid& all_grid);

	// time_direction == "forward"  の場合は正の時間の方向に backtrace する(普通のsemi-Lagrangian と同じ)
	// time_direction == "backward" の場合は逆の時間の方向に backtrace する( MacCormack の2段階目で使う)
//	cell_face calc_backtraced_face(
//		const Grid& all_grid,
//		const MY_FLOAT_TYPE time_step_length,
//		const cell_face original_face,
//		const std::string time_direction,
//		const bool set_zero_normal_velocity_at_boundary,
//		const std::string interpolation_method
//	);
	// density の MacCormack 法でのclamping の処理
	MY_FLOAT_TYPE clamp_substance_density_of_MacCormack(
		const Grid& all_grid,
		MY_FLOAT_TYPE substance_density_after_advect,
		VEC3_TYPE after_backtrace_position,
		const std::string interpolation_method
	);
	MY_FLOAT_TYPE integrate_normal_component_of_psi_on_face(
		const Grid& all_grid,
		cell_face face,
		const int num_gauss_quad_boundary,
		const int num_gauss_quad_bulk,
		const std::string split_method,
		const std::string interpolation_method,
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density
	);
	VEC3_TYPE analytical_integral_of_psi_on_splitted_face(
		const Grid& all_grid,
		cell_face face
	);
	MY_FLOAT_TYPE integrate_normal_component_of_psi_on_face_analytically(
		const Grid& all_grid,
		cell_face face,
		const std::string split_method
	);
	//引数の面を構成する頂点の位置が、グリッドの範囲を超えていたらclampしてグリッド内に収まるように修正する関数
//	cell_face clamp_vertex_position_by_grid(const Grid& all_grid, cell_face face);


}//namespace smoke_simulation

#endif //MOVING_SUBSTANCES_H
