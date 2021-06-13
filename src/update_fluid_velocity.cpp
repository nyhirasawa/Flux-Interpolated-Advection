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

#include "linear_solver.h"
#include "calc_pressure.h"
#include "physical_const.h"
#include "sparse_matrix.h"
#include "utils.h"
#include "define_float_type.h"
#include "calc_cell_volumes.h"
#include "linear_interpolation_2d.h"
#include "calc_psi.h"
#include "calc_backtraced_face.h"
#include "split_face.h"
#include "gauss_quadrature_points.h"

namespace smoke_simulation {
	//advect項の計算
	//速度場を時間 -dt だけバックトレースしてadvect項を計算する
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
	) {
		//移流項の計算
		if (use_flux_advection) {
			std::string interpolation_method_in_calc_psi;
			if(interpolation_method =="linear"){
				interpolation_method_in_calc_psi = "linear";
			}
			else if(interpolation_method=="1Dy_linear"){
	//            interpolation_method_in_calc_psi = "linear";
	//            interpolation_method_in_calc_psi = "1Dy_linear";
				interpolation_method_in_calc_psi = "const";
			}
			else{
				interpolation_method_in_calc_psi = "const";
			}
			advect_velocity_flux_advection(
				all_grid,
				time_step_length,
				use_MacCormack_scheme,
				use_clamping_in_MacCormack_scheme,
				num_gauss_quad_boundary,
				num_gauss_quad_bulk,
				enable_eliminate_zero_velocity_false_diffusion,
				split_method,
				integral_method,
				interpolation_method,
				interpolation_method_in_calc_psi,
				use_integral,
				num_gauss_quadrature_point_for_integrate_density,
				enable_cell_volume_correction
			);
		}
		else {
			advect_velocity_semi_lagrangian(
				all_grid,
				time_step_length,
				interpolation_method
			);
		}
	}

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
	) {
		// 移流前の速度場
		std::vector<MY_FLOAT_TYPE> before_advect_velocity_x((all_grid.Grid_num_x + 1) * all_grid.Grid_num_y);
		std::vector<MY_FLOAT_TYPE> before_advect_velocity_y(all_grid.Grid_num_x * (all_grid.Grid_num_y + 1));
		before_advect_velocity_x = all_grid.velocity_in_voxel_face_x;
		before_advect_velocity_y = all_grid.velocity_in_voxel_face_y;
		// 移流後の速度場
		std::vector<MY_FLOAT_TYPE> velocity_after_advect_x((all_grid.Grid_num_x + 1) * all_grid.Grid_num_y);
		std::vector<MY_FLOAT_TYPE> velocity_after_advect_y(all_grid.Grid_num_x * (all_grid.Grid_num_y + 1));
		velocity_after_advect_x = all_grid.velocity_in_voxel_face_x;
		velocity_after_advect_y = all_grid.velocity_in_voxel_face_y;

		//x成分の移流
		advect_velocity_flux_advection_x(
			all_grid,
			time_step_length,
			use_MacCormack_scheme,
			use_clamping_in_MacCormack_scheme,
			num_gauss_quad_boundary,
			num_gauss_quad_bulk,
			split_method,
			integral_method,
			velocity_after_advect_x,
//			all_grid.velocity_in_voxel_face_x,
			interpolation_method,
			interpolation_method_in_calc_psi,
			use_integral,
			num_gauss_quadrature_point_for_integrate_density,
			enable_cell_volume_correction
		);
/*
		//x成分の移流
		advect_velocity_semi_lagrangian_x(
			all_grid,
			time_step_length,
			interpolation_method,
			velocity_after_advect_x
		);
*/

		//y成分の移流
		advect_velocity_flux_advection_y(
			all_grid,
			time_step_length,
			use_MacCormack_scheme,
			use_clamping_in_MacCormack_scheme,
			num_gauss_quad_boundary,
			num_gauss_quad_bulk,
			split_method,
			integral_method,
			velocity_after_advect_y,
//			all_grid.velocity_in_voxel_face_x,
			interpolation_method,
			interpolation_method_in_calc_psi,
			use_integral,
			num_gauss_quadrature_point_for_integrate_density,
			enable_cell_volume_correction
		);
/*
		//y成分の移流
		advect_velocity_semi_lagrangian_y(
			all_grid,
			time_step_length,
			interpolation_method,
			velocity_after_advect_y
		);
*/
/*
		//計算結果をコピー
		for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= velocity_after_advect_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
			}
		}
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
				all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= velocity_after_advect_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
			}
		}
*/
		////////////////////////////////////////////////////////////
		// 0速度場の拡散補正
		////////////////////////////////////////////////////////////
		if (enable_eliminate_zero_velocity_false_diffusion) {
			/////////////////////////
			// 速度場のx成分の補正量
			/////////////////////////
			std::vector<MY_FLOAT_TYPE> zero_velocity_correction_values_velocity_x((all_grid.Grid_num_x + 1) * all_grid.Grid_num_y);
			////// "バックトレース前の面で定義されたpsiの離散値"を足す
			for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					zero_velocity_correction_values_velocity_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						= before_advect_velocity_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
				}
			}

			for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					//(ix, iy)番目のセルのface を走るループ
					int dx[4], dy[4];
					dx[0] = 0; dx[1] = 1; dx[2] = 1; dx[3] = 0;
					dy[0] = 0; dy[1] = 0; dy[2] = 1; dy[3] = 1;
					MY_FLOAT_TYPE velocity_x_in_the_cell = 0.0;
					MY_FLOAT_TYPE cell_volume = 0.0;
					for (int i_face = 0; i_face < 4; ++i_face) {
						// 考えるcell face を構成するvertex のindex
						std::vector<int>
							vertex_index_1{ ix + dx[i_face], iy + dy[i_face] },
							vertex_index_2{ ix + dx[(i_face + 1) % 4], iy + dy[(i_face + 1) % 4] };
						// バックトレース前の cell face を構成するvertexを定義する
						cell_vertex before_backtrace_vertex_1(vertex_index_1, all_grid._cell_length);
						cell_vertex before_backtrace_vertex_2(vertex_index_2, all_grid._cell_length);
						//// vertex の位置を定義
						// vertex_1 の位置を設定
						before_backtrace_vertex_1._vertex_pos[0] = (ix + dx[i_face] - 0.5) * all_grid._cell_length;
						before_backtrace_vertex_1._vertex_pos[1] = (iy + dy[i_face]) * all_grid._cell_length;
						before_backtrace_vertex_1._vertex_pos[2] = 0.0;
						// vertex_2 の位置を設定
						before_backtrace_vertex_2._vertex_pos[0] = (ix + dx[(i_face + 1) % 4] - 0.5) * all_grid._cell_length;
						before_backtrace_vertex_2._vertex_pos[1] = (iy + dy[(i_face + 1) % 4]) * all_grid._cell_length;
						before_backtrace_vertex_2._vertex_pos[2] = 0.0;

						// バックトレース前の頂点から面を定義
						cell_face before_backtrace_face;
						before_backtrace_face._vertex_list[0] = before_backtrace_vertex_1;
						before_backtrace_face._vertex_list[1] = before_backtrace_vertex_2;
						////// エラーの補正
						////// "バックトレース前の面上で積分した値" と "バックトレース前の面で定義されたpsiの離散値"の差を引いていく
						zero_velocity_correction_values_velocity_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
							-= integrate_normal_component_of_psi_velocity_x(
								all_grid,
								before_backtrace_face,
								num_gauss_quad_boundary,
								num_gauss_quad_bulk,
								split_method,
								interpolation_method,
								use_integral,
								num_gauss_quadrature_point_for_integrate_density
							) / (all_grid._cell_volume);
					}
				}
			}
			// エラーの移流と密度の移流で同じ補間の方法を使う
//	            std::string interpolation_method_of_eccor_advection = interpolation_method;

			std::string interpolation_method_of_eccor_advection;
			if(interpolation_method == "linear" || interpolation_method == "1Dy_linear"){
				interpolation_method_of_eccor_advection = "1Dy_linear";
//				interpolation_method_of_eccor_advection = "linear";
			}
			else if(interpolation_method == "WENO6"||interpolation_method == "WENO6-optimized"){
				interpolation_method_of_eccor_advection = "1Dy_WENO6";
			}
			std::string split_method_of_eccor_advection = "y-CellFaceAligned";
			std::string interpolation_method_in_calc_psi = "const";
			// エラーを移流させる
			advect_velocity_flux_advection_x(
				all_grid,
				time_step_length,
//					set_zero_normal_velocity_at_boundary,
				use_MacCormack_scheme,
				use_clamping_in_MacCormack_scheme,
				num_gauss_quad_boundary,
				num_gauss_quad_bulk,
//				enable_eliminate_zero_velocity_false_diffusion,
				split_method_of_eccor_advection,
				integral_method,
				zero_velocity_correction_values_velocity_x,
				interpolation_method_of_eccor_advection,
				interpolation_method_in_calc_psi,
				use_integral,
				num_gauss_quadrature_point_for_integrate_density,
				/*enable_cell_volume_correction*/false
			);

			/////////////////////////
			// 速度場のy成分の補正量
			/////////////////////////
			std::vector<MY_FLOAT_TYPE> zero_velocity_correction_values_velocity_y(all_grid.Grid_num_x * (all_grid.Grid_num_y + 1));
			////// "バックトレース前の面で定義されたpsiの離散値"を足す
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
					zero_velocity_correction_values_velocity_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						= before_advect_velocity_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
				}
			}

			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
					//(ix, iy)番目のセルのface を走るループ
					int dx[4], dy[4];
					dx[0] = 0; dx[1] = 1; dx[2] = 1; dx[3] = 0;
					dy[0] = 0; dy[1] = 0; dy[2] = 1; dy[3] = 1;
					MY_FLOAT_TYPE velocity_x_in_the_cell = 0.0;
					MY_FLOAT_TYPE cell_volume = 0.0;
					for (int i_face = 0; i_face < 4; ++i_face) {
						// 考えるcell face を構成するvertex のindex
						std::vector<int>
							vertex_index_1{ ix + dx[i_face], iy + dy[i_face] },
							vertex_index_2{ ix + dx[(i_face + 1) % 4], iy + dy[(i_face + 1) % 4] };
						// バックトレース前の cell face を構成するvertexを定義する
						cell_vertex before_backtrace_vertex_1(vertex_index_1, all_grid._cell_length);
						cell_vertex before_backtrace_vertex_2(vertex_index_2, all_grid._cell_length);
						//// vertex の位置を定義
						// vertex_1 の位置を設定
						before_backtrace_vertex_1._vertex_pos[0] = (ix + dx[i_face]) * all_grid._cell_length;
						before_backtrace_vertex_1._vertex_pos[1] = (iy + dy[i_face] - 0.5) * all_grid._cell_length;
						before_backtrace_vertex_1._vertex_pos[2] = 0.0;
						// vertex_2 の位置を設定
						before_backtrace_vertex_2._vertex_pos[0] = (ix + dx[(i_face + 1) % 4]) * all_grid._cell_length;
						before_backtrace_vertex_2._vertex_pos[1] = (iy + dy[(i_face + 1) % 4] - 0.5) * all_grid._cell_length;
						before_backtrace_vertex_2._vertex_pos[2] = 0.0;

						// バックトレース前の頂点から面を定義
						cell_face before_backtrace_face;
						before_backtrace_face._vertex_list[0] = before_backtrace_vertex_1;
						before_backtrace_face._vertex_list[1] = before_backtrace_vertex_2;
						////// エラーの補正
						////// "バックトレース前の面上で積分した値" と "バックトレース前の面で定義されたpsiの離散値"の差を引いていく
						zero_velocity_correction_values_velocity_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
							-= integrate_normal_component_of_psi_velocity_y(
								all_grid,
								before_backtrace_face,
								num_gauss_quad_boundary,
								num_gauss_quad_bulk,
								split_method,
								interpolation_method,
								use_integral,
								num_gauss_quadrature_point_for_integrate_density
							)
							/ (all_grid._cell_volume);
					}
				}
			}
			// エラーの移流と密度の移流で同じ補間の方法を使う
//	            std::string interpolation_method_of_eccor_advection = interpolation_method;
			if(interpolation_method == "linear" || interpolation_method == "1Dy_linear"){
				interpolation_method_of_eccor_advection = "1Dy_linear";
//				interpolation_method_of_eccor_advection = "linear";
			}
			else if(interpolation_method == "WENO6"||interpolation_method == "WENO6-optimized"){
				interpolation_method_of_eccor_advection = "1Dy_WENO6";
			}
			split_method_of_eccor_advection = "y-CellFaceAligned";
			interpolation_method_in_calc_psi = "const";
			// エラーを移流させる
			advect_velocity_flux_advection_y(
				all_grid,
				time_step_length,
//					set_zero_normal_velocity_at_boundary,
				use_MacCormack_scheme,
				use_clamping_in_MacCormack_scheme,
				num_gauss_quad_boundary,
				num_gauss_quad_bulk,
//				enable_eliminate_zero_velocity_false_diffusion,
				split_method_of_eccor_advection,
				integral_method,
				zero_velocity_correction_values_velocity_y,
				interpolation_method_of_eccor_advection,
				interpolation_method_in_calc_psi,
				use_integral,
				num_gauss_quadrature_point_for_integrate_density,
				/*enable_cell_volume_correction*/false
			);
			//計算結果をコピー
			for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						= velocity_after_advect_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
				}
			}
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
					all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						= velocity_after_advect_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
				}
			}
			//velocity_xを補正
			for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+=zero_velocity_correction_values_velocity_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
				}
			}
			//velocity_yを補正
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
					all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+=zero_velocity_correction_values_velocity_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
				}
			}
		}
		// 0速度場の拡散を補正しない場合
		else{
			//計算結果をコピー
			for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						= velocity_after_advect_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
				}
			}
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
					all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						= velocity_after_advect_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
				}
			}
		}

		correct_velocity_x_and_volume_by_pressure_solve(
			all_grid
		);
		correct_velocity_y_and_volume_by_pressure_solve(
			all_grid
		);

//		for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
//			all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(0, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//				= 0.0;
//			all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(all_grid.Grid_num_x, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//				= 0.0;
//		}

	}

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
	) {
		// psi の離散値を計算
		calc_discrete_psi_velocity_x(
			all_grid,
			all_grid.psi_velocity_cell_vertex_x,
			advected_values,
			all_grid.Grid_num_x,
			all_grid.Grid_num_y,
			all_grid._cell_length,
			interpolation_method_in_calc_psi
		);
		//計算結果を格納する変数
		std::vector<MY_FLOAT_TYPE> velocity_x_after_advect((all_grid.Grid_num_x + 1) * all_grid.Grid_num_y);
		//移流後のセル体積を格納する変数
		std::vector<MY_FLOAT_TYPE> cell_volume_after_advect((all_grid.Grid_num_x + 1) * all_grid.Grid_num_y);
		for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				//(ix, iy)番目のセルのface を走るループ
				int dx[4], dy[4];
				dx[0] = 0; dx[1] = 1; dx[2] = 1; dx[3] = 0;
				dy[0] = 0; dy[1] = 0; dy[2] = 1; dy[3] = 1;
				MY_FLOAT_TYPE velocity_x_in_the_cell = 0.0;
				MY_FLOAT_TYPE cell_volume = 0.0;
				for (int i_face = 0; i_face < 4; ++i_face) {
					// 考えるcell face を構成するvertex のindex
					std::vector<int>
						vertex_index_1{ ix + dx[i_face], iy + dy[i_face] },
						vertex_index_2{ ix + dx[(i_face + 1) % 4], iy + dy[(i_face + 1) % 4] };
					// バックトレース前の cell face を構成するvertexを定義する
					cell_vertex before_backtrace_vertex_1(vertex_index_1, all_grid._cell_length);
					cell_vertex before_backtrace_vertex_2(vertex_index_2, all_grid._cell_length);
					//// vertex の位置を定義
					// vertex_1 の位置を設定
					before_backtrace_vertex_1._vertex_pos[0] = (ix + dx[i_face] - 0.5) * all_grid._cell_length;
					before_backtrace_vertex_1._vertex_pos[1] = (iy + dy[i_face]) * all_grid._cell_length;
					before_backtrace_vertex_1._vertex_pos[2] = 0.0;
					// vertex_2 の位置を設定
					before_backtrace_vertex_2._vertex_pos[0] = (ix + dx[(i_face + 1) % 4] - 0.5) * all_grid._cell_length;
					before_backtrace_vertex_2._vertex_pos[1] = (iy + dy[(i_face + 1) % 4]) * all_grid._cell_length;
					before_backtrace_vertex_2._vertex_pos[2] = 0.0;

					// バックトレース前の頂点から面を定義
					cell_face before_backtrace_face;
					before_backtrace_face._vertex_list[0] = before_backtrace_vertex_1;
					before_backtrace_face._vertex_list[1] = before_backtrace_vertex_2;

					// MacCormack を使う場合
					if (use_MacCormack_scheme) {
						////エラーを補正する前のベースとなる面での計算
						// バックトレース後の頂点から面を定義
						cell_face after_backtrace_face
							= calc_backtraced_face(
								all_grid,
								time_step_length,
								before_backtrace_face,
								"forward",
								interpolation_method,
								"velocity_x"
							);
						//バックトレースした先の位置がグリッドの範囲を超えていたら修正する
//						after_backtrace_face = clamp_vertex_position_by_grid(all_grid, after_backtrace_face);
						// バックトレース後の面について psi と 法線の内積を面上で積分した値を計算する
						MY_FLOAT_TYPE integral_of_normal_component_of_psi_on_after_backtrace_face;
						if(integral_method == "gauss"){
							integral_of_normal_component_of_psi_on_after_backtrace_face
								= integrate_normal_component_of_psi_velocity_x(
									all_grid,
									after_backtrace_face,
									num_gauss_quad_boundary,
									num_gauss_quad_bulk,
									split_method,
									interpolation_method,
									use_integral,
									num_gauss_quadrature_point_for_integrate_density
								);
						}
//						if(enable_cell_volume_correction){
//							calc_cell_volume_contribution(
//								after_backtrace_face,
//								cell_volume
//							);
//						}

						//// エラーの計算
						// 元のface(バックトレース前のオリジナルの face)について, psi と 法線の内積を面上で積分した値を計算する
						MY_FLOAT_TYPE integral_of_normal_component_of_psi_on_before_backtrace_face;
						if(integral_method == "gauss"){
							integral_of_normal_component_of_psi_on_before_backtrace_face
							= integrate_normal_component_of_psi_velocity_x(
								all_grid,
								before_backtrace_face,
								num_gauss_quad_boundary,
								num_gauss_quad_bulk,
								split_method,
								interpolation_method,
								use_integral,
								num_gauss_quadrature_point_for_integrate_density
							);
						}
						// バックトレース後の面を逆向きにバックトレースした面の計算(エラーの計算用に使う)
						cell_face backward_trace_face_of_after_backtrace_face
							= calc_backtraced_face(
								all_grid,
								time_step_length,
								after_backtrace_face,
								"backward",
								interpolation_method,
								"velocity_x"
							);

						//バックトレースした先の位置がグリッドの範囲を超えていたら修正する
//						backward_trace_face_of_after_backtrace_face
//							= clamp_vertex_position_by_grid(
//								all_grid,
//								backward_trace_face_of_after_backtrace_face);
						// psi と 法線の内積を面上で積分した値を計算する
						MY_FLOAT_TYPE integral_of_normal_component_of_psi_on_backward_trace_face_of_after_backtrace_face;
						if(integral_method == "gauss"){
							integral_of_normal_component_of_psi_on_backward_trace_face_of_after_backtrace_face
								= integrate_normal_component_of_psi_velocity_x(
									all_grid,
									backward_trace_face_of_after_backtrace_face,
									num_gauss_quad_boundary,
									num_gauss_quad_bulk,
									split_method,
									interpolation_method,
									use_integral,
									num_gauss_quadrature_point_for_integrate_density
								);
						}

						// MacCormack におけるエラーの定義
						MY_FLOAT_TYPE error
							= (integral_of_normal_component_of_psi_on_before_backtrace_face
								- integral_of_normal_component_of_psi_on_backward_trace_face_of_after_backtrace_face)
							/ 2.0;
						// エラーを修正した値
						MY_FLOAT_TYPE error_corrected_integral_of_normal_component_of_psi
							= integral_of_normal_component_of_psi_on_after_backtrace_face + error;
						//clamping
						if (use_clamping_in_MacCormack_scheme) {
							error_corrected_integral_of_normal_component_of_psi
								= clamp_integral_of_normal_component_of_psi_velocity_x_of_MacCormack(
									all_grid,
									error_corrected_integral_of_normal_component_of_psi,
									after_backtrace_face.calc_face_center(),
									interpolation_method);
						}
						// エラーを修正した値を用いて質量への寄与を計算する
						velocity_x_in_the_cell
							+= integral_of_normal_component_of_psi_on_after_backtrace_face
							+ error;
					}
					// semi-Lagrangian を使う場合
					else {
						// バックトレース後の頂点から面を定義
						cell_face after_backtrace_face
							= calc_backtraced_face(
								all_grid,
								time_step_length,
								before_backtrace_face,
								"forward",
								interpolation_method,
								"velocity_x"
							);
						//バックトレースした先の位置がグリッドの範囲を超えていたら修正する
//						after_backtrace_face = clamp_vertex_position_by_grid(all_grid, after_backtrace_face);

						// psi と 法線の内積を面上で積分した値を計算する
						if(integral_method == "gauss"){
							velocity_x_in_the_cell
								+= integrate_normal_component_of_psi_velocity_x(
									all_grid,
									after_backtrace_face,
									num_gauss_quad_boundary,
									num_gauss_quad_bulk,
									split_method,
									interpolation_method,
									use_integral,
									num_gauss_quadrature_point_for_integrate_density
								);
						}
//						else if(integral_method == "analytical"){
//							mass_in_the_cell
//								+= integrate_normal_component_of_psi_on_face_analytically(
//									all_grid,
//									after_backtrace_face,
//									split_method);
//						}
						if(enable_cell_volume_correction){
							calc_cell_volume_contribution(
								after_backtrace_face,
								cell_volume
							);
						}
					}
				}
				//質量の計算結果を格納
				velocity_x_after_advect[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= velocity_x_in_the_cell / (all_grid._cell_volume);

				//体積の計算結果を格納
				if(enable_cell_volume_correction){
					cell_volume_after_advect[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						= cell_volume;
				}
				else{
					cell_volume_after_advect[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						= all_grid._cell_volume;
				}

			}
		}

		//計算結果をコピー
		for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				//体積の計算結果
				all_grid.cell_volume_cell_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= cell_volume_after_advect[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
				//密度の計算結果
				advected_values[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= velocity_x_after_advect[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
//					/ all_grid._cell_volume;
			}
		}

	}

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
	) {
		// psi の離散値を計算
		calc_discrete_psi_velocity_y(
			all_grid,
			all_grid.psi_velocity_cell_center_y,
			advected_values,
			all_grid.Grid_num_x,
			all_grid.Grid_num_y,
			all_grid._cell_length,
			interpolation_method_in_calc_psi
		);
		//計算結果を格納する変数
		std::vector<MY_FLOAT_TYPE> velocity_y_after_advect(all_grid.Grid_num_x * (all_grid.Grid_num_y + 1));
		//移流後のセル体積を格納する変数
		std::vector<MY_FLOAT_TYPE> cell_volume_after_advect(all_grid.Grid_num_x * (all_grid.Grid_num_y + 1));
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
				//(ix, iy)番目のセルのface を走るループ
				int dx[4], dy[4];
				dx[0] = 0; dx[1] = 1; dx[2] = 1; dx[3] = 0;
				dy[0] = 0; dy[1] = 0; dy[2] = 1; dy[3] = 1;
				MY_FLOAT_TYPE velocity_y_in_the_cell = 0.0;
				MY_FLOAT_TYPE cell_volume = 0.0;
				for (int i_face = 0; i_face < 4; ++i_face) {
					// 考えるcell face を構成するvertex のindex
					std::vector<int>
						vertex_index_1{ ix + dx[i_face], iy + dy[i_face] },
						vertex_index_2{ ix + dx[(i_face + 1) % 4], iy + dy[(i_face + 1) % 4] };
					// バックトレース前の cell face を構成するvertexを定義する
					cell_vertex before_backtrace_vertex_1(vertex_index_1, all_grid._cell_length);
					cell_vertex before_backtrace_vertex_2(vertex_index_2, all_grid._cell_length);
					//// vertex の位置を定義
					// vertex_1 の位置を設定
					before_backtrace_vertex_1._vertex_pos[0] = (ix + dx[i_face]) * all_grid._cell_length;
					before_backtrace_vertex_1._vertex_pos[1] = (iy + dy[i_face] - 0.5) * all_grid._cell_length;
					before_backtrace_vertex_1._vertex_pos[2] = 0.0;
					// vertex_2 の位置を設定
					before_backtrace_vertex_2._vertex_pos[0] = (ix + dx[(i_face + 1) % 4]) * all_grid._cell_length;
					before_backtrace_vertex_2._vertex_pos[1] = (iy + dy[(i_face + 1) % 4] - 0.5) * all_grid._cell_length;
					before_backtrace_vertex_2._vertex_pos[2] = 0.0;

					// バックトレース前の頂点から面を定義
					cell_face before_backtrace_face;
					before_backtrace_face._vertex_list[0] = before_backtrace_vertex_1;
					before_backtrace_face._vertex_list[1] = before_backtrace_vertex_2;

					// MacCormack を使う場合
					if (use_MacCormack_scheme) {
						////エラーを補正する前のベースとなる面での計算
						// バックトレース後の頂点から面を定義
						cell_face after_backtrace_face
							= calc_backtraced_face(
								all_grid,
								time_step_length,
								before_backtrace_face,
								"forward",
								interpolation_method,
								"velocity_y"
							);
						//バックトレースした先の位置がグリッドの範囲を超えていたら修正する
//						after_backtrace_face = clamp_vertex_position_by_grid(all_grid, after_backtrace_face);
						// バックトレース後の面について psi と 法線の内積を面上で積分した値を計算する
						MY_FLOAT_TYPE integral_of_normal_component_of_psi_on_after_backtrace_face;
						if(integral_method == "gauss"){
							integral_of_normal_component_of_psi_on_after_backtrace_face
								= integrate_normal_component_of_psi_velocity_y(
									all_grid,
									after_backtrace_face,
									num_gauss_quad_boundary,
									num_gauss_quad_bulk,
									split_method,
									interpolation_method,
									use_integral,
									num_gauss_quadrature_point_for_integrate_density
								);
						}
//						if(enable_cell_volume_correction){
//							calc_cell_volume_contribution(
//								after_backtrace_face,
//								cell_volume
//							);
//						}

						//// エラーの計算
						// 元のface(バックトレース前のオリジナルの face)について, psi と 法線の内積を面上で積分した値を計算する
						MY_FLOAT_TYPE integral_of_normal_component_of_psi_on_before_backtrace_face;
						if(integral_method == "gauss"){
							integral_of_normal_component_of_psi_on_before_backtrace_face
							= integrate_normal_component_of_psi_velocity_y(
								all_grid,
								before_backtrace_face,
								num_gauss_quad_boundary,
								num_gauss_quad_bulk,
								split_method,
								interpolation_method,
								use_integral,
								num_gauss_quadrature_point_for_integrate_density
							);
						}
						// バックトレース後の面を逆向きにバックトレースした面の計算(エラーの計算用に使う)
						cell_face backward_trace_face_of_after_backtrace_face
							= calc_backtraced_face(
								all_grid,
								time_step_length,
								after_backtrace_face,
								"backward",
								interpolation_method,
								"velocity_y"
							);

						//バックトレースした先の位置がグリッドの範囲を超えていたら修正する
//						backward_trace_face_of_after_backtrace_face
//							= clamp_vertex_position_by_grid(
//								all_grid,
//								backward_trace_face_of_after_backtrace_face);
						// psi と 法線の内積を面上で積分した値を計算する
						MY_FLOAT_TYPE integral_of_normal_component_of_psi_on_backward_trace_face_of_after_backtrace_face;
						if(integral_method == "gauss"){
							integral_of_normal_component_of_psi_on_backward_trace_face_of_after_backtrace_face
								= integrate_normal_component_of_psi_velocity_y(
									all_grid,
									backward_trace_face_of_after_backtrace_face,
									num_gauss_quad_boundary,
									num_gauss_quad_bulk,
									split_method,
									interpolation_method,
									use_integral,
									num_gauss_quadrature_point_for_integrate_density
								);
						}

						// MacCormack におけるエラーの定義
						MY_FLOAT_TYPE error
							= (integral_of_normal_component_of_psi_on_before_backtrace_face
								- integral_of_normal_component_of_psi_on_backward_trace_face_of_after_backtrace_face)
							/ 2.0;
						// エラーを修正した値
						MY_FLOAT_TYPE error_corrected_integral_of_normal_component_of_psi
							= integral_of_normal_component_of_psi_on_after_backtrace_face + error;
						//clamping
						if (use_clamping_in_MacCormack_scheme) {
							error_corrected_integral_of_normal_component_of_psi
								= clamp_integral_of_normal_component_of_psi_velocity_y_of_MacCormack(
									all_grid,
									error_corrected_integral_of_normal_component_of_psi,
									after_backtrace_face.calc_face_center(),
									interpolation_method);
						}
						// エラーを修正した値を用いて質量への寄与を計算する
						velocity_y_in_the_cell
							+= integral_of_normal_component_of_psi_on_after_backtrace_face
							+ error;
					}
					// semi-Lagrangian を使う場合
					else {
						// バックトレース後の頂点から面を定義
						cell_face after_backtrace_face
							= calc_backtraced_face(
								all_grid,
								time_step_length,
								before_backtrace_face,
								"forward",
								interpolation_method,
								"velocity_y"
							);
						//バックトレースした先の位置がグリッドの範囲を超えていたら修正する
//						after_backtrace_face = clamp_vertex_position_by_grid(all_grid, after_backtrace_face);

						// psi と 法線の内積を面上で積分した値を計算する
						if(integral_method == "gauss"){
							velocity_y_in_the_cell
								+= integrate_normal_component_of_psi_velocity_y(
									all_grid,
									after_backtrace_face,
									num_gauss_quad_boundary,
									num_gauss_quad_bulk,
									split_method,
									interpolation_method,
									use_integral,
									num_gauss_quadrature_point_for_integrate_density
								);
						}
//						else if(integral_method == "analytical"){
//							mass_in_the_cell
//								+= integrate_normal_component_of_psi_on_face_analytically(
//									all_grid,
//									after_backtrace_face,
//									split_method);
//						}
						if(enable_cell_volume_correction){
							calc_cell_volume_contribution(
								after_backtrace_face,
								cell_volume
							);
						}
					}
				}
				//質量の計算結果を格納
				velocity_y_after_advect[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= velocity_y_in_the_cell / (all_grid._cell_volume);

				//体積の計算結果を格納
				if(enable_cell_volume_correction){
					cell_volume_after_advect[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						= cell_volume;
				}
				else{
					cell_volume_after_advect[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						= all_grid._cell_volume;
				}

			}
		}

		//計算結果をコピー
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
				//体積の計算結果
				all_grid.cell_volume_cell_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= cell_volume_after_advect[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
				//密度の計算結果
				advected_values[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= velocity_y_after_advect[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
//					/ all_grid._cell_volume;
			}
		}

	}

	void correct_velocity_x_and_volume_by_pressure_solve(
		Grid& all_grid
	){
		// 圧力計算
		//グリッドの総数
		const int N = (all_grid.Grid_num_x + 1)
					* all_grid.Grid_num_y;
		//係数行列の計算
		//ノイマン境界条件では係数行列をナイーブに構築するとフルランクにならない(1つランクが落ちる)
		//フルランクにするために圧力の最後の要素を0に固定する。それに合わせて係数行列も最後の行, 列を落とす
//		linear_algebra::sparse_matrix A(N, N);
		linear_algebra::sparse_matrix A(N - 1, N - 1);
//		linear_algebra::sparse_matrix_with_diagonal_element A(N - 1, N - 1);
//		for(int i_xy = 0; i_xy < N; ++i_xy){
		for(int i_xy = 0; i_xy < N - 1; ++i_xy){
			int iy = i_xy % all_grid.Grid_num_y;
			int ix = (i_xy - iy) / all_grid.Grid_num_y;
			//非対角項
			if (ix != 0) {
				int index_0=get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				int index_1=get_voxel_face_index_x(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				if((index_0 != N - 1) && (index_1 != N - 1)){
					A.input_element(index_0, index_1, 1);
				}
			}
			//非対角項
			if (iy != 0) {
				int index_0=get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				int index_1=get_voxel_face_index_x(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y);
				if((index_0 != N - 1) && (index_1 != N - 1)){
					A.input_element(index_0, index_1, 1);
				}
			}
			//対角項
			MY_FLOAT_TYPE diagonal_value = -4;
			if(ix == 0){
				diagonal_value += 1;
			}
			if(iy == 0){
				diagonal_value += 1;
			}
			if(iy == all_grid.Grid_num_y - 1){
				diagonal_value += 1;
			}
			if(ix == all_grid.Grid_num_x){
				diagonal_value += 1;
			}
			A.input_element(
				get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y),
				get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y),
				diagonal_value
			);
			//非対角項
			if (iy != all_grid.Grid_num_y - 1) {
				int index_0=get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				int index_1=get_voxel_face_index_x(ix, iy + 1, all_grid.Grid_num_x, all_grid.Grid_num_y);
				if((index_0 != N - 1) && (index_1 != N - 1)){
					A.input_element(index_0, index_1, 1);
				}
			}
			//非対角項
			if (ix != all_grid.Grid_num_x) {
				int index_0=get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				int index_1=get_voxel_face_index_x(ix + 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				if((index_0 != N - 1) && (index_1 != N - 1)){
					A.input_element(index_0, index_1, 1);
				}
			}
		}
		//連立方程式の右辺のベクトルの計算
//		std::vector<MY_FLOAT_TYPE> b(N);
		std::vector<MY_FLOAT_TYPE> b(N - 1);
//		for(int i_xy = 0; i_xy < N; ++i_xy){
		for(int i_xy = 0; i_xy < N - 1; ++i_xy){
			int iy = i_xy % all_grid.Grid_num_y;
			int ix = (i_xy - iy) / all_grid.Grid_num_y;

			b[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
				= -(all_grid.cell_volume_cell_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
				- all_grid._cell_volume);
		}

		//amgcl ライブラリで解く
		typedef amgcl::backend::builtin<MY_FLOAT_TYPE> Backend;
		typedef amgcl::make_solver<
			// Use AMG as preconditioner:
			amgcl::amg<
				Backend,
				amgcl::coarsening::smoothed_aggregation,
//		        amgcl::relaxation::spai0
				amgcl::relaxation::gauss_seidel
				>,
			// And BiCGStab as iterative solver:
//	    	amgcl::solver::bicgstab<Backend>
			amgcl::solver::cg<Backend>
			> Solver;
//		int n = N;
		int n = N - 1;
//		int n = all_grid.Grid_num_x * all_grid.Grid_num_y;
		Solver::params prm;
	    prm.solver.tol = 1e-1;
		Solver solve( std::tie(n, A.cumulative_num_nonzero_element_in_row, A.column_index, A.element_value) ,prm);
//		Solver solve( std::tie(n, A.cumulative_num_nonzero_element_in_row, A.column_index, A.element_value));
		int    iters;
		MY_FLOAT_TYPE error;
		std::vector<MY_FLOAT_TYPE> pressure(N - 1);
		std::tie(iters, error) = solve(b, pressure);
		//係数行列をフルランクにするために落とした最後の要素を復元する
		pressure.push_back(0.0);
		// 補正前の質量と体積
		std::vector<MY_FLOAT_TYPE> before_correct_velocity_x((all_grid.Grid_num_x + 1) * all_grid.Grid_num_y);
		std::vector<MY_FLOAT_TYPE> before_correct_volume((all_grid.Grid_num_x + 1) * all_grid.Grid_num_y);
		for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				before_correct_velocity_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
//					* all_grid._cell_volume;
				before_correct_volume[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= all_grid.cell_volume_cell_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
			}
		}
		//// 体積の勾配によって質量と体積を拡散させる
		MY_FLOAT_TYPE diffusion_coefficient = 1.0;
		// x 方向の拡散
		for (int ix = 2; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				//フラックス
				MY_FLOAT_TYPE volume_flux
					= -(pressure[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					  - pressure[get_voxel_face_index_x(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					 / all_grid._cell_volume;
				// 面での密度
				MY_FLOAT_TYPE density_at_face
						= (before_correct_velocity_x[get_voxel_face_index_x(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+  before_correct_velocity_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
						/ 2.0;
				// 面での体積
				MY_FLOAT_TYPE volume_at_face
					= (before_correct_volume[get_voxel_face_index_x(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+  before_correct_volume[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					/ 2.0;
				//密度の拡散
				all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					-= diffusion_coefficient * volume_flux * density_at_face;
				all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+= diffusion_coefficient * volume_flux * density_at_face;
				//体積の拡散
				all_grid.cell_volume_cell_face_x[get_voxel_face_index_x(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					-= diffusion_coefficient * volume_flux * volume_at_face;
				all_grid.cell_volume_cell_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+= diffusion_coefficient * volume_flux * volume_at_face;
			}
		}
		// y 方向の拡散
		for (int ix = 1; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 1; iy < all_grid.Grid_num_y; iy++) {
				//フラックス
				MY_FLOAT_TYPE volume_flux
					= -(pressure[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					  - pressure[get_voxel_face_index_x(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					 / all_grid._cell_volume;
				// 面での密度
				MY_FLOAT_TYPE density_at_face
					= (before_correct_velocity_x[get_voxel_face_index_x(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+  before_correct_velocity_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					/ 2.0;
				// 面での体積
				MY_FLOAT_TYPE volume_at_face
					= (before_correct_volume[get_voxel_face_index_x(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+  before_correct_volume[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					/ 2.0;
				//密度の拡散
				all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					-= diffusion_coefficient * volume_flux * density_at_face;
				all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+= diffusion_coefficient * volume_flux * density_at_face;
				//体積の拡散
				all_grid.cell_volume_cell_face_x[get_voxel_face_index_x(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					-= diffusion_coefficient * volume_flux * volume_at_face;
				all_grid.cell_volume_cell_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+= diffusion_coefficient * volume_flux * volume_at_face;
			}
		}
	}

	void correct_velocity_y_and_volume_by_pressure_solve(
		Grid& all_grid
	){
		// 圧力計算
		//グリッドの総数
		const int N = all_grid.Grid_num_x
					* (all_grid.Grid_num_y + 1);
		//係数行列の計算
		//ノイマン境界条件では係数行列をナイーブに構築するとフルランクにならない(1つランクが落ちる)
		//フルランクにするために圧力の最後の要素を0に固定する。それに合わせて係数行列も最後の行, 列を落とす
//		linear_algebra::sparse_matrix A(N, N);
		linear_algebra::sparse_matrix A(N - 1, N - 1);
//		linear_algebra::sparse_matrix_with_diagonal_element A(N - 1, N - 1);
//		for(int i_xy = 0; i_xy < N; ++i_xy){
		for(int i_xy = 0; i_xy < N - 1; ++i_xy){
			int iy = i_xy % (all_grid.Grid_num_y + 1);
			int ix = (i_xy - iy) / (all_grid.Grid_num_y + 1);
			//非対角項
			if (ix != 0) {
				int index_0=get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				int index_1=get_voxel_face_index_y(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				if((index_0 != N - 1) && (index_1 != N - 1)){
					A.input_element(index_0, index_1, 1);
				}
			}
			//非対角項
			if (iy != 0) {
				int index_0=get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				int index_1=get_voxel_face_index_y(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y);
				if((index_0 != N - 1) && (index_1 != N - 1)){
					A.input_element(index_0, index_1, 1);
				}
			}
			//対角項
			MY_FLOAT_TYPE diagonal_value = -4;
			if(ix == 0){
				diagonal_value += 1;
			}
			if(iy == 0){
				diagonal_value += 1;
			}
			if(iy == all_grid.Grid_num_y){
				diagonal_value += 1;
			}
			if(ix == all_grid.Grid_num_x - 1){
				diagonal_value += 1;
			}
			A.input_element(
				get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y),
				get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y),
				diagonal_value
			);
			//非対角項
			if (iy != all_grid.Grid_num_y) {
				int index_0=get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				int index_1=get_voxel_face_index_y(ix, iy + 1, all_grid.Grid_num_x, all_grid.Grid_num_y);
				if((index_0 != N - 1) && (index_1 != N - 1)){
					A.input_element(index_0, index_1, 1);
				}
			}
			//非対角項
			if (ix != all_grid.Grid_num_x - 1) {
				int index_0=get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				int index_1=get_voxel_face_index_y(ix + 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				if((index_0 != N - 1) && (index_1 != N - 1)){
					A.input_element(index_0, index_1, 1);
				}
			}
		}
		//連立方程式の右辺のベクトルの計算
//		std::vector<MY_FLOAT_TYPE> b(N);
		std::vector<MY_FLOAT_TYPE> b(N - 1);
//		for(int i_xy = 0; i_xy < N; ++i_xy){
		for(int i_xy = 0; i_xy < N - 1; ++i_xy){
			int iy = i_xy % (all_grid.Grid_num_y + 1);
			int ix = (i_xy - iy) / (all_grid.Grid_num_y + 1);

			b[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
				= -(all_grid.cell_volume_cell_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
				- all_grid._cell_volume);
		}

		//amgcl ライブラリで解く
		typedef amgcl::backend::builtin<MY_FLOAT_TYPE> Backend;
		typedef amgcl::make_solver<
			// Use AMG as preconditioner:
			amgcl::amg<
				Backend,
				amgcl::coarsening::smoothed_aggregation,
//		        amgcl::relaxation::spai0
				amgcl::relaxation::gauss_seidel
				>,
			// And BiCGStab as iterative solver:
//	    	amgcl::solver::bicgstab<Backend>
			amgcl::solver::cg<Backend>
			> Solver;
//		int n = N;
		int n = N - 1;
//		int n = all_grid.Grid_num_x * all_grid.Grid_num_y;
		Solver::params prm;
	    prm.solver.tol = 1e-1;
		Solver solve( std::tie(n, A.cumulative_num_nonzero_element_in_row, A.column_index, A.element_value) ,prm);
//		Solver solve( std::tie(n, A.cumulative_num_nonzero_element_in_row, A.column_index, A.element_value));
		int    iters;
		MY_FLOAT_TYPE error;
		std::vector<MY_FLOAT_TYPE> pressure(N - 1);
		std::tie(iters, error) = solve(b, pressure);
		//係数行列をフルランクにするために落とした最後の要素を復元する
		pressure.push_back(0.0);
		// 補正前の質量と体積
		std::vector<MY_FLOAT_TYPE> before_correct_velocity_y(all_grid.Grid_num_x * (all_grid.Grid_num_y + 1));
		std::vector<MY_FLOAT_TYPE> before_correct_volume(all_grid.Grid_num_x * (all_grid.Grid_num_y + 1));
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
				before_correct_velocity_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					* all_grid._cell_volume;
				before_correct_volume[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= all_grid.cell_volume_cell_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
			}
		}
		//// 体積の勾配によって質量と体積を拡散させる
		MY_FLOAT_TYPE diffusion_coefficient = 1.0;
		// x 方向の拡散
		for (int ix = 1; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 1; iy < all_grid.Grid_num_y; iy++) {
				//フラックス
				MY_FLOAT_TYPE volume_flux
					= -(pressure[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					  - pressure[get_voxel_face_index_y(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					 / all_grid._cell_volume;
				// 面での密度
				MY_FLOAT_TYPE density_at_face
					= (before_correct_velocity_y[get_voxel_face_index_y(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+  before_correct_velocity_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					/ 2.0;
				// 面での体積
				MY_FLOAT_TYPE volume_at_face
					= (before_correct_volume[get_voxel_face_index_y(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+  before_correct_volume[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					/ 2.0;
				//密度の拡散
				all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					-= diffusion_coefficient * volume_flux * density_at_face;
				all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+= diffusion_coefficient * volume_flux * density_at_face;
				//体積の拡散
				all_grid.cell_volume_cell_face_y[get_voxel_face_index_y(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					-= diffusion_coefficient * volume_flux * volume_at_face;
				all_grid.cell_volume_cell_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+= diffusion_coefficient * volume_flux * volume_at_face;
			}
		}
		// y 方向の拡散
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 2; iy < all_grid.Grid_num_y; iy++) {
				//フラックス
				MY_FLOAT_TYPE volume_flux
					= -(pressure[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					  - pressure[get_voxel_face_index_y(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					 / all_grid._cell_volume;
				// 面での密度
				MY_FLOAT_TYPE density_at_face
					= (before_correct_velocity_y[get_voxel_face_index_y(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+  before_correct_velocity_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					/ 2.0;
				// 面での体積
				MY_FLOAT_TYPE volume_at_face
					= (before_correct_volume[get_voxel_face_index_y(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+  before_correct_volume[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					/ 2.0;
				//密度の拡散
				all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					-= diffusion_coefficient * volume_flux * density_at_face;
				all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+= diffusion_coefficient * volume_flux * density_at_face;
				//体積の拡散
				all_grid.cell_volume_cell_face_y[get_voxel_face_index_y(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					-= diffusion_coefficient * volume_flux * volume_at_face;
				all_grid.cell_volume_cell_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+= diffusion_coefficient * volume_flux * volume_at_face;
			}
		}
	}

	//advect項の計算
	//速度場を時間 -dt だけバックトレースしてadvect項を計算する
	void advect_velocity_semi_lagrangian(
		Grid& all_grid,
		const MY_FLOAT_TYPE time_step_length,
		const std::string interpolation_method
	) {
		//移流後の速度場
		std::vector<MY_FLOAT_TYPE> velocity_after_advect_x((all_grid.Grid_num_x + 1) * all_grid.Grid_num_y);
		std::vector<MY_FLOAT_TYPE> velocity_after_advect_y(all_grid.Grid_num_x * (all_grid.Grid_num_y + 1));
		//x成分の移流
		advect_velocity_semi_lagrangian_x(
			all_grid,
			time_step_length,
			interpolation_method,
			velocity_after_advect_x
		);
		//y成分の移流
		advect_velocity_semi_lagrangian_y(
			all_grid,
			time_step_length,
			interpolation_method,
			velocity_after_advect_y
		);
		//計算結果をコピー
		for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= velocity_after_advect_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
			}
		}
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
				all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= velocity_after_advect_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
			}
		}

/*
		std::vector<MY_FLOAT_TYPE> velocity_after_advect_x((all_grid.Grid_num_x + 1) * all_grid.Grid_num_y);
		std::vector<MY_FLOAT_TYPE> velocity_after_advect_y(all_grid.Grid_num_x * (all_grid.Grid_num_y + 1));

		//velocityのx成分を計算
		for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				MY_FLOAT_TYPE velocity_y;
				if (ix <= 0) {
					velocity_y = (all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy + 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]) / 2.0;
				}
				else if (ix >= all_grid.Grid_num_x) {
					velocity_y = (all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix - 1, iy + 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]) / 2.0;
				}
				else {
					velocity_y = (all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix - 1, iy + 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy + 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]) / 4.0;
				}

				//バックトレース先の位置座標
				MY_FLOAT_TYPE advected_pos_x = ix * all_grid._cell_length - (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] * time_step_length);
				MY_FLOAT_TYPE advected_pos_y = (iy + 0.5) * all_grid._cell_length - (velocity_y * time_step_length);
				velocity_after_advect_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= all_grid.calc_interpolated_velocity_x(
						VEC3_TYPE(advected_pos_x, advected_pos_y, 0.0),
						interpolation_method
					);
			}
		}
		//velocityのy成分を計算
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
				MY_FLOAT_TYPE velocity_x;
				if (iy <= 0) {
					velocity_x = (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix + 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]) / 2.0;
				}
				else if (iy >= all_grid.Grid_num_y) {
					velocity_x = (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix + 1, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]) / 2.0;
				}
				else {
					velocity_x = (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix + 1, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix + 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]) / 4.0;
				}

				//バックトレース先の位置座標
				MY_FLOAT_TYPE advected_pos_x = (ix + 0.5) * all_grid._cell_length - (velocity_x * time_step_length);
				MY_FLOAT_TYPE advected_pos_y = iy * all_grid._cell_length - (all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] * time_step_length);
				velocity_after_advect_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= all_grid.calc_interpolated_velocity_y(
						VEC3_TYPE(advected_pos_x, advected_pos_y, 0.0),
						interpolation_method
					);
			}
		}
		//計算結果をコピー
		for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= velocity_after_advect_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
			}
		}
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
				all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= velocity_after_advect_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
			}
		}
*/
	}

	//advect項の計算
	//速度場を時間 -dt だけバックトレースしてadvect項を計算する
	void advect_velocity_semi_lagrangian_x(
		Grid& all_grid,
		const MY_FLOAT_TYPE time_step_length,
		const std::string interpolation_method,
		std::vector<MY_FLOAT_TYPE> &velocity_after_advect_x
	){
//		std::vector<MY_FLOAT_TYPE> velocity_after_advect_x((all_grid.Grid_num_x + 1) * all_grid.Grid_num_y);

		//velocityのx成分を計算
		for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				MY_FLOAT_TYPE velocity_y;
				if (ix <= 0) {
					velocity_y = (all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy + 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]) / 2.0;
				}
				else if (ix >= all_grid.Grid_num_x) {
					velocity_y = (all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix - 1, iy + 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]) / 2.0;
				}
				else {
					velocity_y = (all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix - 1, iy + 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy + 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]) / 4.0;
				}

				//バックトレース先の位置座標
				MY_FLOAT_TYPE advected_pos_x = ix * all_grid._cell_length - (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] * time_step_length);
				MY_FLOAT_TYPE advected_pos_y = (iy + 0.5) * all_grid._cell_length - (velocity_y * time_step_length);
				velocity_after_advect_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= all_grid.calc_interpolated_velocity_x(
						VEC3_TYPE(advected_pos_x, advected_pos_y, 0.0),
						interpolation_method
					);
			}
		}
/*
		//計算結果をコピー
		for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= velocity_after_advect_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
			}
		}
*/
	}

	//advect項の計算
	//速度場を時間 -dt だけバックトレースしてadvect項を計算する
	void advect_velocity_semi_lagrangian_y(
		Grid& all_grid,
		const MY_FLOAT_TYPE time_step_length,
		const std::string interpolation_method,
		std::vector<MY_FLOAT_TYPE> &velocity_after_advect_y
	){
//		std::vector<MY_FLOAT_TYPE> velocity_after_advect_y(all_grid.Grid_num_x * (all_grid.Grid_num_y + 1));

		//velocityのy成分を計算
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
				MY_FLOAT_TYPE velocity_x;
				if (iy <= 0) {
					velocity_x = (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix + 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]) / 2.0;
				}
				else if (iy >= all_grid.Grid_num_y) {
					velocity_x = (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix + 1, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]) / 2.0;
				}
				else {
					velocity_x = (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix + 1, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+ all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix + 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]) / 4.0;
				}

				//バックトレース先の位置座標
				MY_FLOAT_TYPE advected_pos_x = (ix + 0.5) * all_grid._cell_length - (velocity_x * time_step_length);
				MY_FLOAT_TYPE advected_pos_y = iy * all_grid._cell_length - (all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] * time_step_length);
				velocity_after_advect_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= all_grid.calc_interpolated_velocity_y(
						VEC3_TYPE(advected_pos_x, advected_pos_y, 0.0),
						interpolation_method
					);
			}
		}
/*
		//計算結果をコピー
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
				all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= velocity_after_advect_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
			}
		}
*/
	}

	// velocity の数値がnanを含むかをチェック
	void check_nan_in_velocity(const Grid& all_grid) {
		bool has_nan = false;
		for (int ix = 0; ix < all_grid.Grid_num_x + 1; ++ix) {
			for (int iy = 0; iy < all_grid.Grid_num_y; ++iy) {
				if (std::isnan(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])) {
					has_nan = true;
				}
			}
		}
		for (int ix = 0; ix < all_grid.Grid_num_x; ++ix) {
			for (int iy = 0; iy < all_grid.Grid_num_y + 1; ++iy) {
				if (std::isnan(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])) {
					has_nan = true;
				}
			}
		}
		if (has_nan) {
			std::cout << "!!!!!!! velocity field has nan !!!!!!!" << std::endl;
		}
	}

	//流体の 1 time step
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

	) {
		if(fix_velocity && !enable_cell_volume_correction){
			return;
		}
		if (i_frame == 0) {
			//セルの体積を計算する
			if(enable_cell_volume_correction){
				calc_all_cell_volumes(
					all_grid,
					time_step_length,
					use_MacCormack_scheme,
					interpolation_method
				);
			}
			calc_pressure_gradient_term(
				all_grid,
				enable_cell_volume_correction,
				fix_velocity,
				time_step_length
			);
			if(!fix_velocity){
				advect_velocity(
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
					use_integral,
					num_gauss_quadrature_point_for_integrate_density,
					enable_cell_volume_correction
				);
			}
		}
		else {
			if(!fix_velocity){
				advect_velocity(
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
					use_integral,
					num_gauss_quadrature_point_for_integrate_density,
					enable_cell_volume_correction
				);
			}
			//セルの体積を計算する
			if(enable_cell_volume_correction){
				calc_all_cell_volumes(
					all_grid,
					time_step_length,
					use_MacCormack_scheme,
					interpolation_method
				);
			}
			calc_pressure_gradient_term(
				all_grid,
				enable_cell_volume_correction,
				fix_velocity,
				time_step_length
			);
		}
		//getchar();
	}

	MY_FLOAT_TYPE integrate_normal_component_of_psi_velocity_x(
		const Grid& all_grid,
		cell_face face,
		const int num_gauss_quad_boundary,
		const int num_gauss_quad_bulk,
		const std::string split_method,
		const std::string interpolation_method,
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density
	) {
		//バックトレースした面の切断
		std::vector<cell_face> splitted_faces;
		std::string split_method_from_interpolation;
		if(interpolation_method == "linear" || interpolation_method == "WENO6"){
			split_method_from_interpolation = "y-AxisAligned";
		}
		else if(interpolation_method == "1Dy_linear" || interpolation_method == "1Dy_WENO6"){
			split_method_from_interpolation = "y-CellFaceAligned";
		}
		else{
			split_method_from_interpolation = split_method;
		}

//		splitted_faces = split_face(all_grid, face, split_method);
		splitted_faces = split_face(all_grid, face, split_method_from_interpolation);
		VEC3_TYPE area_weighted_psi_sum(0.0, 0.0, 0.0);
		// バックトレース後の面上で積分した値を足していく
		for (int i_face = 0; i_face < splitted_faces.size(); ++i_face) {
			int num_gauss_quad_points;
			if (face.is_boundary_face) {
				num_gauss_quad_points = num_gauss_quad_boundary;
			}
			else {
				num_gauss_quad_points = num_gauss_quad_bulk;
			}
			std::vector<MY_FLOAT_TYPE> gauss_quadrature_weights = gauss_quadrature_points_1D::get_quadtarure_weights_1D(num_gauss_quad_points);
			std::vector<VEC3_TYPE> gauss_quadrature_positions_on_splitted_face
				= gauss_quadrature_points_1D::calc_quadtarure_positions_on_cell_face(splitted_faces[i_face], num_gauss_quad_points);

			for (int i_quad = 0; i_quad < num_gauss_quad_points; ++i_quad) {
				//face face の中心でのpsiを計算する
				VEC3_TYPE quadrature_position_on_cell_face = gauss_quadrature_positions_on_splitted_face[i_quad];
				VEC3_TYPE psi_density_on_quadrature_position
					= all_grid.calc_psi_velocity_x_by_interpolation(
						quadrature_position_on_cell_face,
						interpolation_method,
						use_integral,
						num_gauss_quadrature_point_for_integrate_density
					);

//				std::cout << "splitted_faces[i_face].calc_face_area(): "<<splitted_faces[i_face].calc_face_area()<<std::endl;
//				std::cout << "gauss_quadrature_weights[i_quad]: "<<gauss_quadrature_weights[i_quad]<<std::endl;
//				std::cout << "splitted_faces[i_face].calc_face_area(): "<<psi_density_on_quadrature_position<<std::endl;
//				getchar();

				area_weighted_psi_sum
					+= splitted_faces[i_face].calc_face_area()
					* gauss_quadrature_weights[i_quad]
					* psi_density_on_quadrature_position;
			}
		}
//		std::cout << "area_weighted_psi_sum.dot(face.calc_normal()): "<<area_weighted_psi_sum.dot(face.calc_normal())<<std::endl;
		return area_weighted_psi_sum.dot(face.calc_normal());
	}

	MY_FLOAT_TYPE integrate_normal_component_of_psi_velocity_y(
		const Grid& all_grid,
		cell_face face,
		const int num_gauss_quad_boundary,
		const int num_gauss_quad_bulk,
		const std::string split_method,
		const std::string interpolation_method,
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density
	) {
		//バックトレースした面の切断
		std::vector<cell_face> splitted_faces;
		std::string split_method_from_interpolation;
		if(interpolation_method == "linear" || interpolation_method == "WENO6"){
			split_method_from_interpolation = "y-CellFaceAligned";
		}
		else if(interpolation_method == "1Dy_linear" || interpolation_method == "1Dy_WENO6"){
			split_method_from_interpolation = "y-AxisAligned";
		}
		else{
			split_method_from_interpolation = split_method;
		}

//		splitted_faces = split_face(all_grid, face, split_method);
		splitted_faces = split_face(all_grid, face, split_method_from_interpolation);
		VEC3_TYPE area_weighted_psi_sum(0.0, 0.0, 0.0);
		// バックトレース後の面上で積分した値を足していく
		for (int i_face = 0; i_face < splitted_faces.size(); ++i_face) {
			int num_gauss_quad_points;
			if (face.is_boundary_face) {
				num_gauss_quad_points = num_gauss_quad_boundary;
			}
			else {
				num_gauss_quad_points = num_gauss_quad_bulk;
			}
			std::vector<MY_FLOAT_TYPE> gauss_quadrature_weights = gauss_quadrature_points_1D::get_quadtarure_weights_1D(num_gauss_quad_points);
			std::vector<VEC3_TYPE> gauss_quadrature_positions_on_splitted_face
				= gauss_quadrature_points_1D::calc_quadtarure_positions_on_cell_face(splitted_faces[i_face], num_gauss_quad_points);

			for (int i_quad = 0; i_quad < num_gauss_quad_points; ++i_quad) {
				//face face の中心でのpsiを計算する
				VEC3_TYPE quadrature_position_on_cell_face = gauss_quadrature_positions_on_splitted_face[i_quad];
				VEC3_TYPE psi_density_on_quadrature_position
					= all_grid.calc_psi_velocity_y_by_interpolation(
						quadrature_position_on_cell_face,
						interpolation_method,
						use_integral,
						num_gauss_quadrature_point_for_integrate_density
					);

				area_weighted_psi_sum
					+= splitted_faces[i_face].calc_face_area()
					* gauss_quadrature_weights[i_quad]
					* psi_density_on_quadrature_position;
			}
		}
//		std::cout << "area_weighted_psi_sum.dot(face.calc_normal()): "<<area_weighted_psi_sum.dot(face.calc_normal())<<std::endl;
		return area_weighted_psi_sum.dot(face.calc_normal());
	}


	// density の MacCormack 法でのclamping の処理
	MY_FLOAT_TYPE clamp_integral_of_normal_component_of_psi_velocity_x_of_MacCormack(
		const Grid& all_grid,
		MY_FLOAT_TYPE integral_of_normal_component_of_psi_after_advect,
		VEC3_TYPE after_backtrace_position,
		const std::string interpolation_method) {
		// interpolation に使った値を記録する変数
		std::vector<MY_FLOAT_TYPE> interpolant_integral_of_normal_component_of_psi_values;
		if (interpolation_method == "linear") {
			// y軸に垂直な壁での積分
			int advected_index_x = floor((after_backtrace_position[0]) / all_grid._cell_length);
			int advected_index_y = floor((after_backtrace_position[1]) / all_grid._cell_length);
			if (advected_index_x < 0) {
				advected_index_x = 0;
			}
			if (advected_index_x >= all_grid.Grid_num_x) {
				advected_index_x = all_grid.Grid_num_x - 1;
			}
			if (advected_index_y <= 0) {
				advected_index_y = 0;
			}
			if (advected_index_y >= all_grid.Grid_num_y) {
				advected_index_y = all_grid.Grid_num_y - 1;
			}
			for (int ix = 0; ix < 2; ++ix) {
				for (int iy = 0; iy < 2; ++iy) {
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						all_grid.psi_velocity_cell_vertex_x[get_voxel_center_index(advected_index_x + ix, advected_index_y + iy, all_grid.Grid_num_x + 1, all_grid.Grid_num_y + 1)]
						* all_grid._cell_length
					);
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						-all_grid.psi_velocity_cell_vertex_x[get_voxel_center_index(advected_index_x + ix, advected_index_y + iy, all_grid.Grid_num_x + 1, all_grid.Grid_num_y + 1)]
						* all_grid._cell_length
					);
				}
			}
		}
		else if(interpolation_method == "1Dy_linear") {
			// y軸に垂直な壁での積分
			int advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
			int advected_index_y = floor((after_backtrace_position[1]) / all_grid._cell_length);
			if (advected_index_x < 0) {
				advected_index_x = 0;
			}
			if (advected_index_x >= all_grid.Grid_num_x) {
				advected_index_x = all_grid.Grid_num_x - 1;
			}
			if (advected_index_y <= 0) {
				advected_index_y = 0;
			}
			if (advected_index_y >= all_grid.Grid_num_y) {
				advected_index_y = all_grid.Grid_num_y - 1;
			}
			for (int iy = 0; iy < 2; ++iy) {
				interpolant_integral_of_normal_component_of_psi_values.push_back(
					all_grid.psi_velocity_cell_vertex_x[get_voxel_center_index(advected_index_x, advected_index_y + iy, all_grid.Grid_num_x + 1, all_grid.Grid_num_y + 1)]
					* all_grid._cell_length
				);
				interpolant_integral_of_normal_component_of_psi_values.push_back(
					-all_grid.psi_velocity_cell_vertex_x[get_voxel_center_index(advected_index_x, advected_index_y + iy, all_grid.Grid_num_x + 1, all_grid.Grid_num_y + 1)]
					* all_grid._cell_length
				);
			}
		}
/*
		else if (interpolation_method == "WENO6") {
			// x軸に垂直な壁での積分
			int advected_index_x = floor((after_backtrace_position[0]) / all_grid._cell_length);
			int advected_index_y = floor((after_backtrace_position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
			for (int ix = 0; ix < 6; ix++) {
				for (int iy = 0; iy < 6; ++iy) {
					int index_x = advected_index_x - 2 + ix;
					int index_y = advected_index_y - 2 + iy;
					if (index_x < 0) {
						index_x = 0;
					}
					if (index_x > all_grid.Grid_num_x) {
						int exceed = index_x - all_grid.Grid_num_x;
						//境界の値が外側までずっと続く場合
						index_x = all_grid.Grid_num_x;
						//境界を境に鏡のように値が反射する場合
						//index_x = all_grid.Grid_num_x- exceed + 1;
					}
					if (index_y < 0) {
						index_y = 0;
					}
					if (index_y > all_grid.Grid_num_y - 1) {
						int exceed = index_y - (all_grid.Grid_num_y - 1);
						//境界の値が外側までずっと続く場合
						index_y = all_grid.Grid_num_y - 1;
						//境界を境に鏡のように値が反射する場合
						//index_y = all_grid.Grid_num_y - 1 - exceed + 1;
					}
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						all_grid.psi_substance_density_cell_face_x[get_voxel_face_index_x(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						* all_grid._cell_length
					);
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						-all_grid.psi_substance_density_cell_face_x[get_voxel_face_index_x(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						* all_grid._cell_length
					);
				}
			}
			// y軸に垂直な壁での積分
			advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
			advected_index_y = floor((after_backtrace_position[1]) / all_grid._cell_length);
			for (int ix = 0; ix < 6; ix++) {
				for (int iy = 0; iy < 6; ++iy) {
					int index_x = advected_index_x - 2 + ix;
					int index_y = advected_index_y - 2 + iy;
					if (index_x < 0) {
						index_x = 0;
					}
					if (index_x > all_grid.Grid_num_x - 1) {
						int exceed = index_x - (all_grid.Grid_num_x - 1);
						//境界の値が外側までずっと続く場合
						index_x = all_grid.Grid_num_x - 1;
						//境界を境に鏡のように値が反射する場合
						//index_x = all_grid.Grid_num_x - 1 - exceed + 1;
					}
					if (index_y < 0) {
						index_y = 0;
					}
					if (index_y > all_grid.Grid_num_y) {
						int exceed = index_y - all_grid.Grid_num_y;
						//境界の値が外側までずっと続く場合
						index_y = all_grid.Grid_num_y;
						//境界を境に鏡のように値が反射する場合
						//index_y = all_grid.Grid_num_y - exceed + 1;
					}
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						* all_grid._cell_length
					);
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						-all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						* all_grid._cell_length
					);
				}
			}
		}
*/
		// clamping の処理
		std::sort(interpolant_integral_of_normal_component_of_psi_values.begin(), interpolant_integral_of_normal_component_of_psi_values.end());
		integral_of_normal_component_of_psi_after_advect = std::min(integral_of_normal_component_of_psi_after_advect, interpolant_integral_of_normal_component_of_psi_values[interpolant_integral_of_normal_component_of_psi_values.size() - 1]);
		integral_of_normal_component_of_psi_after_advect = std::max(integral_of_normal_component_of_psi_after_advect, interpolant_integral_of_normal_component_of_psi_values[0]);
		return integral_of_normal_component_of_psi_after_advect;
	}

	// density の MacCormack 法でのclamping の処理
	MY_FLOAT_TYPE clamp_integral_of_normal_component_of_psi_velocity_y_of_MacCormack(
		const Grid& all_grid,
		MY_FLOAT_TYPE integral_of_normal_component_of_psi_after_advect,
		VEC3_TYPE after_backtrace_position,
		const std::string interpolation_method) {
		// interpolation に使った値を記録する変数
		std::vector<MY_FLOAT_TYPE> interpolant_integral_of_normal_component_of_psi_values;
		if (interpolation_method == "linear") {
			// y軸に垂直な壁での積分
			int advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
			int advected_index_y = floor((after_backtrace_position[1] + 0.5 * all_grid._cell_length) / all_grid._cell_length);
			if (advected_index_x < 0) {
				advected_index_x = 0;
			}
			if (advected_index_x > all_grid.Grid_num_x - 1) {
				advected_index_x = all_grid.Grid_num_x - 1;
			}
			if (advected_index_y <= 0) {
				advected_index_y = 0;
			}
			if (advected_index_y > all_grid.Grid_num_y + 1) {
				advected_index_y = all_grid.Grid_num_y;
			}
			for (int ix = 0; ix < 2; ++ix) {
				for (int iy = 0; iy < 2; ++iy) {
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						all_grid.psi_velocity_cell_center_y[get_voxel_center_index(advected_index_x + ix, advected_index_y + iy, all_grid.Grid_num_x, all_grid.Grid_num_y + 2)]
						* all_grid._cell_length
					);
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						-all_grid.psi_velocity_cell_center_y[get_voxel_center_index(advected_index_x + ix, advected_index_y + iy, all_grid.Grid_num_x, all_grid.Grid_num_y + 2)]
						* all_grid._cell_length
					);
				}
			}
		}
		else if(interpolation_method == "1Dy_linear") {
			// y軸に垂直な壁での積分
			int advected_index_x = floor((after_backtrace_position[0]) / all_grid._cell_length);
			int advected_index_y = floor((after_backtrace_position[1] + 0.5 * all_grid._cell_length) / all_grid._cell_length);
			if (advected_index_x < 0) {
				advected_index_x = 0;
			}
			if (advected_index_x > all_grid.Grid_num_x - 1) {
				advected_index_x = all_grid.Grid_num_x - 1;
			}
			if (advected_index_y <= 0) {
				advected_index_y = 0;
			}
			if (advected_index_y > all_grid.Grid_num_y + 1) {
				advected_index_y = all_grid.Grid_num_y;
			}
			for (int iy = 0; iy < 2; ++iy) {
				interpolant_integral_of_normal_component_of_psi_values.push_back(
					all_grid.psi_velocity_cell_center_y[get_voxel_center_index(advected_index_x, advected_index_y + iy, all_grid.Grid_num_x, all_grid.Grid_num_y + 2)]
					* all_grid._cell_length
				);
				interpolant_integral_of_normal_component_of_psi_values.push_back(
					-all_grid.psi_velocity_cell_center_y[get_voxel_center_index(advected_index_x, advected_index_y + iy, all_grid.Grid_num_x, all_grid.Grid_num_y + 2)]
					* all_grid._cell_length
				);
			}
		}
/*
		else if (interpolation_method == "WENO6") {
			// x軸に垂直な壁での積分
			int advected_index_x = floor((after_backtrace_position[0]) / all_grid._cell_length);
			int advected_index_y = floor((after_backtrace_position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
			for (int ix = 0; ix < 6; ix++) {
				for (int iy = 0; iy < 6; ++iy) {
					int index_x = advected_index_x - 2 + ix;
					int index_y = advected_index_y - 2 + iy;
					if (index_x < 0) {
						index_x = 0;
					}
					if (index_x > all_grid.Grid_num_x) {
						int exceed = index_x - all_grid.Grid_num_x;
						//境界の値が外側までずっと続く場合
						index_x = all_grid.Grid_num_x;
						//境界を境に鏡のように値が反射する場合
						//index_x = all_grid.Grid_num_x- exceed + 1;
					}
					if (index_y < 0) {
						index_y = 0;
					}
					if (index_y > all_grid.Grid_num_y - 1) {
						int exceed = index_y - (all_grid.Grid_num_y - 1);
						//境界の値が外側までずっと続く場合
						index_y = all_grid.Grid_num_y - 1;
						//境界を境に鏡のように値が反射する場合
						//index_y = all_grid.Grid_num_y - 1 - exceed + 1;
					}
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						all_grid.psi_substance_density_cell_face_x[get_voxel_face_index_x(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						* all_grid._cell_length
					);
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						-all_grid.psi_substance_density_cell_face_x[get_voxel_face_index_x(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						* all_grid._cell_length
					);
				}
			}
			// y軸に垂直な壁での積分
			advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
			advected_index_y = floor((after_backtrace_position[1]) / all_grid._cell_length);
			for (int ix = 0; ix < 6; ix++) {
				for (int iy = 0; iy < 6; ++iy) {
					int index_x = advected_index_x - 2 + ix;
					int index_y = advected_index_y - 2 + iy;
					if (index_x < 0) {
						index_x = 0;
					}
					if (index_x > all_grid.Grid_num_x - 1) {
						int exceed = index_x - (all_grid.Grid_num_x - 1);
						//境界の値が外側までずっと続く場合
						index_x = all_grid.Grid_num_x - 1;
						//境界を境に鏡のように値が反射する場合
						//index_x = all_grid.Grid_num_x - 1 - exceed + 1;
					}
					if (index_y < 0) {
						index_y = 0;
					}
					if (index_y > all_grid.Grid_num_y) {
						int exceed = index_y - all_grid.Grid_num_y;
						//境界の値が外側までずっと続く場合
						index_y = all_grid.Grid_num_y;
						//境界を境に鏡のように値が反射する場合
						//index_y = all_grid.Grid_num_y - exceed + 1;
					}
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						* all_grid._cell_length
					);
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						-all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						* all_grid._cell_length
					);
				}
			}
		}
*/
		// clamping の処理
		std::sort(interpolant_integral_of_normal_component_of_psi_values.begin(), interpolant_integral_of_normal_component_of_psi_values.end());
		integral_of_normal_component_of_psi_after_advect = std::min(integral_of_normal_component_of_psi_after_advect, interpolant_integral_of_normal_component_of_psi_values[interpolant_integral_of_normal_component_of_psi_values.size() - 1]);
		integral_of_normal_component_of_psi_after_advect = std::max(integral_of_normal_component_of_psi_after_advect, interpolant_integral_of_normal_component_of_psi_values[0]);
		return integral_of_normal_component_of_psi_after_advect;
	}

}
