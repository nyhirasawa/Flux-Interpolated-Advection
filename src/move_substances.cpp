#include "move_substances.h"

#include <iostream>
#include <algorithm>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/bicgstabl.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

#include "linear_solver.h"
#include "physical_const.h"
#include "utils.h"
#include "cell_vertex.h"
#include "cell_face.h"
#include "gauss_quadrature_points.h"
#include "define_float_type.h"
#include "calc_psi.h"
#include "split_face.h"
#include "calc_backtraced_face.h"
#include "calc_cell_volumes.h"

namespace smoke_simulation {
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
	) {
		std::vector<MY_FLOAT_TYPE> before_advect_density = all_grid.substance_density;
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
			advect_density_flux_advection(
				all_grid,
				time_step_length,
//				set_zero_normal_velocity_at_boundary,
				use_MacCormack_scheme,
				use_clamping_in_MacCormack_scheme,
				num_gauss_quad_boundary,
				num_gauss_quad_bulk,
				enable_eliminate_zero_velocity_false_diffusion,
				split_method,
			    integral_method,
				all_grid.substance_density,
				interpolation_method,
				interpolation_method_in_calc_psi,
				use_integral,
				num_gauss_quadrature_point_for_integrate_density,
				enable_cell_volume_correction
			);
			// セルの体積によって質量密度を補正
			if(enable_cell_volume_correction){
//				correct_mass_and_volume_by_volume_gradient_diffusion(
//					all_grid,
//					10000
//				);
				correct_mass_and_volume_by_pressure_solve(
					all_grid.substance_density,
					all_grid.cell_volume_cell_center,
					all_grid
				);
			}
			if (enable_eliminate_zero_velocity_false_diffusion) {
				// 各セル中心での密度の補正量
	            std::vector<MY_FLOAT_TYPE> zero_velocity_correction_values(all_grid.Grid_num_x * all_grid.Grid_num_y);
				////// "バックトレース前の面で定義されたpsiの離散値"を足す
				for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
					for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
						zero_velocity_correction_values[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
							= before_advect_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
				 	}
				}
				for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
					for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
						//(ix, iy)番目のセルのface を走るループ
						int dx[4], dy[4];
						dx[0] = 0; dx[1] = 1; dx[2] = 1; dx[3] = 0;
						dy[0] = 0; dy[1] = 0; dy[2] = 1; dy[3] = 1;
						for (int i_face = 0; i_face < 4; ++i_face) {
							// 考えるcell face を構成するvertex のindex
							std::vector<int> vertex_index_1{ ix + dx[i_face], iy + dy[i_face] }, vertex_index_2{ ix + dx[(i_face + 1) % 4], iy + dy[(i_face + 1) % 4] };
							// バックトレース前の cell face を構成するvertexを定義する
							cell_vertex before_backtrace_vertex_1(vertex_index_1, all_grid._cell_length);
							cell_vertex before_backtrace_vertex_2(vertex_index_2, all_grid._cell_length);
							// バックトレース前の頂点から面を定義
							cell_face before_backtrace_face;
							before_backtrace_face._vertex_list[0] = before_backtrace_vertex_1;
							before_backtrace_face._vertex_list[1] = before_backtrace_vertex_2;
							////// エラーの補正
							////// "バックトレース前の面上で積分した値" と "バックトレース前の面で定義されたpsiの離散値"の差を引いていく
							zero_velocity_correction_values[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
								-= integrate_normal_component_of_psi_on_face(
									all_grid,
									before_backtrace_face,
									num_gauss_quad_boundary,
									num_gauss_quad_bulk,
									split_method,
									interpolation_method,
									use_integral,
									num_gauss_quadrature_point_for_integrate_density
									) / (all_grid._cell_length * all_grid._cell_length);
						}
					}
				}
				// エラーの移流と密度の移流で同じ補間の方法を使う
//	            std::string interpolation_method_of_eccor_advection = interpolation_method;

	            std::string interpolation_method_of_eccor_advection;
	            if(interpolation_method == "linear" || interpolation_method == "1Dy_linear"){
					interpolation_method_of_eccor_advection = "1Dy_linear";
//					interpolation_method_of_eccor_advection = "linear";
	            }
	            else if(interpolation_method == "WENO6"||interpolation_method == "WENO6-optimized"){
	                interpolation_method_of_eccor_advection = "1Dy_WENO6";
	            }
				std::string split_method_of_eccor_advection = "y-AxisAligned";
				std::string interpolation_method_in_calc_psi = "const";
				// エラーを移流させる
				advect_density_flux_advection(
					all_grid,
					time_step_length,
//					set_zero_normal_velocity_at_boundary,
					use_MacCormack_scheme,
					use_clamping_in_MacCormack_scheme,
					num_gauss_quad_boundary,
					num_gauss_quad_bulk,
					enable_eliminate_zero_velocity_false_diffusion,
					split_method_of_eccor_advection,
				    integral_method,
					zero_velocity_correction_values,
					interpolation_method_of_eccor_advection,
					interpolation_method_in_calc_psi,
					use_integral,
					num_gauss_quadrature_point_for_integrate_density,
					/*enable_cell_volume_correction*/false
				);

				//密度を補正
				for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
					for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
						all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
							+=zero_velocity_correction_values[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
					}
				}
			}
		}
		else {
			if(use_lentines_advection){
				advect_density_Lentine(
					all_grid,
					time_step_length,
					enable_cell_volume_correction
				);
			}
			else{
				advect_density_semi_lagrangian(
					all_grid,
					time_step_length,
					use_MacCormack_scheme,
					use_clamping_in_MacCormack_scheme,
					num_gauss_quad_boundary,
					num_gauss_quad_bulk,
					interpolation_method
				);
			}
		}
		//check_negative_dinsity(all_grid);
	}

	void correct_mass_and_volume_by_volume_gradient_diffusion(
		Grid& all_grid,
		const int iteration_num
	){
		for(int _ = 0; _ < iteration_num; ++_){
			// 補正前の質量と体積
			std::vector<MY_FLOAT_TYPE> before_correct_mass(all_grid.Grid_num_x * all_grid.Grid_num_y);
			std::vector<MY_FLOAT_TYPE> before_correct_volume(all_grid.Grid_num_x * all_grid.Grid_num_y);
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					before_correct_mass[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						= all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						* all_grid._cell_volume;
//						* all_grid.cell_volume_cell_center[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
					before_correct_volume[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						= all_grid.cell_volume_cell_center[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
				}
			}
			//// 体積の勾配によって質量と体積を拡散させる
			MY_FLOAT_TYPE diffusion_coefficient = 0.25;
			// x 方向の拡散
			for (int ix = 1; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					//フラックス
					MY_FLOAT_TYPE volume_flux
						= -(before_correct_volume[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						  - before_correct_volume[get_voxel_center_index(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
						 / all_grid._cell_volume;
					// 面での密度
					MY_FLOAT_TYPE density_at_face
 						= (before_correct_mass[get_voxel_center_index(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
 						+  before_correct_mass[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
 						/ (2.0 * all_grid._cell_volume);
					// 面での体積
					MY_FLOAT_TYPE volume_at_face
						= (before_correct_volume[get_voxel_center_index(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+  before_correct_volume[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
						/ 2.0;
//					if(volume_flux > 0.0){
//						density_at_face = before_correct_mass[get_voxel_center_index(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
//						volume_at_face = before_correct_volume[get_voxel_center_index(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
//					}
//					else{
//						density_at_face = before_correct_mass[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
//						volume_at_face = before_correct_volume[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
//					}
					//密度の拡散
					all_grid.substance_density[get_voxel_center_index(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						-= diffusion_coefficient * volume_flux * density_at_face;
					all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+= diffusion_coefficient * volume_flux * density_at_face;
					//体積の拡散
					all_grid.cell_volume_cell_center[get_voxel_center_index(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						-= diffusion_coefficient * volume_flux * volume_at_face;
					all_grid.cell_volume_cell_center[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+= diffusion_coefficient * volume_flux * volume_at_face;
				}
			}
			// y 方向の拡散
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 1; iy < all_grid.Grid_num_y; iy++) {
					//フラックス
					MY_FLOAT_TYPE volume_flux
						= -(before_correct_volume[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						  - before_correct_volume[get_voxel_center_index(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)])
						 / all_grid._cell_volume;
					// 面での密度
					MY_FLOAT_TYPE density_at_face
 						= (before_correct_mass[get_voxel_center_index(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
 						+  before_correct_mass[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
 						/ (2.0 * all_grid._cell_volume);
					// 面での体積
					MY_FLOAT_TYPE volume_at_face
						= (before_correct_volume[get_voxel_center_index(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+  before_correct_volume[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
						/ 2.0;
//					if(volume_flux > 0.0){
//						density_at_face = before_correct_mass[get_voxel_center_index(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)];
//						volume_at_face = before_correct_volume[get_voxel_center_index(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)];
//					}
//					else{
//						density_at_face = before_correct_mass[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
//						volume_at_face = before_correct_volume[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
//					}
					//密度の拡散
					all_grid.substance_density[get_voxel_center_index(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						-= diffusion_coefficient * volume_flux * density_at_face;
					all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+= diffusion_coefficient * volume_flux * density_at_face;
					//体積の拡散
					all_grid.cell_volume_cell_center[get_voxel_center_index(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						-= diffusion_coefficient * volume_flux * volume_at_face;
					all_grid.cell_volume_cell_center[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+= diffusion_coefficient * volume_flux * volume_at_face;
				}
			}
		}
	}

	void correct_mass_and_volume_by_pressure_solve(
		std::vector<MY_FLOAT_TYPE> &corrected_values,
		const std::vector<MY_FLOAT_TYPE> &cell_volume,
		const Grid& all_grid
	){
		// 圧力計算
		//グリッドの総数
		const int N = all_grid.Grid_num_x
					* all_grid.Grid_num_y;
		//係数行列の計算
		//ノイマン境界条件では係数行列をナイーブに構築するとフルランクにならない(1つランクが落ちる)
		//フルランクにするために圧力の最後の要素を0に固定する。それに合わせて係数行列も最後の行, 列を落とす
		linear_algebra::sparse_matrix A(N - 1, N - 1);

		for(int i_xy = 0; i_xy < N - 1; ++i_xy){
			int iy = i_xy % all_grid.Grid_num_y;
			int ix = (i_xy - iy) / all_grid.Grid_num_y;
			//非対角項
			if (ix != 0) {
				int index_0=get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				int index_1=get_voxel_center_index(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				if((index_0 != N - 1) && (index_1 != N - 1)){
					A.input_element(index_0, index_1, 1);
				}
			}
			//非対角項
			if (iy != 0) {
				int index_0=get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				int index_1=get_voxel_center_index(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y);
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
			if(ix == all_grid.Grid_num_x - 1){
				diagonal_value += 1;
			}
			A.input_element(
				get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y),
				get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y),
				diagonal_value
			);
			//非対角項
			if (iy != all_grid.Grid_num_y - 1) {
				int index_0=get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				int index_1=get_voxel_center_index(ix, iy + 1, all_grid.Grid_num_x, all_grid.Grid_num_y);
				if((index_0 != N - 1) && (index_1 != N - 1)){
					A.input_element(index_0, index_1, 1);
				}
			}
			//非対角項
			if (ix != all_grid.Grid_num_x - 1) {
				int index_0=get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
				int index_1=get_voxel_center_index(ix + 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y);
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

			b[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
				= -(cell_volume[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
				- all_grid._cell_volume);
		}

		//CG法により圧力場を得る
//		std::vector<MY_FLOAT_TYPE> truncated_pressure(N - 1);
//		std::cout<<"linear solve begin"<<std::endl;
		//linear_algebra::conjugate_gradient(A, b, all_grid.pressure, N, 10000, 0.0001);
//		linear_algebra::incomplete_cholesky_conjugate_gradient(A, b, truncated_pressure, N - 1, 10000, 0.0001);
		//gauss seidel法を使う場合(Aはsparse_matrix_with_diagonal_elementにする)
//		gauss_seidel(A, b, truncated_pressure, N - 1, 2000);
		//係数行列をフルランクにするために落とした最後の要素を復元する
//		truncated_pressure.push_back(0.0);
//		all_grid.pressure = truncated_pressure;

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
	    	amgcl::solver::bicgstabl<Backend>
//			amgcl::solver::cg<Backend>
			> Solver;
//		int n = N;
		int n = N - 1;
//		int n = all_grid.Grid_num_x * all_grid.Grid_num_y;

		Solver::params prm;
		prm.solver.maxiter = 2;
//		prm.solver.tol = 1e-1;
		Solver solve( std::tie(n, A.cumulative_num_nonzero_element_in_row, A.column_index, A.element_value) ,prm);
//		Solver solve( std::tie(n, A.cumulative_num_nonzero_element_in_row, A.column_index, A.element_value));

		int    iters;
		MY_FLOAT_TYPE error;
		std::vector<MY_FLOAT_TYPE> pressure(N - 1);
		std::tie(iters, error) = solve(b, pressure);
		//係数行列をフルランクにするために落とした最後の要素を復元する
		pressure.push_back(0.0);

		// 補正前の質量と体積
		std::vector<MY_FLOAT_TYPE> before_correct_mass(all_grid.Grid_num_x * all_grid.Grid_num_y);
		std::vector<MY_FLOAT_TYPE> before_correct_volume(all_grid.Grid_num_x * all_grid.Grid_num_y);
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				before_correct_mass[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= corrected_values[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					* all_grid._cell_volume;
				before_correct_volume[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= cell_volume[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
			}
		}
		//// 体積の勾配によって質量と体積を拡散させる
		MY_FLOAT_TYPE diffusion_coefficient = 1.0;
		// x 方向の拡散
		for (int ix = 1; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				//フラックス
				MY_FLOAT_TYPE volume_flux
					= -(pressure[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					  - pressure[get_voxel_center_index(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					 / all_grid._cell_volume;
				// 面での密度
				MY_FLOAT_TYPE density_at_face
						= (before_correct_mass[get_voxel_center_index(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+  before_correct_mass[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
						/ (2.0 * all_grid._cell_volume);
				// 面での体積
				MY_FLOAT_TYPE volume_at_face
					= (before_correct_volume[get_voxel_center_index(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+  before_correct_volume[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					/ 2.0;
				//密度の拡散
				corrected_values[get_voxel_center_index(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					-= diffusion_coefficient * volume_flux * density_at_face;
				corrected_values[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+= diffusion_coefficient * volume_flux * density_at_face;
				//体積の拡散
//				cell_volume[get_voxel_center_index(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//					-= diffusion_coefficient * volume_flux * volume_at_face;
//				cell_volume[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//					+= diffusion_coefficient * volume_flux * volume_at_face;
			}
		}
		// y 方向の拡散
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 1; iy < all_grid.Grid_num_y; iy++) {
				//フラックス
				MY_FLOAT_TYPE volume_flux
					= -(pressure[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					  - pressure[get_voxel_center_index(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					 / all_grid._cell_volume;
				// 面での密度
				MY_FLOAT_TYPE density_at_face
					= (before_correct_mass[get_voxel_center_index(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+  before_correct_mass[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					/ (2.0 * all_grid._cell_volume);
				// 面での体積
				MY_FLOAT_TYPE volume_at_face
					= (before_correct_volume[get_voxel_center_index(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+  before_correct_volume[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					/ 2.0;
				//密度の拡散
				corrected_values[get_voxel_center_index(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					-= diffusion_coefficient * volume_flux * density_at_face;
				corrected_values[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+= diffusion_coefficient * volume_flux * density_at_face;
				//体積の拡散
//				cell_volume[get_voxel_center_index(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//					-= diffusion_coefficient * volume_flux * volume_at_face;
//				cell_volume[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//					+= diffusion_coefficient * volume_flux * volume_at_face;
			}
		}
	}

	//密度場の移流項をsemi-Lagrangian で計算
	void advect_density_semi_lagrangian(
		Grid& all_grid,
		const MY_FLOAT_TYPE time_step_length,
		const bool use_MacCormack_scheme,
		const bool use_clamping_in_MacCormack_scheme,
		const int num_gauss_quad_boundary,
		const int num_gauss_quad_bulk,
		const std::string interpolation_method
	) {
		MY_FLOAT_TYPE substance_density_after_advect[all_grid.Grid_num_x][all_grid.Grid_num_y];
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				// MacCormack 法を使う場合
				if (use_MacCormack_scheme) {
					//バックトレース前の座標 ( 0.5 は cell vertex から cell center までのオフセット )
					VEC3_TYPE before_backtrace_position((ix + 0.5) * all_grid._cell_length, (iy + 0.5) * all_grid._cell_length, 0.0);
					//バックトレース前の座標での速度の計算
					VEC3_TYPE velocity_at_before_backtrace_position
						= all_grid.calc_interpolated_velocity(
							before_backtrace_position,
							interpolation_method);
					//バックトレース後の座標 ( 0.5 は cell vertex から cell center までのオフセット )
					VEC3_TYPE after_backtrace_position = before_backtrace_position - time_step_length * velocity_at_before_backtrace_position;
					// バックトレース先の密度を補間で求める
					MY_FLOAT_TYPE density_at_advected_position
						= all_grid.calc_substance_density_by_interpolation(
							after_backtrace_position,
							interpolation_method
						);
					////// エラーの計算
					//バックトレース後の座標での速度
					VEC3_TYPE velocity_at_after_backtrace_position
						= all_grid.calc_interpolated_velocity(
							after_backtrace_position,
							interpolation_method);
					//バックトレース後の座標を時間を巻き戻すように逆向きにトレースする
					VEC3_TYPE forwardtrace_position_of_after_backtrace_position = after_backtrace_position + time_step_length * velocity_at_after_backtrace_position;
					MY_FLOAT_TYPE error
						= (all_grid.calc_substance_density_by_interpolation(before_backtrace_position, interpolation_method)
							- all_grid.calc_substance_density_by_interpolation(forwardtrace_position_of_after_backtrace_position, interpolation_method))
						/ 2.0;
					// エラー補正後の値
					substance_density_after_advect[ix][iy] = density_at_advected_position + error;
					// clamping
					if (use_clamping_in_MacCormack_scheme) {
						substance_density_after_advect[ix][iy]
							= clamp_substance_density_of_MacCormack(
								all_grid,
								substance_density_after_advect[ix][iy],
								after_backtrace_position,
								interpolation_method);
					}
				}
				// semi-Lagrangian 法を使う場合
				else {
					MY_FLOAT_TYPE velocity_x = (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] + all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix + 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]) / 2.0;
					MY_FLOAT_TYPE velocity_y = (all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] + all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy + 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]) / 2.0;
					//バックトレース先の座標
					// 0.5 は cell vertex から cell center までのオフセット
					MY_FLOAT_TYPE advected_pos_x = ((ix + 0.5) * all_grid._cell_length) - (velocity_x * time_step_length);
					MY_FLOAT_TYPE advected_pos_y = ((iy + 0.5) * all_grid._cell_length) - (velocity_y * time_step_length);
					//バックトレース先の速度を補間する
					substance_density_after_advect[ix][iy] = all_grid.calc_substance_density_by_interpolation(VEC3_TYPE(advected_pos_x, advected_pos_y, 0.0), interpolation_method);
				}
			}
		}
		//計算結果をコピー
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = substance_density_after_advect[ix][iy];
			}
		}
	}

	static void get_ij_xy( int Gx, int Gy, float dx, const VEC3_TYPE &p,
		int &i, int &j, float &x, float &y) {
		//
		x = (p.x() - 0.5 * dx)/dx;
		y = (p.y() - 0.5 * dx)/dx;
		//
		if( x < 0.0 ) x = 0.0;
		if( x > Gx-1 ) x = Gx-1;
		//
		if( y < 0.0 ) y = 0.0;
		if( y > Gy-1 ) y = Gy-1;
		//
		i = x;
		j = y;
		//
		if( i < 0 ) i = 0;
		if( j < 0 ) j = 0;
		//
		if( i > Gx-2 ) i = Gx-2;
		if( j > Gy-2 ) j = Gy-2;
	}

	VEC3_TYPE advect_pos_linear(
		const VEC3_TYPE pos,
		const std::string trace_mode, // "backtrace" or "forwardtrace",
		const float time_step_length,
		const Grid& all_grid
	) {
		//速度を計算
		VEC3_TYPE velocity = all_grid.calc_interpolated_velocity(
			pos,
			"linear"
		);
		if( trace_mode == "backtrace" ) {
			return pos - time_step_length * velocity;
		}
		else if ( trace_mode == "forwardtrace" ) {
			return pos + time_step_length * velocity;
		}
	}

	void get_bilinear_matrix(
		const int num_grid_cell_x,
		const int num_grid_cell_y,
		const MY_FLOAT_TYPE grid_cell_length,
		const VEC3_TYPE p,
		size_t cell_index_for_interpolation[4],
		float interpolation_weights[4]
	) {
		int i, j, k;
		float x, y, z;
		get_ij_xy(num_grid_cell_x, num_grid_cell_y, grid_cell_length, p, i, j, x, y);
		const size_t idx[] = {
			get_voxel_center_index(i, j, num_grid_cell_x, num_grid_cell_y),
			get_voxel_center_index(i+1, j, num_grid_cell_x, num_grid_cell_y),
			get_voxel_center_index(i+1, j+1, num_grid_cell_x, num_grid_cell_y),
			get_voxel_center_index(i, j+1, num_grid_cell_x, num_grid_cell_y)
		};
		//
		const float xx = x-i;
		const float yy = y-j;
		//
		const float t[] = {
			(1.0f-xx)*(1.0f-yy),
			xx*(1.0f-yy),
			xx*yy,
			(1.0f-xx)*yy
		};
		//
		for( int n=0; n<4; ++n ) {
			cell_index_for_interpolation[n] = idx[n];
			interpolation_weights[n] = t[n];
		}
	}

	//密度場の移流項をsemi-Lagrangian で計算
	void advect_density_Lentine(
		Grid& all_grid,
		const MY_FLOAT_TYPE time_step_length,
		const bool enable_cell_volume_correction
	) {
		const size_t num_grid_cell = all_grid.Grid_num_x * all_grid.Grid_num_y;
		//// weight_sum を計算する
		// weight_sum[i] i番目のセルから出る重みの和
		std::vector<MY_FLOAT_TYPE>weight_sum(num_grid_cell, 0.0);
		for(int i = 0; i < all_grid.Grid_num_x; ++i) {
			for(int j = 0; j < all_grid.Grid_num_y; ++j) {
				VEC3_TYPE advected_pos = advect_pos_linear(
					VEC3_TYPE(i + 0.5, j + 0.5, 0.0) * all_grid._cell_length,
					"backtrace",
					time_step_length,
					all_grid
				);

				size_t cell_index_for_interpolation[4];
				float interpolation_weights[4];
				get_bilinear_matrix(
					all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid._cell_length,
					advected_pos,
					cell_index_for_interpolation, interpolation_weights
				);
				for(int m = 0; m < 4; ++m){
					weight_sum[cell_index_for_interpolation[m]]
						+= interpolation_weights[m];
				}
			}
		}

		std::vector<float> cell_volume(num_grid_cell,0.0f);	// volume[i]: i番目のセルの体積
		std::vector<float> advected_value_results; // tmp[i]: i番目のセルの移流後の結果を格納する
		advected_value_results.resize(num_grid_cell);

		////weight_sum > 1.0 の時の処理
		for(int i = 0; i < all_grid.Grid_num_x; ++i) {
			for(int j = 0; j < all_grid.Grid_num_y; ++j) {
				VEC3_TYPE advected_pos = advect_pos_linear(
					VEC3_TYPE(i + 0.5, j + 0.5, 0.0) * all_grid._cell_length,
					"backtrace",
					time_step_length,
					all_grid
				);

				size_t cell_index_for_interpolation[4];
				float interpolation_weights[4];
				get_bilinear_matrix(
					all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid._cell_length,
					advected_pos,
					cell_index_for_interpolation, interpolation_weights
				);

				float sum = 0.0;
				for(int m = 0; m < 4; ++m){
					const float w = std::max(1.0f, weight_sum[cell_index_for_interpolation[m]]);
					sum += interpolation_weights[m] * all_grid.substance_density[cell_index_for_interpolation[m]] / w;
					cell_volume[get_voxel_center_index(i, j, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						+= interpolation_weights[m] / w;
				}

				advected_value_results[get_voxel_center_index(i, j, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= sum;
			}
		}

		////weight_sum < 1.0 の時の処理
		for(int i = 0; i < all_grid.Grid_num_x; ++i) {
			for(int j = 0; j < all_grid.Grid_num_y; ++j) {
				if( weight_sum[get_voxel_center_index(i, j, all_grid.Grid_num_x, all_grid.Grid_num_y)] < 1.0f ) {
					VEC3_TYPE forward_advected_pos = advect_pos_linear(
						VEC3_TYPE(i + 0.5, j + 0.5, 0.0) * all_grid._cell_length,
						"forwardtrace",
						time_step_length,
						all_grid
					);

					size_t cell_index_for_interpolation[4];
					float interpolation_weights[4];
					get_bilinear_matrix(
						all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid._cell_length,
						forward_advected_pos,
						cell_index_for_interpolation, interpolation_weights
					);

					const float w = 1.0 - weight_sum[get_voxel_center_index(i, j, all_grid.Grid_num_x, all_grid.Grid_num_y)];
					for(int m = 0; m < 4; ++m){
						advected_value_results[cell_index_for_interpolation[m]]
							+= w
							* interpolation_weights[m]
							* all_grid.substance_density[get_voxel_center_index(i, j, all_grid.Grid_num_x, all_grid.Grid_num_y)];
						cell_volume[cell_index_for_interpolation[m]]
							+= w * interpolation_weights[m];
					}
				}
			}
		}
		for(size_t n = 0; n < num_grid_cell; ++n){
			cell_volume[n] *= (all_grid._cell_length * all_grid._cell_length);
		}

		if(enable_cell_volume_correction) {
			correct_mass_and_volume_by_pressure_solve(
				all_grid.substance_density,
				cell_volume,
				all_grid
			);
		}

		//結果を格納
		for(int i = 0; i < all_grid.Grid_num_x; ++i) {
			for(int j = 0; j < all_grid.Grid_num_y; ++j) {
				all_grid.substance_density[get_voxel_center_index(i, j, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= advected_value_results[get_voxel_center_index(i, j, all_grid.Grid_num_x, all_grid.Grid_num_y)];
			}
		}

	}


	// density の MacCormack 法でのclamping の処理
	MY_FLOAT_TYPE clamp_substance_density_of_MacCormack(
		const Grid& all_grid,
		MY_FLOAT_TYPE substance_density_after_advect,
		VEC3_TYPE after_backtrace_position,
		const std::string interpolation_method) {
		// interpolation に使った値を記録する変数
		std::vector<MY_FLOAT_TYPE> interpolant_density_values;
		if (interpolation_method == "linear") {
			int advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
			int advected_index_y = floor((after_backtrace_position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
			if (advected_index_x < 0) {
				advected_index_x = 0;
				after_backtrace_position[0] = (advected_index_x * all_grid._cell_length);
			}
			if (advected_index_x >= all_grid.Grid_num_x - 1) {
				advected_index_x = all_grid.Grid_num_x - 2;
				after_backtrace_position[0] = (advected_index_x * all_grid._cell_length);
			}
			if (advected_index_y < 0) {
				advected_index_y = 0;
				after_backtrace_position[1] = (advected_index_y * all_grid._cell_length);
			}
			if (advected_index_y >= all_grid.Grid_num_y - 1) {
				advected_index_y = all_grid.Grid_num_y - 2;
				after_backtrace_position[1] = (advected_index_y * all_grid._cell_length);
			}

			for (int ix = 0; ix < 2; ++ix) {
				for (int iy = 0; iy < 2; ++iy) {
					interpolant_density_values.push_back(all_grid.substance_density[get_voxel_center_index(advected_index_x + ix, advected_index_y + iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]);
				}
			}
		}
		else if (interpolation_method == "WENO6") {
			int advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
			int advected_index_y = floor((after_backtrace_position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
			for (int ix = 0; ix < 6; ix++) {
				for (int iy = 0; iy < 6; ++iy) {
					int index_x = advected_index_x - 2 + ix;
					int index_y = advected_index_y - 2 + iy;
					if (index_x < 0) {
						index_x = 0;
					}
					if (index_x > all_grid.Grid_num_x - 1) {
						int exceed = index_x - (all_grid.Grid_num_x - 1);
						//境界の値がずっと外側まで続く場合
						//index_x = all_grid.Grid_num_x - 1;
						//境界を境に鏡のように値が反射する場合
						index_x = all_grid.Grid_num_x - 1 - exceed + 1;
					}
					if (index_y < 0) {
						index_y = 0;
					}
					if (index_y > all_grid.Grid_num_y - 1) {
						int exceed = index_y - (all_grid.Grid_num_y - 1);
						//境界の値がずっと外側まで続く場合
						//index_y = all_grid.Grid_num_y - 1;
						//境界を境に鏡のように値が反射する場合
						index_y = all_grid.Grid_num_y - 1 - exceed + 1;
					}
					interpolant_density_values.push_back(all_grid.substance_density[get_voxel_center_index(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)]);
				}
			}
		}
		// clamping の処理
		std::sort(interpolant_density_values.begin(), interpolant_density_values.end());
		substance_density_after_advect = std::min(substance_density_after_advect, interpolant_density_values[interpolant_density_values.size() - 1]);
		substance_density_after_advect = std::max(substance_density_after_advect, interpolant_density_values[0]);
		return substance_density_after_advect;
	}

	MY_FLOAT_TYPE integrate_normal_component_of_psi_on_face(
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
		splitted_faces = split_face(all_grid, face, split_method);

		VEC3_TYPE area_weighted_psi_sum(0.0, 0.0, 0.0);
		// バックトレース後の面上で積分した値を足していく
		for (int i_face = 0; i_face < splitted_faces.size(); ++i_face) {
			int num_gauss_quad_points;
			if (face.is_boundary_face) {
//				num_gauss_quad_points = 11;
				num_gauss_quad_points = num_gauss_quad_boundary;
			}
			else {
//				num_gauss_quad_points = 11;
				num_gauss_quad_points = num_gauss_quad_bulk;
			}
			std::vector<MY_FLOAT_TYPE> gauss_quadrature_weights = gauss_quadrature_points_1D::get_quadtarure_weights_1D(num_gauss_quad_points);
			std::vector<VEC3_TYPE> gauss_quadrature_positions_on_splitted_face
				= gauss_quadrature_points_1D::calc_quadtarure_positions_on_cell_face(splitted_faces[i_face], num_gauss_quad_points);

			for (int i_quad = 0; i_quad < num_gauss_quad_points; ++i_quad) {
				//face face の中心でのpsiを計算する
				VEC3_TYPE quadrature_position_on_cell_face = gauss_quadrature_positions_on_splitted_face[i_quad];
				VEC3_TYPE psi_density_on_quadrature_position
					= all_grid.calc_psi_substance_density_by_interpolation_cell_face(
						quadrature_position_on_cell_face,
						interpolation_method,
						use_integral,
						num_gauss_quadrature_point_for_integrate_density
					);
//				std::cout<<"calc_psi_substance_density_by_interpolation_cell_face() end"<<std::endl;
				// 1次元方向のみの補間を使う場合( 8/30 のノートの方法 )
				//VEC3_TYPE psi_density_on_quadrature_position = all_grid.calc_psi_substance_density_by_1D_interpolation_cell_face(quadrature_position_on_cell_face);

				area_weighted_psi_sum
					+= splitted_faces[i_face].calc_face_area()
					* gauss_quadrature_weights[i_quad]
					* psi_density_on_quadrature_position;

				// face normal の計算
				//VEC3_TYPE face_normal = splitted_faces[i_face].calc_normal();
				//mass_in_the_cell += psi_density_on_cell_face.dot(face_normal) * splitted_faces[i_face].calc_face_area();
			}
			//getchar();
		}
		//考えている面が境界の面かどうか
		//bool is_boundary_cell_face = false;
		//if (  (vertex_1._vertex_index[0] == 0                           && vertex_2._vertex_index[0] == 0                          )
		//	||(vertex_1._vertex_index[0] == all_grid.Grid_num_x && vertex_2._vertex_index[0] == all_grid.Grid_num_x)
		//	||(vertex_1._vertex_index[1] == 0                           && vertex_2._vertex_index[1] == 0                          )
		//	||(vertex_1._vertex_index[1] == all_grid.Grid_num_y && vertex_2._vertex_index[1] == all_grid.Grid_num_y)){
		//	is_boundary_cell_face = true;
		//}
		//if (is_boundary_cell_face) {
		//	VEC3_TYPE psi_at_cell_face;
		//	//x軸に垂直な壁の場合
		//	if (vertex_1._vertex_index[0] == vertex_2._vertex_index[0]) {
		//		psi_at_cell_face = VEC3_TYPE(0.0, 0.0, 0.0);
		//	}
		//	//y軸に垂直な壁の場合
		//	else if (vertex_1._vertex_index[1] == vertex_2._vertex_index[1]) {
		//		const int min_x_index = std::min(vertex_1._vertex_index[0], vertex_2._vertex_index[0]);
		//		psi_at_cell_face = VEC3_TYPE(0.0, all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y(min_x_index, vertex_1._vertex_index[1])], 0.0);
		//	}
		//	mass_in_the_cell += face.calc_face_area() * psi_at_cell_face.dot(face.calc_normal());
		//}
		//else {
		return area_weighted_psi_sum.dot(face.calc_normal());
		//}
	}

	VEC3_TYPE analytical_integral_of_psi_on_splitted_face(
		const Grid& all_grid,
		cell_face face
	) {
		////補間に線型補間をつかう場合
		VEC3_TYPE face_center = face.calc_face_center();
		//補間に使う面のindex
		int advected_index_x[2];
		int advected_index_y[2];
		advected_index_x[0] = floor((face_center[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
		advected_index_x[1] = floor((face_center[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length) + 1;
		advected_index_y[0] = floor((face_center[1]) / all_grid._cell_length);
		advected_index_y[1] = floor((face_center[1]) / all_grid._cell_length) + 1;
		//補間に使う4つの面の位置(をれぞれx座標とy座標)
		MY_FLOAT_TYPE face_center_pos_x[2], face_center_pos_y[2];
		face_center_pos_x[0] = (advected_index_x[0] + 0.5) * all_grid._cell_length;
		face_center_pos_x[1] = (advected_index_x[1] + 0.5) * all_grid._cell_length;
		face_center_pos_y[0] = (advected_index_y[0]) * all_grid._cell_length;
		face_center_pos_y[1] = (advected_index_y[1]) * all_grid._cell_length;
		// 範囲外を参照する場合はクランピングする (これはall_grid.psi_substance_density_cell_face_yにアクセスする際に範囲外を参照するとまずいから行っている
	    // face_center_pos_x, face_center_pos_y 計算にはクランピングする前のindex を使うべき).
		if(advected_index_x[0] < 0){
			advected_index_x[0] = 0;
		}
		if(advected_index_x[0] > all_grid.Grid_num_x - 1){
			advected_index_x[0] = all_grid.Grid_num_x - 1;
		}
		if(advected_index_x[1] < 0){
			advected_index_x[1] = 0;
		}
		if(advected_index_x[1] > all_grid.Grid_num_x - 1){
			advected_index_x[1] = all_grid.Grid_num_x - 1;
		}
		if(advected_index_y[0] < 0){
			advected_index_y[0] = 0;
		}
		if(advected_index_y[0] > all_grid.Grid_num_y){
			advected_index_y[0] = all_grid.Grid_num_y;
		}
		if(advected_index_y[1] < 0){
			advected_index_y[1] = 0;
		}
		if(advected_index_y[1] > all_grid.Grid_num_y){
			advected_index_y[1] = all_grid.Grid_num_y;
		}
		//backtrace face の端点の座標
		VEC3_TYPE backtrace_face_vertex_pos_0 = face._vertex_list[0]._vertex_pos;
		VEC3_TYPE backtrace_face_vertex_pos_1 = face._vertex_list[1]._vertex_pos;
		//補間に使う4つの面でのpsi_y
		MY_FLOAT_TYPE cell_center_psi_y[2][2];
		cell_center_psi_y[0][0] = all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y(advected_index_x[0], advected_index_y[0], all_grid.Grid_num_x, all_grid.Grid_num_y)];
		cell_center_psi_y[1][0] = all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y(advected_index_x[1], advected_index_y[0], all_grid.Grid_num_x, all_grid.Grid_num_y)];
		cell_center_psi_y[0][1] = all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y(advected_index_x[0], advected_index_y[1], all_grid.Grid_num_x, all_grid.Grid_num_y)];
		cell_center_psi_y[1][1] = all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y(advected_index_x[1], advected_index_y[1], all_grid.Grid_num_x, all_grid.Grid_num_y)];
		//psi のy成分の積分
		MY_FLOAT_TYPE integral_of_psi_y;
		MY_FLOAT_TYPE term0 = cell_center_psi_y[1][1] * ( 2.0 * backtrace_face_vertex_pos_0[0] * backtrace_face_vertex_pos_0[1]
												  + backtrace_face_vertex_pos_0[1] * backtrace_face_vertex_pos_1[0]
												  + backtrace_face_vertex_pos_0[0] * backtrace_face_vertex_pos_1[1]
												  + 2.0 * backtrace_face_vertex_pos_1[0] * backtrace_face_vertex_pos_1[1]
												  - 3.0 * backtrace_face_vertex_pos_0[1] * face_center_pos_x[0]
												  - 3.0 * backtrace_face_vertex_pos_1[1] * face_center_pos_x[0]
												  - 3.0 * (backtrace_face_vertex_pos_0[0] + backtrace_face_vertex_pos_1[0] - 2.0 * face_center_pos_x[0]) * face_center_pos_y[0]
												 ) / (6.0 * all_grid._cell_length * all_grid._cell_length);
		MY_FLOAT_TYPE term1 = cell_center_psi_y[0][1] * ( - backtrace_face_vertex_pos_0[0] * (2.0 * backtrace_face_vertex_pos_0[1] + backtrace_face_vertex_pos_1[1])
									  			   - backtrace_face_vertex_pos_1[0] * (backtrace_face_vertex_pos_0[1] + 2.0 * backtrace_face_vertex_pos_1[1])
									  	 		   + 3.0 * (backtrace_face_vertex_pos_0[1] + backtrace_face_vertex_pos_1[1]) * face_center_pos_x[1]
									  	 		   + 3.0 * (backtrace_face_vertex_pos_0[0] + backtrace_face_vertex_pos_1[0] - 2.0 * face_center_pos_x[1]) * face_center_pos_y[0]
								  	 		   	 ) / (6.0 * all_grid._cell_length * all_grid._cell_length);
		MY_FLOAT_TYPE term2 = cell_center_psi_y[1][0] * ( - backtrace_face_vertex_pos_0[0] * (2.0 * backtrace_face_vertex_pos_0[1] + backtrace_face_vertex_pos_1[1])
			  									   - backtrace_face_vertex_pos_1[0] * (backtrace_face_vertex_pos_0[1] + 2.0 * backtrace_face_vertex_pos_1[1])
			  									   + 3.0 * (backtrace_face_vertex_pos_0[1] + backtrace_face_vertex_pos_1[1]) * face_center_pos_x[0]
			  									   + 3.0 * (backtrace_face_vertex_pos_0[0] + backtrace_face_vertex_pos_1[0] - 2.0 * face_center_pos_x[0]) * face_center_pos_y[1]
			  									 ) / (6.0 * all_grid._cell_length * all_grid._cell_length);
		MY_FLOAT_TYPE term3 = cell_center_psi_y[0][0] * ( 2.0 * backtrace_face_vertex_pos_0[0] * backtrace_face_vertex_pos_0[1]
												   + backtrace_face_vertex_pos_0[1] * backtrace_face_vertex_pos_1[0]
												   + backtrace_face_vertex_pos_0[0] * backtrace_face_vertex_pos_1[1]
												   + 2.0 * backtrace_face_vertex_pos_1[0] * backtrace_face_vertex_pos_1[1]
												   - 3.0 * backtrace_face_vertex_pos_0[1] * face_center_pos_x[1]
												   - 3.0 * backtrace_face_vertex_pos_1[1] * face_center_pos_x[1]
											       - 3.0 * (backtrace_face_vertex_pos_0[0] + backtrace_face_vertex_pos_1[0] - 2.0 * face_center_pos_x[1]) * face_center_pos_y[1]
												 ) / (6.0 * all_grid._cell_length * all_grid._cell_length);
		integral_of_psi_y = term0 + term1 + term2 + term3;
		MY_FLOAT_TYPE jacobian = (backtrace_face_vertex_pos_0- backtrace_face_vertex_pos_1).norm();
		return jacobian * VEC3_TYPE(0.0, integral_of_psi_y, 0.0);
	}

	MY_FLOAT_TYPE integrate_normal_component_of_psi_on_face_analytically(
		const Grid& all_grid,
		cell_face face,
		const std::string split_method
	) {
		//バックトレースした面の切断
		std::vector<cell_face> splitted_faces;
		//考えるバックトレース後の面を切断するかどうか
		bool split_this_face = true;
		if (split_this_face) {
			splitted_faces = split_face(all_grid, face, split_method);
		}
		else {
			// 8/22 のノートの方法
			splitted_faces.push_back(face);
		}
		VEC3_TYPE area_weighted_psi_sum(0.0, 0.0, 0.0);
		// バックトレース後の面上で積分した値を足していく
		for (int i_face = 0; i_face < splitted_faces.size(); ++i_face) {
			area_weighted_psi_sum
				+= analytical_integral_of_psi_on_splitted_face(all_grid, splitted_faces[i_face]);
		}
		return area_weighted_psi_sum.dot(face.calc_normal());
	}

	// density の MacCormack 法でのclamping の処理
	MY_FLOAT_TYPE clamp_integral_of_normal_component_of_psi_of_MacCormack(
		const Grid& all_grid,
		MY_FLOAT_TYPE integral_of_normal_component_of_psi_after_advect,
		VEC3_TYPE after_backtrace_position,
		const std::string interpolation_method) {
		// interpolation に使った値を記録する変数
		std::vector<MY_FLOAT_TYPE> interpolant_integral_of_normal_component_of_psi_values;
		if (interpolation_method == "linear") {
			// x軸に垂直な壁での積分
			int advected_index_x = floor((after_backtrace_position[0]) / all_grid._cell_length);
			int advected_index_y = floor((after_backtrace_position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
			if (advected_index_x <= 0) {
				advected_index_x = 0;
			}
			if (advected_index_x >= all_grid.Grid_num_x) {
				advected_index_x = all_grid.Grid_num_x - 1;
			}
			if (advected_index_y < 0) {
				advected_index_y = 0;
			}
			if (advected_index_y >= all_grid.Grid_num_y - 1) {
				advected_index_y = all_grid.Grid_num_y - 1;
			}
			for (int ix = 0; ix < 2; ++ix) {
				for (int iy = 0; iy < 2; ++iy) {
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						all_grid.psi_substance_density_cell_face_x[get_voxel_face_index_x(advected_index_x + ix, advected_index_y + iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						* all_grid._cell_length
					);
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						-all_grid.psi_substance_density_cell_face_x[get_voxel_face_index_x(advected_index_x + ix, advected_index_y + iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						* all_grid._cell_length
					);
				}
			}
			// y軸に垂直な壁での積分
			advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
			advected_index_y = floor((after_backtrace_position[1]) / all_grid._cell_length);
			if (advected_index_x < 0) {
				advected_index_x = 0;
			}
			if (advected_index_x >= all_grid.Grid_num_x - 1) {
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
						all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y(advected_index_x + ix, advected_index_y + iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						* all_grid._cell_length
					);
					interpolant_integral_of_normal_component_of_psi_values.push_back(
						-all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y(advected_index_x + ix, advected_index_y + iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						* all_grid._cell_length
					);
				}
			}
		}
		else if(interpolation_method == "1Dy_linear") {
			// x軸に垂直な壁での積分
			int advected_index_x = floor((after_backtrace_position[0]) / all_grid._cell_length);
			int advected_index_y = floor((after_backtrace_position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
			if (advected_index_x <= 0) {
				advected_index_x = 0;
			}
			if (advected_index_x >= all_grid.Grid_num_x) {
				advected_index_x = all_grid.Grid_num_x - 1;
			}
			if (advected_index_y < 0) {
				advected_index_y = 0;
			}
			if (advected_index_y >= all_grid.Grid_num_y - 1) {
				advected_index_y = all_grid.Grid_num_y - 1;
			}
			for (int iy = 0; iy < 2; ++iy) {
				interpolant_integral_of_normal_component_of_psi_values.push_back(
					all_grid.psi_substance_density_cell_face_x[get_voxel_face_index_x(advected_index_x, advected_index_y + iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					* all_grid._cell_length
				);
				interpolant_integral_of_normal_component_of_psi_values.push_back(
					-all_grid.psi_substance_density_cell_face_x[get_voxel_face_index_x(advected_index_x, advected_index_y + iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					* all_grid._cell_length
				);
			}
			// y軸に垂直な壁での積分
			advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
			advected_index_y = floor((after_backtrace_position[1]) / all_grid._cell_length);
			if (advected_index_x < 0) {
				advected_index_x = 0;
			}
			if (advected_index_x >= all_grid.Grid_num_x - 1) {
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
					all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y(advected_index_x, advected_index_y + iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					* all_grid._cell_length
				);
				interpolant_integral_of_normal_component_of_psi_values.push_back(
					-all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y(advected_index_x, advected_index_y + iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					* all_grid._cell_length
				);
			}
		}
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
		// clamping の処理
		std::sort(interpolant_integral_of_normal_component_of_psi_values.begin(), interpolant_integral_of_normal_component_of_psi_values.end());
		integral_of_normal_component_of_psi_after_advect = std::min(integral_of_normal_component_of_psi_after_advect, interpolant_integral_of_normal_component_of_psi_values[interpolant_integral_of_normal_component_of_psi_values.size() - 1]);
		integral_of_normal_component_of_psi_after_advect = std::max(integral_of_normal_component_of_psi_after_advect, interpolant_integral_of_normal_component_of_psi_values[0]);
		return integral_of_normal_component_of_psi_after_advect;
	}


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
	) {
		// psi の離散値を計算
		calc_psi_on_cell_face_from_density_on_cell_center(
			all_grid,
			all_grid.psi_substance_density_cell_face_x,
			all_grid.psi_substance_density_cell_face_y,
			advected_values,
			all_grid.Grid_num_x,
			all_grid.Grid_num_y,
			all_grid._cell_length,
			interpolation_method_in_calc_psi
		);
		//計算結果を格納する変数
//		std::vector<MY_FLOAT_TYPE> substance_density_after_advect(all_grid.Grid_num_x * all_grid.Grid_num_y);
		std::vector<MY_FLOAT_TYPE> mass_after_advect(all_grid.Grid_num_x * all_grid.Grid_num_y);
		//移流後のセル体積を格納する変数
		std::vector<MY_FLOAT_TYPE> cell_volume_after_advect(all_grid.Grid_num_x * all_grid.Grid_num_y);
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				//(ix, iy)番目のセルのface を走るループ
				int dx[4], dy[4];
				dx[0] = 0; dx[1] = 1; dx[2] = 1; dx[3] = 0;
				dy[0] = 0; dy[1] = 0; dy[2] = 1; dy[3] = 1;
				MY_FLOAT_TYPE mass_in_the_cell = 0.0;
				MY_FLOAT_TYPE cell_volume = 0.0;
				bool is_out_of_range = false;
				for (int i_face = 0; i_face < 4; ++i_face) {
					// 考えるcell face を構成するvertex のindex
					std::vector<int> vertex_index_1{ ix + dx[i_face], iy + dy[i_face] }, vertex_index_2{ ix + dx[(i_face + 1) % 4], iy + dy[(i_face + 1) % 4] };
					// バックトレース前の cell face を構成するvertexを定義する
					cell_vertex before_backtrace_vertex_1(vertex_index_1, all_grid._cell_length);
					cell_vertex before_backtrace_vertex_2(vertex_index_2, all_grid._cell_length);
					// バックトレース前の頂点から面を定義
					cell_face before_backtrace_face;
					before_backtrace_face._vertex_list[0] = before_backtrace_vertex_1;
					before_backtrace_face._vertex_list[1] = before_backtrace_vertex_2;
					// MacCormack を使う場合
					if (use_MacCormack_scheme) {
//						std::cout<<"begin MacCormack (density)"<<std::endl;
//						std::cout << "ix + dx[i_face]: "<<ix + dx[i_face]<<std::endl;
//						std::cout << "iy + dy[i_face]: "<<iy + dy[i_face]<<std::endl;
//						std::cout << "ix + dx[(i_face + 1) % 4]: "<<ix + dx[(i_face + 1) % 4]<<std::endl;
//						std::cout << "iy + dy[(i_face + 1) % 4]: "<<iy + dy[(i_face + 1) % 4]<<std::endl;
						////エラーを補正する前のベースとなる面での計算
						// バックトレース後の頂点から面を定義
						cell_face after_backtrace_face
							= calc_backtraced_face(
								all_grid,
								time_step_length,
								before_backtrace_face,
								"forward",
//								set_zero_normal_velocity_at_boundary,
								interpolation_method,
								"density"
							);
						//バックトレースした先の位置がグリッドの範囲を超えていたら修正する
						after_backtrace_face
							= clamp_vertex_position_by_grid(
								all_grid,
								after_backtrace_face,
								"density"
							);
						// バックトレース後の面について psi と 法線の内積を面上で積分した値を計算する
						MY_FLOAT_TYPE integral_of_normal_component_of_psi_on_after_backtrace_face
							= integrate_normal_component_of_psi_on_face(
								all_grid,
								after_backtrace_face,
								num_gauss_quad_boundary,
								num_gauss_quad_bulk,
								split_method,
								interpolation_method,
								use_integral,
								num_gauss_quadrature_point_for_integrate_density
							);
						////////////////////////////////////////
						////////// エラーの計算
						////////////////////////////////////////
						// 元のface(バックトレース前のオリジナルの face)について, psi と 法線の内積を面上で積分した値を計算する
						MY_FLOAT_TYPE integral_of_normal_component_of_psi_on_before_backtrace_face
							= integrate_normal_component_of_psi_on_face(
								all_grid,
								before_backtrace_face,
								num_gauss_quad_boundary,
								num_gauss_quad_bulk,
								split_method,
								interpolation_method,
								use_integral,
								num_gauss_quadrature_point_for_integrate_density
							);
						// バックトレース後の面を逆向きにバックトレースした面の計算(エラーの計算用に使う)
						cell_face backward_trace_face_of_after_backtrace_face
							= calc_backtraced_face(
								all_grid,
								time_step_length,
								after_backtrace_face,
								"backward",
//								set_zero_normal_velocity_at_boundary,
								interpolation_method,
								"density"
							);
						//バックトレースした先の位置がグリッドの範囲を超えていたら修正する
						backward_trace_face_of_after_backtrace_face
							= clamp_vertex_position_by_grid(
								all_grid,
								backward_trace_face_of_after_backtrace_face,
								"density"
							);
						// psi と 法線の内積を面上で積分した値を計算する
						MY_FLOAT_TYPE integral_of_normal_component_of_psi_on_backward_trace_face_of_after_backtrace_face
							= integrate_normal_component_of_psi_on_face(
								all_grid,
								backward_trace_face_of_after_backtrace_face,
								num_gauss_quad_boundary,
								num_gauss_quad_bulk,
								split_method,
								interpolation_method,
								use_integral,
								num_gauss_quadrature_point_for_integrate_density
							);
						// MacCormack におけるエラーの定義
						MY_FLOAT_TYPE error
							= (integral_of_normal_component_of_psi_on_before_backtrace_face
								- integral_of_normal_component_of_psi_on_backward_trace_face_of_after_backtrace_face)
							/ 2.0;
						// エラーを修正した値
						MY_FLOAT_TYPE error_corrected_integral_of_normal_component_of_psi = integral_of_normal_component_of_psi_on_after_backtrace_face + error;
						//clamping
						if (use_clamping_in_MacCormack_scheme) {
							error_corrected_integral_of_normal_component_of_psi
								= clamp_integral_of_normal_component_of_psi_of_MacCormack(
									all_grid,
									error_corrected_integral_of_normal_component_of_psi,
									after_backtrace_face.calc_face_center(),
									interpolation_method);
						}
//						std::cout<<"end MacCormack (density) use_clamping_in_MacCormack_scheme"<<std::endl;
						// エラーを修正した値を用いて質量への寄与を計算する
						mass_in_the_cell
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
//								set_zero_normal_velocity_at_boundary,
								interpolation_method,
								"density"
							);
						//バックトレースした先の位置がグリッドの範囲を超えていたら修正する
						after_backtrace_face
							= clamp_vertex_position_by_grid(
								all_grid,
								after_backtrace_face,
								"density"
							);
						// psi と 法線の内積を面上で積分した値を計算する
						if(integral_method == "gauss"){
							mass_in_the_cell
								+= integrate_normal_component_of_psi_on_face(
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
						else if(integral_method == "analytical"){
							mass_in_the_cell
								+= integrate_normal_component_of_psi_on_face_analytically(
									all_grid,
									after_backtrace_face,
									split_method);
						}
						if(enable_cell_volume_correction){
							calc_cell_volume_contribution(
								after_backtrace_face,
								cell_volume
							);
						}
					}
					////------------------------------
				}
				//質量の計算結果を格納
				mass_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= mass_in_the_cell;
				//体積の計算結果を格納
				if(enable_cell_volume_correction){
//					if(abs(cell_volume) < 0.001 * all_grid._cell_volume){
//						cell_volume_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//							= all_grid._cell_volume;
//					}
//					else{
						cell_volume_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
							= cell_volume;
//					}
				}
				else{
					cell_volume_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						= all_grid._cell_volume;
				}
/*
				if(abs(cell_volume) < 0.000000001){
					substance_density_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = mass_in_the_cell / all_grid._cell_volume;
//					substance_density_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = mass_in_the_cell / 0.000000001;
//					substance_density_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 0.0;
				}
				else{
					substance_density_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = mass_in_the_cell / cell_volume;
				}
*/
			}
		}

		//全セルの体積の和がグリッド全体の体積に等しくなるように再スケーリング
/*
		MY_FLOAT_TYPE sum_cell_volumes = 0.0;
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				sum_cell_volumes += cell_volume_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
			}
		}
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				cell_volume_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					*= ((all_grid._cell_length * all_grid.Grid_num_x * all_grid._cell_length * all_grid.Grid_num_y) / sum_cell_volumes);
			}
		}
*/

		//psi の補正
/*
		if(enable_cell_volume_correction) {
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					advected_values[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
						*=(all_grid.cell_volume_cell_center[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
							/ all_grid._cell_volume);
				}
			}
		}
*/
		//計算結果をコピー
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				//体積の計算結果
				all_grid.cell_volume_cell_center[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= cell_volume_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
//				all_grid.cell_volume_cell_center[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//					*= (cell_volume_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] / all_grid._cell_volume);

				//密度の計算結果
//				advected_values[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//					= mass_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//					/ cell_volume_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
				advected_values[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					= mass_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					/ all_grid._cell_volume;

//				advected_values[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//					*=(all_grid.cell_volume_cell_center[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//						/ all_grid._cell_volume);

//				advected_values[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = substance_density_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
//				all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = substance_density_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
			}
		}

	}

	// psi から計算した密度が正しいかのチェック
	void check_psi_density(Grid& all_grid) {
	}

	void check_negative_dinsity(Grid& all_grid) {
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			if (abs(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, 0, all_grid.Grid_num_x, all_grid.Grid_num_y)]) > 0.0000000001) {
				std::cout << "boundary velocity_in_voxel_face_y: " << all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, 0, all_grid.Grid_num_x, all_grid.Grid_num_y)] << std::endl;
				getchar();
			}
		}
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			if (abs(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, all_grid.Grid_num_y, all_grid.Grid_num_x, all_grid.Grid_num_y)]) > 0.0000000001) {
				std::cout << "boundary velocity_in_voxel_face_y: " << all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, all_grid.Grid_num_y, all_grid.Grid_num_x, all_grid.Grid_num_y)] << std::endl;
				getchar();
			}
		}
		for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
			if (abs(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(0, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]) > 0.0000000001) {
				std::cout << "boundary velocity_in_voxel_face_x: " << all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(0, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] << std::endl;
				getchar();
			}
		}
		for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
			if (abs(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(all_grid.Grid_num_x, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]) > 0.0000000001) {
				std::cout << "boundary velocity_in_voxel_face_x: " << all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(all_grid.Grid_num_x, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] << std::endl;
				getchar();
			}
		}
	}
}// namespace smoke_simulation
