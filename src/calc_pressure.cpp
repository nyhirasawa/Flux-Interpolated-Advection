#include "calc_pressure.h"

#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

#include "sparse_matrix.h"
#include "utils.h"

namespace smoke_simulation{
    //Poisson eq. を解くことによる圧力の計算
	void calc_pressure(
		Grid& all_grid,
		const bool enable_cell_volume_correction,
		const bool fix_velocity,
		const MY_FLOAT_TYPE time_step_length
	) {

		//グリッドの総数
		const int N = all_grid.Grid_num_x
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
			//b[get_voxel_center_index(ix, iy)] = (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix + 1, iy)] - all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy)]
			//	+ all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy + 1)] - all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy)])
			//	* all_grid._cell_length
			//	/ (smoke_simulation::physical_const::kDt);
			b[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
				= 0.0;

			if(!fix_velocity){
				b[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+= (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix + 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					-  all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					+  all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy + 1, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					-  all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					* all_grid._cell_length;
			}

//			if(enable_cell_volume_correction){
//				b[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//					+= +(all_grid._cell_volume - all_grid.cell_volume_cell_center[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
//					;
//					/ all_grid._cell_length;
//					* all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//					/ (all_grid._cell_length * all_grid._cell_length);
//			}

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
        	amgcl::solver::cg<Backend>
		    > Solver;
//		int n = N;
		int n = N - 1;
//		int n = all_grid.Grid_num_x * all_grid.Grid_num_y;

//		Solver::params prm;
//	    prm.solver.tol = 1e-4;
//		Solver solve( std::tie(n, A.cumulative_num_nonzero_element_in_row, A.column_index, A.element_value) ,prm);
		Solver solve( std::tie(n, A.cumulative_num_nonzero_element_in_row, A.column_index, A.element_value));

//		std::cout<<"getchar()"<<std::endl;
//		getchar();
		int    iters;
		MY_FLOAT_TYPE error;
//		std::cout<<"getchar()"<<std::endl;
//		getchar();
		///////////////
//		std::tie(iters, error) = solve(b, all_grid.pressure);

		std::vector<MY_FLOAT_TYPE> truncated_pressure(N - 1);
		std::tie(iters, error) = solve(b, truncated_pressure);
		//係数行列をフルランクにするために落とした最後の要素を復元する
		truncated_pressure.push_back(0.0);
		all_grid.pressure = truncated_pressure;

//		std::cout<<"iters: "<<iters<<std::endl;
//		std::cout<<"error: "<<error<<std::endl;

//		std::cout<<"linear solve end"<<std::endl;
	}

	//pressure gradient termの計算
	void calc_pressure_gradient_term(
		Grid& all_grid,
		const bool enable_cell_volume_correction,
		const bool fix_velocity,
		const MY_FLOAT_TYPE time_step_length
	) {
		//圧力の計算
		calc_pressure(
			all_grid,
			enable_cell_volume_correction,
			fix_velocity,
			time_step_length
		);
		//圧力場に境界条件をセット
//		set_boundary_pressure(all_grid);
		//pressure gradeint term によって速度場を更新
		for (int ix = 1; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
				//all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy)] -= (all_grid.pressure[get_voxel_center_index(ix, iy)] - all_grid.pressure[get_voxel_center_index(ix - 1, iy)])
				//	* smoke_simulation::physical_const::kDt
				//	/ (all_grid._cell_length);
				all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					-= (all_grid.pressure[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] - all_grid.pressure[get_voxel_center_index(ix - 1, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					/ (all_grid._cell_length);
			}
		}
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
			for (int iy = 1; iy < all_grid.Grid_num_y; iy++) {
				//all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy)] -= (all_grid.pressure[get_voxel_center_index(ix, iy)] - all_grid.pressure[get_voxel_center_index(ix, iy - 1)])
				//	* smoke_simulation::physical_const::kDt
				//	/ (all_grid._cell_length);
				all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
					-= (all_grid.pressure[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] - all_grid.pressure[get_voxel_center_index(ix, iy - 1, all_grid.Grid_num_x, all_grid.Grid_num_y)])
					/ (all_grid._cell_length);
			}
		}
	}
}
