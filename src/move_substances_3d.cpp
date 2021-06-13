#include "move_substances_3d.h"

#include <iostream>
#include <chrono>
#include <thread>
#include <unordered_map>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

#include "linear_solver.h"
#include "physical_const.h"
#include "utils.h"
#include "cell_vertex_3D.h"
#include "cell_face_3D.h"
#include "split_face_3D.h"
#include "triangle_3D.h"
#include "polygon_vector.h"
#include "define_float_type.h"
#include "gauss_quadrature_points.h"
#include "sginterp3.h"
#include "calc_psi.h"
#include "calc_psi_3D.h"
#include "calc_backtraced_face_3D.h"
#include "backtrace_and_calc_all_quadrature_point_list.h"
#include "calc_density_in_cell_from_quadrature_point_list.h"
#include "clamping_in_MacCormack_scheme.h"
#include "backtrace_and_calc_all_quadrature_point_list_for_integrate_density.h"
#include "calc_density_in_cell_from_quadrature_point_list_for_integrate_density.h"

namespace smoke_simulation{
    //上の4ステップをまとめただけの関数(substance densityの1時間ステップ分の更新に相当)
    void move_substances_3D(
        Grid_3D& all_grid,
        int i_frame,
    	const MY_FLOAT_TYPE time_step_length,
    	const bool use_flux_advection,
//    	const bool set_zero_normal_velocity_at_boundary,
    	const bool use_MacCormack_scheme,
    	const bool use_clamping_in_MacCormack_scheme,
    	const int num_gauss_quad_boundary,
    	const int num_gauss_quad_bulk,
    	const bool enable_eliminate_zero_velocity_false_diffusion,
        const std::string split_method,
        const std::string integral_method,
        const int num_gauss_quadrature_points,
        const std::string interpolation_method,
        bool use_integral,
	    const int num_gauss_quadrature_point_for_integrate_density,
		const bool enable_cell_volume_correction,
		const MY_FLOAT_TYPE minimal_area
    ){
        if (use_flux_advection) {
            std::vector<MY_FLOAT_TYPE> before_advect_density = all_grid.substance_density;
            ////psiの離散値の計算
            std::string interpolation_method_in_calc_psi;
            if(interpolation_method =="linear"){
                interpolation_method_in_calc_psi = "linear";
            }
            else if(interpolation_method=="1Dy_linear"){
    //            interpolation_method_in_calc_psi = "linear";
    //            interpolation_method_in_calc_psi = "1Dy_linear";
                interpolation_method_in_calc_psi = "const";
                // piecewise-constant の場合, psiの計算に補間を使っても積分で厳密に評価してもは同じなので
                // use_integral は false にしておく(補間を使うほうが計算が速い & piecewise-constantの積分を実装してないから)
                use_integral = false;
            }
            else if(interpolation_method=="WENO6"){
                interpolation_method_in_calc_psi = "WENO6";
//                interpolation_method_in_calc_psi = "const";
            }
            else if(interpolation_method=="WENO6-optimized"){
                interpolation_method_in_calc_psi = "WENO6-optimized";
//                interpolation_method_in_calc_psi = "const";
            }
//            else{
//                interpolation_method_in_calc_psi = "const";
//            }
            advect_density_flux_advection_3D_std_thread(
        		all_grid,
        		time_step_length,
//        		set_zero_normal_velocity_at_boundary,
        		use_MacCormack_scheme,
        		use_clamping_in_MacCormack_scheme,
        		num_gauss_quad_boundary,
        		num_gauss_quad_bulk,
        		enable_eliminate_zero_velocity_false_diffusion,
                split_method,
                integral_method,
                num_gauss_quadrature_points,
                all_grid.substance_density,
                interpolation_method,
                interpolation_method_in_calc_psi,
                false,
                use_integral,
        	    num_gauss_quadrature_point_for_integrate_density,
                enable_cell_volume_correction,
        		minimal_area
            );
            // セルの体積によって質量密度を補正
            if(enable_cell_volume_correction){
                correct_mass_and_volume_by_pressure_solve_3D(
                    all_grid.substance_density,
                    all_grid.cell_volume_cell_center,
                    all_grid
                );
            }
            ////// エラーの補正
            ////// "バックトレース前の面上で積分した値" と "バックトレース前の面で定義されたpsiの離散値"の差を引いていく
            if (enable_eliminate_zero_velocity_false_diffusion) {
                // 各セル中心での密度の補正量
                std::vector<MY_FLOAT_TYPE> zero_velocity_correction_values(all_grid.Grid_num_x * all_grid.Grid_num_y * all_grid.Grid_num_z);
                for(int ix = 0; ix < all_grid.Grid_num_x; ++ix){
                    for(int iy = 0; iy < all_grid.Grid_num_y; ++iy){
                        for(int iz = 0; iz < all_grid.Grid_num_z; ++iz){
                            zero_velocity_correction_values[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                =before_advect_density[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                        }
                    }
                }
                ////-------------------------------------
                //// 各セル中心での密度の補正量を計算する
                ////-------------------------------------
                // zero velocityで 移流させたときの密度場を求める
                // ゼロ速度場で移流を行っているのでバックトレース前と後のfaceは等しい。よって psi の計算には補間を使っても積分で評価しても同じ
                // なので use_integral は true でも false でも良いが、計算量は補間を使う false の方が少ないのでuse_integral は false にしておく
                // (あと use_integral を trueにするとアーティファクトがでる?(要確認))
                advect_density_flux_advection_3D_std_thread(
            		all_grid,
            		time_step_length,
//            		set_zero_normal_velocity_at_boundary,
            		use_MacCormack_scheme,
            		use_clamping_in_MacCormack_scheme,
            		num_gauss_quad_boundary,
            		num_gauss_quad_bulk,
            		enable_eliminate_zero_velocity_false_diffusion,
                    split_method,
                    integral_method,
                    num_gauss_quadrature_points,
                    zero_velocity_correction_values,
                    interpolation_method,
                    interpolation_method_in_calc_psi,
                    true,
                    /*use_integral*/false,
            	    num_gauss_quadrature_point_for_integrate_density,
                    /*enable_cell_volume_correction*/false,
            		minimal_area
                );
                // 密度場の離散地との差を求める. これが移流前のエラー
                for(int ix = 0; ix < all_grid.Grid_num_x; ++ix){
                    for(int iy = 0; iy < all_grid.Grid_num_y; ++iy){
                        for(int iz = 0; iz < all_grid.Grid_num_z; ++iz){
                            zero_velocity_correction_values[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                *= -1.0;
                            zero_velocity_correction_values[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                +=before_advect_density[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
    //                            =all_grid.substance_density[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                        }
                    }
                }
                ////-------------------------------------
                //// エラーの移流
                ////-------------------------------------
    /*
                // エラーの移流と密度の移流で同じ補間の方法を使う
                std::string interpolation_method_of_error_advection = interpolation_method;
                std::string split_method_of_error_advection = split_method;
    //            std::string interpolation_method_in_calc_psi = interpolation_method;
                std::string interpolation_method_in_calc_psi = "const";
    */
                // エラーの移流に piecewise-const(一次元補間) を使う
                std::string interpolation_method_of_error_advection;
                std::string interpolation_method_in_calc_psi;
                if(interpolation_method == "linear" || interpolation_method == "1Dy_linear"){
                    interpolation_method_of_error_advection = "1Dy_linear";
                    interpolation_method_in_calc_psi = "const";
                }
                else if(interpolation_method == "WENO6"||interpolation_method == "WENO6-optimized"){
//                    interpolation_method_of_error_advection = "1Dy_WENO6";
//                    interpolation_method_in_calc_psi = "1Dy_WENO6";
                    interpolation_method_of_error_advection = "1Dy_linear";
                    interpolation_method_in_calc_psi = "const";
                }
                std::string split_method_of_error_advection = "xz-AxisAligned";

                // flux advection をエラーの移流に使う場合
                // エラーの移流には1次元補間を使うが、1次元補間ではpsiの計算に補間を使っても積分で厳密に計算しても結果は同じである。
                // なので use_integral は true でも false でも良いが、計算量は補間を使う false の方が少ないのでuse_integral は false にしておく
                // (あと use_integral を trueにするとアーティファクトがでる?(要確認))
                advect_density_flux_advection_3D_std_thread(
                    all_grid,
                    time_step_length,
//                    set_zero_normal_velocity_at_boundary,
                    use_MacCormack_scheme,
                    use_clamping_in_MacCormack_scheme,
                    num_gauss_quad_boundary,
                    num_gauss_quad_bulk,
                    enable_eliminate_zero_velocity_false_diffusion,
                    split_method_of_error_advection,
                    integral_method,
                    num_gauss_quadrature_points,
                    zero_velocity_correction_values,
                    interpolation_method_of_error_advection,
                    interpolation_method_in_calc_psi,
                    false,
                    /*use_integral=*/false,
            	    num_gauss_quadrature_point_for_integrate_density,
                    /*enable_cell_volume_correction*/false,
            		minimal_area
                );

                //密度を補正
                for(int ix = 0; ix < all_grid.Grid_num_x; ++ix){
                    for(int iy = 0; iy < all_grid.Grid_num_y; ++iy){
                        for(int iz = 0; iz < all_grid.Grid_num_z; ++iz){
                            all_grid.substance_density[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                +=zero_velocity_correction_values[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                        }
                    }
                }
            }
    	}
        // semi-Lagrangian を密度の移流に使う場合
    	else {
    		advect_density_semi_lagrangian_3D(
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

    // セルの体積によって質量密度を補正
    void correct_mass_and_volume_by_pressure_solve_3D(
        std::vector<MY_FLOAT_TYPE> &corrected_values,
        const std::vector<MY_FLOAT_TYPE> &cell_volume,
        const Grid_3D& all_grid
    ){
        //グリッドの総数
        const int N=all_grid.Grid_num_x
                   *all_grid.Grid_num_y
                   *all_grid.Grid_num_z;
        //係数行列の計算
    //    linear_algebra::sparse_matrix A(N,N);
        linear_algebra::sparse_matrix A(N - 1, N - 1);
    //    linear_algebra::sparse_matrix_with_diagonal_element A(N,N);
    //    for(size_t i_xyz = 0; i_xyz < N; ++i_xyz){
        for(size_t i_xyz = 0; i_xyz < N - 1; ++i_xyz){
            int iz = i_xyz % all_grid.Grid_num_z;
            int iy = ((i_xyz - iz) / all_grid.Grid_num_z) % all_grid.Grid_num_y;
            int ix = (i_xyz - iz - iy * all_grid.Grid_num_z) / (all_grid.Grid_num_y * all_grid.Grid_num_z);
            //非対角項
            if(ix != 0){
                int index_0 = get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z);
                int index_1 = get_voxel_center_index_3D(ix-1, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z);
                if((index_0 != N - 1) && (index_1 != N - 1)){
                    A.input_element(index_0, index_1, 1);
                }
            }
            //非対角項
            if(iy != 0){
                int index_0 = get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z);
                int index_1 = get_voxel_center_index_3D(ix, iy-1, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z);
                if((index_0 != N - 1) && (index_1 != N - 1)){
                    A.input_element(index_0, index_1, 1);
                }
            }
            //非対角項
            if(iz != 0){
                int index_0 = get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z);
                int index_1 = get_voxel_center_index_3D(ix, iy, iz-1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z);
                if((index_0 != N - 1) && (index_1 != N - 1)){
                    A.input_element(index_0, index_1, 1);
                }
            }
            //対角項
            MY_FLOAT_TYPE diagonal_value = -6;
            if(ix == 0){
                diagonal_value += 1;
            }
            if(iy == 0){
                diagonal_value += 1;
            }
            if(iz == 0){
                diagonal_value += 1;
            }
            if(iz == all_grid.Grid_num_z - 1){
                diagonal_value += 1;
            }
            if(iy == all_grid.Grid_num_y - 1){
                diagonal_value += 1;
            }
            if(ix == all_grid.Grid_num_x - 1){
                diagonal_value += 1;
            }
            A.input_element(
                get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z),
                get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z),
                diagonal_value
            );
            //非対角項
            if(iz != all_grid.Grid_num_z-1){
                int index_0 = get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z);
                int index_1 = get_voxel_center_index_3D(ix, iy, iz+1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z);
                if((index_0 != N - 1) && (index_1 != N - 1)){
                    A.input_element(index_0, index_1, 1);
                }
            }
            //非対角項
            if(iy != all_grid.Grid_num_y-1){
                int index_0 = get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z);
                int index_1 = get_voxel_center_index_3D(ix, iy+1, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z);
                if((index_0 != N - 1) && (index_1 != N - 1)){
                    A.input_element(index_0, index_1, 1);
                }
            }
            //非対角項
            if(ix != all_grid.Grid_num_x-1){
                int index_0 = get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z);
                int index_1 = get_voxel_center_index_3D(ix+1, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z);
                if((index_0 != N - 1) && (index_1 != N - 1)){
                    A.input_element(index_0, index_1, 1);
                }
            }
        }

		//連立方程式の右辺のベクトルの計算
//		std::vector<MY_FLOAT_TYPE> b(N);
		std::vector<MY_FLOAT_TYPE> b(N - 1);
//		for(int i_xy = 0; i_xy < N; ++i_xy){
        for(size_t i_xyz = 0; i_xyz < N - 1; ++i_xyz){
            int iz = i_xyz % all_grid.Grid_num_z;
            int iy = ((i_xyz - iz) / all_grid.Grid_num_z) % all_grid.Grid_num_y;
            int ix = (i_xyz - iz - iy * all_grid.Grid_num_z) / (all_grid.Grid_num_y * all_grid.Grid_num_z);

            b[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
				= -(cell_volume[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
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
		std::vector<MY_FLOAT_TYPE> before_correct_mass(all_grid.Grid_num_x * all_grid.Grid_num_y * all_grid.Grid_num_z);
		std::vector<MY_FLOAT_TYPE> before_correct_volume(all_grid.Grid_num_x * all_grid.Grid_num_y * all_grid.Grid_num_z);
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
    				before_correct_mass[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					= corrected_values[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					* all_grid._cell_volume;
    				before_correct_volume[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					= cell_volume[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                }
			}
		}
		//// 体積の勾配によって質量と体積を拡散させる
		MY_FLOAT_TYPE diffusion_coefficient = 1.0;
        // x 方向の拡散
		for (int ix = 1; ix < all_grid.Grid_num_x; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
    				//フラックス
    				MY_FLOAT_TYPE volume_flux
    					= -(pressure[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					  - pressure[get_voxel_center_index_3D(ix - 1, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])
    					 / all_grid._cell_volume;
    				// 面での密度
    				MY_FLOAT_TYPE density_at_face
    						= (before_correct_mass[get_voxel_center_index_3D(ix - 1, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    						+  before_correct_mass[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])
    						/ (2.0 * all_grid._cell_volume);
    				// 面での体積
    				MY_FLOAT_TYPE volume_at_face
    					= (before_correct_volume[get_voxel_center_index_3D(ix - 1, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					+  before_correct_volume[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])
    					/ 2.0;
    				//密度の拡散
    				corrected_values[get_voxel_center_index_3D(ix - 1, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					-= diffusion_coefficient * volume_flux * density_at_face;
    				corrected_values[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					+= diffusion_coefficient * volume_flux * density_at_face;
    				//体積の拡散
//    				cell_volume[get_voxel_center_index_3D(ix - 1, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
//    					-= diffusion_coefficient * volume_flux * volume_at_face;
//    				cell_volume[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
//    					+= diffusion_coefficient * volume_flux * volume_at_face;
                }
			}
		}
        // y 方向の拡散
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
            for (int iy = 1; iy < all_grid.Grid_num_y; iy++) {
                for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
    				//フラックス
    				MY_FLOAT_TYPE volume_flux
    					= -(pressure[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					  - pressure[get_voxel_center_index_3D(ix, iy - 1, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])
    					 / all_grid._cell_volume;
    				// 面での密度
    				MY_FLOAT_TYPE density_at_face
    						= (before_correct_mass[get_voxel_center_index_3D(ix, iy - 1, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    						+  before_correct_mass[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])
    						/ (2.0 * all_grid._cell_volume);
    				// 面での体積
    				MY_FLOAT_TYPE volume_at_face
    					= (before_correct_volume[get_voxel_center_index_3D(ix, iy - 1, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					+  before_correct_volume[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])
    					/ 2.0;
    				//密度の拡散
    				corrected_values[get_voxel_center_index_3D(ix, iy - 1, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					-= diffusion_coefficient * volume_flux * density_at_face;
    				corrected_values[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					+= diffusion_coefficient * volume_flux * density_at_face;
    				//体積の拡散
//    				cell_volume[get_voxel_center_index_3D(ix, iy - 1, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
//    					-= diffusion_coefficient * volume_flux * volume_at_face;
//    				cell_volume[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
//    					+= diffusion_coefficient * volume_flux * volume_at_face;
                }
			}
		}
        // z 方向の拡散
		for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                for (int iz = 1; iz < all_grid.Grid_num_z; iz++) {
    				//フラックス
    				MY_FLOAT_TYPE volume_flux
    					= -(pressure[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					  - pressure[get_voxel_center_index_3D(ix, iy, iz - 1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])
    					 / all_grid._cell_volume;
    				// 面での密度
    				MY_FLOAT_TYPE density_at_face
    						= (before_correct_mass[get_voxel_center_index_3D(ix, iy, iz - 1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    						+  before_correct_mass[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])
    						/ (2.0 * all_grid._cell_volume);
    				// 面での体積
    				MY_FLOAT_TYPE volume_at_face
    					= (before_correct_volume[get_voxel_center_index_3D(ix, iy, iz - 1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					+  before_correct_volume[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])
    					/ 2.0;
    				//密度の拡散
    				corrected_values[get_voxel_center_index_3D(ix, iy, iz - 1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					-= diffusion_coefficient * volume_flux * density_at_face;
    				corrected_values[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
    					+= diffusion_coefficient * volume_flux * density_at_face;
    				//体積の拡散
//    				cell_volume[get_voxel_center_index_3D(ix, iy, iz - 1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
//    					-= diffusion_coefficient * volume_flux * volume_at_face;
//    				cell_volume[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
//    					+= diffusion_coefficient * volume_flux * volume_at_face;
                }
			}
		}
    }


//流体の速度場によってsubstanceが運ばれる項(流体のadvect項に相当)
//速度場を時間 -dt だけバックトレースして計算する
void advect_density_semi_lagrangian_3D(
    Grid_3D& all_grid,
    const MY_FLOAT_TYPE time_step_length,
    const bool use_MacCormack_scheme,
    const bool use_clamping_in_MacCormack_scheme,
    const int num_gauss_quad_boundary,
    const int num_gauss_quad_bulk,
    const std::string interpolation_method
){
    std::vector<MY_FLOAT_TYPE> substance_density_after_advect(all_grid.Grid_num_x * all_grid.Grid_num_y * all_grid.Grid_num_z);
    for(int ix=0;ix<all_grid.Grid_num_x;ix++){
        for(int iy=0;iy<all_grid.Grid_num_y;iy++){
            for(int iz=0;iz<all_grid.Grid_num_z;iz++){
                // MacCormack 法を使う場合
                if (use_MacCormack_scheme) {
//                    std::cout << "calc MacCormack"<<std::endl;
                    //バックトレース前の座標 ( 0.5 は cell vertex から cell center までのオフセット )
                    VEC3_TYPE before_backtrace_position((ix + 0.5) * all_grid._cell_length, (iy + 0.5) * all_grid._cell_length, (iz + 0.5) * all_grid._cell_length);
                    //バックトレース前の座標での速度の計算
                    VEC3_TYPE velocity_at_before_backtrace_position = all_grid.calc_interpolated_velocity(before_backtrace_position, interpolation_method);
                    //バックトレース後の座標 ( 0.5 は cell vertex から cell center までのオフセット )
                    VEC3_TYPE after_backtrace_position = before_backtrace_position - time_step_length * velocity_at_before_backtrace_position;
                    // バックトレース先の密度を補間で求める
                    MY_FLOAT_TYPE density_at_advected_position = all_grid.calc_substance_density_by_interpolation(after_backtrace_position, interpolation_method);
                    ////// エラーの計算
                    //バックトレース後の座標での速度
                    VEC3_TYPE velocity_at_after_backtrace_position = all_grid.calc_interpolated_velocity(after_backtrace_position, interpolation_method);
                    //バックトレース後の座標を時間を巻き戻すように逆向きにトレースする
                    VEC3_TYPE forwardtrace_position_of_after_backtrace_position = after_backtrace_position + time_step_length * velocity_at_after_backtrace_position;
                    MY_FLOAT_TYPE error
                        = (all_grid.calc_substance_density_by_interpolation(before_backtrace_position, interpolation_method)
                            - all_grid.calc_substance_density_by_interpolation(forwardtrace_position_of_after_backtrace_position, interpolation_method))
                        / 2.0;
                    // エラー補正後の値
                    substance_density_after_advect[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = density_at_advected_position + error;
                    // clamping
                    if (use_clamping_in_MacCormack_scheme) {
                        substance_density_after_advect[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        = clamp_substance_density_of_MacCormack(all_grid, substance_density_after_advect[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)], after_backtrace_position, interpolation_method);
                    }
//                    std::cout << "getchar()"<<std::endl;
//                    getchar();
                }
                // semi-Lagrangian 法を使う場合
                else {
                    //バックトレース前の座標 ( 0.5 は cell vertex から cell center までのオフセット )
                    VEC3_TYPE before_backtrace_position((ix + 0.5) * all_grid._cell_length, (iy + 0.5) * all_grid._cell_length, (iz + 0.5) * all_grid._cell_length);
                    //バックトレース前の座標での速度の計算
                    VEC3_TYPE velocity_at_before_backtrace_position = all_grid.calc_interpolated_velocity(before_backtrace_position, interpolation_method);
                    //バックトレース後の座標 ( 0.5 は cell vertex から cell center までのオフセット )
                    VEC3_TYPE after_backtrace_position = before_backtrace_position - time_step_length * velocity_at_before_backtrace_position;
                    // バックトレース先の密度を補間で求める
                    MY_FLOAT_TYPE density_at_advected_position = all_grid.calc_substance_density_by_interpolation(after_backtrace_position, interpolation_method);
/*
                    MY_FLOAT_TYPE velocity_x
                        =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    MY_FLOAT_TYPE velocity_y
                        =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    MY_FLOAT_TYPE velocity_z
                        =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz+1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    //バックトレース先の座標
                    MY_FLOAT_TYPE advected_x=ix-((velocity_x*time_step_length)/all_grid._cell_length);
                    MY_FLOAT_TYPE advected_y=iy-((velocity_y*time_step_length)/all_grid._cell_length);
                    MY_FLOAT_TYPE advected_z=iz-((velocity_z*time_step_length)/all_grid._cell_length);
*/
                    //バックトレース先の座標のindex
                    //バックトレース先の速度を補間する
                    substance_density_after_advect[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = density_at_advected_position;
//                    substance_density_after_advect[ix][iy][iz]=all_grid.calc_substance_density_by_interpolation(VEC3_TYPE(advected_x, advected_y, advected_z));
//                    substance_density_after_advect[ix][iy][iz]=linear_interpolation_substances_3D(advected_x, advected_y, advected_z, all_grid);
//                    substance_density_after_advect[ix][iy][iz]=monotonic_cubic_substances_3D(advected_x, advected_y, advected_z, all_grid);
                }
            }
        }
    }
    //計算結果をコピー
    for(int ix=0;ix<all_grid.Grid_num_x;ix++){
        for(int iy=0;iy<all_grid.Grid_num_y;iy++){
            for(int iz=0;iz<all_grid.Grid_num_z;iz++){
                all_grid.substance_density[get_voxel_center_index_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                    =substance_density_after_advect[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
            }
        }
    }
}

//密度場の移流項を flux advection で計算
void advect_density_flux_advection_3D_std_thread(
    Grid_3D& all_grid,
    const MY_FLOAT_TYPE time_step_length,
//    const bool set_zero_normal_velocity_at_boundary,
    const bool use_MacCormack_scheme,
    const bool use_clamping_in_MacCormack_scheme,
    const int num_gauss_quad_boundary,
    const int num_gauss_quad_bulk,
    const bool enable_eliminate_zero_velocity_false_diffusion,
    const std::string split_method,
    const std::string integral_method,
    const int num_gauss_quadrature_points,
    std::vector<MY_FLOAT_TYPE> &advected_values,
    const std::string interpolation_method,
    const std::string interpolation_method_in_calc_psi,
    const bool use_zero_velocity_for_backtrace,
    const bool use_integral,
    const int num_gauss_quadrature_point_for_integrate_density,
    const bool enable_cell_volume_correction,
    const MY_FLOAT_TYPE minimal_area
){
    const size_t N_xyzface = all_grid.Grid_num_x * all_grid.Grid_num_y * all_grid.Grid_num_z * 6;
    ////時間計測開始
    auto start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
    ////-------------------------------------
    ////cell face での psiを計算する
    ////-------------------------------------
    calc_psi_on_cell_face_from_density_on_cell_center_3D(
        all_grid,
        all_grid.psi_substance_density_cell_face_x,
        all_grid.psi_substance_density_cell_face_y,
        all_grid.psi_substance_density_cell_face_z,
        advected_values,
        all_grid.Grid_num_x,
        all_grid.Grid_num_y,
        all_grid.Grid_num_z,
        all_grid._cell_length,
        interpolation_method_in_calc_psi,
        num_gauss_quadrature_point_for_integrate_density,
        use_integral
    );
    ////時間計測終了
    auto end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
    auto duration_time = end_time - start_time;        // 要した時間を計算
    auto duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
    // 要した時間をミリ秒（1/1000秒）に変換して表示
    std::cout <<"calc_psi_substance_density_cell_face_3D(): "<< duration_time_msec << " milli sec \n";
    //計算結果を格納する変数
    std::vector<MY_FLOAT_TYPE> substance_density_after_advect(all_grid.Grid_num_x * all_grid.Grid_num_y * all_grid.Grid_num_z);
    for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
        for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
            for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
                substance_density_after_advect[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
            }
        }
    }
    //セル体積の計算結果を格納する変数
    std::vector<MY_FLOAT_TYPE> cell_volume_after_advect(all_grid.Grid_num_x * all_grid.Grid_num_y * all_grid.Grid_num_z);
    for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
        for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
            for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
                cell_volume_after_advect[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
            }
        }
    }

    //psiの計算に積分を使う場合
    if(use_integral){
        ////全てのグリッドセルの面をバックトレースし、ポリゴンへ切断する
        ////時間計測開始
        start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
        // 求積点のリスト (density からの寄与)
        quadrature_point_vector quadrature_point_list_interpolate_density;
        // 求積点のリスト (psi からの寄与)
        quadrature_point_vector quadrature_point_list_interpolate_psi;
        // セル体積を計算するための求積点のリスト
        quadrature_point_vector quadrature_point_list_for_calc_cell_volume;
        //全ての面のバックトレースを行い、求積点のリストを作る(並列化したバージョン)
        parallelize_backtrace_and_calc_all_quadrature_point_list_for_integrate_density(
            backtrace_and_calc_all_quadrature_point_list_for_integrate_density,
            -1,
            N_xyzface,
            quadrature_point_list_interpolate_density,
            quadrature_point_list_interpolate_psi,
            quadrature_point_list_for_calc_cell_volume,
            use_MacCormack_scheme,
            all_grid,
            time_step_length,
//            set_zero_normal_velocity_at_boundary,
            split_method,
//            false,
            integral_method,
            num_gauss_quadrature_points,
            interpolation_method,
            use_zero_velocity_for_backtrace,
    	    num_gauss_quadrature_point_for_integrate_density,
            enable_cell_volume_correction,
            "density",
    		minimal_area
        );
        ////時間計測終了
        end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
        duration_time = end_time - start_time;        // 要した時間を計算
        duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
        // 要した時間をミリ秒（1/1000秒）に変換して表示
        std::cout <<"backtrace and split cell_face_3D: "<< duration_time_msec << " milli sec \n";

        //全ポリゴンから質量への寄与を計算
        ////時間計測開始
        start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
        ////
        ////全求積点から質量密度を求める
        ////
        //並列化バージョン
        parallelize_calc_density_in_cell_from_quadrature_point_list_for_integrate_density (
            calc_density_in_cell_from_quadrature_point_list_for_integrate_density,
            quadrature_point_list_interpolate_density,
            quadrature_point_list_interpolate_psi,
            all_grid,
            substance_density_after_advect,
            interpolation_method,
            advected_values
        );
        ////時間計測終了
        end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
        duration_time = end_time - start_time;        // 要した時間を計算
        duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
        // 要した時間をミリ秒（1/1000秒）に変換して表示
        std::cout <<"calc mass from all_splitted_faces: "<< duration_time_msec << " milli sec \n";

        if(enable_cell_volume_correction){
            ////
            ////全求積点からセル体積を求める
            ////
            ////時間計測開始
            start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
            //並列化バージョン
            parallelize_calc_cell_volumes_from_quadrature_point_list (
                calc_cell_volumes_from_quadrature_point_list,
                quadrature_point_list_for_calc_cell_volume,
                all_grid,
                cell_volume_after_advect
            );
            ////時間計測終了
            end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
            duration_time = end_time - start_time;        // 要した時間を計算
            duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
            // 要した時間をミリ秒（1/1000秒）に変換して表示
            std::cout <<"calc cell volumes: "<< duration_time_msec << " milli sec \n";
        }
    }
    //psiの計算に補間を使う場合
    else{
        ////全てのグリッドセルの面をバックトレースし、ポリゴンへ切断する
        ////時間計測開始
        start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
        // 求積点のリスト
        quadrature_point_vector quadrature_point_list;
        //全ての面のバックトレースを行い、求積点のリストを作る(並列化したバージョン)
        parallelize_backtrace_and_calc_all_quadrature_point_list(
            backtrace_and_calc_all_quadrature_point_list,
            -1,
            N_xyzface,
            quadrature_point_list,
            use_MacCormack_scheme,
            all_grid,
            time_step_length,
//            set_zero_normal_velocity_at_boundary,
            split_method,
//            false,
            integral_method,
            num_gauss_quadrature_points,
            interpolation_method,
            use_zero_velocity_for_backtrace,
    		minimal_area
        );

        ////時間計測終了
        end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
        duration_time = end_time - start_time;        // 要した時間を計算
        duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
        // 要した時間をミリ秒（1/1000秒）に変換して表示
        std::cout <<"backtrace and split cell_face_3D: "<< duration_time_msec << " milli sec \n";

        //全ポリゴンから質量への寄与を計算
        ////時間計測開始
        start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
        ////
        ////全求積点から質量密度を求める
        ////
        //並列化バージョン
        parallelize_calc_density_in_cell_from_quadrature_point_list (
            calc_density_in_cell_from_quadrature_point_list,
            quadrature_point_list,
            all_grid,
            substance_density_after_advect,
            interpolation_method
        );
        ////時間計測終了
        end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
        duration_time = end_time - start_time;        // 要した時間を計算
        duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
        // 要した時間をミリ秒（1/1000秒）に変換して表示
        std::cout <<"calc mass from all_splitted_faces: "<< duration_time_msec << " milli sec \n";

        if(enable_cell_volume_correction){
            ////
            ////全求積点からセル体積を求める
            ////
            ////時間計測開始
            start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
            //並列化バージョン
            parallelize_calc_cell_volumes_from_quadrature_point_list (
                calc_cell_volumes_from_quadrature_point_list,
                quadrature_point_list,
                all_grid,
                cell_volume_after_advect
            );
            ////時間計測終了
            end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
            duration_time = end_time - start_time;        // 要した時間を計算
            duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
            // 要した時間をミリ秒（1/1000秒）に変換して表示
            std::cout <<"calc cell volumes: "<< duration_time_msec << " milli sec \n";
        }
    }

    ////時間計測開始
    start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
    //計算結果をコピー
    for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
        for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
            for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
                advected_values[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                    = substance_density_after_advect[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
            }
        }
    }
    if(enable_cell_volume_correction){
        //計算結果をコピー
        for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
                    all_grid.cell_volume_cell_center[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        = cell_volume_after_advect[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                }
            }
        }
    }

    ////時間計測終了
    end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
    duration_time = end_time - start_time;        // 要した時間を計算
    duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
    // 要した時間をミリ秒（1/1000秒）に変換して表示
    std::cout <<"copy results(all_grid.substance_density): "<< duration_time_msec << " milli sec \n";
}

//------------------------------------------------------------------------------
// 並列化関連
//------------------------------------------------------------------------------

//backtrace_and_split_cell_faces()をstd::thread で並列に実行する関数
void parallelize_backtrace_and_calc_all_quadrature_point_list(
    auto func,
    const int dim,
    size_t size,
    quadrature_point_vector &quadrature_point_list,
    bool use_MacCormack_scheme,
    const Grid_3D &all_grid,
    MY_FLOAT_TYPE time_step_length,
//    bool set_zero_normal_velocity_at_boundary,
    std::string split_method,
//    bool is_zero_velocity_false_diffusion_correction_mode,
    const std::string integral_method,
    const int num_gauss_quadrature_points,
    const std::string interpolation_method,
    const bool use_zero_velocity_for_backtrace,
    const MY_FLOAT_TYPE minimal_area
) {
    const unsigned num_threads (std::thread::hardware_concurrency());
    const size_t chunk_size (std::ceil(size/(MY_FLOAT_TYPE)num_threads));

    std::vector<std::thread> threads;
    int i_thread = 0;
    if( size <= num_threads ) {
        for( size_t n=0; n<size; ++n ) {
            threads.push_back(std::thread(
            func,
            dim,
            std::ref(quadrature_point_list),
            n,
            n+1,
            use_MacCormack_scheme,
            std::ref(all_grid),
            time_step_length,
//            set_zero_normal_velocity_at_boundary,
            split_method,
//            is_zero_velocity_false_diffusion_correction_mode,
            integral_method,
            num_gauss_quadrature_points,
            i_thread,
            interpolation_method,
            use_zero_velocity_for_backtrace,
    		minimal_area
        ));
        ++i_thread;
        }
    } else if( num_threads*chunk_size != size ) {
        const size_t end_n (size+num_threads-num_threads*chunk_size);
        const size_t end_m (num_threads*chunk_size-size);
        size_t end (0);
        for( size_t n=0; n<end_n; ++n ) {
            const size_t start (end);
            end = start+chunk_size;
            threads.push_back(std::thread(
                func,
                dim,
                std::ref(quadrature_point_list),
                start,
                end,
                use_MacCormack_scheme,
                std::ref(all_grid),
                time_step_length,
//                set_zero_normal_velocity_at_boundary,
                split_method,
//                is_zero_velocity_false_diffusion_correction_mode,
                integral_method,
                num_gauss_quadrature_points,
                i_thread,
                interpolation_method,
                use_zero_velocity_for_backtrace,
        		minimal_area
        ));
        ++i_thread;
        }
        for( size_t m=0; m<end_m; ++m ) {
            const size_t start (end);
            end = start+chunk_size-1;
            threads.push_back(std::thread(
                func,
                dim,
                std::ref(quadrature_point_list),
                start,
                end,
                use_MacCormack_scheme,
                std::ref(all_grid),
                time_step_length,
//                set_zero_normal_velocity_at_boundary,
                split_method,
//                is_zero_velocity_false_diffusion_correction_mode,
                integral_method,
                num_gauss_quadrature_points,
                i_thread,
                interpolation_method,
                use_zero_velocity_for_backtrace,
        		minimal_area
        ));
        ++i_thread;
        }
    } else {
        for( size_t n=0; n<num_threads; ++n ) {
            const size_t start (n*chunk_size);
            threads.push_back(std::thread(
                func,
                dim,
                std::ref(quadrature_point_list),
                start,
                start+chunk_size,
                use_MacCormack_scheme,
                std::ref(all_grid),
                time_step_length,
//                set_zero_normal_velocity_at_boundary,
                split_method,
//                is_zero_velocity_false_diffusion_correction_mode,
                integral_method,
                num_gauss_quadrature_points,
                i_thread,
                interpolation_method,
                use_zero_velocity_for_backtrace,
        		minimal_area
        ));
        ++i_thread;
        }
    }
    for( auto& t : threads ) t.join();
}

//calc_density_in_cell_from_quadrature_point_list()をstd::thread で並列に実行する関数
void parallelize_calc_density_in_cell_from_quadrature_point_list (
    auto func,
    const quadrature_point_vector &quadrature_point_list,
    const Grid_3D &all_grid,
    std::vector<MY_FLOAT_TYPE> &substance_density_after_advect,
    const std::string interpolation_method
) {
    const unsigned num_threads (std::thread::hardware_concurrency());

    std::vector<std::thread> threads;
    for(int i_thread = 0; i_thread < num_threads; ++i_thread){
        threads.push_back(std::thread(
            func,
            std::ref(quadrature_point_list),
//            0,
//            quadrature_point_list._quadtarure_position[i_thread].size(),
            std::ref(all_grid),
            std::ref(substance_density_after_advect),
            i_thread,
            interpolation_method
        ));
    }
    for( auto& t : threads ) t.join();
}

//calc_cell_volumes_from_quadrature_point_list()をstd::thread で並列に実行する関数
void parallelize_calc_cell_volumes_from_quadrature_point_list (
    auto func,
    const quadrature_point_vector &quadrature_point_list,
    const Grid_3D &all_grid,
    std::vector<MY_FLOAT_TYPE> &cell_volume_after_advect
) {
    const unsigned num_threads (std::thread::hardware_concurrency());

    std::vector<std::thread> threads;
    for(int i_thread = 0; i_thread < num_threads; ++i_thread){
        threads.push_back(std::thread(
            func,
            std::ref(quadrature_point_list),
//            0,
//            quadrature_point_list._quadtarure_position[i_thread].size(),
            std::ref(all_grid),
            std::ref(cell_volume_after_advect),
            i_thread
        ));
    }
    for( auto& t : threads ) t.join();
}


//backtrace_and_split_cell_faces()をstd::thread で並列に実行する関数
void parallelize_backtrace_and_calc_all_quadrature_point_list_for_integrate_density(
    auto func,
    const int dim,
    size_t size,
    quadrature_point_vector &quadrature_point_list_for_interpolate_density,
    quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
    quadrature_point_vector &quadrature_point_list_for_calc_cell_volume,
    bool use_MacCormack_scheme,
    const Grid_3D &all_grid,
    MY_FLOAT_TYPE time_step_length,
//    bool set_zero_normal_velocity_at_boundary,
    std::string split_method,
//    bool is_zero_velocity_false_diffusion_correction_mode,
    const std::string integral_method,
    const int num_gauss_quadrature_points,
    const std::string interpolation_method,
    const bool use_zero_velocity_for_backtrace,
    const int num_gauss_quadrature_point_for_integrate_density,
    const bool enable_cell_volume_correction,
    const std::string backtrace_usage,
    const MY_FLOAT_TYPE minimal_area
) {
    const unsigned num_threads (std::thread::hardware_concurrency());
    const size_t chunk_size (std::ceil(size/(MY_FLOAT_TYPE)num_threads));

    std::vector<std::thread> threads;
    int i_thread = 0;
    if( size <= num_threads ) {
        for( size_t n=0; n<size; ++n ) {
            threads.push_back(std::thread(
            func,
            dim,
            std::ref(quadrature_point_list_for_interpolate_density),
            std::ref(quadrature_point_list_for_interpolate_psi),
            std::ref(quadrature_point_list_for_calc_cell_volume),
            n,
            n+1,
            use_MacCormack_scheme,
            std::ref(all_grid),
            time_step_length,
//            set_zero_normal_velocity_at_boundary,
            split_method,
//            is_zero_velocity_false_diffusion_correction_mode,
            integral_method,
            num_gauss_quadrature_points,
            i_thread,
            interpolation_method,
            use_zero_velocity_for_backtrace,
    	    num_gauss_quadrature_point_for_integrate_density,
            enable_cell_volume_correction,
            backtrace_usage,
    		minimal_area
        ));
        ++i_thread;
        }
    } else if( num_threads*chunk_size != size ) {
        const size_t end_n (size+num_threads-num_threads*chunk_size);
        const size_t end_m (num_threads*chunk_size-size);
        size_t end (0);
        for( size_t n=0; n<end_n; ++n ) {
            const size_t start (end);
            end = start+chunk_size;
            threads.push_back(std::thread(
                func,
                dim,
                std::ref(quadrature_point_list_for_interpolate_density),
                std::ref(quadrature_point_list_for_interpolate_psi),
                std::ref(quadrature_point_list_for_calc_cell_volume),
                start,
                end,
                use_MacCormack_scheme,
                std::ref(all_grid),
                time_step_length,
//                set_zero_normal_velocity_at_boundary,
                split_method,
//                is_zero_velocity_false_diffusion_correction_mode,
                integral_method,
                num_gauss_quadrature_points,
                i_thread,
                interpolation_method,
                use_zero_velocity_for_backtrace,
        	    num_gauss_quadrature_point_for_integrate_density,
                enable_cell_volume_correction,
                backtrace_usage,
        		minimal_area
        ));
        ++i_thread;
        }
        for( size_t m=0; m<end_m; ++m ) {
            const size_t start (end);
            end = start+chunk_size-1;
            threads.push_back(std::thread(
                func,
                dim,
                std::ref(quadrature_point_list_for_interpolate_density),
                std::ref(quadrature_point_list_for_interpolate_psi),
                std::ref(quadrature_point_list_for_calc_cell_volume),
                start,
                end,
                use_MacCormack_scheme,
                std::ref(all_grid),
                time_step_length,
//                set_zero_normal_velocity_at_boundary,
                split_method,
//                is_zero_velocity_false_diffusion_correction_mode,
                integral_method,
                num_gauss_quadrature_points,
                i_thread,
                interpolation_method,
                use_zero_velocity_for_backtrace,
        	    num_gauss_quadrature_point_for_integrate_density,
                enable_cell_volume_correction,
                backtrace_usage,
        		minimal_area
        ));
        ++i_thread;
        }
    } else {
        for( size_t n=0; n<num_threads; ++n ) {
            const size_t start (n*chunk_size);
            threads.push_back(std::thread(
                func,
                dim,
                std::ref(quadrature_point_list_for_interpolate_density),
                std::ref(quadrature_point_list_for_interpolate_psi),
                std::ref(quadrature_point_list_for_calc_cell_volume),
                start,
                start+chunk_size,
                use_MacCormack_scheme,
                std::ref(all_grid),
                time_step_length,
//                set_zero_normal_velocity_at_boundary,
                split_method,
//                is_zero_velocity_false_diffusion_correction_mode,
                integral_method,
                num_gauss_quadrature_points,
                i_thread,
                interpolation_method,
                use_zero_velocity_for_backtrace,
        	    num_gauss_quadrature_point_for_integrate_density,
                enable_cell_volume_correction,
                backtrace_usage,
        		minimal_area
        ));
        ++i_thread;
        }
    }
    for( auto& t : threads ) t.join();
}

//calc_density_in_cell_from_quadrature_point_list()をstd::thread で並列に実行する関数
void parallelize_calc_density_in_cell_from_quadrature_point_list_for_integrate_density (
    auto func,
    const quadrature_point_vector &quadrature_point_list_for_interpolate_density,
    const quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
    const Grid_3D &all_grid,
    std::vector<MY_FLOAT_TYPE> &substance_density_after_advect,
    const std::string interpolation_method,
    const std::vector<MY_FLOAT_TYPE> &advected_values
) {
    const unsigned num_threads (std::thread::hardware_concurrency());

    std::vector<std::thread> threads;
    for(int i_thread = 0; i_thread < num_threads; ++i_thread){
        threads.push_back(std::thread(
            func,
            std::ref(quadrature_point_list_for_interpolate_density),
            std::ref(quadrature_point_list_for_interpolate_psi),
            std::ref(all_grid),
            std::ref(substance_density_after_advect),
            i_thread,
            interpolation_method,
            advected_values
        ));
    }
    for( auto& t : threads ) t.join();
}


}
