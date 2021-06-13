#include "calc_pressure_3d.h"

#include <chrono>//時間計測用
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

#include "linear_solver.h"
#include "utils.h"

namespace smoke_simulation{
    //Poisson eq. を解くことによる圧力の計算
    void calc_pressure_3D(Grid_3D& all_grid){
        //グリッドの総数
        const int N = all_grid.Grid_num_x
                    * all_grid.Grid_num_y
                    * all_grid.Grid_num_z;
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
    //    std::vector<MY_FLOAT_TYPE> b(N);
        std::vector<MY_FLOAT_TYPE> b(N - 1);

    //    for(size_t i_xyz = 0; i_xyz < N; ++i_xyz){
        for(size_t i_xyz = 0; i_xyz < N - 1; ++i_xyz){
            int iz = i_xyz % all_grid.Grid_num_z;
            int iy = ((i_xyz - iz) / all_grid.Grid_num_z) % all_grid.Grid_num_y;
            int ix = (i_xyz - iz - iy * all_grid.Grid_num_z) / (all_grid.Grid_num_y * all_grid.Grid_num_z);
    //                b[get_voxel_center_index_3D(ix, iy, iz)]
    //                   =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz)]-all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]
    //                    +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz)]-all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]
    //                    +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz+1)]-all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)])
    //                   *all_grid._cell_length
    //                   /(time_step_length);
                    b[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                         -all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                         -all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz+1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                         -all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])
                         *all_grid._cell_length;
    //            }
    //        }
        }
    /*
        //時間計測用
        std::chrono::system_clock::time_point  start, end;
        start = std::chrono::system_clock::now(); // 時間計測開始
    */
        //CG法により圧力場を得る
    //    linear_algebra::conjugate_gradient(A, b, all_grid.pressure, N, 10000, 0.0001);
    //    linear_algebra::incomplete_cholesky_conjugate_gradient(A, b, all_grid.pressure, N, 10000, 0.0001);
        //amgcl ライブラリで解く
        typedef amgcl::backend::builtin<MY_FLOAT_TYPE> Backend;
        typedef amgcl::make_solver<
            // Use AMG as preconditioner:
            amgcl::amg<
                Backend,
                amgcl::coarsening::smoothed_aggregation,
    //            amgcl::relaxation::spai0
                amgcl::relaxation::gauss_seidel
                >,
            // And BiCGStab as iterative solver:
    //        amgcl::solver::bicgstab<Backend>
            amgcl::solver::cg<Backend>
            > Solver;
    //    int n=all_grid.Grid_num_x*all_grid.Grid_num_y*all_grid.Grid_num_z;
    //    int n = N;
        int n = N - 1;
        Solver::params prm;
        prm.solver.tol = 1e-4;
    //    Solver solve( std::tie(n, A.cumulative_num_nonzero_element_in_row, A.column_index, A.element_value) ,prm);
        Solver solve( std::tie(n, A.cumulative_num_nonzero_element_in_row, A.column_index, A.element_value));
        int    iters;
        MY_FLOAT_TYPE error;
    //    std::tie(iters, error) = solve(b, all_grid.pressure);

        std::vector<MY_FLOAT_TYPE> truncated_pressure(N - 1);
        std::tie(iters, error) = solve(b, truncated_pressure);
        //係数行列をフルランクにするために落とした最後の要素を復元する
        truncated_pressure.push_back(0.0);
        all_grid.pressure = truncated_pressure;

        //gauss seidel法を使う場合(Aはsparse_matrix_with_diagonal_elementにする)
    //    gauss_seidel(A,b,all_grid.pressure,N,200);
    /*
        end = std::chrono::system_clock::now();  // 時間計測終了
        std::ofstream writing_file;
        writing_file.open("length_of_time_ICCG_64_3d.dat", std::ios::app);
        writing_file << (std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()) << std::endl;
        writing_file.close();
    */
    }

    //pressure gradient termの計算
    void calc_pressure_gradient_term_3D(Grid_3D& all_grid){
        ////時間計測開始
        auto start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
        //圧力の計算
        calc_pressure_3D(all_grid);
        ////時間計測終了
        auto end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
        auto duration_time = end_time - start_time;        // 要した時間を計算
        auto duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
        // 要した時間をミリ秒（1/1000秒）に変換して表示
        std::cout <<"calc_pressure_3D: "<< duration_time_msec << " milli sec \n";

        //pressure gradeint term によって速度場を更新
        for(int ix=1;ix<all_grid.Grid_num_x;ix++){
            for(int iy=0;iy<all_grid.Grid_num_y;iy++){
                for(int iz=0;iz<all_grid.Grid_num_z;iz++){
    //                all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]-=(all_grid.pressure[get_voxel_center_index_3D(ix, iy, iz)]-all_grid.pressure[get_voxel_center_index_3D(ix-1, iy, iz)])
    //                *time_step_length
    //                /(all_grid._cell_length);
                    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                    -=(all_grid.pressure[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]-
                       all_grid.pressure[get_voxel_center_index_3D(ix-1, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])
                        /(all_grid._cell_length);
                }
            }
        }
        for(int ix=0;ix<all_grid.Grid_num_x;ix++){
            for(int iy=1;iy<all_grid.Grid_num_y;iy++){
                for(int iz=0;iz<all_grid.Grid_num_z;iz++){
    //                all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]-=(all_grid.pressure[get_voxel_center_index_3D(ix, iy, iz)]-all_grid.pressure[get_voxel_center_index_3D(ix, iy-1, iz)])
    //                *time_step_length
    //                /(all_grid._cell_length);
                    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        -=(all_grid.pressure[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                          -all_grid.pressure[get_voxel_center_index_3D(ix, iy-1, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])
                         /(all_grid._cell_length);
                }
            }
        }
        for(int ix=0;ix<all_grid.Grid_num_x;ix++){
            for(int iy=0;iy<all_grid.Grid_num_y;iy++){
                for(int iz=1;iz<all_grid.Grid_num_z;iz++){
    //                all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)]-=(all_grid.pressure[get_voxel_center_index_3D(ix, iy, iz)]-all_grid.pressure[get_voxel_center_index_3D(ix, iy, iz-1)])
    //                *time_step_length
    //                /(all_grid._cell_length);
                    all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        -=(all_grid.pressure[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                          -all_grid.pressure[get_voxel_center_index_3D(ix, iy, iz-1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])
                         /(all_grid._cell_length);
                }
            }
        }
        //壁の圧力によって系の外に逃げ出した運動量を計算して足し合わせる
        //momentum_x のロス
		for (int iy = 0; iy < all_grid.Grid_num_y; ++iy) {
			for (int iz = 0; iz < all_grid.Grid_num_z; ++iz) {
				all_grid.total_loss_of_momentum_x
					-= all_grid.pressure[smoke_simulation::get_voxel_center_index_3D(0, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
					* all_grid._cell_length * all_grid._cell_length;
				all_grid.total_loss_of_momentum_x
					+= all_grid.pressure[smoke_simulation::get_voxel_center_index_3D(all_grid.Grid_num_x - 1, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
					* all_grid._cell_length * all_grid._cell_length;
			}
		}
        //momentum_y のロス
		for (int ix = 0; ix < all_grid.Grid_num_x; ++ix) {
			for (int iz = 0; iz < all_grid.Grid_num_z; ++iz) {
				all_grid.total_loss_of_momentum_y
					-= all_grid.pressure[smoke_simulation::get_voxel_center_index_3D(ix, 0, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
					* all_grid._cell_length * all_grid._cell_length;
				all_grid.total_loss_of_momentum_y
					+= all_grid.pressure[smoke_simulation::get_voxel_center_index_3D(ix, all_grid.Grid_num_y - 1, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
					* all_grid._cell_length * all_grid._cell_length;
			}
		}
        //momentum_z のロス
		for (int ix = 0; ix < all_grid.Grid_num_x; ++ix) {
			for (int iy = 0; iy < all_grid.Grid_num_y; ++iy) {
				all_grid.total_loss_of_momentum_z
					-= all_grid.pressure[smoke_simulation::get_voxel_center_index_3D(ix, iy, 0, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
					* all_grid._cell_length * all_grid._cell_length;
				all_grid.total_loss_of_momentum_z
					+= all_grid.pressure[smoke_simulation::get_voxel_center_index_3D(ix, iy, all_grid.Grid_num_z - 1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
					* all_grid._cell_length * all_grid._cell_length;
			}
		}
    }
} // namespace smoke_simulation
