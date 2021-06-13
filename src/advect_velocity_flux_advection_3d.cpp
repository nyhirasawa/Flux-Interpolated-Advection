#include "advect_velocity_flux_advection_3d.h"

#include <chrono>//時間計測用

#include "calc_psi_3D.h"
#include "utils.h"
#include "backtrace_and_calc_all_quadrature_point_list_for_integrate_density.h"
#include "calc_velocity_from_quadrature_point_list_use_integral.h"
#include "backtrace_and_calc_all_quadrature_point_list.h"
#include "move_substances_3d.h"
#include "calc_velocity_from_quadrature_point_list.h"
#include "parallelize_functions.h"

namespace smoke_simulation{
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
                std::ref(velocity_after_advect),
                i_thread,
                interpolation_method
            ));
        }
        for( auto& t : threads ) t.join();
    }

    //calc_density_in_cell_from_quadrature_point_list()をstd::thread で並列に実行する関数
    void parallelize_calc_velocity_from_quadrature_point_list (
        auto func,
        const int dim,
        const quadrature_point_vector &quadrature_point_list,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &velocity_after_advect,
        const std::string interpolation_method
    ) {
        const unsigned num_threads (std::thread::hardware_concurrency());

        std::vector<std::thread> threads;
        for(int i_thread = 0; i_thread < num_threads; ++i_thread){
            threads.push_back(std::thread(
                func,
                dim,
                std::ref(quadrature_point_list),
    //            0,
    //            quadrature_point_list._quadtarure_position[i_thread].size(),
                std::ref(all_grid),
                std::ref(velocity_after_advect),
                i_thread,
                interpolation_method
            ));
        }
        for( auto& t : threads ) t.join();
    }

    //backtrace_and_split_cell_faces()をstd::thread で並列に実行する関数
    void parallelize_backtrace_and_calc_all_quadrature_point_list_for_integrate_velocity(
        auto func,
        const int dim,
        size_t size,
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
//                    set_zero_normal_velocity_at_boundary,
                    split_method,
//                    is_zero_velocity_false_diffusion_correction_mode,
                    integral_method,
                    num_gauss_quadrature_points,
                    i_thread,
                    interpolation_method,
                    use_zero_velocity_for_backtrace,
                    num_gauss_quadrature_point_for_integrate_density,
                    enable_cell_volume_correction,
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
//                    set_zero_normal_velocity_at_boundary,
                    split_method,
//                    is_zero_velocity_false_diffusion_correction_mode,
                    integral_method,
                    num_gauss_quadrature_points,
                    i_thread,
                    interpolation_method,
                    use_zero_velocity_for_backtrace,
                    num_gauss_quadrature_point_for_integrate_density,
                    enable_cell_volume_correction,
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
//                    set_zero_normal_velocity_at_boundary,
                    split_method,
//                    is_zero_velocity_false_diffusion_correction_mode,
                    integral_method,
                    num_gauss_quadrature_points,
                    i_thread,
                    interpolation_method,
                    use_zero_velocity_for_backtrace,
                    num_gauss_quadrature_point_for_integrate_density,
                    enable_cell_volume_correction,
                    minimal_area
            ));
            ++i_thread;
            }
        }
        for( auto& t : threads ) t.join();
    }

    //calc_velocity_from_quadrature_point_list_use_integral()をstd::thread で並列に実行する関数
    void parallelize_calc_velocity_from_quadrature_point_list_use_integral (
        auto func,
        const int dim,
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
                dim,
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
    ) {
        const int dx = (dim == 0);
        const int dy = (dim == 1);
        const int dz = (dim == 2);
        const size_t N_xyzface = (all_grid.Grid_num_x + dx) * (all_grid.Grid_num_y + dy) * (all_grid.Grid_num_z + dz) * 6;
        ////時間計測開始
        auto start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
        ////-------------------------------------
        ////cell face での psiを計算する
        ////-------------------------------------
        if(dim == 0){
            calc_discrete_psi_velocity_3D(
                dim,
                all_grid,
                all_grid.psi_velocity_x_cell_face,
                advected_values,
                all_grid.Grid_num_x,
                all_grid.Grid_num_y,
                all_grid.Grid_num_z,
                all_grid._cell_length,
                interpolation_method_in_calc_psi,
                num_gauss_quadrature_point_for_integrate_density,
                use_integral
            );
        }
        else if (dim == 1){
            calc_discrete_psi_velocity_3D(
                dim,
                all_grid,
                all_grid.psi_velocity_y_cell_face,
                advected_values,
                all_grid.Grid_num_x,
                all_grid.Grid_num_y,
                all_grid.Grid_num_z,
                all_grid._cell_length,
                interpolation_method_in_calc_psi,
                num_gauss_quadrature_point_for_integrate_density,
                use_integral
            );
        }
        else if (dim == 2){
            calc_discrete_psi_velocity_3D(
                dim,
                all_grid,
                all_grid.psi_velocity_z_cell_face,
                advected_values,
                all_grid.Grid_num_x,
                all_grid.Grid_num_y,
                all_grid.Grid_num_z,
                all_grid._cell_length,
                interpolation_method_in_calc_psi,
                num_gauss_quadrature_point_for_integrate_density,
                use_integral
            );
        }
        ////時間計測終了
        auto end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
        auto duration_time = end_time - start_time;        // 要した時間を計算
        auto duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
        // 要した時間をミリ秒（1/1000秒）に変換して表示
        std::cout <<"calc_psi_substance_density_cell_face_3D(): "<< duration_time_msec << " milli sec \n";
        //計算結果を格納する変数
        std::vector<MY_FLOAT_TYPE> velocity_after_advect((all_grid.Grid_num_x + dx) * (all_grid.Grid_num_y + dy) * (all_grid.Grid_num_z + dz));

        for (int ix = 0; ix < all_grid.Grid_num_x + dx; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y + dy; iy++) {
                for (int iz = 0; iz < all_grid.Grid_num_z + dz; iz++) {
                    velocity_after_advect[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x + dx, all_grid.Grid_num_y + dy, all_grid.Grid_num_z + dz)] = 0.0;
                }
            }
        }
    /*
        //セル体積の計算結果を格納する変数
        std::vector<MY_FLOAT_TYPE> cell_volume_after_advect(all_grid.Grid_num_x * all_grid.Grid_num_y * all_grid.Grid_num_z);
        for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
                    cell_volume_after_advect[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
                }
            }
        }
    */
        //psiの計算に積分を使う場合
        if(use_integral){
            ////全てのグリッドセルの面をバックトレースし、ポリゴンへ切断する
            ////時間計測開始
            start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
            // 求積点のリスト (velocity からの寄与)
            quadrature_point_vector quadrature_point_list_interpolate_velocity;
            // 求積点のリスト (psi からの寄与)
            quadrature_point_vector quadrature_point_list_interpolate_psi;
            // セル体積を計算するための求積点のリスト
            quadrature_point_vector quadrature_point_list_for_calc_cell_volume;
            //全ての面のバックトレースを行い、求積点のリストを作る(並列化したバージョン)
    //        std::string split_method_velocity_x;
    //        if(interpolation_method == "linear" || interpolation_method == "WENO6" || interpolation_method == "WENO6-optimized"){
    //            split_method_velocity_x = "xz-AxisAligned";
    //        }
    //        else if(interpolation_method == "1Dy_linear" || interpolation_method == "1Dy_WENO6"){
    //            split_method_velocity_x = "xz-CellFaceAligned";
    //        }
    //        else{
    //            std::cout << "interpolation method is not correct" <<std::endl;
    //        }
            parallelize_backtrace_and_calc_all_quadrature_point_list_for_integrate_velocity(
                backtrace_and_calc_all_quadrature_point_list_for_integrate_velocity,
                dim,
                N_xyzface,
                quadrature_point_list_interpolate_velocity,
                quadrature_point_list_interpolate_psi,
                quadrature_point_list_for_calc_cell_volume,
                use_MacCormack_scheme,
                all_grid,
                time_step_length,
//                /*set_zero_normal_velocity_at_boundary*/false,
    //            split_method,
                split_method_velocity,
//                false,
                integral_method,
                num_gauss_quadrature_points_on_one_triangle,
                interpolation_method,
                use_zero_velocity_for_backtrace,
        	    num_gauss_quadrature_point_for_integrate_density,
                enable_cell_volume_correction,
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
            parallelize_calc_velocity_from_quadrature_point_list_use_integral (
                calc_velocity_from_quadrature_point_list_use_integral,
                dim,
                quadrature_point_list_interpolate_velocity,
                quadrature_point_list_interpolate_psi,
                all_grid,
                velocity_after_advect,
                interpolation_method,
                advected_values
            );

            ////時間計測終了
            end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
            duration_time = end_time - start_time;        // 要した時間を計算
            duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
            // 要した時間をミリ秒（1/1000秒）に変換して表示
            std::cout <<"calc mass from all_splitted_faces: "<< duration_time_msec << " milli sec \n";
    /*
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
    */
        }
        //psiの計算に補間を使う場合
        else{
            ////全てのグリッドセルの面をバックトレースし、ポリゴンへ切断する
            ////時間計測開始
            start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
            // 求積点のリスト
            quadrature_point_vector quadrature_point_list;
            //全ての面のバックトレースを行い、求積点のリストを作る(並列化したバージョン)
    //        std::string split_method_velocity_x;
    //        if(interpolation_method == "linear" || interpolation_method == "WENO6" || interpolation_method == "WENO6-optimized"){
    //            split_method_velocity_x = "xz-AxisAligned";
    //        }
    //        else if(interpolation_method == "1Dy_linear" || interpolation_method == "1Dy_WENO6"){
    //            split_method_velocity_x = "xz-CellFaceAligned";
    //        }
    //        else{
    //            std::cout << "interpolation method is not correct" <<std::endl;
    //        }
            std::cout<<"parallelize_backtrace_and_calc_all_quadrature_point_list() begin velocity x"<<std::endl;
            parallelize_backtrace_and_calc_all_quadrature_point_list(
                backtrace_and_calc_all_quadrature_point_list,
                dim,
                N_xyzface,
                quadrature_point_list,
                use_MacCormack_scheme,
                all_grid,
                time_step_length,
//                /*set_zero_normal_velocity_at_boundary = */false,
    //            split_method,
                split_method_velocity,
//                false,
                integral_method,
                num_gauss_quadrature_points_on_one_triangle,
                interpolation_method,
                use_zero_velocity_for_backtrace,
                minimal_area
            );
            std::cout<<"parallelize_backtrace_and_calc_all_quadrature_point_list() end velocity x"<<std::endl;

            ////時間計測終了
            end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
            duration_time = end_time - start_time;        // 要した時間を計算
            duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
            // 要した時間をミリ秒（1/1000秒）に変換して表示
            std::cout <<"backtrace and split cell_face_3D velocity x: "<< duration_time_msec << " milli sec \n";

            //全ポリゴンから質量への寄与を計算
            ////時間計測開始
            start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
            ////
            ////全求積点から質量密度を求める
            ////
            //並列化バージョン
            std::cout<<"parallelize_calc_velocity_from_quadrature_point_list() begin velocity x"<<std::endl;
            parallelize_calc_velocity_from_quadrature_point_list (
                calc_velocity_from_quadrature_point_list,
                dim,
                quadrature_point_list,
                all_grid,
                velocity_after_advect,
                interpolation_method
            );
            std::cout<<"parallelize_calc_velocity_from_quadrature_point_list() end velocity x"<<std::endl;

            ////時間計測終了
            end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
            duration_time = end_time - start_time;        // 要した時間を計算
            duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
            // 要した時間をミリ秒（1/1000秒）に変換して表示
            std::cout <<"calc mass from all_splitted_faces: "<< duration_time_msec << " milli sec \n";
    /*
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
    */
        }

        ////時間計測開始
        start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
        //計算結果をコピー
        for (int ix = 0; ix < all_grid.Grid_num_x + dx; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y + dy; iy++) {
                for (int iz = 0; iz < all_grid.Grid_num_z + dz; iz++) {
                    advected_values[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x + dx, all_grid.Grid_num_y + dy, all_grid.Grid_num_z + dz)]
                        = velocity_after_advect[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x + dx, all_grid.Grid_num_y + dy, all_grid.Grid_num_z + dz)];
                }
            }
        }
    /*
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
    */
        ////時間計測終了
        end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
        duration_time = end_time - start_time;        // 要した時間を計算
        duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
        // 要した時間をミリ秒（1/1000秒）に変換して表示
        std::cout <<"copy results(all_grid.substance_density): "<< duration_time_msec << " milli sec \n";
    }

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
    ) {
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
        //// x成分の移流における切断方法の指定
        std::string split_method_velocity_x;
        if(interpolation_method == "linear" || interpolation_method == "WENO6" || interpolation_method == "WENO6-optimized"){
            split_method_velocity_x = "x-AxisAligned_z-CellFaceAligned";
        }
        else if(interpolation_method == "1Dy_linear" || interpolation_method == "1Dy_WENO6"){
            split_method_velocity_x = "x-CellFaceAligned_z-AxisAligned";
        }
        //x成分の移流
        advect_velocity_flux_advection_3D_std_thread(
            0,
            all_grid,
            time_step_length,
            use_MacCormack_scheme,
            use_clamping_in_MacCormack_scheme,
            num_gauss_quad_boundary,
            num_gauss_quad_bulk,
    //		const bool enable_eliminate_zero_velocity_false_diffusion,
            split_method_velocity_x,
            integral_method,
            num_gauss_quadrature_points_on_one_triangle,
            after_advect_values,
            interpolation_method,
            interpolation_method_in_calc_psi,
            /*use_zero_velocity_for_backtrace*/false,
            use_integral,
            num_gauss_quadrature_point_for_integrate_density,
            enable_cell_volume_correction,
            minimal_area
        );
        ////////////////////////////////////////////////////////////
        // 0速度場の拡散補正
        ////////////////////////////////////////////////////////////
        if (enable_eliminate_zero_velocity_false_diffusion) {
            /////////////////////////
            // 速度場のx成分の補正量
            /////////////////////////
            std::vector<MY_FLOAT_TYPE> zero_velocity_correction_values_velocity_x((all_grid.Grid_num_x + 1) * all_grid.Grid_num_y * all_grid.Grid_num_z);
            ////// "バックトレース前の面で定義されたpsiの離散値"を足す
            for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
                for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                    for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
                        zero_velocity_correction_values_velocity_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                            = before_advect_values[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
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
            advect_velocity_flux_advection_3D_std_thread(
                0,
                all_grid,
                time_step_length,
                use_MacCormack_scheme,
                use_clamping_in_MacCormack_scheme,
                num_gauss_quad_boundary,
                num_gauss_quad_bulk,
    //            enable_eliminate_zero_velocity_false_diffusion,
                split_method_velocity_x,
                integral_method,
                num_gauss_quadrature_points_on_one_triangle,
                zero_velocity_correction_values_velocity_x,
                interpolation_method,
                interpolation_method_in_calc_psi,
                /*use_zero_velocity_for_backtrace*/true,
                /*use_integral*/false,
                num_gauss_quadrature_point_for_integrate_density,
                /*enable_cell_volume_correction*/false,
                minimal_area
            );

            // 速度場の離散値との差を求める. これが移流前のエラー
            for(int ix = 0; ix < all_grid.Grid_num_x + 1; ++ix){
                for(int iy = 0; iy < all_grid.Grid_num_y; ++iy){
                    for(int iz = 0; iz < all_grid.Grid_num_z; ++iz){
                        zero_velocity_correction_values_velocity_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                            *= -1.0;
                        zero_velocity_correction_values_velocity_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                            +=before_advect_values[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                    }
                }
            }
            ////-------------------------------------
            //// エラーの移流
            ////-------------------------------------
            // エラーの移流に piecewise-const(一次元補間) を使う
            std::string interpolation_method_of_error_advection;
            if(interpolation_method == "linear" || interpolation_method == "1Dy_linear"){
                interpolation_method_of_error_advection = "1Dy_linear";
            }
            else if(interpolation_method == "WENO6"||interpolation_method == "WENO6-optimized"||interpolation_method == "1Dy_WENO6"){
    //                    interpolation_method_of_error_advection = "1Dy_WENO6";
                interpolation_method_of_error_advection = "1Dy_linear";
            }
    //        std::string split_method_of_error_advection = "xz-AxisAligned";
            std::string split_method_of_error_advection = "x-CellFaceAligned_z-AxisAligned";
            std::string interpolation_method_in_calc_psi = "const";

            // flux advection をエラーの移流に使う場合
            // エラーの移流には1次元補間を使うが、1次元補間ではpsiの計算に補間を使っても積分で厳密に計算しても結果は同じである。
            // なので use_integral は true でも false でも良いが、計算量は補間を使う false の方が少ないのでuse_integral は false にしておく
            // (あと use_integral を trueにするとアーティファクトがでる?(要確認))
            advect_velocity_flux_advection_3D_std_thread(
                0,
                all_grid,
                time_step_length,
                use_MacCormack_scheme,
                use_clamping_in_MacCormack_scheme,
                num_gauss_quad_boundary,
                num_gauss_quad_bulk,
    //            enable_eliminate_zero_velocity_false_diffusion,
                split_method_of_error_advection,
                integral_method,
                num_gauss_quadrature_points_on_one_triangle,
                zero_velocity_correction_values_velocity_x,
                interpolation_method_of_error_advection,
                interpolation_method_in_calc_psi,
                /*use_zero_velocity_for_backtrace*/false,
                /*use_integral=*/false,
                num_gauss_quadrature_point_for_integrate_density,
                /*enable_cell_volume_correction*/false,
                minimal_area
            );
            // 補正量を結果に足す
            for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
                for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                    for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
                        after_advect_values[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        += zero_velocity_correction_values_velocity_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                    }
                }
            }
        }
    }

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
    ) {
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
    //        interpolation_method_in_calc_psi = "const";
        }
        else if(interpolation_method=="WENO6-optimized"){
            interpolation_method_in_calc_psi = "WENO6-optimized";
    //        interpolation_method_in_calc_psi = "const";
        }
        //// y成分の移流における切断方法の指定
        std::string split_method_velocity_y;
        if(interpolation_method == "linear" || interpolation_method == "WENO6" || interpolation_method == "WENO6-optimized"){
             split_method_velocity_y = "xz-CellFaceAligned";
        }
        else if(interpolation_method == "1Dy_linear" || interpolation_method == "1Dy_WENO6"){
            split_method_velocity_y = "xz-AxisAligned";
        }
        //y成分の移流
        advect_velocity_flux_advection_3D_std_thread(
            1,
            all_grid,
            time_step_length,
            use_MacCormack_scheme,
            use_clamping_in_MacCormack_scheme,
            num_gauss_quad_boundary,
            num_gauss_quad_bulk,
    //		const bool enable_eliminate_zero_velocity_false_diffusion,
            split_method_velocity_y,
            integral_method,
            num_gauss_quadrature_points_on_one_triangle,
            after_advect_values,
    //			all_grid.velocity_in_voxel_face_x,
            interpolation_method,
            interpolation_method_in_calc_psi,
            /*use_zero_velocity_for_backtrace*/false,
            use_integral,
            num_gauss_quadrature_point_for_integrate_density,
            enable_cell_volume_correction,
            minimal_area
        );

    /*
        //y成分の移流
        advect_velocity_semi_lagrangian_y_3D(
            all_grid,
            time_step_length,
            interpolation_method,
            velocity_after_advect_y
        );
    */
        ////////////////////////////////////////////////////////////
        // 0速度場の拡散補正
        ////////////////////////////////////////////////////////////
        if (enable_eliminate_zero_velocity_false_diffusion) {
            /////////////////////////
            // 速度場のy成分の補正量
            /////////////////////////
            std::vector<MY_FLOAT_TYPE> zero_velocity_correction_values_velocity_y(all_grid.Grid_num_x * (all_grid.Grid_num_y + 1) * all_grid.Grid_num_z);
            ////// "バックトレース前の面で定義されたpsiの離散値"を足す
            for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
                for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
                    for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
                        zero_velocity_correction_values_velocity_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                            = before_advect_values[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
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
            advect_velocity_flux_advection_3D_std_thread(
                1,
                all_grid,
                time_step_length,
                use_MacCormack_scheme,
                use_clamping_in_MacCormack_scheme,
                num_gauss_quad_boundary,
                num_gauss_quad_bulk,
    //            enable_eliminate_zero_velocity_false_diffusion,
                split_method_velocity_y,
                integral_method,
                num_gauss_quadrature_points_on_one_triangle,
                zero_velocity_correction_values_velocity_y,
                interpolation_method,
                interpolation_method_in_calc_psi,
                /*use_zero_velocity_for_backtrace*/true,
                /*use_integral*/false,
                num_gauss_quadrature_point_for_integrate_density,
                /*enable_cell_volume_correction*/false,
                minimal_area
            );

            // 速度場の離散値との差を求める. これが移流前のエラー
            for(int ix = 0; ix < all_grid.Grid_num_x; ++ix){
                for(int iy = 0; iy < all_grid.Grid_num_y + 1; ++iy){
                    for(int iz = 0; iz < all_grid.Grid_num_z; ++iz){
                        zero_velocity_correction_values_velocity_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                            *= -1.0;
                        zero_velocity_correction_values_velocity_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                            +=before_advect_values[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                    }
                }
            }
            ////-------------------------------------
            //// エラーの移流
            ////-------------------------------------
            // エラーの移流に piecewise-const(一次元補間) を使う
            std::string interpolation_method_of_error_advection;
            if(interpolation_method == "linear" || interpolation_method == "1Dy_linear"){
                interpolation_method_of_error_advection = "1Dy_linear";
            }
            else if(interpolation_method == "WENO6"||interpolation_method == "WENO6-optimized"||interpolation_method == "1Dy_WENO6"){
    //                    interpolation_method_of_error_advection = "1Dy_WENO6";
                interpolation_method_of_error_advection = "1Dy_linear";
            }
            std::string split_method_of_error_advection = "xz-AxisAligned";
    //        const std::string split_method_of_error_advection = "xz-CellFaceAligned";
            std::string interpolation_method_in_calc_psi = "const";

            // flux advection をエラーの移流に使う場合
            // エラーの移流には1次元補間を使うが、1次元補間ではpsiの計算に補間を使っても積分で厳密に計算しても結果は同じである。
            // なので use_integral は true でも false でも良いが、計算量は補間を使う false の方が少ないのでuse_integral は false にしておく
            // (あと use_integral を trueにするとアーティファクトがでる?(要確認))
            advect_velocity_flux_advection_3D_std_thread(
                1,
                all_grid,
                time_step_length,
                use_MacCormack_scheme,
                use_clamping_in_MacCormack_scheme,
                num_gauss_quad_boundary,
                num_gauss_quad_bulk,
    //            enable_eliminate_zero_velocity_false_diffusion,
                split_method_of_error_advection,
                integral_method,
                num_gauss_quadrature_points_on_one_triangle,
                zero_velocity_correction_values_velocity_y,
                interpolation_method_of_error_advection,
                interpolation_method_in_calc_psi,
                /*use_zero_velocity_for_backtrace*/false,
                /*use_integral=*/false,
                num_gauss_quadrature_point_for_integrate_density,
                /*enable_cell_volume_correction*/false,
                minimal_area
            );
            // 補正量を結果に足す
            for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
                for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
                    for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
                        after_advect_values[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        += zero_velocity_correction_values_velocity_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                    }
                }
            }
        }
    }

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
    ) {
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
        //// z成分の移流における切断方法の指定
        std::string split_method_velocity_z;
        if(interpolation_method == "linear" || interpolation_method == "WENO6" || interpolation_method == "WENO6-optimized"){
            split_method_velocity_z = "x-CellFaceAligned_z-AxisAligned";
        }
        else if(interpolation_method == "1Dy_linear" || interpolation_method == "1Dy_WENO6"){
            split_method_velocity_z = "x-AxisAligned_z-CellFaceAligned";
        }
        //z成分の移流
        advect_velocity_flux_advection_3D_std_thread(
            2,
            all_grid,
            time_step_length,
            use_MacCormack_scheme,
            use_clamping_in_MacCormack_scheme,
            num_gauss_quad_boundary,
            num_gauss_quad_bulk,
            //		const bool enable_eliminate_zero_velocity_false_diffusion,
            split_method_velocity_z,
            integral_method,
            num_gauss_quadrature_points_on_one_triangle,
            after_advect_values,
            interpolation_method,
            interpolation_method_in_calc_psi,
            /*use_zero_velocity_for_backtrace*/false,
            use_integral,
            num_gauss_quadrature_point_for_integrate_density,
            enable_cell_volume_correction,
            minimal_area
        );
        ////////////////////////////////////////////////////////////
        // 0速度場の拡散補正
        ////////////////////////////////////////////////////////////
        if (enable_eliminate_zero_velocity_false_diffusion) {
            /////////////////////////
            // 速度場のz成分の補正量
            /////////////////////////
            std::vector<MY_FLOAT_TYPE> zero_velocity_correction_values_velocity_z(all_grid.Grid_num_x * all_grid.Grid_num_y * (all_grid.Grid_num_z + 1));
            ////// "バックトレース前の面で定義されたpsiの離散値"を足す
            for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
                for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                    for (int iz = 0; iz < all_grid.Grid_num_z + 1; iz++) {
                        zero_velocity_correction_values_velocity_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                            = before_advect_values[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
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
            advect_velocity_flux_advection_3D_std_thread(
                2,
                all_grid,
                time_step_length,
                use_MacCormack_scheme,
                use_clamping_in_MacCormack_scheme,
                num_gauss_quad_boundary,
                num_gauss_quad_bulk,
    //            enable_eliminate_zero_velocity_false_diffusion,
                split_method_velocity_z,
                integral_method,
                num_gauss_quadrature_points_on_one_triangle,
                zero_velocity_correction_values_velocity_z,
                interpolation_method,
                interpolation_method_in_calc_psi,
                /*use_zero_velocity_for_backtrace*/true,
                /*use_integral*/false,
                num_gauss_quadrature_point_for_integrate_density,
                /*enable_cell_volume_correction*/false,
                minimal_area
            );

            // 速度場の離散値との差を求める. これが移流前のエラー
            for(int ix = 0; ix < all_grid.Grid_num_x; ++ix){
                for(int iy = 0; iy < all_grid.Grid_num_y; ++iy){
                    for(int iz = 0; iz < all_grid.Grid_num_z + 1; ++iz){
                        zero_velocity_correction_values_velocity_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                            *= -1.0;
                        zero_velocity_correction_values_velocity_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                            +=before_advect_values[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                    }
                }
            }
            ////-------------------------------------
            //// エラーの移流
            ////-------------------------------------
            // エラーの移流に piecewise-const(一次元補間) を使う
            std::string interpolation_method_of_error_advection;
            if(interpolation_method == "linear" || interpolation_method == "1Dy_linear"){
                interpolation_method_of_error_advection = "1Dy_linear";
            }
            else if(interpolation_method == "WENO6"||interpolation_method == "WENO6-optimized"||interpolation_method == "1Dy_WENO6"){
    //                    interpolation_method_of_error_advection = "1Dy_WENO6";
                interpolation_method_of_error_advection = "1Dy_linear";
            }
    //        std::string split_method_of_error_advection = "xz-AxisAligned";
            std::string split_method_of_error_advection = "x-AxisAligned_z-CellFaceAligned";
            std::string interpolation_method_in_calc_psi = "const";

            // flux advection をエラーの移流に使う場合
            // エラーの移流には1次元補間を使うが、1次元補間ではpsiの計算に補間を使っても積分で厳密に計算しても結果は同じである。
            // なので use_integral は true でも false でも良いが、計算量は補間を使う false の方が少ないのでuse_integral は false にしておく
            // (あと use_integral を trueにするとアーティファクトがでる?(要確認))
            advect_velocity_flux_advection_3D_std_thread(
                2,
                all_grid,
                time_step_length,
                use_MacCormack_scheme,
                use_clamping_in_MacCormack_scheme,
                num_gauss_quad_boundary,
                num_gauss_quad_bulk,
    //            enable_eliminate_zero_velocity_false_diffusion,
                split_method_of_error_advection,
                integral_method,
                num_gauss_quadrature_points_on_one_triangle,
                zero_velocity_correction_values_velocity_z,
                interpolation_method_of_error_advection,
                interpolation_method_in_calc_psi,
                /*use_zero_velocity_for_backtrace*/false,
                /*use_integral=*/false,
                num_gauss_quadrature_point_for_integrate_density,
                /*enable_cell_volume_correction*/false,
                minimal_area
            );
            // 補正量を結果に足す
            for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
                for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                    for (int iz = 0; iz < all_grid.Grid_num_z + 1; iz++) {
                        after_advect_values[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        += zero_velocity_correction_values_velocity_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                    }
                }
            }
        }
    }

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
    ) {
        //////////////////////////////////////////////////
        // x 成分の移流
        //////////////////////////////////////////////////
        // 移流前の速度場
        std::vector<MY_FLOAT_TYPE> before_advect_velocity_x((all_grid.Grid_num_x + 1) * all_grid.Grid_num_y * all_grid.Grid_num_z);
        before_advect_velocity_x = all_grid.velocity_in_voxel_face_x;
        // 移流後の速度場
        std::vector<MY_FLOAT_TYPE> velocity_after_advect_x((all_grid.Grid_num_x + 1) * all_grid.Grid_num_y * all_grid.Grid_num_z);
        velocity_after_advect_x = all_grid.velocity_in_voxel_face_x;
        //移流の処理
        advect_velocity_x_flux_advection_3D_std_thread(
            before_advect_velocity_x,
            velocity_after_advect_x,
            all_grid,
            time_step_length,
            use_MacCormack_scheme,
            use_clamping_in_MacCormack_scheme,
            num_gauss_quad_boundary,
            num_gauss_quad_bulk,
            enable_eliminate_zero_velocity_false_diffusion,
    //    const std::string split_method,
            integral_method,
            num_gauss_quadrature_points_on_one_triangle,
            interpolation_method,
    //    const std::string interpolation_method_in_calc_psi,
            use_integral,
            num_gauss_quadrature_point_for_integrate_density,
            enable_cell_volume_correction,
            minimal_area
        );
        //////////////////////////////////////////////////
        // y 成分の移流
        //////////////////////////////////////////////////
        // 移流前の速度場
        std::vector<MY_FLOAT_TYPE> before_advect_velocity_y(all_grid.Grid_num_x * (all_grid.Grid_num_y + 1) * all_grid.Grid_num_z);
        before_advect_velocity_y = all_grid.velocity_in_voxel_face_y;
        // 移流後の速度場
        std::vector<MY_FLOAT_TYPE> velocity_after_advect_y(all_grid.Grid_num_x * (all_grid.Grid_num_y + 1) * all_grid.Grid_num_z);
        velocity_after_advect_y = all_grid.velocity_in_voxel_face_y;
        //移流の処理
        advect_velocity_y_flux_advection_3D_std_thread(
            before_advect_velocity_y,
            velocity_after_advect_y,
            all_grid,
            time_step_length,
            use_MacCormack_scheme,
            use_clamping_in_MacCormack_scheme,
            num_gauss_quad_boundary,
            num_gauss_quad_bulk,
            enable_eliminate_zero_velocity_false_diffusion,
    //    const std::string split_method,
            integral_method,
            num_gauss_quadrature_points_on_one_triangle,
            interpolation_method,
    //    const std::string interpolation_method_in_calc_psi,
            use_integral,
            num_gauss_quadrature_point_for_integrate_density,
            enable_cell_volume_correction,
            minimal_area
        );
        //////////////////////////////////////////////////
        // z 成分の移流
        //////////////////////////////////////////////////
        // 移流前の速度場
        std::vector<MY_FLOAT_TYPE> before_advect_velocity_z(all_grid.Grid_num_x * all_grid.Grid_num_y * (all_grid.Grid_num_z + 1));
        before_advect_velocity_z = all_grid.velocity_in_voxel_face_z;
        // 移流後の速度場
        std::vector<MY_FLOAT_TYPE> velocity_after_advect_z(all_grid.Grid_num_x * all_grid.Grid_num_y * (all_grid.Grid_num_z + 1));
        velocity_after_advect_z = all_grid.velocity_in_voxel_face_z;
        //移流の処理
        advect_velocity_z_flux_advection_3D_std_thread(
            before_advect_velocity_z,
            velocity_after_advect_z,
            all_grid,
            time_step_length,
            use_MacCormack_scheme,
            use_clamping_in_MacCormack_scheme,
            num_gauss_quad_boundary,
            num_gauss_quad_bulk,
            enable_eliminate_zero_velocity_false_diffusion,
    //    const std::string split_method,
            integral_method,
            num_gauss_quadrature_points_on_one_triangle,
            interpolation_method,
    //    const std::string interpolation_method_in_calc_psi,
            use_integral,
            num_gauss_quadrature_point_for_integrate_density,
            enable_cell_volume_correction,
            minimal_area
        );

        //計算結果をコピー
        for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
                    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        = velocity_after_advect_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                }
            }
        }
        for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
                for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
                    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        = velocity_after_advect_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                }
            }
        }
        for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                for (int iz = 0; iz < all_grid.Grid_num_z + 1; iz++) {
                    all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        = velocity_after_advect_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                }
            }
        }
    }
} // namespace smoke_simulation
