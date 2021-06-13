#include "calc_density_in_cell_from_quadrature_point_list.h"

#include "sginterp3.h"
#include "utils.h"

namespace smoke_simulation{
    void calc_density_in_cell_from_quadrature_point_list(
        const quadrature_point_vector &quadrature_point_list,
//      size_t i_start,
//      size_t i_end,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &substance_density_after_advect,
        const int i_thread,
        const std::string interpolation_method
    ){
        //最適化されたWENO6 を使う場合
        if(interpolation_method == "WENO6-optimized"){
            //最適化オプションを設定
            int backend (0);
//          backend |= sginterp3::backend::multithread;
    		backend |= sginterp3::backend::simd;
    		backend |= sginterp3::backend::presort;
            //補間の手法を指定
    		sginterp3::filter scheme;
    		scheme = sginterp3::filter::WENO6;

            //求積点の位置のリストをsginterp3::point3 のstd::vector に直す
            const int num_quadrature_points = quadrature_point_list._quadtarure_position[i_thread].size();
            std::vector<sginterp3::point3<MY_FLOAT_TYPE>> quadrature_point_position_list(num_quadrature_points);

            //---------------------------------------------------
            //補間の計算の部分
            //---------------------------------------------------
            //// psiのx成分を補間
            // psi はセルの面で定義されているので補完する位置をセルの長さの1/2だけずらす
            for(int i_point = 0; i_point< num_quadrature_points; ++i_point){
                quadrature_point_position_list[i_point].x = quadrature_point_list._quadtarure_position[i_thread][i_point][2] - 0.5 * all_grid._cell_length;
//              quadrature_point_position_list[i_point].x = quadrature_point_list._quadtarure_position[i_thread][i_point][2];
                quadrature_point_position_list[i_point].y = quadrature_point_list._quadtarure_position[i_thread][i_point][1] - 0.5 * all_grid._cell_length;
//              quadrature_point_position_list[i_point].y = quadrature_point_list._quadtarure_position[i_thread][i_point][1];
//              quadrature_point_position_list[i_point].z = quadrature_point_list._quadtarure_position[i_thread][i_point][0] - 0.5 * all_grid._cell_length;
                quadrature_point_position_list[i_point].z = quadrature_point_list._quadtarure_position[i_thread][i_point][0];
            }
            auto interpolation_result_list_of_psi_density_x = sginterp3::interpolate(
        		quadrature_point_position_list, all_grid.Grid_num_z, all_grid.Grid_num_y, all_grid.Grid_num_x + 1, all_grid._cell_length, all_grid.psi_substance_density_cell_face_x.data(), scheme, backend
            );
            ////psiのy成分を補間
            // psi はセルの面で定義されているので補完する位置をセルの長さの1/2だけずらす
            for(int i_point = 0; i_point< num_quadrature_points; ++i_point){
                quadrature_point_position_list[i_point].x = quadrature_point_list._quadtarure_position[i_thread][i_point][2] - 0.5 * all_grid._cell_length;
    //            quadrature_point_position_list[i_point].x = quadrature_point_list._quadtarure_position[i_thread][i_point][2];
    //            quadrature_point_position_list[i_point].y = quadrature_point_list._quadtarure_position[i_thread][i_point][1] - 0.5 * all_grid._cell_length;
                quadrature_point_position_list[i_point].y = quadrature_point_list._quadtarure_position[i_thread][i_point][1];
                quadrature_point_position_list[i_point].z = quadrature_point_list._quadtarure_position[i_thread][i_point][0] - 0.5 * all_grid._cell_length;
    //            quadrature_point_position_list[i_point].z = quadrature_point_list._quadtarure_position[i_thread][i_point][0];
            }
            auto interpolation_result_list_of_psi_density_y = sginterp3::interpolate(
        		quadrature_point_position_list, all_grid.Grid_num_z, all_grid.Grid_num_y + 1, all_grid.Grid_num_x, all_grid._cell_length, all_grid.psi_substance_density_cell_face_y.data(), scheme, backend
            );
            ////psiのz成分を補間
            for(int i_point = 0; i_point< num_quadrature_points; ++i_point){
    //            quadrature_point_position_list[i_point].x = quadrature_point_list._quadtarure_position[i_thread][i_point][2] - 0.5 * all_grid._cell_length;
                quadrature_point_position_list[i_point].x = quadrature_point_list._quadtarure_position[i_thread][i_point][2];
                quadrature_point_position_list[i_point].y = quadrature_point_list._quadtarure_position[i_thread][i_point][1] - 0.5 * all_grid._cell_length;
    //            quadrature_point_position_list[i_point].y = quadrature_point_list._quadtarure_position[i_thread][i_point][1];
                quadrature_point_position_list[i_point].z = quadrature_point_list._quadtarure_position[i_thread][i_point][0] - 0.5 * all_grid._cell_length;
    //            quadrature_point_position_list[i_point].z = quadrature_point_list._quadtarure_position[i_thread][i_point][0];
            }
            auto interpolation_result_list_of_psi_density_z = sginterp3::interpolate(
        		quadrature_point_position_list, all_grid.Grid_num_z + 1, all_grid.Grid_num_y, all_grid.Grid_num_x, all_grid._cell_length, all_grid.psi_substance_density_cell_face_z.data(), scheme, backend
            );
            //これ以降OK
            //補間で求めたpsi から質量密度を計算
            for(int i_point = 0; i_point< num_quadrature_points; ++i_point){
                substance_density_after_advect[
                    get_voxel_center_index_3D(
                        quadrature_point_list._included_cell_index[i_thread][i_point][0],
                        quadrature_point_list._included_cell_index[i_thread][i_point][1],
                        quadrature_point_list._included_cell_index[i_thread][i_point][2],
                        all_grid.Grid_num_x,
                        all_grid.Grid_num_y,
                        all_grid.Grid_num_z)]
    //                +=quadrature_point_list._quadtarure_weight[i_thread][i_point]
    //                * VEC3_TYPE(
    //                    interpolation_result_list_of_psi_density_x[i_point],
    //                    interpolation_result_list_of_psi_density_y[i_point],
    //                    interpolation_result_list_of_psi_density_z[i_point]
    //                ).dot(quadrature_point_list._normal[i_thread][i_point])
    //                / all_grid._cell_volume;
                    +=VEC3_TYPE(
                        interpolation_result_list_of_psi_density_x[i_point],
                        interpolation_result_list_of_psi_density_y[i_point],
                        interpolation_result_list_of_psi_density_z[i_point]
                    ).dot(quadrature_point_list._weighted_area_normal[i_thread][i_point])
                    / all_grid._cell_volume;
            }
        }
        else{
            const int num_quadrature_points=quadrature_point_list._quadtarure_position[i_thread].size();
            //求積点を走査するループ
            for(size_t i_point = 0; i_point < num_quadrature_points; ++i_point){
                ////補間に依ってpsiを求める
                // x, y 方向の補間を使う場合( 8/22 のノートの方法 )
                std::string interpolation_method_for_calc_interpolated_psi;
//                if(interpolation_method == "linear"){
//                    interpolation_method_for_calc_interpolated_psi = "linear_integral";
//                }
    //            if(interpolation_method == "1Dy_linear"){
    //                interpolation_method_for_calc_interpolated_psi = "1Dy_linear_integral";
    //            }
//                else{
                    interpolation_method_for_calc_interpolated_psi = interpolation_method;
//                }
                VEC3_TYPE psi_density_on_barycenter_of_i_tri
                    = all_grid.calc_interpolated_psi_density(quadrature_point_list._quadtarure_position[i_thread][i_point], interpolation_method_for_calc_interpolated_psi);
                // 1次元方向のみの補間を使う場合( 8/30 のノートの方法 )
                //VEC3_TYPE psi_density_on_barycenter_of_i_tri = all_grid.calc_psi_substance_density_by_1D_interpolation_cell_face(quadrature_position_on_cell_face);
                substance_density_after_advect[
                    get_voxel_center_index_3D(
                        quadrature_point_list._included_cell_index[i_thread][i_point][0],
                        quadrature_point_list._included_cell_index[i_thread][i_point][1],
                        quadrature_point_list._included_cell_index[i_thread][i_point][2],
                        all_grid.Grid_num_x,
                        all_grid.Grid_num_y,
                        all_grid.Grid_num_z)]
    //                +=quadrature_point_list._quadtarure_weight[i_thread][i_point]
    //                * psi_density_on_barycenter_of_i_tri.dot(quadrature_point_list._normal[i_thread][i_point])
    //                / all_grid._cell_volume;
                    +=psi_density_on_barycenter_of_i_tri.dot(quadrature_point_list._weighted_area_normal[i_thread][i_point])
                    / all_grid._cell_volume;
            }
        }
    }
    // i_thread 番目のスレッドで計算する求積点のリスト quadrature_point_list_of_i_thread から 移流後のセル体積 cell_volume_after_advect への寄与を計算する
    void calc_cell_volumes_from_quadrature_point_list(
        const quadrature_point_vector &quadrature_point_list,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &cell_volume_after_advect,
        const int i_thread
    ){
        const int num_quadrature_points=quadrature_point_list._quadtarure_position[i_thread].size();
        //求積点を走査するループ
        for(size_t i_point = 0; i_point < num_quadrature_points; ++i_point){
            // セル体積への寄与を計算
            cell_volume_after_advect[
                get_voxel_center_index_3D(
                    quadrature_point_list._included_cell_index[i_thread][i_point][0],
                    quadrature_point_list._included_cell_index[i_thread][i_point][1],
                    quadrature_point_list._included_cell_index[i_thread][i_point][2],
                    all_grid.Grid_num_x,
                    all_grid.Grid_num_y,
                    all_grid.Grid_num_z)]
                +=quadrature_point_list._quadtarure_position[i_thread][i_point][1] * quadrature_point_list._weighted_area_normal[i_thread][i_point][1];

        }
    }

}
