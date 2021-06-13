#include "calc_density_in_cell_from_quadrature_point_list_for_integrate_density.h"
#include "linear_interpolation_1d.h"
#include "linear_interpolation_3d.h"
#include "WENO6_interpolation_1d.h"
#include "WENO6_interpolation_3d.h"

#include "sginterp3.h"
#include "utils.h"

namespace smoke_simulation{
    // i_thread 番目のスレッドで計算する求積点のリスト quadrature_point_list_of_i_thread から 移流後の質量密度場substance_density_after_advectへの寄与を計算する
    // psiの計算を補間でなく積分によって行うパターン
	// quadrature_point_list_interpolate_density には積分するdensityからの寄与
	// quadrature_point_list_interpolate_psi にはpsiの値からの寄与
    void calc_density_in_cell_from_quadrature_point_list_for_integrate_density(
        const quadrature_point_vector &quadrature_point_list_for_interpolate_density,
        const quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &substance_density_after_advect,
        const int i_thread,
        const std::string interpolation_method,
        const std::vector<MY_FLOAT_TYPE> &advected_values
    ){
        //最適化されたWENO6 を使う場合
        if(interpolation_method == "WENO6-optimized"){
            //最適化オプションを設定
            int backend (0);
    //        backend |= sginterp3::backend::multithread;
    		backend |= sginterp3::backend::simd;
    		backend |= sginterp3::backend::presort;
            //補間の手法を指定
    		sginterp3::filter scheme;
    		scheme = sginterp3::filter::WENO6;
            ////求積点の位置のリストをsginterp3::point3 のstd::vector に直す
            // density の求積点
            const int num_quadrature_points_interpolate_density = quadrature_point_list_for_interpolate_density._quadtarure_position[i_thread].size();
            std::vector<sginterp3::point3<MY_FLOAT_TYPE>> quadrature_point_position_list_interpolate_density(num_quadrature_points_interpolate_density);
            // psi の求積点
            const int num_quadrature_points_interpolate_psi = quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread].size();
            std::vector<sginterp3::point3<MY_FLOAT_TYPE>> quadrature_point_position_list_interpolate_psi(num_quadrature_points_interpolate_psi);
            //---------------------------------------------------
            //補間の計算の部分
            //---------------------------------------------------

            ////////////////////////////////////////
            // densityからの寄与
            ////////////////////////////////////////
            //// densityを補間
            // density はセルの中心で定義されているので補完する位置をセルの長さの1/2だけずらす
            for(int i_point = 0; i_point< num_quadrature_points_interpolate_density; ++i_point){
                quadrature_point_position_list_interpolate_density[i_point].x = quadrature_point_list_for_interpolate_density._quadtarure_position[i_thread][i_point][2] - 0.5 * all_grid._cell_length;
//                quadrature_point_position_list_interpolate_density[i_point].x = quadrature_point_list_for_interpolate_density._quadtarure_position[i_thread][i_point][2];
                quadrature_point_position_list_interpolate_density[i_point].y = quadrature_point_list_for_interpolate_density._quadtarure_position[i_thread][i_point][1] - 0.5 * all_grid._cell_length;
//                quadrature_point_position_list_interpolate_density[i_point].y = quadrature_point_list_for_interpolate_density._quadtarure_position[i_thread][i_point][1];
                quadrature_point_position_list_interpolate_density[i_point].z = quadrature_point_list_for_interpolate_density._quadtarure_position[i_thread][i_point][0] - 0.5 * all_grid._cell_length;
//                quadrature_point_position_list_interpolate_density[i_point].z = quadrature_point_list_for_interpolate_density._quadtarure_position[i_thread][i_point][0];
            }
            auto interpolation_result_list_of_density = sginterp3::interpolate(
        		quadrature_point_position_list_interpolate_density, all_grid.Grid_num_z, all_grid.Grid_num_y, all_grid.Grid_num_x, all_grid._cell_length, all_grid.substance_density.data(), scheme, backend
            );

            //これ以降OK
            //補間で求めたpsi から質量密度を計算
            for(int i_point = 0; i_point< num_quadrature_points_interpolate_density; ++i_point){
                substance_density_after_advect[
                    get_voxel_center_index_3D(
                        quadrature_point_list_for_interpolate_density._included_cell_index[i_thread][i_point][0],
                        quadrature_point_list_for_interpolate_density._included_cell_index[i_thread][i_point][1],
                        quadrature_point_list_for_interpolate_density._included_cell_index[i_thread][i_point][2],
                        all_grid.Grid_num_x,
                        all_grid.Grid_num_y,
                        all_grid.Grid_num_z)]
                    +=interpolation_result_list_of_density[i_point]
                    * quadrature_point_list_for_interpolate_density._weighted_area_normal[i_thread][i_point][1]
                    / all_grid._cell_volume;
            }
            ////////////////////////////////////////
            // psiからの寄与
            ////////////////////////////////////////
            ////psiのy成分を補間
            // psi はセルの面で定義されているので補完する位置をセルの長さの1/2だけずらす
            for(int i_point = 0; i_point< num_quadrature_points_interpolate_psi; ++i_point){
                quadrature_point_position_list_interpolate_psi[i_point].x = quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread][i_point][2] - 0.5 * all_grid._cell_length;
    //            quadrature_point_position_list[i_point].x = quadrature_point_list._quadtarure_position[i_thread][i_point][2];
    //            quadrature_point_position_list[i_point].y = quadrature_point_list._quadtarure_position[i_thread][i_point][1] - 0.5 * all_grid._cell_length;
                quadrature_point_position_list_interpolate_psi[i_point].y = quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread][i_point][1];
                quadrature_point_position_list_interpolate_psi[i_point].z = quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread][i_point][0] - 0.5 * all_grid._cell_length;
    //            quadrature_point_position_list[i_point].z = quadrature_point_list._quadtarure_position[i_thread][i_point][0];
            }
            auto interpolation_result_list_of_psi_density_y = sginterp3::interpolate(
        		quadrature_point_position_list_interpolate_psi, all_grid.Grid_num_z, all_grid.Grid_num_y + 1, all_grid.Grid_num_x, all_grid._cell_length, all_grid.psi_substance_density_cell_face_y.data(), scheme, backend
            );
            //これ以降OK
            //補間で求めたpsi から質量密度を計算
            for(int i_point = 0; i_point< num_quadrature_points_interpolate_psi; ++i_point){
                substance_density_after_advect[
                    get_voxel_center_index_3D(
                        quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][0],
                        quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][1],
                        quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][2],
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
                        0.0,
                        interpolation_result_list_of_psi_density_y[i_point],
                        0.0
                    ).dot(quadrature_point_list_for_interpolate_psi._weighted_area_normal[i_thread][i_point])
                    / all_grid._cell_volume;
            }
        }
        else{
            ////////////////////////////////////////
            // densityからの寄与
            ////////////////////////////////////////
            const int num_quadrature_points_for_interpolate_density = quadrature_point_list_for_interpolate_density._quadtarure_position[i_thread].size();
            //求積点を走査するループ
            for(size_t i_point = 0; i_point < num_quadrature_points_for_interpolate_density; ++i_point){
                ////補間に依ってdensityを求める
//                MY_FLOAT_TYPE density_on_barycenter_of_i_tri
//                    = linear_interpolation_3D_cell_center_values(
//                        quadrature_point_list_for_interpolate_density._quadtarure_position[i_thread][i_point],
//                        all_grid.substance_density,
//                        all_grid
//                );

                MY_FLOAT_TYPE density_on_barycenter_of_i_tri;
                if(interpolation_method == "linear"){
                    density_on_barycenter_of_i_tri
                        = linear_interpolation_3D_cell_center_values(
                            quadrature_point_list_for_interpolate_density._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
                else if(interpolation_method == "1Dy_linear"){
                    density_on_barycenter_of_i_tri
                        = linear_interpolation_1Dy_cell_center_values(
                            quadrature_point_list_for_interpolate_density._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
                else if(interpolation_method == "WENO6"){
                    density_on_barycenter_of_i_tri
                        = WENO6_interpolation_3D_cell_center_values(
                            quadrature_point_list_for_interpolate_density._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
                else if(interpolation_method == "1Dy_WENO6"){
                    density_on_barycenter_of_i_tri
                        = WENO6_interpolation_1Dy_cell_center_values(
                            quadrature_point_list_for_interpolate_density._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }

                substance_density_after_advect[
                    get_voxel_center_index_3D(
                        quadrature_point_list_for_interpolate_density._included_cell_index[i_thread][i_point][0],
                        quadrature_point_list_for_interpolate_density._included_cell_index[i_thread][i_point][1],
                        quadrature_point_list_for_interpolate_density._included_cell_index[i_thread][i_point][2],
                        all_grid.Grid_num_x,
                        all_grid.Grid_num_y,
                        all_grid.Grid_num_z)
                ]   +=density_on_barycenter_of_i_tri
                    * quadrature_point_list_for_interpolate_density._weighted_area_normal[i_thread][i_point][1]
                    / all_grid._cell_volume;
            }
            ////////////////////////////////////////
            // psiからの寄与
            ////////////////////////////////////////
            const int num_quadrature_points_for_interpolate_psi=quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread].size();
            //求積点を走査するループ
            for(size_t i_point = 0; i_point < num_quadrature_points_for_interpolate_psi; ++i_point){
                ////補間に依ってpsiを求める
                VEC3_TYPE psi_density_on_barycenter_of_i_tri
                    = all_grid.calc_interpolated_psi_density(quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread][i_point], interpolation_method);

                substance_density_after_advect[
                    get_voxel_center_index_3D(
                        quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][0],
                        quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][1],
                        quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][2],
                        all_grid.Grid_num_x,
                        all_grid.Grid_num_y,
                        all_grid.Grid_num_z)]
                    +=psi_density_on_barycenter_of_i_tri.dot(quadrature_point_list_for_interpolate_psi._weighted_area_normal[i_thread][i_point])
                    / all_grid._cell_volume;
            }
        }
    }

}
