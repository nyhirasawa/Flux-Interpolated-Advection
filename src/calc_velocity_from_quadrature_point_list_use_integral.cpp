#include "calc_velocity_from_quadrature_point_list_use_integral.h"

#include "sginterp3.h"
#include "utils.h"
#include "linear_interpolation_1d.h"
#include "linear_interpolation_3d.h"
#include "WENO6_interpolation_1d.h"
#include "WENO6_interpolation_3d.h"

namespace smoke_simulation{
    // i_thread 番目のスレッドで計算する求積点のリスト quadrature_point_list_of_i_thread から 移流後の質量密度場 velocity_x_after_advectへの寄与を計算する
    void calc_velocity_from_quadrature_point_list_use_integral(
        const int dim,
        const quadrature_point_vector &quadrature_point_list_for_interpolate_velocity,
        const quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &velocity_after_advect,
        const int i_thread,
        const std::string interpolation_method,
        const std::vector<MY_FLOAT_TYPE> &advected_values
    ){
        const int dx = (dim == 0);
        const int dy = (dim == 1);
        const int dz = (dim == 2);
        ////////////////////////////////////////
        // densityからの寄与
        ////////////////////////////////////////
        const int num_quadrature_points_for_interpolate_velocity = quadrature_point_list_for_interpolate_velocity._quadtarure_position[i_thread].size();
        //求積点を走査するループ
        for(size_t i_point = 0; i_point < num_quadrature_points_for_interpolate_velocity; ++i_point){
            ////補間に依ってdensityを求める
            MY_FLOAT_TYPE velocity_on_barycenter_of_i_tri;
            if(interpolation_method == "linear"){
                if(dim == 0){
                    velocity_on_barycenter_of_i_tri
                        = linear_interpolation_3D_cell_face_x_values(
                            quadrature_point_list_for_interpolate_velocity._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
                else if(dim == 1){
                    velocity_on_barycenter_of_i_tri
                        = linear_interpolation_3D_cell_face_y_values(
                            quadrature_point_list_for_interpolate_velocity._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
                else if(dim == 2){
                    velocity_on_barycenter_of_i_tri
                        = linear_interpolation_3D_cell_face_z_values(
                            quadrature_point_list_for_interpolate_velocity._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
            }
            else if(interpolation_method == "1Dy_linear"){
                if(dim == 0){
                    velocity_on_barycenter_of_i_tri
                        = linear_interpolation_1Dy_cell_face_x_values(
                            quadrature_point_list_for_interpolate_velocity._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
                else if(dim == 1){
                    velocity_on_barycenter_of_i_tri
                        = linear_interpolation_1Dy_cell_face_y_values(
                            quadrature_point_list_for_interpolate_velocity._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
                else if(dim == 2){
                    velocity_on_barycenter_of_i_tri
                        = linear_interpolation_1Dy_cell_face_z_values(
                            quadrature_point_list_for_interpolate_velocity._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
            }
            else if(interpolation_method == "WENO6"){
                if(dim == 0){
                    velocity_on_barycenter_of_i_tri
                        = WENO6_interpolation_3D_cell_face_x_values(
                            quadrature_point_list_for_interpolate_velocity._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
                else if(dim == 1){
                    velocity_on_barycenter_of_i_tri
                        = WENO6_interpolation_3D_cell_face_y_values(
                            quadrature_point_list_for_interpolate_velocity._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
                else if(dim == 2){
                    velocity_on_barycenter_of_i_tri
                        = WENO6_interpolation_3D_cell_face_z_values(
                            quadrature_point_list_for_interpolate_velocity._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
            }
            else if(interpolation_method == "1Dy_WENO6"){
                if(dim == 0){
                    velocity_on_barycenter_of_i_tri
                        = WENO6_interpolation_1Dy_cell_face_x_values(
                            quadrature_point_list_for_interpolate_velocity._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
                else if(dim == 1){
                    velocity_on_barycenter_of_i_tri
                        = WENO6_interpolation_1Dy_cell_face_y_values(
                            quadrature_point_list_for_interpolate_velocity._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
                else if(dim == 2){
                    velocity_on_barycenter_of_i_tri
                        = WENO6_interpolation_1Dy_cell_face_z_values(
                            quadrature_point_list_for_interpolate_velocity._quadtarure_position[i_thread][i_point],
                            advected_values,
                            all_grid
                    );
                }
            }
            velocity_after_advect[
                get_voxel_center_index_3D(
                    quadrature_point_list_for_interpolate_velocity._included_cell_index[i_thread][i_point][0],
                    quadrature_point_list_for_interpolate_velocity._included_cell_index[i_thread][i_point][1],
                    quadrature_point_list_for_interpolate_velocity._included_cell_index[i_thread][i_point][2],
                    all_grid.Grid_num_x + dx,
                    all_grid.Grid_num_y + dy,
                    all_grid.Grid_num_z + dz)
            ]   +=velocity_on_barycenter_of_i_tri
                * quadrature_point_list_for_interpolate_velocity._weighted_area_normal[i_thread][i_point][1]
                / all_grid._cell_volume;
        }
        ////////////////////////////////////////
        // psiからの寄与
        ////////////////////////////////////////
        const int num_quadrature_points_for_interpolate_psi=quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread].size();
        //求積点を走査するループ
        for(size_t i_point = 0; i_point < num_quadrature_points_for_interpolate_psi; ++i_point){
            ////補間に依ってpsiを求める
            VEC3_TYPE psi_velocity_on_barycenter_of_i_tri;
            if(dim == 0){
                psi_velocity_on_barycenter_of_i_tri
                    = all_grid.calc_interpolated_psi_velocity_x(quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread][i_point], interpolation_method);
            }
            else if(dim == 1){
                psi_velocity_on_barycenter_of_i_tri
                    = all_grid.calc_interpolated_psi_velocity_y(quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread][i_point], interpolation_method);
            }
            else if(dim == 2){
                psi_velocity_on_barycenter_of_i_tri
                    = all_grid.calc_interpolated_psi_velocity_z(quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread][i_point], interpolation_method);
            }

            velocity_after_advect[
                get_voxel_center_index_3D(
                    quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][0],
                    quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][1],
                    quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][2],
                    all_grid.Grid_num_x + dx,
                    all_grid.Grid_num_y + dy,
                    all_grid.Grid_num_z + dz)]
                +=psi_velocity_on_barycenter_of_i_tri.dot(quadrature_point_list_for_interpolate_psi._weighted_area_normal[i_thread][i_point])
                / all_grid._cell_volume;
        }
    }

    void calc_velocity_x_from_quadrature_point_list_use_integral(
        const quadrature_point_vector &quadrature_point_list_for_interpolate_velocity_x,
        const quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &velocity_x_after_advect,
        const int i_thread,
        const std::string interpolation_method,
        const std::vector<MY_FLOAT_TYPE> &advected_values
    ){
        ////////////////////////////////////////
        // densityからの寄与
        ////////////////////////////////////////
        const int num_quadrature_points_for_interpolate_velocity_x = quadrature_point_list_for_interpolate_velocity_x._quadtarure_position[i_thread].size();
        //求積点を走査するループ
        for(size_t i_point = 0; i_point < num_quadrature_points_for_interpolate_velocity_x; ++i_point){
            ////補間に依ってdensityを求める
            MY_FLOAT_TYPE velocity_x_on_barycenter_of_i_tri;
            if(interpolation_method == "linear"){
                velocity_x_on_barycenter_of_i_tri
                    = linear_interpolation_3D_cell_face_x_values(
                        quadrature_point_list_for_interpolate_velocity_x._quadtarure_position[i_thread][i_point],
                        advected_values,
                        all_grid
                );
            }
            else if(interpolation_method == "1Dy_linear"){
                velocity_x_on_barycenter_of_i_tri
                    = linear_interpolation_1Dy_cell_face_x_values(
                        quadrature_point_list_for_interpolate_velocity_x._quadtarure_position[i_thread][i_point],
                        advected_values,
                        all_grid
                );
            }
            else if(interpolation_method == "WENO6"){
                velocity_x_on_barycenter_of_i_tri
                    = WENO6_interpolation_3D_cell_face_x_values(
                        quadrature_point_list_for_interpolate_velocity_x._quadtarure_position[i_thread][i_point],
                        advected_values,
                        all_grid
                );
            }
            else if(interpolation_method == "1Dy_WENO6"){
                velocity_x_on_barycenter_of_i_tri
                    = WENO6_interpolation_1Dy_cell_face_x_values(
                        quadrature_point_list_for_interpolate_velocity_x._quadtarure_position[i_thread][i_point],
                        advected_values,
                        all_grid
                );
            }

            velocity_x_after_advect[
                get_voxel_face_index_x_3D(
                    quadrature_point_list_for_interpolate_velocity_x._included_cell_index[i_thread][i_point][0],
                    quadrature_point_list_for_interpolate_velocity_x._included_cell_index[i_thread][i_point][1],
                    quadrature_point_list_for_interpolate_velocity_x._included_cell_index[i_thread][i_point][2],
                    all_grid.Grid_num_x,
                    all_grid.Grid_num_y,
                    all_grid.Grid_num_z)
            ]   +=velocity_x_on_barycenter_of_i_tri
                * quadrature_point_list_for_interpolate_velocity_x._weighted_area_normal[i_thread][i_point][1]
                / all_grid._cell_volume;
        }
        ////////////////////////////////////////
        // psiからの寄与
        ////////////////////////////////////////
        const int num_quadrature_points_for_interpolate_psi=quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread].size();
        //求積点を走査するループ
        for(size_t i_point = 0; i_point < num_quadrature_points_for_interpolate_psi; ++i_point){
            ////補間に依ってpsiを求める
            VEC3_TYPE psi_velocity_x_on_barycenter_of_i_tri
                = all_grid.calc_interpolated_psi_velocity_x(quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread][i_point], interpolation_method);

            velocity_x_after_advect[
                get_voxel_face_index_x_3D(
                    quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][0],
                    quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][1],
                    quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][2],
                    all_grid.Grid_num_x,
                    all_grid.Grid_num_y,
                    all_grid.Grid_num_z)]
                +=psi_velocity_x_on_barycenter_of_i_tri.dot(quadrature_point_list_for_interpolate_psi._weighted_area_normal[i_thread][i_point])
                / all_grid._cell_volume;
        }
    }

    void calc_velocity_y_from_quadrature_point_list_use_integral(
        const quadrature_point_vector &quadrature_point_list_for_interpolate_velocity_y,
        const quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &velocity_y_after_advect,
        const int i_thread,
        const std::string interpolation_method,
        const std::vector<MY_FLOAT_TYPE> &advected_values
    ){
        ////////////////////////////////////////
        // densityからの寄与
        ////////////////////////////////////////
        const int num_quadrature_points_for_interpolate_velocity_y = quadrature_point_list_for_interpolate_velocity_y._quadtarure_position[i_thread].size();
        //求積点を走査するループ
        for(size_t i_point = 0; i_point < num_quadrature_points_for_interpolate_velocity_y; ++i_point){
            ////補間に依ってdensityを求める
            MY_FLOAT_TYPE velocity_y_on_barycenter_of_i_tri;
            if(interpolation_method == "linear"){
                velocity_y_on_barycenter_of_i_tri
                    = linear_interpolation_3D_cell_face_y_values(
                        quadrature_point_list_for_interpolate_velocity_y._quadtarure_position[i_thread][i_point],
                        advected_values,
                        all_grid
                );
            }
            else if(interpolation_method == "1Dy_linear"){
                velocity_y_on_barycenter_of_i_tri
                    = linear_interpolation_1Dy_cell_face_y_values(
                        quadrature_point_list_for_interpolate_velocity_y._quadtarure_position[i_thread][i_point],
                        advected_values,
                        all_grid
                );
            }
            else if(interpolation_method == "WENO6"){
                velocity_y_on_barycenter_of_i_tri
                    = WENO6_interpolation_3D_cell_face_y_values(
                        quadrature_point_list_for_interpolate_velocity_y._quadtarure_position[i_thread][i_point],
                        advected_values,
                        all_grid
                );
            }
            else if(interpolation_method == "1Dy_WENO6"){
                velocity_y_on_barycenter_of_i_tri
                    = WENO6_interpolation_1Dy_cell_face_y_values(
                        quadrature_point_list_for_interpolate_velocity_y._quadtarure_position[i_thread][i_point],
                        advected_values,
                        all_grid
                );
            }

            velocity_y_after_advect[
                get_voxel_face_index_y_3D(
                    quadrature_point_list_for_interpolate_velocity_y._included_cell_index[i_thread][i_point][0],
                    quadrature_point_list_for_interpolate_velocity_y._included_cell_index[i_thread][i_point][1],
                    quadrature_point_list_for_interpolate_velocity_y._included_cell_index[i_thread][i_point][2],
                    all_grid.Grid_num_x,
                    all_grid.Grid_num_y,
                    all_grid.Grid_num_z)
            ]   +=velocity_y_on_barycenter_of_i_tri
                * quadrature_point_list_for_interpolate_velocity_y._weighted_area_normal[i_thread][i_point][1]
                / all_grid._cell_volume;
        }
        ////////////////////////////////////////
        // psiからの寄与
        ////////////////////////////////////////
        const int num_quadrature_points_for_interpolate_psi=quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread].size();
        //求積点を走査するループ
        for(size_t i_point = 0; i_point < num_quadrature_points_for_interpolate_psi; ++i_point){
            ////補間に依ってpsiを求める
            VEC3_TYPE psi_velocity_y_on_barycenter_of_i_tri
                = all_grid.calc_interpolated_psi_velocity_y(quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread][i_point], interpolation_method);

            velocity_y_after_advect[
                get_voxel_face_index_y_3D(
                    quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][0],
                    quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][1],
                    quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][2],
                    all_grid.Grid_num_x,
                    all_grid.Grid_num_y,
                    all_grid.Grid_num_z)]
                +=psi_velocity_y_on_barycenter_of_i_tri.dot(quadrature_point_list_for_interpolate_psi._weighted_area_normal[i_thread][i_point])
                / all_grid._cell_volume;
        }
    }

    void calc_velocity_z_from_quadrature_point_list_use_integral(
        const quadrature_point_vector &quadrature_point_list_for_interpolate_velocity_z,
        const quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &velocity_z_after_advect,
        const int i_thread,
        const std::string interpolation_method,
        const std::vector<MY_FLOAT_TYPE> &advected_values
    ){
        ////////////////////////////////////////
        // densityからの寄与
        ////////////////////////////////////////
        const int num_quadrature_points_for_interpolate_velocity_z = quadrature_point_list_for_interpolate_velocity_z._quadtarure_position[i_thread].size();
        //求積点を走査するループ
        for(size_t i_point = 0; i_point < num_quadrature_points_for_interpolate_velocity_z; ++i_point){
            ////補間に依ってdensityを求める
            MY_FLOAT_TYPE velocity_z_on_barycenter_of_i_tri;
            if(interpolation_method == "linear"){
                velocity_z_on_barycenter_of_i_tri
                    = linear_interpolation_3D_cell_face_z_values(
                        quadrature_point_list_for_interpolate_velocity_z._quadtarure_position[i_thread][i_point],
                        advected_values,
                        all_grid
                );
            }
            else if(interpolation_method == "1Dy_linear"){
                velocity_z_on_barycenter_of_i_tri
                    = linear_interpolation_1Dy_cell_face_z_values(
                        quadrature_point_list_for_interpolate_velocity_z._quadtarure_position[i_thread][i_point],
                        advected_values,
                        all_grid
                );
            }
            else if(interpolation_method == "WENO6"){
                velocity_z_on_barycenter_of_i_tri
                    = WENO6_interpolation_3D_cell_face_z_values(
                        quadrature_point_list_for_interpolate_velocity_z._quadtarure_position[i_thread][i_point],
                        advected_values,
                        all_grid
                );
            }
            else if(interpolation_method == "1Dy_WENO6"){
                velocity_z_on_barycenter_of_i_tri
                    = WENO6_interpolation_1Dy_cell_face_z_values(
                        quadrature_point_list_for_interpolate_velocity_z._quadtarure_position[i_thread][i_point],
                        advected_values,
                        all_grid
                );
            }

            velocity_z_after_advect[
                get_voxel_face_index_z_3D(
                    quadrature_point_list_for_interpolate_velocity_z._included_cell_index[i_thread][i_point][0],
                    quadrature_point_list_for_interpolate_velocity_z._included_cell_index[i_thread][i_point][1],
                    quadrature_point_list_for_interpolate_velocity_z._included_cell_index[i_thread][i_point][2],
                    all_grid.Grid_num_x,
                    all_grid.Grid_num_y,
                    all_grid.Grid_num_z)
            ]   +=velocity_z_on_barycenter_of_i_tri
                * quadrature_point_list_for_interpolate_velocity_z._weighted_area_normal[i_thread][i_point][1]
                / all_grid._cell_volume;
        }
        ////////////////////////////////////////
        // psiからの寄与
        ////////////////////////////////////////
        const int num_quadrature_points_for_interpolate_psi=quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread].size();
        //求積点を走査するループ
        for(size_t i_point = 0; i_point < num_quadrature_points_for_interpolate_psi; ++i_point){
            ////補間に依ってpsiを求める
            VEC3_TYPE psi_velocity_z_on_barycenter_of_i_tri
                = all_grid.calc_interpolated_psi_velocity_z(quadrature_point_list_for_interpolate_psi._quadtarure_position[i_thread][i_point], interpolation_method);

            velocity_z_after_advect[
                get_voxel_face_index_z_3D(
                    quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][0],
                    quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][1],
                    quadrature_point_list_for_interpolate_psi._included_cell_index[i_thread][i_point][2],
                    all_grid.Grid_num_x,
                    all_grid.Grid_num_y,
                    all_grid.Grid_num_z)]
                +=psi_velocity_z_on_barycenter_of_i_tri.dot(quadrature_point_list_for_interpolate_psi._weighted_area_normal[i_thread][i_point])
                / all_grid._cell_volume;
        }
    }


}// namespace smoke_simulation
