#include "calc_cell_volumes.h"

#include "calc_backtraced_face.h"
#include "utils.h"

namespace smoke_simulation{
    void calc_all_cell_volumes(
        Grid &all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const bool use_MacCormack_scheme,
        const std::string interpolation_method_for_velocity
    ){
        std::vector<MY_FLOAT_TYPE> cell_volumes_after_advect(all_grid.Grid_num_x * all_grid.Grid_num_y);

        for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                //(ix, iy)番目のセルのface を走るループ
                int dx[4], dy[4];
                dx[0] = 0; dx[1] = 1; dx[2] = 1; dx[3] = 0;
                dy[0] = 0; dy[1] = 0; dy[2] = 1; dy[3] = 1;
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
/*
                        ////エラーを補正する前のベースとなる面での計算
                        // バックトレース後の頂点から面を定義
                        cell_face after_backtrace_face
                            = calc_backtraced_face(
                                all_grid,
                                time_step_length,
                                before_backtrace_face,
                                "forward",
                                set_zero_normal_velocity_at_boundary,
                                interpolation_method);
                        //バックトレースした先の位置がグリッドの範囲を超えていたら修正する
                        after_backtrace_face
                            = clamp_vertex_position_by_grid(
                                all_grid,
                                after_backtrace_face);
                        // バックトレース後の面について psi と 法線の内積を面上で積分した値を計算する
                        MY_FLOAT_TYPE integral_of_normal_component_of_psi_on_after_backtrace_face
                            = integrate_normal_component_of_psi_on_face(
                                all_grid,
                                after_backtrace_face,
                                num_gauss_quad_boundary,
                                num_gauss_quad_bulk,
                                split_method,
                                interpolation_method);

                        //// エラーの計算
                        // 元のface(バックトレース前のオリジナルの face)について, psi と 法線の内積を面上で積分した値を計算する
                        MY_FLOAT_TYPE integral_of_normal_component_of_psi_on_before_backtrace_face
                            = integrate_normal_component_of_psi_on_face(
                                all_grid,
                                before_backtrace_face,
                                num_gauss_quad_boundary,
                                num_gauss_quad_bulk,
                                split_method,
                                interpolation_method);
                        // バックトレース後の面を逆向きにバックトレースした面の計算(エラーの計算用に使う)
                        cell_face backward_trace_face_of_after_backtrace_face
                            = calc_backtraced_face(
                                all_grid,
                                time_step_length,
                                after_backtrace_face,
                                "backward",
                                set_zero_normal_velocity_at_boundary,
                                interpolation_method);
                        //バックトレースした先の位置がグリッドの範囲を超えていたら修正する
                        backward_trace_face_of_after_backtrace_face
                            = clamp_vertex_position_by_grid(
                                all_grid,
                                backward_trace_face_of_after_backtrace_face);
                        // psi と 法線の内積を面上で積分した値を計算する
                        MY_FLOAT_TYPE integral_of_normal_component_of_psi_on_backward_trace_face_of_after_backtrace_face
                            = integrate_normal_component_of_psi_on_face(
                                all_grid,
                                backward_trace_face_of_after_backtrace_face,
                                num_gauss_quad_boundary,
                                num_gauss_quad_bulk,
                                split_method,
                                interpolation_method);
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
                        // エラーを修正した値を用いて質量への寄与を計算する
                        mass_in_the_cell
                            += integral_of_normal_component_of_psi_on_after_backtrace_face
                            + error;
*/
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
//                                set_zero_normal_velocity_at_boundary,
                                interpolation_method_for_velocity,
								"density"
                            );
                        //バックトレースした先の位置がグリッドの範囲を超えていたら修正する
                        after_backtrace_face = clamp_vertex_position_by_grid(
                            all_grid,
                            after_backtrace_face,
                            "density"
                        );
                        // バックトレース後のセルの体積を計算する
                        calc_cell_volume_contribution(
                            after_backtrace_face,
                            cell_volume
                        );
//                        std::cout<<"getchar()"<<std::endl;
//                        getchar();
                    }
                }
                cell_volumes_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = cell_volume;
            }
        }
        //計算結果をコピー
        for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                all_grid.cell_volume_cell_center[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
                    = cell_volumes_after_advect[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
            }
        }
    }

    void calc_cell_volume_contribution(
        const cell_face &face,
        MY_FLOAT_TYPE &cell_volume
    ){
        VEC3_TYPE center_position_of_cell_face = face.calc_face_center();
        VEC3_TYPE normal_of_cell_face = face.calc_normal();
        MY_FLOAT_TYPE face_area = face.calc_face_area();

        cell_volume += face_area * center_position_of_cell_face[1] * normal_of_cell_face[1];
    }
}
