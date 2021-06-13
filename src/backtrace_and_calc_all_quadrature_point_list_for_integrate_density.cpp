#include "backtrace_and_calc_all_quadrature_point_list_for_integrate_density.h"

#include "calc_backtraced_face_3D.h"

namespace smoke_simulation{
    //バックトレースを行った後求積点のリストを作る関数
    void backtrace_and_calc_all_quadrature_point_list_for_integrate_density(
        const int dim,
        quadrature_point_vector &quadrature_point_list_interpolate_density,
        quadrature_point_vector &quadrature_point_list_interpolate_psi,
        quadrature_point_vector &quadrature_point_list_for_calc_cell_volume,
        size_t i_start,
        size_t i_end,
        bool use_MacCormack_scheme,
        const Grid_3D &all_grid,
        MY_FLOAT_TYPE time_step_length,
//        bool set_zero_normal_velocity_at_boundary,
        std::string split_method,
//        bool is_zero_velocity_false_diffusion_correction_mode,
        const std::string integral_method,
        const int num_gauss_quadrature_points,
        int i_thread,
        const std::string interpolation_method,
        const bool use_zero_velocity_for_backtrace,
	    const int num_gauss_quadrature_point_for_integrate_density,
        const bool enable_cell_volume_correction,
        const std::string backtrace_usage,
		const MY_FLOAT_TYPE minimal_area
    ){
        //グリッドの原点の位置
        VEC3_TYPE origin_pos_of_grid;
        if(backtrace_usage == "velocity_x"){
            origin_pos_of_grid = VEC3_TYPE(-0.5 * all_grid._cell_length, 0.0, 0.0);
        }
        else if(backtrace_usage == "velocity_y"){
            origin_pos_of_grid = VEC3_TYPE(0.0, -0.5 * all_grid._cell_length, 0.0);
        }
        else if(backtrace_usage == "velocity_z"){
            origin_pos_of_grid = VEC3_TYPE(0.0, 0.0, -0.5 * all_grid._cell_length);
        }
        else if(backtrace_usage == "density"){
            origin_pos_of_grid = VEC3_TYPE(0.0, 0.0, 0.0);
        }

        for(size_t i_xyzface = i_start; i_xyzface < i_end; ++i_xyzface){
            int i_face;
            int iz;
            int iy;
            int ix;
            if(backtrace_usage == "velocity_x"){
                i_face = i_xyzface % 6;
                iz = ((i_xyzface - i_face) / 6) % all_grid.Grid_num_z;
                iy = (((i_xyzface - i_face) / 6) / all_grid.Grid_num_z) % all_grid.Grid_num_y;
                ix = (i_xyzface - i_face - 6 * iz - 6 * all_grid.Grid_num_z * iy) / (6 * all_grid.Grid_num_z * all_grid.Grid_num_y);
            }
            else if(backtrace_usage == "velocity_y"){
                i_face = i_xyzface % 6;
                iz = ((i_xyzface - i_face) / 6) % all_grid.Grid_num_z;
                iy = (((i_xyzface - i_face) / 6) / all_grid.Grid_num_z) % (all_grid.Grid_num_y + 1);
                ix = (i_xyzface - i_face - 6 * iz - 6 * all_grid.Grid_num_z * iy) / (6 * all_grid.Grid_num_z * (all_grid.Grid_num_y + 1));
            }
            else if(backtrace_usage == "velocity_z"){
                i_face = i_xyzface % 6;
                iz = ((i_xyzface - i_face) / 6) % (all_grid.Grid_num_z + 1);
                iy = (((i_xyzface - i_face) / 6) / (all_grid.Grid_num_z + 1)) % all_grid.Grid_num_y;
                ix = (i_xyzface - i_face - 6 * iz - 6 * (all_grid.Grid_num_z + 1) * iy) / (6 * (all_grid.Grid_num_z + 1) * all_grid.Grid_num_y);
            }
            else if(backtrace_usage == "density"){
                i_face = i_xyzface % 6;
                iz = ((i_xyzface - i_face) / 6) % all_grid.Grid_num_z;
                iy = (((i_xyzface - i_face) / 6) / all_grid.Grid_num_z) % all_grid.Grid_num_y;
                ix = (i_xyzface - i_face - 6 * iz - 6 * all_grid.Grid_num_z * iy) / (6 * all_grid.Grid_num_z * all_grid.Grid_num_y);
            }
            // 考えるcell face を構成するvertex のindex
            // OK
            Eigen::Vector3i vertex_index_1, vertex_index_2, vertex_index_3, vertex_index_4;
            if(i_face == 0){
                vertex_index_1[0] = ix    ; vertex_index_1[1] = iy    ; vertex_index_1[2] = iz;
                vertex_index_2[0] = ix    ; vertex_index_2[1] = iy + 1; vertex_index_2[2] = iz;
                vertex_index_3[0] = ix + 1; vertex_index_3[1] = iy + 1; vertex_index_3[2] = iz;
                vertex_index_4[0] = ix + 1; vertex_index_4[1] = iy    ; vertex_index_4[2] = iz;
            }
            else if(i_face == 1){
                vertex_index_1[0] = ix + 1; vertex_index_1[1] = iy    ; vertex_index_1[2] = iz    ;
                vertex_index_2[0] = ix + 1; vertex_index_2[1] = iy + 1; vertex_index_2[2] = iz    ;
                vertex_index_3[0] = ix + 1; vertex_index_3[1] = iy + 1; vertex_index_3[2] = iz + 1;
                vertex_index_4[0] = ix + 1; vertex_index_4[1] = iy    ; vertex_index_4[2] = iz + 1;
            }
            else if(i_face == 2){
                vertex_index_1[0] = ix    ; vertex_index_1[1] = iy + 1; vertex_index_1[2] = iz;
                vertex_index_2[0] = ix    ; vertex_index_2[1] = iy + 1; vertex_index_2[2] = iz + 1;
                vertex_index_3[0] = ix + 1; vertex_index_3[1] = iy + 1; vertex_index_3[2] = iz + 1;
                vertex_index_4[0] = ix + 1; vertex_index_4[1] = iy + 1; vertex_index_4[2] = iz;
            }
            else if(i_face == 3){
                vertex_index_1[0] = ix    ; vertex_index_1[1] = iy    ; vertex_index_1[2] = iz    ;
                vertex_index_2[0] = ix    ; vertex_index_2[1] = iy    ; vertex_index_2[2] = iz + 1;
                vertex_index_3[0] = ix    ; vertex_index_3[1] = iy + 1; vertex_index_3[2] = iz + 1;
                vertex_index_4[0] = ix    ; vertex_index_4[1] = iy + 1; vertex_index_4[2] = iz    ;
            }
            else if(i_face == 4){
                vertex_index_1[0] = ix    ; vertex_index_1[1] = iy    ; vertex_index_1[2] = iz;
                vertex_index_2[0] = ix + 1; vertex_index_2[1] = iy    ; vertex_index_2[2] = iz;
                vertex_index_3[0] = ix + 1; vertex_index_3[1] = iy    ; vertex_index_3[2] = iz + 1;
                vertex_index_4[0] = ix    ; vertex_index_4[1] = iy    ; vertex_index_4[2] = iz + 1;
            }
            else if(i_face == 5){
                vertex_index_1[0] = ix    ; vertex_index_1[1] = iy    ; vertex_index_1[2] = iz + 1;
                vertex_index_2[0] = ix + 1; vertex_index_2[1] = iy    ; vertex_index_2[2] = iz + 1;
                vertex_index_3[0] = ix + 1; vertex_index_3[1] = iy + 1; vertex_index_3[2] = iz + 1;
                vertex_index_4[0] = ix    ; vertex_index_4[1] = iy + 1; vertex_index_4[2] = iz + 1;
            }
            // バックトレース前の cell face を構成するvertexを定義する
            cell_vertex_3D before_backtrace_vertex_1(vertex_index_1, all_grid._cell_length);
            cell_vertex_3D before_backtrace_vertex_2(vertex_index_2, all_grid._cell_length);
            cell_vertex_3D before_backtrace_vertex_3(vertex_index_3, all_grid._cell_length);
            cell_vertex_3D before_backtrace_vertex_4(vertex_index_4, all_grid._cell_length);
            if(backtrace_usage == "velocity_x"){
                before_backtrace_vertex_1._vertex_pos[0] -= 0.5 *all_grid._cell_length;
                before_backtrace_vertex_2._vertex_pos[0] -= 0.5 *all_grid._cell_length;
                before_backtrace_vertex_3._vertex_pos[0] -= 0.5 *all_grid._cell_length;
                before_backtrace_vertex_4._vertex_pos[0] -= 0.5 *all_grid._cell_length;
            }
            else if(backtrace_usage == "velocity_y"){
                before_backtrace_vertex_1._vertex_pos[1] -= 0.5 *all_grid._cell_length;
                before_backtrace_vertex_2._vertex_pos[1] -= 0.5 *all_grid._cell_length;
                before_backtrace_vertex_3._vertex_pos[1] -= 0.5 *all_grid._cell_length;
                before_backtrace_vertex_4._vertex_pos[1] -= 0.5 *all_grid._cell_length;
            }
            else if(backtrace_usage == "velocity_z"){
                before_backtrace_vertex_1._vertex_pos[2] -= 0.5 *all_grid._cell_length;
                before_backtrace_vertex_2._vertex_pos[2] -= 0.5 *all_grid._cell_length;
                before_backtrace_vertex_3._vertex_pos[2] -= 0.5 *all_grid._cell_length;
                before_backtrace_vertex_4._vertex_pos[2] -= 0.5 *all_grid._cell_length;
            }
            else if(backtrace_usage == "density"){
            }

            // バックトレース前の頂点から面を定義
            cell_face_3D before_backtrace_face;
            before_backtrace_face._vertex_list[0] = before_backtrace_vertex_1;
            before_backtrace_face._vertex_list[1] = before_backtrace_vertex_2;
            before_backtrace_face._vertex_list[2] = before_backtrace_vertex_3;
            before_backtrace_face._vertex_list[3] = before_backtrace_vertex_4;

            // MacCormack を使う場合
            if(use_MacCormack_scheme){
                ////移流しない時の求積点からの寄与
                //求積点を追加
                split_face_and_add_quadrature_points_3D_for_integrate_density_ref(
                    quadrature_point_list_interpolate_density,
                    quadrature_point_list_interpolate_psi,
                    quadrature_point_list_for_calc_cell_volume,
                    all_grid,
                    before_backtrace_face,
                    split_method,
                    i_thread,
//                    negative_weight,
                    ix,
                    iy,
                    iz,
                    integral_method,
                    num_gauss_quadrature_points,
                    0.5,
            	    num_gauss_quadrature_point_for_integrate_density,
                    enable_cell_volume_correction,
                    origin_pos_of_grid,
            		minimal_area
                );
                // バックトレース後の頂点から面を定義
                cell_face_3D after_backtrace_face;
                //速度場0の拡散の補正を計算する際はバックトレースしない
                if(use_zero_velocity_for_backtrace){
                    after_backtrace_face = before_backtrace_face;
                }
                //バックトレースの処理
                else{
                    calc_backtraced_face_3D(
                        dim,
                        all_grid,
                        time_step_length,
                        before_backtrace_face,
                        after_backtrace_face,
                        "forward",
//                        set_zero_normal_velocity_at_boundary,
                        interpolation_method,
                        true,
                        before_backtrace_face._vertex_list[0]._vertex_index,
                        before_backtrace_face._vertex_list[1]._vertex_index,
                        before_backtrace_face._vertex_list[2]._vertex_index,
                        before_backtrace_face._vertex_list[3]._vertex_index
                    );
                    //バックトレースした先の位置がグリッドの範囲を超えていたら修正する
                    clamp_vertex_position_by_grid_3D(dim, all_grid, after_backtrace_face);
                }
                //求積点を追加
                split_face_and_add_quadrature_points_3D_for_integrate_density_ref(
                    quadrature_point_list_interpolate_density,
                    quadrature_point_list_interpolate_psi,
                    quadrature_point_list_for_calc_cell_volume,
                    all_grid,
                    after_backtrace_face,
                    split_method,
                    i_thread,
//                    negative_weight,
                    ix,
                    iy,
                    iz,
                    integral_method,
                    num_gauss_quadrature_points,
                    1.0,
            	    num_gauss_quadrature_point_for_integrate_density,
                    enable_cell_volume_correction,
                    origin_pos_of_grid,
            		minimal_area
                );
                //// 2回目のバックトレース
                // バックトレース後の面を逆向きにバックトレースした面の計算(エラーの計算用に使う)
                cell_face_3D backward_trace_face_of_after_backtrace_face;
                if(use_zero_velocity_for_backtrace){
                    backward_trace_face_of_after_backtrace_face = after_backtrace_face;
                }
                //バックトレースの処理
                else{
                    calc_backtraced_face_3D(
                        dim,
                        all_grid,
                        time_step_length,
                        after_backtrace_face,
                        backward_trace_face_of_after_backtrace_face,
                        "backward",
//                        set_zero_normal_velocity_at_boundary,
                        interpolation_method,
                        true,
                        before_backtrace_face._vertex_list[0]._vertex_index,
                        before_backtrace_face._vertex_list[1]._vertex_index,
                        before_backtrace_face._vertex_list[2]._vertex_index,
                        before_backtrace_face._vertex_list[3]._vertex_index
                    );
                    //バックトレースした先の位置がグリッドの範囲を超えていたら修正する
                    clamp_vertex_position_by_grid_3D(dim, all_grid, backward_trace_face_of_after_backtrace_face);
                }
                //求積点を追加
                split_face_and_add_quadrature_points_3D_for_integrate_density_ref(
                    quadrature_point_list_interpolate_density,
                    quadrature_point_list_interpolate_psi,
                    quadrature_point_list_for_calc_cell_volume,
                    all_grid,
                    backward_trace_face_of_after_backtrace_face,
                    split_method,
                    i_thread,
//                    negative_weight,
                    ix,
                    iy,
                    iz,
                    integral_method,
                    num_gauss_quadrature_points,
                    -0.5,
            	    num_gauss_quadrature_point_for_integrate_density,
                    enable_cell_volume_correction,
                    origin_pos_of_grid,
            		minimal_area
                );
            }
            // MacCormack を使わない(semi-Lagrangian を使う)場合
            else{
                // バックトレース後の頂点から面を定義
                cell_face_3D after_backtrace_face;
                //速度場0の拡散の補正を計算する際はバックトレースしない
                if(use_zero_velocity_for_backtrace){
                    after_backtrace_face = before_backtrace_face;
                }
                else{
                    calc_backtraced_face_3D(
                        dim,
                        all_grid,
                        time_step_length,
                        before_backtrace_face,
                        after_backtrace_face,
                        "forward",
//                        set_zero_normal_velocity_at_boundary,
                        interpolation_method,
                        true,
                        before_backtrace_face._vertex_list[0]._vertex_index,
                        before_backtrace_face._vertex_list[1]._vertex_index,
                        before_backtrace_face._vertex_list[2]._vertex_index,
                        before_backtrace_face._vertex_list[3]._vertex_index
                    );
                    //バックトレースした先の位置がグリッドの範囲を超えていたら修正する
                    clamp_vertex_position_by_grid_3D(dim, all_grid, after_backtrace_face);
                }
                split_face_and_add_quadrature_points_3D_for_integrate_density_ref(
                    quadrature_point_list_interpolate_density,
                    quadrature_point_list_interpolate_psi,
                    quadrature_point_list_for_calc_cell_volume,
                    all_grid,
                    after_backtrace_face,
                    split_method,
                    i_thread,
//                    negative_weight,
                    ix,
                    iy,
                    iz,
                    integral_method,
                    num_gauss_quadrature_points,
                    1.0,
            	    num_gauss_quadrature_point_for_integrate_density,
                    enable_cell_volume_correction,
                    origin_pos_of_grid,
            		minimal_area
                );
            }
        }
    }

    //バックトレースを行った後求積点のリストを作る関数
    void backtrace_and_calc_all_quadrature_point_list_for_integrate_velocity(
        const int dim,
        quadrature_point_vector &quadrature_point_list_interpolate_velocity,
        quadrature_point_vector &quadrature_point_list_interpolate_psi,
        quadrature_point_vector &quadrature_point_list_for_calc_cell_volume,
        size_t i_start,
        size_t i_end,
        bool use_MacCormack_scheme,
        const Grid_3D &all_grid,
        MY_FLOAT_TYPE time_step_length,
//        bool set_zero_normal_velocity_at_boundary,
        std::string split_method,
//        bool is_zero_velocity_false_diffusion_correction_mode,
        const std::string integral_method,
        const int num_gauss_quadrature_points,
        int i_thread,
        const std::string interpolation_method,
        const bool use_zero_velocity_for_backtrace,
	    const int num_gauss_quadrature_point_for_integrate_density,
        const bool enable_cell_volume_correction,
		const MY_FLOAT_TYPE minimal_area
    ){
        const int dx = (dim == 0);
        const int dy = (dim == 1);
        const int dz = (dim == 2);

        //グリッドの原点の位置
        VEC3_TYPE origin_pos_of_grid;
        origin_pos_of_grid
            = VEC3_TYPE(
                - dx * 0.5 * all_grid._cell_length,
                - dy * 0.5 * all_grid._cell_length,
                - dz * 0.5 * all_grid._cell_length
            );

        for(size_t i_xyzface = i_start; i_xyzface < i_end; ++i_xyzface){
            int i_face;
            int iz;
            int iy;
            int ix;
            i_face = i_xyzface % 6;
            iz = ((i_xyzface - i_face) / 6) % (all_grid.Grid_num_z + dz);
            iy = (((i_xyzface - i_face) / 6) / (all_grid.Grid_num_z + dz)) % (all_grid.Grid_num_y + dy);
            ix = (i_xyzface - i_face - 6 * iz - 6 * (all_grid.Grid_num_z + dz) * iy) / (6 * (all_grid.Grid_num_z + dz) * (all_grid.Grid_num_y + dy));
            // 考えるcell face を構成するvertex のindex
            // OK
            Eigen::Vector3i vertex_index_1, vertex_index_2, vertex_index_3, vertex_index_4;
            if(i_face == 0){
                vertex_index_1[0] = ix    ; vertex_index_1[1] = iy    ; vertex_index_1[2] = iz;
                vertex_index_2[0] = ix    ; vertex_index_2[1] = iy + 1; vertex_index_2[2] = iz;
                vertex_index_3[0] = ix + 1; vertex_index_3[1] = iy + 1; vertex_index_3[2] = iz;
                vertex_index_4[0] = ix + 1; vertex_index_4[1] = iy    ; vertex_index_4[2] = iz;
            }
            else if(i_face == 1){
                vertex_index_1[0] = ix + 1; vertex_index_1[1] = iy    ; vertex_index_1[2] = iz    ;
                vertex_index_2[0] = ix + 1; vertex_index_2[1] = iy + 1; vertex_index_2[2] = iz    ;
                vertex_index_3[0] = ix + 1; vertex_index_3[1] = iy + 1; vertex_index_3[2] = iz + 1;
                vertex_index_4[0] = ix + 1; vertex_index_4[1] = iy    ; vertex_index_4[2] = iz + 1;
            }
            else if(i_face == 2){
                vertex_index_1[0] = ix    ; vertex_index_1[1] = iy + 1; vertex_index_1[2] = iz;
                vertex_index_2[0] = ix    ; vertex_index_2[1] = iy + 1; vertex_index_2[2] = iz + 1;
                vertex_index_3[0] = ix + 1; vertex_index_3[1] = iy + 1; vertex_index_3[2] = iz + 1;
                vertex_index_4[0] = ix + 1; vertex_index_4[1] = iy + 1; vertex_index_4[2] = iz;
            }
            else if(i_face == 3){
                vertex_index_1[0] = ix    ; vertex_index_1[1] = iy    ; vertex_index_1[2] = iz    ;
                vertex_index_2[0] = ix    ; vertex_index_2[1] = iy    ; vertex_index_2[2] = iz + 1;
                vertex_index_3[0] = ix    ; vertex_index_3[1] = iy + 1; vertex_index_3[2] = iz + 1;
                vertex_index_4[0] = ix    ; vertex_index_4[1] = iy + 1; vertex_index_4[2] = iz    ;
            }
            else if(i_face == 4){
                vertex_index_1[0] = ix    ; vertex_index_1[1] = iy    ; vertex_index_1[2] = iz;
                vertex_index_2[0] = ix + 1; vertex_index_2[1] = iy    ; vertex_index_2[2] = iz;
                vertex_index_3[0] = ix + 1; vertex_index_3[1] = iy    ; vertex_index_3[2] = iz + 1;
                vertex_index_4[0] = ix    ; vertex_index_4[1] = iy    ; vertex_index_4[2] = iz + 1;
            }
            else if(i_face == 5){
                vertex_index_1[0] = ix    ; vertex_index_1[1] = iy    ; vertex_index_1[2] = iz + 1;
                vertex_index_2[0] = ix + 1; vertex_index_2[1] = iy    ; vertex_index_2[2] = iz + 1;
                vertex_index_3[0] = ix + 1; vertex_index_3[1] = iy + 1; vertex_index_3[2] = iz + 1;
                vertex_index_4[0] = ix    ; vertex_index_4[1] = iy + 1; vertex_index_4[2] = iz + 1;
            }
            // バックトレース前の cell face を構成するvertexを定義する
            cell_vertex_3D before_backtrace_vertex_1(vertex_index_1, all_grid._cell_length);
            cell_vertex_3D before_backtrace_vertex_2(vertex_index_2, all_grid._cell_length);
            cell_vertex_3D before_backtrace_vertex_3(vertex_index_3, all_grid._cell_length);
            cell_vertex_3D before_backtrace_vertex_4(vertex_index_4, all_grid._cell_length);

            before_backtrace_vertex_1._vertex_pos[0] -= dx * 0.5 * all_grid._cell_length;
            before_backtrace_vertex_2._vertex_pos[0] -= dx * 0.5 * all_grid._cell_length;
            before_backtrace_vertex_3._vertex_pos[0] -= dx * 0.5 * all_grid._cell_length;
            before_backtrace_vertex_4._vertex_pos[0] -= dx * 0.5 * all_grid._cell_length;

            before_backtrace_vertex_1._vertex_pos[1] -= dy * 0.5 * all_grid._cell_length;
            before_backtrace_vertex_2._vertex_pos[1] -= dy * 0.5 * all_grid._cell_length;
            before_backtrace_vertex_3._vertex_pos[1] -= dy * 0.5 * all_grid._cell_length;
            before_backtrace_vertex_4._vertex_pos[1] -= dy * 0.5 * all_grid._cell_length;

            before_backtrace_vertex_1._vertex_pos[2] -= dz * 0.5 * all_grid._cell_length;
            before_backtrace_vertex_2._vertex_pos[2] -= dz * 0.5 * all_grid._cell_length;
            before_backtrace_vertex_3._vertex_pos[2] -= dz * 0.5 * all_grid._cell_length;
            before_backtrace_vertex_4._vertex_pos[2] -= dz * 0.5 * all_grid._cell_length;

            // バックトレース前の頂点から面を定義
            cell_face_3D before_backtrace_face;
            before_backtrace_face._vertex_list[0] = before_backtrace_vertex_1;
            before_backtrace_face._vertex_list[1] = before_backtrace_vertex_2;
            before_backtrace_face._vertex_list[2] = before_backtrace_vertex_3;
            before_backtrace_face._vertex_list[3] = before_backtrace_vertex_4;

            // MacCormack を使う場合
            if(use_MacCormack_scheme){
                ////移流しない時の求積点からの寄与
                //求積点を追加
                split_face_and_add_quadrature_points_3D_for_integrate_density_ref(
                    quadrature_point_list_interpolate_velocity,
                    quadrature_point_list_interpolate_psi,
                    quadrature_point_list_for_calc_cell_volume,
                    all_grid,
                    before_backtrace_face,
                    split_method,
                    i_thread,
//                    negative_weight,
                    ix,
                    iy,
                    iz,
                    integral_method,
                    num_gauss_quadrature_points,
                    0.5,
            	    num_gauss_quadrature_point_for_integrate_density,
                    enable_cell_volume_correction,
                    origin_pos_of_grid,
            		minimal_area
                );
                // バックトレース後の頂点から面を定義
                cell_face_3D after_backtrace_face;
                //速度場0の拡散の補正を計算する際はバックトレースしない
                if(use_zero_velocity_for_backtrace){
                    after_backtrace_face = before_backtrace_face;
                }
                //バックトレースの処理
                else{
                    calc_backtraced_face_3D(
                        dim,
                        all_grid,
                        time_step_length,
                        before_backtrace_face,
                        after_backtrace_face,
                        "forward",
//                        set_zero_normal_velocity_at_boundary,
                        interpolation_method,
                        true,
                        before_backtrace_face._vertex_list[0]._vertex_index,
                        before_backtrace_face._vertex_list[1]._vertex_index,
                        before_backtrace_face._vertex_list[2]._vertex_index,
                        before_backtrace_face._vertex_list[3]._vertex_index
                    );
                    //バックトレースした先の位置がグリッドの範囲を超えていたら修正する
                    clamp_vertex_position_by_grid_3D(dim, all_grid, after_backtrace_face);
                }
                //求積点を追加
                split_face_and_add_quadrature_points_3D_for_integrate_density_ref(
                    quadrature_point_list_interpolate_velocity,
                    quadrature_point_list_interpolate_psi,
                    quadrature_point_list_for_calc_cell_volume,
                    all_grid,
                    after_backtrace_face,
                    split_method,
                    i_thread,
//                    negative_weight,
                    ix,
                    iy,
                    iz,
                    integral_method,
                    num_gauss_quadrature_points,
                    1.0,
            	    num_gauss_quadrature_point_for_integrate_density,
                    enable_cell_volume_correction,
                    origin_pos_of_grid,
            		minimal_area
                );
                //// 2回目のバックトレース
                // バックトレース後の面を逆向きにバックトレースした面の計算(エラーの計算用に使う)
                cell_face_3D backward_trace_face_of_after_backtrace_face;
                if(use_zero_velocity_for_backtrace){
                    backward_trace_face_of_after_backtrace_face = after_backtrace_face;
                }
                //バックトレースの処理
                else{
                    calc_backtraced_face_3D(
                        dim,
                        all_grid,
                        time_step_length,
                        after_backtrace_face,
                        backward_trace_face_of_after_backtrace_face,
                        "backward",
//                        set_zero_normal_velocity_at_boundary,
                        interpolation_method,
                        true,
                        before_backtrace_face._vertex_list[0]._vertex_index,
                        before_backtrace_face._vertex_list[1]._vertex_index,
                        before_backtrace_face._vertex_list[2]._vertex_index,
                        before_backtrace_face._vertex_list[3]._vertex_index
                    );
                    //バックトレースした先の位置がグリッドの範囲を超えていたら修正する
                    clamp_vertex_position_by_grid_3D(dim, all_grid, backward_trace_face_of_after_backtrace_face);
                }
                //求積点を追加
                split_face_and_add_quadrature_points_3D_for_integrate_density_ref(
                    quadrature_point_list_interpolate_velocity,
                    quadrature_point_list_interpolate_psi,
                    quadrature_point_list_for_calc_cell_volume,
                    all_grid,
                    backward_trace_face_of_after_backtrace_face,
                    split_method,
                    i_thread,
//                    negative_weight,
                    ix,
                    iy,
                    iz,
                    integral_method,
                    num_gauss_quadrature_points,
                    -0.5,
            	    num_gauss_quadrature_point_for_integrate_density,
                    enable_cell_volume_correction,
                    origin_pos_of_grid,
            		minimal_area
                );
            }
            // MacCormack を使わない(semi-Lagrangian を使う)場合
            else{
                // バックトレース後の頂点から面を定義
                cell_face_3D after_backtrace_face;
                //速度場0の拡散の補正を計算する際はバックトレースしない
                if(use_zero_velocity_for_backtrace){
                    after_backtrace_face = before_backtrace_face;
                }
                else{
                    calc_backtraced_face_3D(
                        dim,
                        all_grid,
                        time_step_length,
                        before_backtrace_face,
                        after_backtrace_face,
                        "forward",
//                        set_zero_normal_velocity_at_boundary,
                        interpolation_method,
                        true,
                        before_backtrace_face._vertex_list[0]._vertex_index,
                        before_backtrace_face._vertex_list[1]._vertex_index,
                        before_backtrace_face._vertex_list[2]._vertex_index,
                        before_backtrace_face._vertex_list[3]._vertex_index
                    );
                    //バックトレースした先の位置がグリッドの範囲を超えていたら修正する
                    clamp_vertex_position_by_grid_3D(dim, all_grid, after_backtrace_face);
                }
                split_face_and_add_quadrature_points_3D_for_integrate_density_ref(
                    quadrature_point_list_interpolate_velocity,
                    quadrature_point_list_interpolate_psi,
                    quadrature_point_list_for_calc_cell_volume,
                    all_grid,
                    after_backtrace_face,
                    split_method,
                    i_thread,
//                    negative_weight,
                    ix,
                    iy,
                    iz,
                    integral_method,
                    num_gauss_quadrature_points,
                    1.0,
            	    num_gauss_quadrature_point_for_integrate_density,
                    enable_cell_volume_correction,
                    origin_pos_of_grid,
            		minimal_area
                );
            }
        }

    }
} //namespace smoke_simulation
