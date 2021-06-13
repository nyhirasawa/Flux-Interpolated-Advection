#include "calc_backtraced_face.h"

#include "cell_face.h"
#include "define_float_type.h"

namespace smoke_simulation{

    // time_direction == "forward"  の場合は正の時間の方向に backtrace する(普通のsemi-Lagrangian と同じ)
    // time_direction == "backward" の場合は逆の時間の方向に backtrace する( MacCormack の2段階目で使う)
    // backtrace_usage == "density" or "velocity_x" or "velocity_y"
    cell_face calc_backtraced_face(
        const Grid& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const cell_face original_face,
        const std::string time_direction,
//        const bool set_zero_normal_velocity_at_boundary,
        const std::string interpolation_method,
        const std::string backtrace_usage
    ) {
        //cell vertex の速度を計算する
        VEC3_TYPE vertex_velocity_1;
        VEC3_TYPE vertex_velocity_2;

        if(backtrace_usage == "velocity_x") {
            // 0番目の頂点の速度を計算
            //壁を含んだセルの速度は0にする
            if(original_face._vertex_list[0]._vertex_index[0] <= 1 || original_face._vertex_list[0]._vertex_index[0] >= all_grid.Grid_num_x){
                vertex_velocity_1 = VEC3_TYPE(0.0, 0.0, 0.0);
            }
            else{
                vertex_velocity_1
                    = all_grid.calc_interpolated_velocity(
                        original_face._vertex_list[0]._vertex_pos,
                        interpolation_method);
            }
            // 1番目の頂点の速度を計算
            //壁を含んだセルの速度は0にする
            if(original_face._vertex_list[1]._vertex_index[0] <= 1 || original_face._vertex_list[1]._vertex_index[0] >= all_grid.Grid_num_x){
                vertex_velocity_2 = VEC3_TYPE(0.0, 0.0, 0.0);
            }
            else{
                vertex_velocity_2
                    = all_grid.calc_interpolated_velocity(
                        original_face._vertex_list[1]._vertex_pos,
                        interpolation_method);
            }
        }
        else if(backtrace_usage == "velocity_y") {
            // 0番目の頂点の速度を計算
            //壁を含んだセルの速度は0にする
            if(original_face._vertex_list[0]._vertex_index[1] <= 1 || original_face._vertex_list[0]._vertex_index[1] >= all_grid.Grid_num_y){
                vertex_velocity_1 = VEC3_TYPE(0.0, 0.0, 0.0);
            }
            else{
                vertex_velocity_1
                    = all_grid.calc_interpolated_velocity(
                        original_face._vertex_list[0]._vertex_pos,
                        interpolation_method);
            }
            // 1番目の頂点の速度を計算
            //壁を含んだセルの速度は0にする
            if(original_face._vertex_list[1]._vertex_index[1] <= 1 || original_face._vertex_list[1]._vertex_index[1] >= all_grid.Grid_num_y){
                vertex_velocity_2 = VEC3_TYPE(0.0, 0.0, 0.0);
            }
            else{
                vertex_velocity_2
                    = all_grid.calc_interpolated_velocity(
                        original_face._vertex_list[1]._vertex_pos,
                        interpolation_method);
            }
/*
            vertex_velocity_1
                = all_grid.calc_interpolated_velocity(
                    original_face._vertex_list[0]._vertex_pos,
                    interpolation_method);
            vertex_velocity_2
                = all_grid.calc_interpolated_velocity(
                    original_face._vertex_list[1]._vertex_pos,
                    interpolation_method);
*/
        }
        else{
            vertex_velocity_1
                = all_grid.calc_interpolated_velocity(
                    original_face._vertex_list[0]._vertex_pos,
                    interpolation_method);
            vertex_velocity_2
                = all_grid.calc_interpolated_velocity(
                    original_face._vertex_list[1]._vertex_pos,
                    interpolation_method);
        }

        // バックトレース後の cell face を構成するvertexを定義する
        cell_vertex after_backtrace_vertex_1 = original_face._vertex_list[0];
        cell_vertex after_backtrace_vertex_2 = original_face._vertex_list[1];
        if (time_direction == "forward") {
            after_backtrace_vertex_1._vertex_pos -= vertex_velocity_1 * time_step_length;
            after_backtrace_vertex_2._vertex_pos -= vertex_velocity_2 * time_step_length;
        }
        else if (time_direction == "backward") {
            after_backtrace_vertex_1._vertex_pos += vertex_velocity_1 * time_step_length;
            after_backtrace_vertex_2._vertex_pos += vertex_velocity_2 * time_step_length;
        }
        // バックトレース後の頂点から面を定義
        cell_face after_backtrace_face;
        after_backtrace_face._vertex_list[0] = after_backtrace_vertex_1;
        after_backtrace_face._vertex_list[1] = after_backtrace_vertex_2;
        return after_backtrace_face;
    }


    //引数の面を構成する頂点の位置が、グリッドの範囲を超えていたらclampしてグリッド内に収まるように修正する関数
    cell_face clamp_vertex_position_by_grid(
        const Grid& all_grid,
        cell_face face,
        const std::string backtrace_usage
    ) {
        if(backtrace_usage == "velocity_x"){
            if (face._vertex_list[0]._vertex_pos[0] < all_grid.min_pos_x - 0.5 *all_grid._cell_length) {
                face._vertex_list[0]._vertex_pos[0] = all_grid.min_pos_x - 0.5 *all_grid._cell_length;
            }
            if (face._vertex_list[0]._vertex_pos[0] > all_grid.min_pos_x + all_grid._cell_length * all_grid.Grid_num_x + 0.5 *all_grid._cell_length) {
                face._vertex_list[0]._vertex_pos[0] = all_grid.min_pos_x + all_grid._cell_length * all_grid.Grid_num_x + 0.5 *all_grid._cell_length;
            }
            if (face._vertex_list[1]._vertex_pos[0] < all_grid.min_pos_x - 0.5 *all_grid._cell_length) {
                face._vertex_list[1]._vertex_pos[0] = all_grid.min_pos_x - 0.5 *all_grid._cell_length;
            }
            if (face._vertex_list[1]._vertex_pos[0] > all_grid.min_pos_x + all_grid._cell_length * all_grid.Grid_num_x + 0.5 *all_grid._cell_length) {
                face._vertex_list[1]._vertex_pos[0] = all_grid.min_pos_x + all_grid._cell_length * all_grid.Grid_num_x + 0.5 *all_grid._cell_length;
            }
        }
        else{
            if (face._vertex_list[0]._vertex_pos[0] < all_grid.min_pos_x) {
                face._vertex_list[0]._vertex_pos[0] = all_grid.min_pos_x;
            }
            if (face._vertex_list[0]._vertex_pos[0] > all_grid.min_pos_x + all_grid._cell_length * all_grid.Grid_num_x) {
                face._vertex_list[0]._vertex_pos[0] = all_grid.min_pos_x + all_grid._cell_length * all_grid.Grid_num_x;
            }
            if (face._vertex_list[1]._vertex_pos[0] < all_grid.min_pos_x) {
                face._vertex_list[1]._vertex_pos[0] = all_grid.min_pos_x;
            }
            if (face._vertex_list[1]._vertex_pos[0] > all_grid.min_pos_x + all_grid._cell_length * all_grid.Grid_num_x) {
                face._vertex_list[1]._vertex_pos[0] = all_grid.min_pos_x + all_grid._cell_length * all_grid.Grid_num_x;
            }
        }
        if(backtrace_usage == "velocity_y"){
            if (face._vertex_list[0]._vertex_pos[1] < all_grid.min_pos_y - 0.5 *all_grid._cell_length) {
                face._vertex_list[0]._vertex_pos[1] = all_grid.min_pos_y - 0.5 *all_grid._cell_length;
            }
            if (face._vertex_list[0]._vertex_pos[1] > all_grid.min_pos_y + all_grid._cell_length * all_grid.Grid_num_y + 0.5 *all_grid._cell_length) {
                face._vertex_list[0]._vertex_pos[1] = all_grid.min_pos_y + all_grid._cell_length * all_grid.Grid_num_y + 0.5 *all_grid._cell_length;
            }
            if (face._vertex_list[1]._vertex_pos[1] < all_grid.min_pos_y - 0.5 *all_grid._cell_length) {
                face._vertex_list[1]._vertex_pos[1] = all_grid.min_pos_y - 0.5 *all_grid._cell_length;
            }
            if (face._vertex_list[1]._vertex_pos[1] > all_grid.min_pos_y + all_grid._cell_length * all_grid.Grid_num_y + 0.5 *all_grid._cell_length) {
                face._vertex_list[1]._vertex_pos[1] = all_grid.min_pos_y + all_grid._cell_length * all_grid.Grid_num_y + 0.5 *all_grid._cell_length;
            }
        }
        else{
            if (face._vertex_list[0]._vertex_pos[1] < all_grid.min_pos_y) {
                face._vertex_list[0]._vertex_pos[1] = all_grid.min_pos_y;
            }
            if (face._vertex_list[0]._vertex_pos[1] > all_grid.min_pos_y + all_grid._cell_length * all_grid.Grid_num_y) {
                face._vertex_list[0]._vertex_pos[1] = all_grid.min_pos_y + all_grid._cell_length * all_grid.Grid_num_y;
            }
            if (face._vertex_list[1]._vertex_pos[1] < all_grid.min_pos_y) {
                face._vertex_list[1]._vertex_pos[1] = all_grid.min_pos_y;
            }
            if (face._vertex_list[1]._vertex_pos[1] > all_grid.min_pos_y + all_grid._cell_length * all_grid.Grid_num_y) {
                face._vertex_list[1]._vertex_pos[1] = all_grid.min_pos_y + all_grid._cell_length * all_grid.Grid_num_y;
            }
        }

        return face;
    }

}
