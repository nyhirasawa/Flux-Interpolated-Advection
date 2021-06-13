#include "calc_backtraced_face_3D.h"

#include "grid_3d.h"

namespace smoke_simulation{
    // time_direction == "forward"  の場合は正の時間の方向に backtrace する(普通のsemi-Lagrangian と同じ)
    // time_direction == "backward" の場合は逆の時間の方向に backtrace する( MacCormack の2段階目で使う)
    // use_interpolated_velocity は MacCormack の二回目のバックトレースの時は絶対 true にする
    void calc_backtraced_face_3D(
        const int dim,
        const Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const cell_face_3D original_face,
        cell_face_3D &after_backtrace_face,
        const std::string time_direction,
//        const bool set_zero_normal_velocity_at_boundary,
        const std::string interpolation_method,
        const bool use_interpolated_velocity,
        const Eigen::Vector3i cell_vertex_index_0,
        const Eigen::Vector3i cell_vertex_index_1,
        const Eigen::Vector3i cell_vertex_index_2,
        const Eigen::Vector3i cell_vertex_index_3
    ) {
        //cell vertex の速度を計算する
        std::vector<VEC3_TYPE> vertex_velocity(4);

        //境界の面かどうか
        bool is_boundary_vertex[4];
        for(int i_vert = 0; i_vert < 4; ++i_vert){
            is_boundary_vertex[i_vert] = false;
        }
        if(dim == 0){
            if(original_face._vertex_list[0]._vertex_index[0] <= 1 ||original_face._vertex_list[0]._vertex_index[0] >= all_grid.Grid_num_x){
//            if(cell_vertex_index_0[0] <= 1 || cell_vertex_index_0[0] >= all_grid.Grid_num_x){
                is_boundary_vertex[0] = true;
            }
            if(original_face._vertex_list[1]._vertex_index[0] <= 1 ||original_face._vertex_list[1]._vertex_index[0] >= all_grid.Grid_num_x){
//            if(cell_vertex_index_1[0] <= 1 || cell_vertex_index_1[0] >= all_grid.Grid_num_x){
                is_boundary_vertex[1] = true;
            }
            if(original_face._vertex_list[2]._vertex_index[0] <= 1 ||original_face._vertex_list[2]._vertex_index[0] >= all_grid.Grid_num_x){
//            if(cell_vertex_index_2[0] <= 1 || cell_vertex_index_2[0] >= all_grid.Grid_num_x){
                is_boundary_vertex[2] = true;
            }
            if(original_face._vertex_list[3]._vertex_index[0] <= 1 ||original_face._vertex_list[3]._vertex_index[0] >= all_grid.Grid_num_x){
//            if(cell_vertex_index_3[0] <= 1 || cell_vertex_index_3[0] >= all_grid.Grid_num_x){
                is_boundary_vertex[3] = true;
            }
        }
        else if(dim == 1){
            if(original_face._vertex_list[0]._vertex_index[1] <= 1 ||original_face._vertex_list[0]._vertex_index[1] >= all_grid.Grid_num_y){
//            if(cell_vertex_index_0[1] <= 1 || cell_vertex_index_0[1] >= all_grid.Grid_num_y){
                is_boundary_vertex[0] = true;
            }
            if(original_face._vertex_list[1]._vertex_index[1] <= 1 ||original_face._vertex_list[1]._vertex_index[1] >= all_grid.Grid_num_y){
//            if(cell_vertex_index_1[1] <= 1 || cell_vertex_index_1[1] >= all_grid.Grid_num_y){
                is_boundary_vertex[1] = true;
            }
            if(original_face._vertex_list[2]._vertex_index[1] <= 1 ||original_face._vertex_list[2]._vertex_index[1] >= all_grid.Grid_num_y){
//            if(cell_vertex_index_2[1] <= 1 || cell_vertex_index_2[1] >= all_grid.Grid_num_y){
                is_boundary_vertex[2] = true;
            }
            if(original_face._vertex_list[3]._vertex_index[1] <= 1 ||original_face._vertex_list[3]._vertex_index[1] >= all_grid.Grid_num_y){
//            if(cell_vertex_index_3[1] <= 1 || cell_vertex_index_3[1] >= all_grid.Grid_num_y){
                is_boundary_vertex[3] = true;
            }
        }
        else if(dim == 2){
            if(original_face._vertex_list[0]._vertex_index[2] <= 1 ||original_face._vertex_list[0]._vertex_index[2] >= all_grid.Grid_num_z){
//            if(cell_vertex_index_0[2] <= 1 || cell_vertex_index_0[2] >= all_grid.Grid_num_z){
                is_boundary_vertex[0] = true;
            }
            if(original_face._vertex_list[1]._vertex_index[2] <= 1 ||original_face._vertex_list[1]._vertex_index[2] >= all_grid.Grid_num_z){
//            if(cell_vertex_index_1[2] <= 1 || cell_vertex_index_1[2] >= all_grid.Grid_num_z){
                is_boundary_vertex[1] = true;
            }
            if(original_face._vertex_list[2]._vertex_index[2] <= 1 ||original_face._vertex_list[2]._vertex_index[2] >= all_grid.Grid_num_z){
//            if(cell_vertex_index_2[2] <= 1 || cell_vertex_index_2[2] >= all_grid.Grid_num_z){
                is_boundary_vertex[2] = true;
            }
            if(original_face._vertex_list[3]._vertex_index[2] <= 1 ||original_face._vertex_list[3]._vertex_index[2] >= all_grid.Grid_num_z){
//            if(cell_vertex_index_3[2] <= 1 || cell_vertex_index_3[2] >= all_grid.Grid_num_z){
                is_boundary_vertex[3] = true;
            }
        }

        //補間によって4頂点の速度を求める場合
        if(use_interpolated_velocity){
            std::string interpolation_method_velocity;
            interpolation_method_velocity = interpolation_method;
            for(int i_vert = 0; i_vert < 4; ++i_vert){
                vertex_velocity[i_vert] = all_grid.calc_interpolated_velocity(original_face._vertex_list[i_vert]._vertex_pos, interpolation_method_velocity);
            }
        }
        //周りのセル面から頂点の速度を求める場合
        else{
            vertex_velocity[0] = all_grid.calc_cell_vertex_velocity(cell_vertex_index_0);
            vertex_velocity[1] = all_grid.calc_cell_vertex_velocity(cell_vertex_index_1);
            vertex_velocity[2] = all_grid.calc_cell_vertex_velocity(cell_vertex_index_2);
            vertex_velocity[3] = all_grid.calc_cell_vertex_velocity(cell_vertex_index_3);
        }
        // バックトレース後の cell face を構成するvertexを定義する
        std::vector<cell_vertex_3D> after_backtrace_vertex(4);
        for(int i_vert=0;i_vert<4; ++i_vert){
            after_backtrace_vertex[i_vert] = original_face._vertex_list[i_vert];
        }
        // バックトレースの処理
        if (time_direction == "forward") {
            for(int i_vert=0;i_vert<4; ++i_vert){
                //壁を含むセルは移流しない
                if(!is_boundary_vertex[i_vert]){
                    after_backtrace_vertex[i_vert]._vertex_pos -= vertex_velocity[i_vert] * time_step_length;
                }
            }
        }
        else if (time_direction == "backward") {
            for(int i_vert=0;i_vert<4; ++i_vert){
                //壁を含むセルは移流しない
                if(!is_boundary_vertex[i_vert]){
                    after_backtrace_vertex[i_vert]._vertex_pos += vertex_velocity[i_vert] * time_step_length;
                }
            }
        }
        // バックトレース後の頂点から面を定義
        for(int i_vert = 0; i_vert < 4; ++i_vert){
            after_backtrace_face._vertex_list[i_vert] = after_backtrace_vertex[i_vert];
        }
    }

    //引数の面を構成する頂点の位置が、グリッドの範囲を超えていたらclampしてグリッド内に収まるように修正する関数
    void clamp_vertex_position_by_grid_3D(
        const int dim,
        const Grid_3D& all_grid,
        cell_face_3D& face
    ) {
        const int dx = (dim == 0);
        const int dy = (dim == 1);
        const int dz = (dim == 2);
        for(int i_vert = 0; i_vert < 4; ++i_vert){
            if (face._vertex_list[i_vert]._vertex_pos[0] < all_grid.min_pos_x - dx * 0.5 * all_grid._cell_length) {
                face._vertex_list[i_vert]._vertex_pos[0] = all_grid.min_pos_x - dx * 0.5 * all_grid._cell_length;
            }
            if (face._vertex_list[i_vert]._vertex_pos[0] > all_grid.min_pos_x + all_grid._cell_length * all_grid.Grid_num_x + dx * 0.5 * all_grid._cell_length) {
                face._vertex_list[i_vert]._vertex_pos[0] = all_grid.min_pos_x + all_grid._cell_length * all_grid.Grid_num_x + dx * 0.5 * all_grid._cell_length;
            }
            if (face._vertex_list[i_vert]._vertex_pos[1] < all_grid.min_pos_y - dy * 0.5 * all_grid._cell_length) {
                face._vertex_list[i_vert]._vertex_pos[1] = all_grid.min_pos_y - dy * 0.5 * all_grid._cell_length;
            }
            if (face._vertex_list[i_vert]._vertex_pos[1] > all_grid.min_pos_y + all_grid._cell_length * all_grid.Grid_num_y + dy * 0.5 * all_grid._cell_length) {
                face._vertex_list[i_vert]._vertex_pos[1] = all_grid.min_pos_y + all_grid._cell_length * all_grid.Grid_num_y + dy * 0.5 * all_grid._cell_length;
            }
            if (face._vertex_list[i_vert]._vertex_pos[2] < all_grid.min_pos_z - dz * 0.5 * all_grid._cell_length) {
                face._vertex_list[i_vert]._vertex_pos[2] = all_grid.min_pos_z - dz * 0.5 * all_grid._cell_length;
            }
            if (face._vertex_list[i_vert]._vertex_pos[2] > all_grid.min_pos_z + all_grid._cell_length * all_grid.Grid_num_z + dz * 0.5 * all_grid._cell_length) {
                face._vertex_list[i_vert]._vertex_pos[2] = all_grid.min_pos_z + all_grid._cell_length * all_grid.Grid_num_z + dz * 0.5 * all_grid._cell_length;
            }
        }
    }
}//namespace smoke_simulation
