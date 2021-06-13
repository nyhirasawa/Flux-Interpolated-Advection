#pragma once

#include "grid_3d.h"
#include "cell_face_3D.h"

namespace smoke_simulation{
    // use_interpolated_velocity は MacCormack の二回目のバックトレースの時は絶対 true にする
    void calc_backtraced_face_3D(
        const int dim,
        const Grid_3D& all_grid,
    	const MY_FLOAT_TYPE time_step_length,
    	const cell_face_3D original_face,
        cell_face_3D &after_backtrace_face,
    	const std::string time_direction,
//    	const bool set_zero_normal_velocity_at_boundary,
        const std::string interpolation_method,
        const bool use_interpolated_velocity,
        const Eigen::Vector3i cell_vertex_index_0,
        const Eigen::Vector3i cell_vertex_index_1,
        const Eigen::Vector3i cell_vertex_index_2,
        const Eigen::Vector3i cell_vertex_index_3
    );

    //引数の面を構成する頂点の位置が、グリッドの範囲を超えていたらclampしてグリッド内に収まるように修正する関数
    void clamp_vertex_position_by_grid_3D(
        const int dim,
        const Grid_3D& all_grid,
        cell_face_3D& face
    );
}//namespace smoke_simulation
