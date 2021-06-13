#pragma once

#include "cell_face.h"
#include "grid.h"

namespace smoke_simulation{
    // time_direction == "forward"  の場合は正の時間の方向に backtrace する(普通のsemi-Lagrangian と同じ)
    // time_direction == "backward" の場合は逆の時間の方向に backtrace する( MacCormack の2段階目で使う)
    cell_face calc_backtraced_face(
        const Grid& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const cell_face original_face,
        const std::string time_direction,
//        const bool set_zero_normal_velocity_at_boundary,
        const std::string interpolation_method,
        const std::string backtrace_usage
    );

    //引数の面を構成する頂点の位置が、グリッドの範囲を超えていたらclampしてグリッド内に収まるように修正する関数
    cell_face clamp_vertex_position_by_grid(
        const Grid& all_grid,
        cell_face face,
        const std::string backtrace_usage
    );
}
