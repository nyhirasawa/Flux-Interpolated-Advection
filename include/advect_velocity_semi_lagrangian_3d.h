#pragma once

#include "grid_3d.h"

namespace smoke_simulation{
    //advect項の計算 (semi-Lagrangian)
    void advect_velocity_semi_lagrangian_3D(
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const std::string interpolation_method
    );
    //advect項の計算
    //速度場を時間 -dt だけバックトレースしてadvect項を計算する
    void advect_velocity_semi_lagrangian_x_3D(
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const std::string interpolation_method,
        std::vector<MY_FLOAT_TYPE> &velocity_after_advect_x
    );
    //advect項の計算
    //速度場を時間 -dt だけバックトレースしてadvect項を計算する
    void advect_velocity_semi_lagrangian_y_3D(
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const std::string interpolation_method,
        std::vector<MY_FLOAT_TYPE> &velocity_after_advect_y
    );
    //advect項の計算
    //速度場を時間 -dt だけバックトレースしてadvect項を計算する
    void advect_velocity_semi_lagrangian_z_3D(
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const std::string interpolation_method,
        std::vector<MY_FLOAT_TYPE> &velocity_after_advect_z
    );
} // namespace smoke_simulation
