#pragma once

#include "define_float_type.h"
#include "grid_3d.h"
#include "utils.h"

namespace smoke_simulation{
    // density の MacCormack 法でのclamping の処理
    MY_FLOAT_TYPE clamp_substance_density_of_MacCormack(
        const Grid_3D& all_grid,
        MY_FLOAT_TYPE substance_density_after_advect,
        VEC3_TYPE after_backtrace_position,
        const std::string interpolation_method
    );
    // density の MacCormack 法でのclamping の処理
    MY_FLOAT_TYPE clamp_integral_of_normal_component_of_psi_of_MacCormack_3D(
        const Grid_3D& all_grid,
        MY_FLOAT_TYPE integral_of_normal_component_of_psi_after_advect,
        VEC3_TYPE after_backtrace_position,
        const std::string interpolation_method
    );

}
