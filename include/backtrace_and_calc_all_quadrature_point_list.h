#pragma once

#include "gauss_quadrature_points.h"
#include "grid_3d.h"
#include "split_face_3D.h"

namespace smoke_simulation{
    //バックトレースを行った後求積点のリストを作る関数
    void backtrace_and_calc_all_quadrature_point_list(
        const int dim,
        quadrature_point_vector &quadrature_point_list,
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
		const MY_FLOAT_TYPE minimal_area
    );
}//namespace smoke_simulation
