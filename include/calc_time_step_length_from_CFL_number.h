#pragma once

#include "grid.h"
#include "grid_3d.h"
#include "grid_1d.h"
#include "physical_const.h"
#include "utils.h"

#include "define_float_type.h"

namespace smoke_simulation {
	MY_FLOAT_TYPE calc_time_step_length_from_CFL_number_1D(
		const Grid_1D& all_grid,
		const MY_FLOAT_TYPE now_time,
		const int i_movie_frame,
		bool& write_frame_to_movie_file,
		const int CFL_number,
		const MY_FLOAT_TYPE movie_fps
	);
	MY_FLOAT_TYPE calc_time_step_length_from_CFL_number(
		const Grid& all_grid,
		const MY_FLOAT_TYPE now_time,
		const int i_movie_frame,
		bool& write_frame_to_movie_file,
		const int CFL_number,
		const MY_FLOAT_TYPE movie_fps
	);
	MY_FLOAT_TYPE calc_time_step_length_from_fixed_time_step_length(
		const MY_FLOAT_TYPE now_time,
		const int i_movie_frame,
		bool& write_frame_to_movie_file,
		const MY_FLOAT_TYPE movie_fps);
	MY_FLOAT_TYPE calc_time_step_length_from_CFL_number_3D(
		const Grid_3D& all_grid,
		const MY_FLOAT_TYPE now_time,
		const int i_movie_frame,
		bool& write_frame_to_movie_file,
		const int CFL_number,
		const MY_FLOAT_TYPE movie_fps);
}//namespace smoke_simulation
