#pragma once

#include "grid.h"

namespace smoke_simulation{
    //Poisson eq. を解くことによる圧力の計算
	void calc_pressure(
		Grid& all_grid,
		const bool enable_cell_volume_correction,
		const bool fix_velocity,
		const MY_FLOAT_TYPE time_step_length
	);
	//pressure gradient termの計算
	void calc_pressure_gradient_term(
		Grid& all_grid,
		const bool enable_cell_volume_correction,
		const bool fix_velocity,
		const MY_FLOAT_TYPE time_step_length
	);
}
