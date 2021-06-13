#ifndef MOVING_SUBSTANCES_1D_H
#define MOVING_SUBSTANCES_1D_H

#include <vector>
#include "grid_1d.h"

namespace smoke_simulation {
    void move_substances_1D(
        Grid_1D& all_grid,
        int i_frame,
        const MY_FLOAT_TYPE time_step_length,
        const bool use_flux_advection,
//		const bool set_zero_normal_velocity_at_boundary,
        const bool use_MacCormack_scheme,
        const bool use_clamping_in_MacCormack_scheme,
        const int num_gauss_quad_boundary,
        const int num_gauss_quad_bulk,
        const bool enable_eliminate_zero_velocity_false_diffusion,
        const std::string split_method,
        const std::string integral_method,
        const std::string interpolation_method,
        const bool use_integral,
        const int num_gauss_quadrature_point_for_integrate_density,
        const bool enable_cell_volume_correction,
        const bool use_lentines_advection
    );

    //密度場の移流項を flux advection で計算
	void advect_density_flux_advection_1D(
		Grid_1D& all_grid,
		const MY_FLOAT_TYPE time_step_length,
		const int num_gauss_quad_boundary,
		const int num_gauss_quad_bulk,
		const bool enable_eliminate_zero_velocity_false_diffusion,
		const std::string split_method,
		const std::string integral_method,
		std::vector<MY_FLOAT_TYPE> &advected_values,
		const std::string interpolation_method,
	    const std::string interpolation_method_in_calc_psi,
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density
	);
}

#endif // MOVING_SUBSTANCES_1D_H
