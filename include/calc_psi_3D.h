#pragma once

#include "physical_const.h"
#include "utils.h"
#include "grid_3d.h"
#include "define_float_type.h"

namespace smoke_simulation {
    void calc_psi_on_cell_face_from_density_on_cell_center_3D(
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &psi_on_cell_face_x,
        std::vector<MY_FLOAT_TYPE> &psi_on_cell_face_y,
        std::vector<MY_FLOAT_TYPE> &psi_on_cell_face_z,
        const std::vector<MY_FLOAT_TYPE> &density_on_cell_center,
        const int Grid_num_x,
        const int Grid_num_y,
        const int Grid_num_z,
        const MY_FLOAT_TYPE cell_length,
        const std::string interpolation_method,
        const int num_gauss_quadrature_point_for_integrate_density,
        const bool use_integral
    );

    // dim = 0(x), 1(y), 2(z)
    void calc_discrete_psi_velocity_3D(
        const int dim,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &psi_velocity,
        const std::vector<MY_FLOAT_TYPE> &velocity,
        const int Grid_num_x,
        const int Grid_num_y,
        const int Grid_num_z,
        const MY_FLOAT_TYPE cell_length,
        const std::string interpolation_method,
        const int num_gauss_quadrature_point_for_integrate_density,
        const bool use_integral
    );
} //namespace smoke_simulation
