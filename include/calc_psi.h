#pragma once

#include "physical_const.h"
#include "utils.h"
#include "grid.h"
#include "define_float_type.h"

namespace smoke_simulation {
    void calc_discrete_psi_velocity_x(
        const Grid &all_grid,
        std::vector<MY_FLOAT_TYPE> &psi_velocity_x_on_cell_vertex_y,
        const std::vector<MY_FLOAT_TYPE> &velocity_x,
        const int Grid_num_x,
        const int Grid_num_y,
        const MY_FLOAT_TYPE cell_length,
        const std::string interpolation_method
    );

    void calc_discrete_psi_velocity_y(
        const Grid &all_grid,
        std::vector<MY_FLOAT_TYPE> &psi_velocity_cell_center_y,
        const std::vector<MY_FLOAT_TYPE> &velocity_y,
        const int Grid_num_x,
        const int Grid_num_y,
        const MY_FLOAT_TYPE cell_length,
        const std::string interpolation_method
    );

    void calc_psi_on_cell_face_from_density_on_cell_center(
        const Grid &all_grid,
        std::vector<MY_FLOAT_TYPE> &psi_on_cell_face_x,
        std::vector<MY_FLOAT_TYPE> &psi_on_cell_face_y,
        const std::vector<MY_FLOAT_TYPE> &density_on_cell_center,
        const int Grid_num_x,
        const int Grid_num_y,
        const MY_FLOAT_TYPE cell_length,
        const std::string interpolation_method
    );
}//namespace smoke_simulation
