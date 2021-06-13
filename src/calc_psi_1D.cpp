#include "calc_psi_1D.h"

#include "linear_interpolation_1d.h"

namespace smoke_simulation {

    void calc_psi_on_cell_face_from_density_on_cell_center_1D(
        const Grid_1D &all_grid,
        std::vector<MY_FLOAT_TYPE> &psi_on_cell_face_y,
        const std::vector<MY_FLOAT_TYPE> &density_on_cell_center,
        const int Grid_num_y,
        const MY_FLOAT_TYPE cell_length,
        const std::string interpolation_method
    ) {
//        if(interpolation_method == "const" || interpolation_method == "1Dy_linear"){
            //y成分の計算
            psi_on_cell_face_y[0] = 0.0;

            for (int iy = 1; iy < Grid_num_y + 1; ++iy) {
                psi_on_cell_face_y[iy]
                    = psi_on_cell_face_y[iy - 1]
                    + density_on_cell_center[iy - 1] * cell_length;
            }
//        }
/*
        else if(interpolation_method == "linear"){
            //y成分の計算
            psi_on_cell_face_y[0] = 0.0;
            for (int iy = 1; iy < Grid_num_y + 1; ++iy) {
                MY_FLOAT_TYPE interpolated_density_0, interpolated_density_1;
                interpolated_density_0
                    = linear_interpolation_1D_cell_center_values(
                        ((iy - 1) + 0.25) * cell_length,
                        density_on_cell_center,
                        all_grid
                    );
                interpolated_density_1
                    = linear_interpolation_1D_cell_center_values(
                        ((iy - 1) + 0.75) * cell_length,
                        density_on_cell_center,
                        all_grid
                    );
                psi_on_cell_face_y[iy]
                    = psi_on_cell_face_y[iy - 1]
                    + interpolated_density_0 * (cell_length / 2.0)
                    + interpolated_density_1 * (cell_length / 2.0);
            }
        }
*/
    }
} //namespace smoke_simulation
