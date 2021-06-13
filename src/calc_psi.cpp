#include "calc_psi.h"

#include "linear_interpolation_3d.h"
#include "define_float_type.h"
#include "linear_interpolation_1d.h"
#include "WENO6_interpolation_1d.h"
#include "gauss_quadrature_points.h"
#include "linear_interpolation_2d.h"

namespace smoke_simulation {
    void calc_discrete_psi_velocity_x(
        const Grid &all_grid,
        std::vector<MY_FLOAT_TYPE> &psi_velocity_x_on_cell_vertex_y,
        const std::vector<MY_FLOAT_TYPE> &velocity_x,
        const int Grid_num_x,
        const int Grid_num_y,
        const MY_FLOAT_TYPE cell_length,
        const std::string interpolation_method
    ){
            if(interpolation_method == "const" || interpolation_method == "1Dy_linear"){
                //y成分の計算
                for (int ix = 0; ix < Grid_num_x + 1; ++ix) {
                    psi_velocity_x_on_cell_vertex_y[get_voxel_center_index(ix, 0, Grid_num_x + 1, Grid_num_y + 1)] = 0.0;
                }
                for (int iy = 1; iy < Grid_num_y + 1; ++iy) {
                    for (int ix = 0; ix < Grid_num_x + 1; ++ix) {
                        psi_velocity_x_on_cell_vertex_y[get_voxel_center_index(ix, iy, Grid_num_x + 1, Grid_num_y + 1)]
                            = psi_velocity_x_on_cell_vertex_y[get_voxel_center_index(ix, iy - 1, Grid_num_x + 1, Grid_num_y + 1)]
                            + velocity_x[get_voxel_face_index_x(ix, iy - 1, Grid_num_x, Grid_num_y)] * cell_length;
                    }
                }
            }
            else if(interpolation_method == "linear"){
                //y成分の計算
                for (int ix = 0; ix < Grid_num_x + 1; ++ix) {
                    psi_velocity_x_on_cell_vertex_y[get_voxel_center_index(ix, 0, Grid_num_x + 1, Grid_num_y + 1)] = 0.0;
                }
                for (int iy = 1; iy < Grid_num_y + 1; ++iy) {
                    for (int ix = 0; ix < Grid_num_x + 1; ++ix) {
                        MY_FLOAT_TYPE interpolated_velocity_x_0, interpolated_velocity_x_1;
                        interpolated_velocity_x_0
                            = linear_interpolation_1Dy_cell_face_x_values(
                                VEC3_TYPE(ix * cell_length, ((iy - 1) + 0.25) * cell_length, 0.0),
                                velocity_x,
                                all_grid
                            );
                        interpolated_velocity_x_1
                            = linear_interpolation_1Dy_cell_face_x_values(
                                VEC3_TYPE(ix * cell_length, ((iy - 1) + 0.75) * cell_length, 0.0),
                                velocity_x,
                                all_grid
                            );

                        psi_velocity_x_on_cell_vertex_y[get_voxel_center_index(ix, iy, Grid_num_x + 1, Grid_num_y + 1)]
                            = psi_velocity_x_on_cell_vertex_y[get_voxel_center_index(ix, iy - 1, Grid_num_x + 1, Grid_num_y + 1)]
                            + interpolated_velocity_x_0 * (cell_length / 2.0)
                            + interpolated_velocity_x_1 * (cell_length / 2.0);
                    }
                }
            }
//            else if(interpolation_method == "WENO6" || interpolation_method == "WENO6-optimized"){
//            }
    }

    void calc_discrete_psi_velocity_y(
        const Grid &all_grid,
        std::vector<MY_FLOAT_TYPE> &psi_velocity_cell_center_y,
        const std::vector<MY_FLOAT_TYPE> &velocity_y,
        const int Grid_num_x,
        const int Grid_num_y,
        const MY_FLOAT_TYPE cell_length,
        const std::string interpolation_method
    ){
            if(interpolation_method == "const" || interpolation_method == "1Dy_linear"){
                //y成分の計算
                for (int ix = 0; ix < Grid_num_x; ++ix) {
                    psi_velocity_cell_center_y[get_voxel_center_index(ix, 0, Grid_num_x, Grid_num_y + 2)] = 0.0;
                }
                for (int iy = 1; iy < Grid_num_y + 2; ++iy) {
                    for (int ix = 0; ix < Grid_num_x; ++ix) {
                        psi_velocity_cell_center_y[get_voxel_center_index(ix, iy, Grid_num_x, Grid_num_y + 2)]
                            = psi_velocity_cell_center_y[get_voxel_center_index(ix, iy - 1, Grid_num_x, Grid_num_y + 2)]
                            + velocity_y[get_voxel_face_index_y(ix, iy - 1, Grid_num_x, Grid_num_y)] * cell_length;
                    }
                }
            }
            else if(interpolation_method == "linear"){
                //y成分の計算
                for (int ix = 0; ix < Grid_num_x; ++ix) {
                    psi_velocity_cell_center_y[get_voxel_center_index(ix, 0, Grid_num_x, Grid_num_y + 2)] = 0.0;
                }
                for (int iy = 1; iy < Grid_num_y + 2; ++iy) {
                    for (int ix = 0; ix < Grid_num_x; ++ix) {
                        MY_FLOAT_TYPE interpolated_velocity_y_0, interpolated_velocity_y_1;
                        interpolated_velocity_y_0
                            = linear_interpolation_1Dy_cell_face_y_values(
                                VEC3_TYPE((ix + 0.5) * cell_length, ((iy - 1) - 0.25) * cell_length, 0.0),
                                velocity_y,
                                all_grid
                            );
                        interpolated_velocity_y_1
                            = linear_interpolation_1Dy_cell_face_y_values(
                                VEC3_TYPE((ix + 0.5) * cell_length, ((iy - 1) + 0.25) * cell_length, 0.0),
                                velocity_y,
                                all_grid
                            );
                        psi_velocity_cell_center_y[get_voxel_center_index(ix, iy, Grid_num_x, Grid_num_y + 2)]
                            = psi_velocity_cell_center_y[get_voxel_center_index(ix, iy - 1, Grid_num_x, Grid_num_y + 2)]
                            + interpolated_velocity_y_0 * (cell_length / 2.0)
                            + interpolated_velocity_y_1 * (cell_length / 2.0);
                    }
                }
            }
//            else if(interpolation_method == "WENO6" || interpolation_method == "WENO6-optimized"){
//            }
    }


    void calc_psi_on_cell_face_from_density_on_cell_center(
        const Grid &all_grid,
        std::vector<MY_FLOAT_TYPE> &psi_on_cell_face_x,
        std::vector<MY_FLOAT_TYPE> &psi_on_cell_face_y,
        const std::vector<MY_FLOAT_TYPE> &density_on_cell_center,
        const int Grid_num_x,
        const int Grid_num_y,
        const MY_FLOAT_TYPE cell_length,
        const std::string interpolation_method
    ){
        if (physical_const::kPsi_definition == "x") {
            //x成分の計算
            for (int iy = 0; iy < Grid_num_y; ++iy) {
                psi_on_cell_face_x[get_voxel_face_index_x(0, iy, Grid_num_x, Grid_num_y)] = 0.0;
            }
            for (int ix = 1; ix < Grid_num_x + 1; ++ix) {
                for (int iy = 0; iy < Grid_num_y; ++iy) {
                    psi_on_cell_face_x[get_voxel_face_index_x(ix, iy, Grid_num_x, Grid_num_y)]
                        = psi_on_cell_face_x[get_voxel_face_index_x(ix - 1, iy, Grid_num_x, Grid_num_y)]
                        + density_on_cell_center[get_voxel_center_index(ix - 1, iy, Grid_num_x, Grid_num_y)] * cell_length;
                }
            }
            //y成分の計算
            for (int ix = 0; ix < Grid_num_x; ++ix) {
                psi_on_cell_face_y[get_voxel_face_index_y(ix, 0, Grid_num_x, Grid_num_y)] = 0.0;
            }
            for (int iy = 1; iy < Grid_num_y + 1; ++iy) {
                for (int ix = 0; ix < Grid_num_x; ++ix) {
                    psi_on_cell_face_y[get_voxel_face_index_y(ix, iy, Grid_num_x, Grid_num_y)] = 0.0;
                }
            }
        }
        else if (physical_const::kPsi_definition == "y") {
            if(interpolation_method == "const" || interpolation_method == "1Dy_linear"){
                //x成分の計算
                for (int iy = 0; iy < Grid_num_y; ++iy) {
                    psi_on_cell_face_x[get_voxel_face_index_x(0, iy, Grid_num_x, Grid_num_y)] = 0.0;
                }
                for (int ix = 1; ix < Grid_num_x + 1; ++ix) {
                    for (int iy = 0; iy < Grid_num_y; ++iy) {
                        psi_on_cell_face_x[get_voxel_face_index_x(ix, iy, Grid_num_x, Grid_num_y)] = 0.0;
                    }
                }
                //y成分の計算
                for (int ix = 0; ix < Grid_num_x; ++ix) {
                    psi_on_cell_face_y[get_voxel_face_index_y(ix, 0, Grid_num_x, Grid_num_y)] = 0.0;
                }
                for (int iy = 1; iy < Grid_num_y + 1; ++iy) {
                    for (int ix = 0; ix < Grid_num_x; ++ix) {
                        psi_on_cell_face_y[get_voxel_face_index_y(ix, iy, Grid_num_x, Grid_num_y)]
                            = psi_on_cell_face_y[get_voxel_face_index_y(ix, iy - 1, Grid_num_x, Grid_num_y)]
                            + density_on_cell_center[get_voxel_center_index(ix, iy - 1, Grid_num_x, Grid_num_y)] * cell_length;
                    }
                }
            }
            else if(interpolation_method == "linear"){
                //x成分の計算
                for (int ix = 0; ix < Grid_num_x + 1; ++ix) {
                    for (int iy = 0; iy < Grid_num_y; ++iy) {
                        psi_on_cell_face_x[get_voxel_face_index_x(ix, iy, Grid_num_x, Grid_num_y)] = 0.0;
                    }
                }
                //y成分の計算
                for (int ix = 0; ix < Grid_num_x; ++ix) {
                    psi_on_cell_face_y[get_voxel_face_index_y(ix, 0, Grid_num_x, Grid_num_y)] = 0.0;
                }
                for (int iy = 1; iy < Grid_num_y + 1; ++iy) {
                    for (int ix = 0; ix < Grid_num_x; ++ix) {
                        MY_FLOAT_TYPE interpolated_density_0, interpolated_density_1;
                        interpolated_density_0
                            = all_grid.interpolate_cell_center_defined_values_trilinear(VEC3_TYPE((ix + 0.5) * cell_length, ((iy - 1) + 0.25) * cell_length, 0.0), density_on_cell_center);
                        interpolated_density_1
                            = all_grid.interpolate_cell_center_defined_values_trilinear(VEC3_TYPE((ix + 0.5) * cell_length, ((iy - 1) + 0.75) * cell_length, 0.0), density_on_cell_center);
                        psi_on_cell_face_y[get_voxel_face_index_y(ix, iy, Grid_num_x, Grid_num_y)]
                            = psi_on_cell_face_y[get_voxel_face_index_y(ix, iy - 1, Grid_num_x, Grid_num_y)]
                            + interpolated_density_0 * (cell_length / 2.0)
                            + interpolated_density_1 * (cell_length / 2.0);
                    }
                }
            }
//            else if(interpolation_method == "WENO6" || interpolation_method == "WENO6-optimized"){
//            }
        }
        else if (physical_const::kPsi_definition == "xy") {
            //x成分の計算
            for (int iy = 0; iy < Grid_num_y; ++iy) {
                psi_on_cell_face_x[get_voxel_face_index_x(0, iy, Grid_num_x, Grid_num_y)] = 0.0;
            }
            for (int ix = 1; ix < Grid_num_x + 1; ++ix) {
                for (int iy = 0; iy < Grid_num_y; ++iy) {
                    psi_on_cell_face_x[get_voxel_face_index_x(ix, iy, Grid_num_x, Grid_num_y)]
                        = psi_on_cell_face_x[get_voxel_face_index_x(ix - 1, iy, Grid_num_x, Grid_num_y)]
                        + density_on_cell_center[get_voxel_center_index(ix - 1, iy, Grid_num_x, Grid_num_y)] * cell_length;
                }
            }
            //y成分の計算
            for (int ix = 0; ix < Grid_num_x; ++ix) {
                psi_on_cell_face_y[get_voxel_face_index_y(ix, 0, Grid_num_x, Grid_num_y)] = 0.0;
            }
            for (int iy = 1; iy < Grid_num_y + 1; ++iy) {
                for (int ix = 0; ix < Grid_num_x; ++ix) {
                    psi_on_cell_face_y[get_voxel_face_index_y(ix, iy, Grid_num_x, Grid_num_y)]
                        = psi_on_cell_face_y[get_voxel_face_index_y(ix, iy - 1, Grid_num_x, Grid_num_y)]
                        + density_on_cell_center[get_voxel_center_index(ix, iy - 1, Grid_num_x, Grid_num_y)] * cell_length;
                }
            }
        }
    }

}//namespace smoke_simulation
