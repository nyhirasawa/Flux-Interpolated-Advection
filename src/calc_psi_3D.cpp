#include "calc_psi_3D.h"

#include "linear_interpolation_3d.h"
#include "define_float_type.h"
#include "linear_interpolation_1d.h"
#include "WENO6_interpolation_1d.h"
#include "gauss_quadrature_points.h"

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
    ){
        if (physical_const::kPsi_definition == "x") {
            //x成分の計算
            for (int iy = 0; iy < Grid_num_y; ++iy) {
                for (int iz = 0; iz < Grid_num_z; ++iz) {
                    psi_on_cell_face_x[get_voxel_face_index_x_3D(0, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                }
            }
            for (int ix = 1; ix < Grid_num_x + 1; ++ix) {
                for (int iy = 0; iy < Grid_num_y; ++iy) {
                    for (int iz = 0; iz < Grid_num_z; ++iz) {
                        psi_on_cell_face_x[get_voxel_face_index_x_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                        = psi_on_cell_face_x[get_voxel_face_index_x_3D(ix - 1, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                            + density_on_cell_center[get_voxel_center_index_3D(ix - 1, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)] * cell_length;
                    }
                }
            }
            //y成分の計算
            for (int ix = 0; ix < Grid_num_x; ++ix) {
                for (int iy = 0; iy < Grid_num_y + 1; ++iy) {
                    for (int iz = 0; iz < Grid_num_z; ++iz) {
                        psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                    }
                }
            }
            // z成分の計算
            for (int iy = 0; iy < Grid_num_y; ++iy) {
                for (int ix = 0; ix < Grid_num_x; ++ix) {
                    for (int iz = 0; iz < Grid_num_z + 1; ++iz) {
                        psi_on_cell_face_z[get_voxel_face_index_z_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                    }
                }
            }
        }
        else if (physical_const::kPsi_definition == "y") {
            if(use_integral){
                if(interpolation_method == "const" || interpolation_method == "1Dy_linear"){
                    //x成分の計算
                    for (int ix = 0; ix < Grid_num_x + 1; ++ix) {
                        for (int iy = 0; iy < Grid_num_y; ++iy) {
                            for (int iz = 0; iz < Grid_num_z; ++iz) {
                                psi_on_cell_face_x[get_voxel_face_index_x_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                            }
                        }
                    }
                    //y成分の計算
                    for (int ix = 0; ix < Grid_num_x; ++ix) {
                        for (int iz = 0; iz < Grid_num_z; ++iz) {
                            psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, 0, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                        }
                    }
                    for (int iy = 1; iy < Grid_num_y + 1; ++iy) {
                        for (int ix = 0; ix < Grid_num_x; ++ix) {
                            for (int iz = 0; iz < Grid_num_z; ++iz) {
                                psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                                    = psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                                    + density_on_cell_center[get_voxel_center_index_3D(ix, iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)] * cell_length;
                            }
                        }
                    }
                    //z成分の計算
                    for (int ix = 0; ix < Grid_num_x; ++ix) {
                        for (int iy = 0; iy < Grid_num_y; ++iy) {
                            for (int iz = 0; iz < Grid_num_z + 1; ++iz) {
                                psi_on_cell_face_z[get_voxel_face_index_z_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                            }
                        }
                    }
                }
                else if(interpolation_method == "linear"){
                    //x成分の計算
                    for (int ix = 0; ix < Grid_num_x + 1; ++ix) {
                        for (int iy = 0; iy < Grid_num_y; ++iy) {
                            for (int iz = 0; iz < Grid_num_z; ++iz) {
                                psi_on_cell_face_x[get_voxel_face_index_x_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                            }
                        }
                    }
                    //y成分の計算
                    for (int ix = 0; ix < Grid_num_x; ++ix) {
                        for (int iz = 0; iz < Grid_num_z; ++iz) {
                            psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, 0, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                        }
                    }
                    for (int iy = 1; iy < Grid_num_y + 1; ++iy) {
                        for (int ix = 0; ix < Grid_num_x; ++ix) {
                            for (int iz = 0; iz < Grid_num_z; ++iz) {

                                MY_FLOAT_TYPE interpolated_density_0, interpolated_density_1;
                                   interpolated_density_0
                                       = linear_interpolation_3D_cell_center_values(
                                           VEC3_TYPE((ix + 0.5) * cell_length, ((iy - 1) + 0.25) * cell_length, (iz + 0.5) * cell_length),
                                           density_on_cell_center,
                                           all_grid
                                       );
    //                                   = all_grid.interpolate_cell_center_defined_values_trilinear(VEC3_TYPE((ix + 0.5) * cell_length, ((iy - 1) + 0.25) * cell_length, (iz + 0.5) * cell_length), density_on_cell_center);
                                   interpolated_density_1
                                       = linear_interpolation_3D_cell_center_values(
                                           VEC3_TYPE((ix + 0.5) * cell_length, ((iy - 1) + 0.75) * cell_length, (iz + 0.5) * cell_length),
                                           density_on_cell_center,
                                           all_grid
                                       );
    //                                   = all_grid.interpolate_cell_center_defined_values_trilinear(VEC3_TYPE((ix + 0.5) * cell_length, ((iy - 1) + 0.75) * cell_length, (iz + 0.5) * cell_length), density_on_cell_center);

                                psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                                    = psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                                    + interpolated_density_0 * (cell_length / 2.0)
                                    + interpolated_density_1 * (cell_length / 2.0);

    /*
                                MY_FLOAT_TYPE interpolated_density_0, interpolated_density_1,
                                       interpolated_density_2, interpolated_density_3,
                                       interpolated_density_4, interpolated_density_5,
                                       interpolated_density_6, interpolated_density_7;
                                   interpolated_density_0
                                       = all_grid.interpolate_cell_center_defined_values_trilinear(VEC3_TYPE((ix + 0.25) * cell_length, ((iy - 1) + 0.25) * cell_length, (iz + 0.25) * cell_length), density_on_cell_center);
                                   interpolated_density_1
                                       = all_grid.interpolate_cell_center_defined_values_trilinear(VEC3_TYPE((ix + 0.25) * cell_length, ((iy - 1) + 0.75) * cell_length, (iz + 0.25) * cell_length), density_on_cell_center);
                                   interpolated_density_2
                                       = all_grid.interpolate_cell_center_defined_values_trilinear(VEC3_TYPE((ix + 0.75) * cell_length, ((iy - 1) + 0.25) * cell_length, (iz + 0.25) * cell_length), density_on_cell_center);
                                   interpolated_density_3
                                       = all_grid.interpolate_cell_center_defined_values_trilinear(VEC3_TYPE((ix + 0.75) * cell_length, ((iy - 1) + 0.75) * cell_length, (iz + 0.25) * cell_length), density_on_cell_center);
                                   interpolated_density_4
                                       = all_grid.interpolate_cell_center_defined_values_trilinear(VEC3_TYPE((ix + 0.25) * cell_length, ((iy - 1) + 0.25) * cell_length, (iz + 0.75) * cell_length), density_on_cell_center);
                                   interpolated_density_5
                                       = all_grid.interpolate_cell_center_defined_values_trilinear(VEC3_TYPE((ix + 0.25) * cell_length, ((iy - 1) + 0.75) * cell_length, (iz + 0.75) * cell_length), density_on_cell_center);
                                   interpolated_density_6
                                       = all_grid.interpolate_cell_center_defined_values_trilinear(VEC3_TYPE((ix + 0.75) * cell_length, ((iy - 1) + 0.25) * cell_length, (iz + 0.75) * cell_length), density_on_cell_center);
                                   interpolated_density_7
                                       = all_grid.interpolate_cell_center_defined_values_trilinear(VEC3_TYPE((ix + 0.75) * cell_length, ((iy - 1) + 0.75) * cell_length, (iz + 0.75) * cell_length), density_on_cell_center);

                                psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                                    = psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                                    + interpolated_density_0 * (cell_length / (2.0 * 4.0))
                                    + interpolated_density_1 * (cell_length / (2.0 * 4.0))
                                    + interpolated_density_2 * (cell_length / (2.0 * 4.0))
                                    + interpolated_density_3 * (cell_length / (2.0 * 4.0))
                                    + interpolated_density_4 * (cell_length / (2.0 * 4.0))
                                    + interpolated_density_5 * (cell_length / (2.0 * 4.0))
                                    + interpolated_density_6 * (cell_length / (2.0 * 4.0))
                                    + interpolated_density_7 * (cell_length / (2.0 * 4.0));
    */
                            }
                        }
                    }
                    //z成分の計算
                    for (int ix = 0; ix < Grid_num_x; ++ix) {
                        for (int iy = 0; iy < Grid_num_y; ++iy) {
                            for (int iz = 0; iz < Grid_num_z + 1; ++iz) {
                                psi_on_cell_face_z[get_voxel_face_index_z_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                            }
                        }
                    }
                }
                else if(interpolation_method == "WENO6" || interpolation_method == "WENO6-optimized" || interpolation_method == "1Dy_WENO6"){
                    //x成分の計算
                    for (int ix = 0; ix < Grid_num_x + 1; ++ix) {
                        for (int iy = 0; iy < Grid_num_y; ++iy) {
                            for (int iz = 0; iz < Grid_num_z; ++iz) {
                                psi_on_cell_face_x[get_voxel_face_index_x_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                            }
                        }
                    }
                    //y成分の計算
                    for (int ix = 0; ix < Grid_num_x; ++ix) {
                        for (int iz = 0; iz < Grid_num_z; ++iz) {
                            psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, 0, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                        }
                    }
                    for (int iy = 1; iy < Grid_num_y + 1; ++iy) {
                        for (int ix = 0; ix < Grid_num_x; ++ix) {
                            for (int iz = 0; iz < Grid_num_z; ++iz) {
                                psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                                    = psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)];
                                //求積点の位置
                                std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral;
                                //求積の重み
                                std::vector<MY_FLOAT_TYPE> quadrature_weight_list_1D_density_integral;
                                quadrature_weight_list_1D_density_integral
                                    = gauss_quadrature_points_1D::get_quadtarure_weights_1D(
                                		num_gauss_quadrature_point_for_integrate_density
                                    );
                                //// positionがセルの中心より上にある場合はセル中心より下の1点と上の1点の
                                //// 2つの求積点でセルでの積分ができる
                                //セル中心より下の点
                                //求積点の位置を計算
                                quadrature_position_list_1D_density_integral
                                    = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                                        (iy - 1.0) * all_grid._cell_length,
                                		(iy - 0.5) * all_grid._cell_length,
                                		num_gauss_quadrature_point_for_integrate_density
                                    );
                                //積分計算
                                for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                                    //求積点の位置
                                    VEC3_TYPE quadrature_point_position(
                                        (ix + 0.5) * all_grid._cell_length,
                                        quadrature_position_list_1D_density_integral[i_point][0],
                                        (iz + 0.5) * all_grid._cell_length
                                    );
                                    //求積点での密度
                                    MY_FLOAT_TYPE denisty_at_quadrature_point
                                        = WENO6_interpolation_1Dy_cell_center_values(
                                            quadrature_point_position,
                                            density_on_cell_center,
                                            all_grid
                                        );
                                    psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                                        += denisty_at_quadrature_point
                                        * (all_grid._cell_length / 2.0)
                                        * quadrature_weight_list_1D_density_integral[i_point];
                                }
                                //セル中心より上の点
                                //求積点の位置を計算
                                quadrature_position_list_1D_density_integral
                                    = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                                        (iy - 0.5) * all_grid._cell_length,
                                        (iy      ) * all_grid._cell_length,
                                		num_gauss_quadrature_point_for_integrate_density
                                );
                                //積分計算
                                for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                                    VEC3_TYPE quadrature_point_position(
                                        (ix + 0.5) * all_grid._cell_length,
                                        quadrature_position_list_1D_density_integral[i_point][0],
                                        (iz + 0.5) * all_grid._cell_length
                                    );
                                    MY_FLOAT_TYPE denisty_at_quadrature_point
                                        = WENO6_interpolation_1Dy_cell_center_values(
                                            quadrature_point_position,
                                            density_on_cell_center,
                                            all_grid
                                        );
                                    psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                                        += denisty_at_quadrature_point
                                        * (all_grid._cell_length / 2.0)
                                        * quadrature_weight_list_1D_density_integral[i_point];
                                }
                            }
                        }
                    }
                    //z成分の計算
                    for (int ix = 0; ix < Grid_num_x; ++ix) {
                        for (int iy = 0; iy < Grid_num_y; ++iy) {
                            for (int iz = 0; iz < Grid_num_z + 1; ++iz) {
                                psi_on_cell_face_z[get_voxel_face_index_z_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                            }
                        }
                    }
                }
            }
            else{
                //x成分の計算
                for (int ix = 0; ix < Grid_num_x + 1; ++ix) {
                    for (int iy = 0; iy < Grid_num_y; ++iy) {
                        for (int iz = 0; iz < Grid_num_z; ++iz) {
                            psi_on_cell_face_x[get_voxel_face_index_x_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                        }
                    }
                }
                //y成分の計算
                for (int ix = 0; ix < Grid_num_x; ++ix) {
                    for (int iz = 0; iz < Grid_num_z; ++iz) {
                        psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, 0, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                    }
                }
                for (int iy = 1; iy < Grid_num_y + 1; ++iy) {
                    for (int ix = 0; ix < Grid_num_x; ++ix) {
                        for (int iz = 0; iz < Grid_num_z; ++iz) {
                            psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                                = psi_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                                + density_on_cell_center[get_voxel_center_index_3D(ix, iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)] * cell_length;
                        }
                    }
                }
                //z成分の計算
                for (int ix = 0; ix < Grid_num_x; ++ix) {
                    for (int iy = 0; iy < Grid_num_y; ++iy) {
                        for (int iz = 0; iz < Grid_num_z + 1; ++iz) {
                            psi_on_cell_face_z[get_voxel_face_index_z_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                        }
                    }
                }
            }
        }
    }

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
    ) {
        const int dx = (dim == 0);
        const int dy = (dim == 1);
        const int dz = (dim == 2);
        if(use_integral){
            if(interpolation_method == "const" || interpolation_method == "1Dy_linear"){
                //y成分の計算
                for (int ix = 0; ix < Grid_num_x + dx; ++ix) {
                    for (int iz = 0; iz < Grid_num_z + dz; ++iz) {
                        psi_velocity[get_voxel_center_index_3D(ix, 0, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)] = 0.0;
                    }
                }
    /*
                if(dim == 0){
                    MY_FLOAT_TYPE increment;
                    for (int iy = 1; iy < Grid_num_y + 1 + dy; ++iy) {
                        for (int ix = 0; ix < Grid_num_x + dx; ++ix) {
                            for (int iz = 0; iz < Grid_num_z + dz; ++iz) {
                                if(ix == 0 || ix == Grid_num_x + dx - 1){
                                    increment = 0;
                                }
                                else{
                                    increment = psi_velocity[get_voxel_center_index_3D(ix, iy - 1, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)]
                                              + velocity[get_voxel_center_index_3D(ix, iy - 1, iz, Grid_num_x + dx, Grid_num_y + dy, Grid_num_z + dz)] * cell_length;
                                }
                                psi_velocity[get_voxel_center_index_3D(ix, iy, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)]
                                    = increment;
                            }
                        }
                    }
                }
                else if(dim == 1){
                    MY_FLOAT_TYPE increment;
                    for (int iy = 1; iy < Grid_num_y + 1 + dy; ++iy) {
                        for (int ix = 0; ix < Grid_num_x + dx; ++ix) {
                            for (int iz = 0; iz < Grid_num_z + dz; ++iz) {
                                if(iy == 0 || iy == Grid_num_y + dy){
                                    increment = 0;
                                }
                                else{
                                    increment = psi_velocity[get_voxel_center_index_3D(ix, iy - 1, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)]
                                              + velocity[get_voxel_center_index_3D(ix, iy - 1, iz, Grid_num_x + dx, Grid_num_y + dy, Grid_num_z + dz)] * cell_length;
                                }
                                psi_velocity[get_voxel_center_index_3D(ix, iy, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)]
                                    = increment;
                            }
                        }
                    }
                }
                else if(dim == 2){
                    MY_FLOAT_TYPE increment;
                    for (int iy = 1; iy < Grid_num_y + 1 + dy; ++iy) {
                        for (int ix = 0; ix < Grid_num_x + dx; ++ix) {
                            for (int iz = 0; iz < Grid_num_z + dz; ++iz) {
                                if(iz == 0 || iz == Grid_num_z + dz - 1){
                                    increment = 0;
                                }
                                else{
                                    increment = psi_velocity[get_voxel_center_index_3D(ix, iy - 1, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)]
                                              + velocity[get_voxel_center_index_3D(ix, iy - 1, iz, Grid_num_x + dx, Grid_num_y + dy, Grid_num_z + dz)] * cell_length;
                                }
                                psi_velocity[get_voxel_center_index_3D(ix, iy, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)]
                                    = increment;
                            }
                        }
                    }
                }
                else{
    */
                    for (int iy = 1; iy < Grid_num_y + 1 + dy; ++iy) {
                        for (int ix = 0; ix < Grid_num_x + dx; ++ix) {
                            for (int iz = 0; iz < Grid_num_z + dz; ++iz) {
                                psi_velocity[get_voxel_center_index_3D(ix, iy, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)]
                                    = psi_velocity[get_voxel_center_index_3D(ix, iy - 1, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)]
                                    + velocity[get_voxel_center_index_3D(ix, iy - 1, iz, Grid_num_x + dx, Grid_num_y + dy, Grid_num_z + dz)] * cell_length;
                            }
                        }
                    }
    //            }
            }
            else if(interpolation_method == "linear"){
                //y成分の計算
                for (int ix = 0; ix < Grid_num_x + dx; ++ix) {
                    for (int iz = 0; iz < Grid_num_z + dz; ++iz) {
                        psi_velocity[get_voxel_center_index_3D(ix, 0, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)] = 0.0;
                    }
                }
                for (int iy = 1; iy < Grid_num_y + 1 + dy; ++iy) {
                    for (int ix = 0; ix < Grid_num_x + dx; ++ix) {
                        for (int iz = 0; iz < Grid_num_z + dz; ++iz) {
                            MY_FLOAT_TYPE interpolated_density_0, interpolated_density_1;
                               interpolated_density_0
                                   = linear_interpolation_3D_cell_face_values(
                                       dim,
                                       VEC3_TYPE((ix + 0.5 - 0.5 * dx) * cell_length, ((iy - 1) + 0.25 - 0.5 * dy) * cell_length, (iz + 0.5 - 0.5 * dz) * cell_length),
                                       velocity,
                                       all_grid
                                   );
                               interpolated_density_1
                                   = linear_interpolation_3D_cell_face_values(
                                       dim,
                                       VEC3_TYPE((ix + 0.5 - 0.5 * dx) * cell_length, ((iy - 1) + 0.75 - 0.5 * dy) * cell_length, (iz + 0.5 - 0.5 * dz) * cell_length),
                                       velocity,
                                       all_grid
                                   );
                            psi_velocity[get_voxel_center_index_3D(ix, iy, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)]
                                = psi_velocity[get_voxel_center_index_3D(ix, iy - 1, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)]
                                + interpolated_density_0 * (cell_length / 2.0)
                                + interpolated_density_1 * (cell_length / 2.0);
                        }
                    }
                }
            }
    /*
            else if(interpolation_method == "WENO6" || interpolation_method == "WENO6-optimized"){
                //y成分の計算
                for (int ix = 0; ix < Grid_num_x; ++ix) {
                    for (int iz = 0; iz < Grid_num_z; ++iz) {
                        psi_velocity_x_on_cell_face_y[get_voxel_face_index_y_3D(ix, 0, iz, Grid_num_x, Grid_num_y, Grid_num_z)] = 0.0;
                    }
                }
                for (int iy = 1; iy < Grid_num_y + 1; ++iy) {
                    for (int ix = 0; ix < Grid_num_x; ++ix) {
                        for (int iz = 0; iz < Grid_num_z; ++iz) {
                            psi_velocity_x_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                                = psi_velocity_x_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)];
                            //求積点の位置
                            std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral;
                            //求積の重み
                            std::vector<MY_FLOAT_TYPE> quadrature_weight_list_1D_density_integral;
                            quadrature_weight_list_1D_density_integral
                                = gauss_quadrature_points_1D::get_quadtarure_weights_1D(
                                    num_gauss_quadrature_point_for_integrate_density
                                );
                            //// positionがセルの中心より上にある場合はセル中心より下の1点と上の1点の
                            //// 2つの求積点でセルでの積分ができる
                            //セル中心より下の点
                            //求積点の位置を計算
                            quadrature_position_list_1D_density_integral
                                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                                    (iy - 1.0) * all_grid._cell_length,
                                    (iy - 0.5) * all_grid._cell_length,
                                    num_gauss_quadrature_point_for_integrate_density
                                );
                            //積分計算
                            for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                                //求積点の位置
                                VEC3_TYPE quadrature_point_position(
                                    (ix + 0.5) * all_grid._cell_length,
                                    quadrature_position_list_1D_density_integral[i_point][0],
                                    (iz + 0.5) * all_grid._cell_length
                                );
                                //求積点での密度
                                MY_FLOAT_TYPE denisty_at_quadrature_point
                                    = WENO6_interpolation_1Dy_cell_center_values(
                                        quadrature_point_position,
                                        velocity_x,
                                        all_grid
                                    );
                                psi_velocity_x_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                                    += denisty_at_quadrature_point
                                    * (all_grid._cell_length / 2.0)
                                    * quadrature_weight_list_1D_density_integral[i_point];
                            }
                            //セル中心より上の点
                            //求積点の位置を計算
                            quadrature_position_list_1D_density_integral
                                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                                    (iy - 0.5) * all_grid._cell_length,
                                    (iy      ) * all_grid._cell_length,
                                    num_gauss_quadrature_point_for_integrate_density
                            );
                            //積分計算
                            for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                                VEC3_TYPE quadrature_point_position(
                                    (ix + 0.5) * all_grid._cell_length,
                                    quadrature_position_list_1D_density_integral[i_point][0],
                                    (iz + 0.5) * all_grid._cell_length
                                );
                                MY_FLOAT_TYPE denisty_at_quadrature_point
                                    = WENO6_interpolation_1Dy_cell_center_values(
                                        quadrature_point_position,
                                        velocity_x,
                                        all_grid
                                    );
                                psi_velocity_x_on_cell_face_y[get_voxel_face_index_y_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                                    += denisty_at_quadrature_point
                                    * (all_grid._cell_length / 2.0)
                                    * quadrature_weight_list_1D_density_integral[i_point];
                            }
                        }
                    }
                }
            }
        }
    */
        }
        else{
            //y成分の計算
            for (int ix = 0; ix < Grid_num_x + dx; ++ix) {
                for (int iz = 0; iz < Grid_num_z + dz; ++iz) {
                    psi_velocity[get_voxel_center_index_3D(ix, 0, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)] = 0.0;
                }
            }
            for (int iy = 1; iy < Grid_num_y + 1 + dy; ++iy) {
                for (int ix = 0; ix < Grid_num_x + dx; ++ix) {
                    for (int iz = 0; iz < Grid_num_z + dz; ++iz) {
                        psi_velocity[get_voxel_center_index_3D(ix, iy, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)]
                            = psi_velocity[get_voxel_center_index_3D(ix, iy - 1, iz, Grid_num_x + dx, Grid_num_y + 1 + dy, Grid_num_z + dz)]
                            + velocity[get_voxel_center_index_3D(ix, iy - 1, iz, Grid_num_x + dx, Grid_num_y + dy, Grid_num_z + dz)] * cell_length;
                    }
                }
            }
        }
    }
} //namespace smoke_simulation
