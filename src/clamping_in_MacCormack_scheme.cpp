#include "clamping_in_MacCormack_scheme.h"

namespace smoke_simulation{

// density の MacCormack 法でのclamping の処理
MY_FLOAT_TYPE clamp_substance_density_of_MacCormack(
    const Grid_3D& all_grid,
    MY_FLOAT_TYPE substance_density_after_advect,
    VEC3_TYPE after_backtrace_position,
    const std::string interpolation_method
) {
    // interpolation に使った値を記録する変数
    std::vector<MY_FLOAT_TYPE> interpolant_substance_density_values;
    if (interpolation_method == "linear") {
        // x軸に垂直な壁での積分
        int advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        int advected_index_y = floor((after_backtrace_position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        int advected_index_z = floor((after_backtrace_position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        if (advected_index_x <= 0) {
            advected_index_x = 0;
        }
        if (advected_index_x >= all_grid.Grid_num_x - 1) {
            advected_index_x = all_grid.Grid_num_x - 2;
        }
        if (advected_index_y < 0) {
            advected_index_y = 0;
        }
        if (advected_index_y >= all_grid.Grid_num_y - 1) {
            advected_index_y = all_grid.Grid_num_y - 2;
        }
        if (advected_index_z < 0) {
            advected_index_z = 0;
        }
        if (advected_index_z >= all_grid.Grid_num_z - 1) {
            advected_index_z = all_grid.Grid_num_z - 2;
        }
        for (int ix = 0; ix < 2; ++ix) {
            for (int iy = 0; iy < 2; ++iy) {
                for (int iz = 0; iz < 2; ++iz) {
                    interpolant_substance_density_values.push_back(
                        all_grid.substance_density[get_voxel_center_index_3D(advected_index_x + ix, advected_index_y + iy, advected_index_z + iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                    );
                }
            }
        }
    }
    else if (interpolation_method == "WENO4") {
        // x軸に垂直な壁での積分
        int advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        int advected_index_y = floor((after_backtrace_position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        int advected_index_z = floor((after_backtrace_position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        for (int ix = 0; ix < 4; ix++) {
            for (int iy = 0; iy < 4; ++iy) {
                for (int iz = 0; iz < 4; ++iz) {
                    int index_x = advected_index_x - 1 + ix;
                    int index_y = advected_index_y - 1 + iy;
                    int index_z = advected_index_z - 1 + iz;
                    if (index_x < 0) {
                        index_x = 0;
                    }
                    if (index_x > all_grid.Grid_num_x - 1) {
                        int exceed = index_x - (all_grid.Grid_num_x - 1);
                        //境界の値が外側までずっと続く場合
                        index_x = all_grid.Grid_num_x - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_x = all_grid.Grid_num_x - 1 - exceed + 1;
                    }
                    if (index_y < 0) {
                        index_y = 0;
                    }
                    if (index_y > all_grid.Grid_num_y - 1) {
                        int exceed = index_y - (all_grid.Grid_num_y - 1);
                        //境界の値が外側までずっと続く場合
                        index_y = all_grid.Grid_num_y - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_y = all_grid.Grid_num_y - 1 - exceed + 1;
                    }
                    if (index_z < 0) {
                        index_z = 0;
                    }
                    if (index_z > all_grid.Grid_num_z - 1) {
                        int exceed = index_z - (all_grid.Grid_num_z - 1);
                        //境界の値が外側までずっと続く場合
                        index_z = all_grid.Grid_num_z - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_z = all_grid.Grid_num_z - 1 - exceed + 1;
                    }
                    interpolant_substance_density_values.push_back(
                        all_grid.substance_density[get_voxel_center_index_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                    );
                }
            }
        }
    }
    else if (interpolation_method == "WENO6") {
        // x軸に垂直な壁での積分
        int advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        int advected_index_y = floor((after_backtrace_position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        int advected_index_z = floor((after_backtrace_position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        for (int ix = 0; ix < 6; ix++) {
            for (int iy = 0; iy < 6; ++iy) {
                for (int iz = 0; iz < 6; ++iz) {
                    int index_x = advected_index_x - 2 + ix;
                    int index_y = advected_index_y - 2 + iy;
                    int index_z = advected_index_z - 2 + iz;
                    if (index_x < 0) {
                        index_x = 0;
                    }
                    if (index_x > all_grid.Grid_num_x - 1) {
                        int exceed = index_x - (all_grid.Grid_num_x - 1);
                        //境界の値が外側までずっと続く場合
                        index_x = all_grid.Grid_num_x - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_x = all_grid.Grid_num_x - 1 - exceed + 1;
                    }
                    if (index_y < 0) {
                        index_y = 0;
                    }
                    if (index_y > all_grid.Grid_num_y - 1) {
                        int exceed = index_y - (all_grid.Grid_num_y - 1);
                        //境界の値が外側までずっと続く場合
                        index_y = all_grid.Grid_num_y - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_y = all_grid.Grid_num_y - 1 - exceed + 1;
                    }
                    if (index_z < 0) {
                        index_z = 0;
                    }
                    if (index_z > all_grid.Grid_num_z - 1) {
                        int exceed = index_z - (all_grid.Grid_num_z - 1);
                        //境界の値が外側までずっと続く場合
                        index_z = all_grid.Grid_num_z - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_z = all_grid.Grid_num_z - 1 - exceed + 1;
                    }
                    interpolant_substance_density_values.push_back(
                        all_grid.substance_density[get_voxel_center_index_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                    );
                }
            }
        }
    }
    // clamping の処理
    std::sort(interpolant_substance_density_values.begin(), interpolant_substance_density_values.end());
    substance_density_after_advect = std::min(substance_density_after_advect, interpolant_substance_density_values[interpolant_substance_density_values.size() - 1]);
    substance_density_after_advect = std::max(substance_density_after_advect, interpolant_substance_density_values[0]);
    return substance_density_after_advect;
}
// density の MacCormack 法でのclamping の処理
MY_FLOAT_TYPE clamp_integral_of_normal_component_of_psi_of_MacCormack_3D(
    const Grid_3D& all_grid,
    MY_FLOAT_TYPE integral_of_normal_component_of_psi_after_advect,
    VEC3_TYPE after_backtrace_position,
    const std::string interpolation_method
) {
    // interpolation に使った値を記録する変数
    std::vector<MY_FLOAT_TYPE> interpolant_integral_of_normal_component_of_psi_values;
    if (interpolation_method == "linear") {
        // x軸に垂直な壁での積分
        int advected_index_x = floor((after_backtrace_position[0]) / all_grid._cell_length);
        int advected_index_y = floor((after_backtrace_position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        int advected_index_z = floor((after_backtrace_position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        if (advected_index_x <= 0) {
            advected_index_x = 0;
        }
        if (advected_index_x >= all_grid.Grid_num_x) {
            advected_index_x = all_grid.Grid_num_x - 1;
        }
        if (advected_index_y < 0) {
            advected_index_y = 0;
        }
        if (advected_index_y >= all_grid.Grid_num_y - 1) {
            advected_index_y = all_grid.Grid_num_y - 2;
        }
        if (advected_index_z < 0) {
            advected_index_z = 0;
        }
        if (advected_index_z >= all_grid.Grid_num_z - 1) {
            advected_index_z = all_grid.Grid_num_z - 2;
        }
        for (int ix = 0; ix < 2; ++ix) {
            for (int iy = 0; iy < 2; ++iy) {
                for (int iz = 0; iz < 2; ++iz) {
                    interpolant_integral_of_normal_component_of_psi_values.push_back(
                        all_grid.psi_substance_density_cell_face_x[get_voxel_face_index_x_3D(advected_index_x + ix, advected_index_y + iy, advected_index_z + iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        * all_grid._cell_length
                    );
                    interpolant_integral_of_normal_component_of_psi_values.push_back(
                        - all_grid.psi_substance_density_cell_face_x[get_voxel_face_index_x_3D(advected_index_x + ix, advected_index_y + iy, advected_index_z + iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        * all_grid._cell_length
                    );
                }
            }
        }
        // y軸に垂直な壁での積分
        advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        advected_index_y = floor((after_backtrace_position[1]) / all_grid._cell_length);
        advected_index_z = floor((after_backtrace_position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        if (advected_index_x < 0) {
            advected_index_x = 0;
        }
        if (advected_index_x >= all_grid.Grid_num_x - 1) {
            advected_index_x = all_grid.Grid_num_x - 2;
        }
        if (advected_index_y <= 0) {
            advected_index_y = 0;
        }
        if (advected_index_y >= all_grid.Grid_num_y) {
            advected_index_y = all_grid.Grid_num_y - 1;
        }
        if (advected_index_z < 0) {
            advected_index_z = 0;
        }
        if (advected_index_z >= all_grid.Grid_num_z - 1) {
            advected_index_z = all_grid.Grid_num_z - 2;
        }
        for (int ix = 0; ix < 2; ++ix) {
            for (int iy = 0; iy < 2; ++iy) {
                for (int iz = 0; iz < 2; ++iz) {
                    interpolant_integral_of_normal_component_of_psi_values.push_back(
                        all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y_3D(advected_index_x + ix, advected_index_y + iy, advected_index_z + iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        * all_grid._cell_length
                    );
                    interpolant_integral_of_normal_component_of_psi_values.push_back(
                        -all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y_3D(advected_index_x + ix, advected_index_y + iy, advected_index_z + iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        * all_grid._cell_length
                    );
                }
            }
        }
        // z軸に垂直な壁での積分
        advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        advected_index_y = floor((after_backtrace_position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        advected_index_z = floor((after_backtrace_position[2]) / all_grid._cell_length);
        if (advected_index_x < 0) {
            advected_index_x = 0;
        }
        if (advected_index_x >= all_grid.Grid_num_x - 1) {
            advected_index_x = all_grid.Grid_num_x - 2;
        }
        if (advected_index_y <= 0) {
            advected_index_y = 0;
        }
        if (advected_index_y >= all_grid.Grid_num_y - 1) {
            advected_index_y = all_grid.Grid_num_y - 2;
        }
        if (advected_index_z < 0) {
            advected_index_z = 0;
        }
        if (advected_index_z >= all_grid.Grid_num_z) {
            advected_index_z = all_grid.Grid_num_z - 1;
        }
        for (int ix = 0; ix < 2; ++ix) {
            for (int iy = 0; iy < 2; ++iy) {
                for (int iz = 0; iz < 2; ++iz) {
                    interpolant_integral_of_normal_component_of_psi_values.push_back(
                        all_grid.psi_substance_density_cell_face_z[get_voxel_face_index_z_3D(advected_index_x + ix, advected_index_y + iy, advected_index_z + iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        * all_grid._cell_length
                    );
                    interpolant_integral_of_normal_component_of_psi_values.push_back(
                        -all_grid.psi_substance_density_cell_face_z[get_voxel_face_index_z_3D(advected_index_x + ix, advected_index_y + iy, advected_index_z + iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        * all_grid._cell_length
                    );
                }
            }
        }

    }
    else if (interpolation_method == "WENO6") {
        // x軸に垂直な壁での積分
        int advected_index_x = floor((after_backtrace_position[0]) / all_grid._cell_length);
        int advected_index_y = floor((after_backtrace_position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        int advected_index_z = floor((after_backtrace_position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        for (int ix = 0; ix < 6; ix++) {
            for (int iy = 0; iy < 6; ++iy) {
                for (int iz = 0; iz < 6; ++iz) {
                    int index_x = advected_index_x - 2 + ix;
                    int index_y = advected_index_y - 2 + iy;
                    int index_z = advected_index_z - 2 + iz;
                    if (index_x < 0) {
                        index_x = 0;
                    }
                    if (index_x > all_grid.Grid_num_x) {
                        int exceed = index_x - all_grid.Grid_num_x;
                        //境界の値が外側までずっと続く場合
                        index_x = all_grid.Grid_num_x;
                        //境界を境に鏡のように値が反射する場合
                        //index_x = all_grid.Grid_num_x- exceed + 1;
                    }
                    if (index_y < 0) {
                        index_y = 0;
                    }
                    if (index_y > all_grid.Grid_num_y - 1) {
                        int exceed = index_y - (all_grid.Grid_num_y - 1);
                        //境界の値が外側までずっと続く場合
                        index_y = all_grid.Grid_num_y - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_y = all_grid.Grid_num_y - 1 - exceed + 1;
                    }
                    if (index_z < 0) {
                        index_z = 0;
                    }
                    if (index_z > all_grid.Grid_num_z - 1) {
                        int exceed = index_z - (all_grid.Grid_num_z - 1);
                        //境界の値が外側までずっと続く場合
                        index_z = all_grid.Grid_num_z - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_z = all_grid.Grid_num_z - 1 - exceed + 1;
                    }
                    interpolant_integral_of_normal_component_of_psi_values.push_back(
                        all_grid.psi_substance_density_cell_face_x[get_voxel_face_index_x_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        * all_grid._cell_length
                    );
                    interpolant_integral_of_normal_component_of_psi_values.push_back(
                        -all_grid.psi_substance_density_cell_face_x[get_voxel_face_index_x_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        * all_grid._cell_length
                    );
                }
            }
        }
        // y軸に垂直な壁での積分
        advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        advected_index_y = floor((after_backtrace_position[1]) / all_grid._cell_length);
        advected_index_z = floor((after_backtrace_position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        for (int ix = 0; ix < 6; ix++) {
            for (int iy = 0; iy < 6; ++iy) {
                for (int iz = 0; iz < 6; ++iz) {
                    int index_x = advected_index_x - 2 + ix;
                    int index_y = advected_index_y - 2 + iy;
                    int index_z = advected_index_z - 2 + iz;
                    if (index_x < 0) {
                        index_x = 0;
                    }
                    if (index_x > all_grid.Grid_num_x - 1) {
                        int exceed = index_x - (all_grid.Grid_num_x - 1);
                        //境界の値が外側までずっと続く場合
                        index_x = all_grid.Grid_num_x - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_x = all_grid.Grid_num_x - 1 - exceed + 1;
                    }
                    if (index_y < 0) {
                        index_y = 0;
                    }
                    if (index_y > all_grid.Grid_num_y) {
                        int exceed = index_y - all_grid.Grid_num_y;
                        //境界の値が外側までずっと続く場合
                        index_y = all_grid.Grid_num_y;
                        //境界を境に鏡のように値が反射する場合
                        //index_y = all_grid.Grid_num_y - exceed + 1;
                    }
                    if (index_z < 0) {
                        index_z = 0;
                    }
                    if (index_z > all_grid.Grid_num_z - 1) {
                        int exceed = index_z - (all_grid.Grid_num_z - 1);
                        //境界の値が外側までずっと続く場合
                        index_z = all_grid.Grid_num_z - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_z = all_grid.Grid_num_z - 1 - exceed + 1;
                    }
                    interpolant_integral_of_normal_component_of_psi_values.push_back(
                        all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        * all_grid._cell_length
                    );
                    interpolant_integral_of_normal_component_of_psi_values.push_back(
                        -all_grid.psi_substance_density_cell_face_y[get_voxel_face_index_y_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        * all_grid._cell_length
                    );
                }
            }
        }
        // z軸に垂直な壁での積分
        advected_index_x = floor((after_backtrace_position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        advected_index_y = floor((after_backtrace_position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        advected_index_z = floor((after_backtrace_position[2]) / all_grid._cell_length);
        for (int ix = 0; ix < 6; ix++) {
            for (int iy = 0; iy < 6; ++iy) {
                for (int iz = 0; iz < 6; ++iz) {
                    int index_x = advected_index_x - 2 + ix;
                    int index_y = advected_index_y - 2 + iy;
                    int index_z = advected_index_z - 2 + iz;
                    if (index_x < 0) {
                        index_x = 0;
                    }
                    if (index_x > all_grid.Grid_num_x - 1) {
                        int exceed = index_x - (all_grid.Grid_num_x - 1);
                        //境界の値が外側までずっと続く場合
                        index_x = all_grid.Grid_num_x - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_x = all_grid.Grid_num_x - 1 - exceed + 1;
                    }
                    if (index_y < 0) {
                        index_y = 0;
                    }
                    if (index_y > all_grid.Grid_num_y - 1) {
                        int exceed = index_y - (all_grid.Grid_num_y - 1);
                        //境界の値が外側までずっと続く場合
                        index_y = all_grid.Grid_num_y - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_y = all_grid.Grid_num_y - 1 - exceed + 1;
                    }
                    if (index_z < 0) {
                        index_z = 0;
                    }
                    if (index_z > all_grid.Grid_num_z) {
                        int exceed = index_z - all_grid.Grid_num_z;
                        //境界の値が外側までずっと続く場合
                        index_z = all_grid.Grid_num_z;
                        //境界を境に鏡のように値が反射する場合
                        //index_z = all_grid.Grid_num_z - exceed + 1;
                    }
                    interpolant_integral_of_normal_component_of_psi_values.push_back(
                        all_grid.psi_substance_density_cell_face_z[get_voxel_face_index_z_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        * all_grid._cell_length
                    );
                    interpolant_integral_of_normal_component_of_psi_values.push_back(
                        -all_grid.psi_substance_density_cell_face_z[get_voxel_face_index_z_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        * all_grid._cell_length
                    );
                }
            }
        }
    }
    // clamping の処理
    std::sort(interpolant_integral_of_normal_component_of_psi_values.begin(), interpolant_integral_of_normal_component_of_psi_values.end());
    integral_of_normal_component_of_psi_after_advect = std::min(integral_of_normal_component_of_psi_after_advect, interpolant_integral_of_normal_component_of_psi_values[interpolant_integral_of_normal_component_of_psi_values.size() - 1]);
    integral_of_normal_component_of_psi_after_advect = std::max(integral_of_normal_component_of_psi_after_advect, interpolant_integral_of_normal_component_of_psi_values[0]);
    return integral_of_normal_component_of_psi_after_advect;
}
}//namespace smoke_simulation
