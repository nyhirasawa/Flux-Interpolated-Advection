#include "linear_interpolation_3d.h"

#include <vector>
#include "utils.h"

namespace smoke_simulation{
    MY_FLOAT_TYPE linear_interpolation_3D(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_3D &all_grid,
        const VEC3_TYPE origin_pos,
        const int grid_num_x,
        const int grid_num_y,
        const int grid_num_z,
        const MY_FLOAT_TYPE cell_length
    ){
        int advected_index_x = floor((position[0] - origin_pos[0]) / cell_length);
        int advected_index_y = floor((position[1] - origin_pos[1]) / cell_length);
        int advected_index_z = floor((position[2] - origin_pos[2]) / cell_length);

        /////linear interpolation の処理
        //z方向のinterpolationの結果を格納するための変数
        MY_FLOAT_TYPE z_interpolation_values_of_substance_density[2][2];
        //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
        for (int ix = 0; ix < 2; ix++) {
            for(int iy = 0; iy < 2; ++iy){
                //補間に使う離散値をセット
                MY_FLOAT_TYPE psi_val[2];
                for (int iz = 0; iz < 2; ++iz) {
                    int index_x = advected_index_x + ix;
                    int index_y = advected_index_y + iy;
                    int index_z = advected_index_z + iz;
                    // グリッドの外側を参照しようとしたときの処理(x方向)
                    if (index_x < 0) {
                        index_x = 0;
                    }
                    if (index_x > grid_num_x - 1) {
                        int exceed = index_x - (grid_num_x - 1);
                        //境界の値が外側までずっと続く場合
                        index_x = grid_num_x - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_x = grid_num_x - 1 - exceed + 1;
                    }
                    // グリッドの外側を参照しようとしたときの処理(y方向)
                    if (index_y < 0) {
                        index_y = 0;
                    }
                    if (index_y > grid_num_y - 1) {
                        int exceed = index_y - (grid_num_y - 1);
                        //境界の値が外側までずっと続く場合
                        index_y = grid_num_y - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_y = grid_num_y - 1 - exceed + 1;
                    }
                    // グリッドの外側を参照しようとしたときの処理(y方向)
                    if (index_z < 0) {
                        index_z = 0;
                    }
                    if (index_z > grid_num_z - 1) {
                        int exceed = index_z - (grid_num_z - 1);
                        //境界の値が外側までずっと続く場合
                        index_z = grid_num_z - 1;
                        //境界を境に鏡のように値が反射する場合
                        //index_z = grid_num_z - 1 - exceed + 1;
                    }
                    psi_val[iz] = cell_center_values[get_voxel_center_index_3D(index_x, index_y, index_z, grid_num_x, grid_num_y, grid_num_z)];
                }
                //interpolation の係数
                MY_FLOAT_TYPE c0, c1;
                c0 = (position[2] - origin_pos[2]) / cell_length - advected_index_z;
                c1 = 1.0 - c0;
                z_interpolation_values_of_substance_density[ix][iy] = c1 * psi_val[0] + c0 * psi_val[1];
            }
        }

        ////ここからy方向の補間
        //yz方向のinterpolationの結果を格納するための変数
        MY_FLOAT_TYPE yz_interpolation_values_of_substance_density[2];
        for(int ix = 0; ix < 2; ++ix){
            //interpolation の係数
            MY_FLOAT_TYPE b0, b1;
            b0 = (position[1] - origin_pos[1]) / cell_length - advected_index_y;
            b1 = 1.0 - b0;
            yz_interpolation_values_of_substance_density[ix]
                = b1 * z_interpolation_values_of_substance_density[ix][0] + b0 * z_interpolation_values_of_substance_density[ix][1];
        }

        ////ここからx方向の補間
        //interpolation の係数
        MY_FLOAT_TYPE a0, a1;
        a0 = (position[0] - origin_pos[0]) / cell_length - advected_index_x;
        a1 = 1.0 - a0;
        //結果
        return a1 * yz_interpolation_values_of_substance_density[0] + a0 * yz_interpolation_values_of_substance_density[1];
    }

MY_FLOAT_TYPE linear_interpolation_3D_cell_center_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_center_values,
    const Grid_3D &all_grid
){
    return linear_interpolation_3D(
        position,
        cell_center_values,
        all_grid,
        VEC3_TYPE(0.5 * all_grid._cell_length, 0.5 * all_grid._cell_length, 0.5 * all_grid._cell_length),
        all_grid.Grid_num_x,
        all_grid.Grid_num_y,
        all_grid.Grid_num_z,
        all_grid._cell_length
    );

/*
    int advected_index_x = floor((position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_z = floor((position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length);

    /////linear interpolation の処理
    //z方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE z_interpolation_values_of_substance_density[2][2];
    //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    for (int ix = 0; ix < 2; ix++) {
        for(int iy = 0; iy < 2; ++iy){
            //補間に使う離散値をセット
            MY_FLOAT_TYPE psi_val[2];
            for (int iz = 0; iz < 2; ++iz) {
                int index_x = advected_index_x + ix;
                int index_y = advected_index_y + iy;
                int index_z = advected_index_z + iz;
                // グリッドの外側を参照しようとしたときの処理(x方向)
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
                // グリッドの外側を参照しようとしたときの処理(y方向)
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
                // グリッドの外側を参照しようとしたときの処理(y方向)
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
                psi_val[iz] = cell_center_values[get_voxel_center_index_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
            }
            //interpolation の係数
            MY_FLOAT_TYPE c0, c1;
            c0 = (position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_z;
            c1 = 1.0 - c0;
            z_interpolation_values_of_substance_density[ix][iy] = c1 * psi_val[0] + c0 * psi_val[1];
        }
    }

    ////ここからy方向の補間
    //yz方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE yz_interpolation_values_of_substance_density[2];
    for(int ix = 0; ix < 2; ++ix){
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_y;
        b1 = 1.0 - b0;
        yz_interpolation_values_of_substance_density[ix]
            = b1 * z_interpolation_values_of_substance_density[ix][0] + b0 * z_interpolation_values_of_substance_density[ix][1];
    }

    ////ここからx方向の補間
    //interpolation の係数
    MY_FLOAT_TYPE a0, a1;
    a0 = (position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_x;
    a1 = 1.0 - a0;
    //結果
    return a1 * yz_interpolation_values_of_substance_density[0] + a0 * yz_interpolation_values_of_substance_density[1];
*/
}

MY_FLOAT_TYPE linear_interpolation_3D_cell_face_values(
    const int dim,
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
    const Grid_3D &all_grid
){
    if(dim == 0){
        return linear_interpolation_3D_cell_face_x_values(
            position,
            cell_face_x_values,
            all_grid
        );
    }
    else if(dim == 1){
        return linear_interpolation_3D_cell_face_y_values(
            position,
            cell_face_x_values,
            all_grid
        );
    }
    else if(dim == 2){
        return linear_interpolation_3D_cell_face_z_values(
            position,
            cell_face_x_values,
            all_grid
        );
    }
}

MY_FLOAT_TYPE linear_interpolation_3D_cell_face_x_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
    const Grid_3D &all_grid
){
    return linear_interpolation_3D(
        position,
        cell_face_x_values,
        all_grid,
        VEC3_TYPE(0.0, 0.5 * all_grid._cell_length, 0.5 * all_grid._cell_length),
        all_grid.Grid_num_x + 1,
        all_grid.Grid_num_y,
        all_grid.Grid_num_z,
        all_grid._cell_length
    );
/*
    int advected_index_x = floor((position[0]) / all_grid._cell_length);
    int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_z = floor((position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length);

    /////linear interpolation の処理
    //z方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE z_interpolation_values_of_psi_substance_density[2][2];
    //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    for (int ix = 0; ix < 2; ix++) {
        for(int iy = 0; iy < 2; ++iy){
            //補間に使う離散値をセット
            MY_FLOAT_TYPE psi_val[2];
            for (int iz = 0; iz < 2; ++iz) {
                int index_x = advected_index_x + ix;
                int index_y = advected_index_y + iy;
                int index_z = advected_index_z + iz;
                // グリッドの外側を参照しようとしたときの処理(x方向)
                if (index_x < 0) {
                    index_x = 0;
                }
                if (index_x > all_grid.Grid_num_x) {
                    int exceed = index_x - all_grid.Grid_num_x;
                    //境界の値が外側までずっと続く場合
                    index_x = all_grid.Grid_num_x;
                    //境界を境に鏡のように値が反射する場合
                    //index_x = all_grid.Grid_num_x - exceed + 1;
                }
                // グリッドの外側を参照しようとしたときの処理(y方向)
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
                // グリッドの外側を参照しようとしたときの処理(y方向)
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
                psi_val[iz] = cell_face_x_values[get_voxel_face_index_x_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
            }
            //interpolation の係数
            MY_FLOAT_TYPE c0, c1;
            c0 = (position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_z;
            c1 = 1.0 - c0;
            z_interpolation_values_of_psi_substance_density[ix][iy] = c1 * psi_val[0] + c0 * psi_val[1];
        }
    }

    ////ここからy方向の補間
    //yz方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE yz_interpolation_values_of_psi_substance_density[2];
    for(int ix = 0; ix < 2; ++ix){
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_y;
        b1 = 1.0 - b0;
        yz_interpolation_values_of_psi_substance_density[ix]
            = b1 * z_interpolation_values_of_psi_substance_density[ix][0] + b0 * z_interpolation_values_of_psi_substance_density[ix][1];
    }

    ////ここからx方向の補間
    //interpolation の係数
    MY_FLOAT_TYPE a0, a1;
    a0 = (position[0]) / all_grid._cell_length - advected_index_x;
    a1 = 1.0 - a0;
    //結果
    return a1 * yz_interpolation_values_of_psi_substance_density[0] + a0 * yz_interpolation_values_of_psi_substance_density[1];
*/
}

MY_FLOAT_TYPE linear_interpolation_3D_cell_face_y_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_y_values,
    const Grid_3D &all_grid
){
    return linear_interpolation_3D(
        position,
        cell_face_y_values,
        all_grid,
        VEC3_TYPE(0.5 * all_grid._cell_length, 0.0, 0.5 * all_grid._cell_length),
        all_grid.Grid_num_x,
        all_grid.Grid_num_y + 1,
        all_grid.Grid_num_z,
        all_grid._cell_length
    );

/*
    int advected_index_x = floor((position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_y = floor((position[1]) / all_grid._cell_length);
    int advected_index_z = floor((position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length);

    /////linear interpolation の処理
    //z方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE z_interpolation_values_of_psi_substance_density[2][2];
    //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    for (int ix = 0; ix < 2; ix++) {
        for(int iy = 0; iy < 2; ++iy){
            //補間に使う離散値をセット
            MY_FLOAT_TYPE psi_val[2];
            for (int iz = 0; iz < 2; ++iz) {
                int index_x = advected_index_x + ix;
                int index_y = advected_index_y + iy;
                int index_z = advected_index_z + iz;
                // グリッドの外側を参照しようとしたときの処理(x方向)
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
                // グリッドの外側を参照しようとしたときの処理(y方向)
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
                // グリッドの外側を参照しようとしたときの処理(y方向)
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
                psi_val[iz] = cell_face_y_values[get_voxel_face_index_y_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
            }
            //interpolation の係数
            MY_FLOAT_TYPE c0, c1;
            c0 = (position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_z;
            c1 = 1.0 - c0;
            z_interpolation_values_of_psi_substance_density[ix][iy] = c1 * psi_val[0] + c0 * psi_val[1];
        }
    }

    ////ここからy方向の補間
    //yz方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE yz_interpolation_values_of_psi_substance_density[2];
    for(int ix = 0; ix < 2; ++ix){
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1]) / all_grid._cell_length - advected_index_y;
        b1 = 1.0 - b0;
        yz_interpolation_values_of_psi_substance_density[ix]
            = b1 * z_interpolation_values_of_psi_substance_density[ix][0] + b0 * z_interpolation_values_of_psi_substance_density[ix][1];
    }

    ////ここからx方向の補間
    //interpolation の係数
    MY_FLOAT_TYPE a0, a1;
    a0 = (position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_x;
    a1 = 1.0 - a0;
    //結果
    return a1 * yz_interpolation_values_of_psi_substance_density[0] + a0 * yz_interpolation_values_of_psi_substance_density[1];
*/
}

MY_FLOAT_TYPE linear_interpolation_3D_cell_face_z_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_z_values,
    const Grid_3D &all_grid
){
    return linear_interpolation_3D(
        position,
        cell_face_z_values,
        all_grid,
        VEC3_TYPE(0.5 * all_grid._cell_length, 0.5 * all_grid._cell_length, 0.0),
        all_grid.Grid_num_x,
        all_grid.Grid_num_y,
        all_grid.Grid_num_z + 1,
        all_grid._cell_length
    );
/*
    int advected_index_x = floor((position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_z = floor((position[2]) / all_grid._cell_length);

    /////linear interpolation の処理
    //z方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE z_interpolation_values_of_psi_substance_density[2][2];
    //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    for (int ix = 0; ix < 2; ix++) {
        for(int iy = 0; iy < 2; ++iy){
            //補間に使う離散値をセット
            MY_FLOAT_TYPE psi_val[2];
            for (int iz = 0; iz < 2; ++iz) {
                int index_x = advected_index_x + ix;
                int index_y = advected_index_y + iy;
                int index_z = advected_index_z + iz;
                // グリッドの外側を参照しようとしたときの処理(x方向)
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
                // グリッドの外側を参照しようとしたときの処理(y方向)
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
                // グリッドの外側を参照しようとしたときの処理(y方向)
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
                psi_val[iz] = cell_face_z_values[get_voxel_face_index_z_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
            }
            //interpolation の係数
            MY_FLOAT_TYPE c0, c1;
            c0 = (position[2]) / all_grid._cell_length - advected_index_z;
            c1 = 1.0 - c0;
            z_interpolation_values_of_psi_substance_density[ix][iy] = c1 * psi_val[0] + c0 * psi_val[1];
        }
    }

    ////ここからy方向の補間
    //yz方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE yz_interpolation_values_of_psi_substance_density[2];
    for(int ix = 0; ix < 2; ++ix){
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_y;
        b1 = 1.0 - b0;
        yz_interpolation_values_of_psi_substance_density[ix]
            = b1 * z_interpolation_values_of_psi_substance_density[ix][0] + b0 * z_interpolation_values_of_psi_substance_density[ix][1];
    }

    ////ここからx方向の補間
    //interpolation の係数
    MY_FLOAT_TYPE a0, a1;
    a0 = (position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_x;
    a1 = 1.0 - a0;
    //結果
    return a1 * yz_interpolation_values_of_psi_substance_density[0] + a0 * yz_interpolation_values_of_psi_substance_density[1];
*/
}

MY_FLOAT_TYPE linear_interpolation_3D_psi_velocity_x(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_center_values,
    const Grid_3D &all_grid
) {
    return linear_interpolation_3D(
        position,
        cell_center_values,
        all_grid,
        VEC3_TYPE(0.0, 0.0, 0.5 * all_grid._cell_length),
        all_grid.Grid_num_x + 1,
        all_grid.Grid_num_y + 1,
        all_grid.Grid_num_z,
        all_grid._cell_length
    );
}

MY_FLOAT_TYPE linear_interpolation_3D_psi_velocity_y(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_center_values,
    const Grid_3D &all_grid
) {
    return linear_interpolation_3D(
        position,
        cell_center_values,
        all_grid,
        VEC3_TYPE(0.5 * all_grid._cell_length, -0.5 * all_grid._cell_length, 0.5 * all_grid._cell_length),
        all_grid.Grid_num_x,
        all_grid.Grid_num_y + 2,
        all_grid.Grid_num_z,
        all_grid._cell_length
    );
}

MY_FLOAT_TYPE linear_interpolation_3D_psi_velocity_z(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_center_values,
    const Grid_3D &all_grid
) {
    return linear_interpolation_3D(
        position,
        cell_center_values,
        all_grid,
        VEC3_TYPE(0.5 * all_grid._cell_length, 0.0, 0.0),
        all_grid.Grid_num_x,
        all_grid.Grid_num_y + 1,
        all_grid.Grid_num_z + 1,
        all_grid._cell_length
    );
}

}//namespace smoke_simulation
