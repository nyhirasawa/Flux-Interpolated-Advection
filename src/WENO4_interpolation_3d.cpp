#include "linear_interpolation_3d.h"

#include <vector>
#include "utils.h"

namespace smoke_simulation{
MY_FLOAT_TYPE WENO4_interpolation_3D_cell_center_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_center_values,
    const Grid_3D &all_grid
){
    int advected_index_x = floor((position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_z = floor((position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length);

    /////linear interpolation の処理
    //z方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE z_interpolation_values_of_psi_substance_density[4][4];
    //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    for (int ix = 0; ix < 4; ix++) {
        for(int iy = 0; iy < 4; ++iy){
            //補間に使う離散値をセット
            MY_FLOAT_TYPE psi_val[4];
            for (int iz = 0; iz < 4; ++iz) {
                int index_x = advected_index_x - 1 + ix;
                int index_y = advected_index_y - 1 + iy;
                int index_z = advected_index_z - 1 + iz;
                // グリッドの外側を参照しようとしたときの処理(x方向)
                if (index_x < 0) {
                    index_x = 0;
                }
                if (index_x > all_grid.Grid_num_x - 1) {
                    int exceed = index_x - (all_grid.Grid_num_x - 1);
                    //境界の値が外側までずっと続く場合
                    index_x = all_grid.Grid_num_x - 1;
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
                psi_val[iz] = cell_center_values[get_voxel_center_index_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
            }
            // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
            MY_FLOAT_TYPE z_distance_from_face_to_position[4];
            for (int i = 0; i < 4; ++i) {
                z_distance_from_face_to_position[i] = position[2] - (advected_index_z - 1 + i + 0.5) * all_grid._cell_length;
            }

            // candidate interpolants の計算
            MY_FLOAT_TYPE candidate_interpolants[2];
            candidate_interpolants[0]
                = psi_val[1]
                + ((psi_val[2] - psi_val[0]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
                + ((psi_val[2] - 2.0 * psi_val[1] + psi_val[0]) * (z_distance_from_face_to_position[1]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);
            candidate_interpolants[1]
                = psi_val[1]
                + ((-psi_val[3] + 4.0 * psi_val[2] - 3.0 * psi_val[1]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
                + ((psi_val[3] - 2.0 * psi_val[2] + psi_val[1]) * (z_distance_from_face_to_position[1]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);

            //ideal weights の計算
            MY_FLOAT_TYPE ideal_weights[2];
            ideal_weights[0]
                = -z_distance_from_face_to_position[3] / (3.0 * all_grid._cell_length);
            ideal_weights[1]
                = z_distance_from_face_to_position[0] / (3.0 * all_grid._cell_length);

            // smoothness indicatior の計算
            MY_FLOAT_TYPE smoothness_indicators[2];
            smoothness_indicators[0]
                = (26.0 * psi_val[2] * psi_val[0] - 52.0 * psi_val[1] * psi_val[0] - 76.0 * psi_val[2] * psi_val[1] + 25.0 * psi_val[2] * psi_val[2] + 64.0 * psi_val[1] * psi_val[1] + 13.0 * psi_val[0] * psi_val[0]) / 12.0;
            smoothness_indicators[1]
                = (26.0 * psi_val[3] * psi_val[1] - 52.0 * psi_val[3] * psi_val[2] - 76.0 * psi_val[2] * psi_val[1] + 25.0 * psi_val[1] * psi_val[1] + 64.0 * psi_val[2] * psi_val[2] + 13.0 * psi_val[3] * psi_val[3]) / 12.0;

            //weight を作るための1時的な変数
            MY_FLOAT_TYPE alpha[2];
            for (int i = 0; i < 2; ++i) {
                MY_FLOAT_TYPE eps = 0.000001;
                alpha[i] = ideal_weights[i] / ((eps + smoothness_indicators[i]) * (eps + smoothness_indicators[i]));
            }
            //weights の計算
            MY_FLOAT_TYPE weights[2];
            for (int i = 0; i < 2; ++i) {
                weights[i] = alpha[i] / (alpha[0] + alpha[1]);
            }
            z_interpolation_values_of_psi_substance_density[ix][iy] = weights[0] * candidate_interpolants[0] + weights[1] * candidate_interpolants[1];
        }
    }

    ////ここからy方向の補間
    //yz方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE yz_interpolation_values_of_psi_substance_density[4];
    for(int ix = 0; ix < 4; ++ix){
        // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
        MY_FLOAT_TYPE y_distance_from_face_to_position[4];
        for (int iy = 0; iy < 4; ++iy) {
            // (advected_index_x - 2 + ix + 0.5) の 0.5 は考えているfaceがy軸に垂直なものを考えているため、面の端から中心までのオフセット分の距離
            y_distance_from_face_to_position[iy] = position[1] - (advected_index_y - 1 + iy + 0.5) * all_grid._cell_length;
        }
        // candidate interpolants の計算
        MY_FLOAT_TYPE y_candidate_interpolants[2];
        y_candidate_interpolants[0]
            = z_interpolation_values_of_psi_substance_density[ix][1]
            + ((z_interpolation_values_of_psi_substance_density[ix][2] - z_interpolation_values_of_psi_substance_density[ix][0]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
            + ((z_interpolation_values_of_psi_substance_density[ix][2] - 2.0 * z_interpolation_values_of_psi_substance_density[ix][1] + z_interpolation_values_of_psi_substance_density[ix][0]) * (y_distance_from_face_to_position[1]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);
        y_candidate_interpolants[1]
            = z_interpolation_values_of_psi_substance_density[ix][1]
            + ((-z_interpolation_values_of_psi_substance_density[ix][3] + 4.0 * z_interpolation_values_of_psi_substance_density[ix][2] - 3.0 * z_interpolation_values_of_psi_substance_density[ix][1]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
            + ((z_interpolation_values_of_psi_substance_density[ix][3] - 2.0 * z_interpolation_values_of_psi_substance_density[ix][2] + z_interpolation_values_of_psi_substance_density[ix][1]) * (y_distance_from_face_to_position[1]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);

        //ideal weights の計算
        MY_FLOAT_TYPE y_ideal_weights[2];
        y_ideal_weights[0]
            = -y_distance_from_face_to_position[3] / (3.0 * all_grid._cell_length);
        y_ideal_weights[1]
            = y_distance_from_face_to_position[0] / (3.0 * all_grid._cell_length);

        // smoothness indicatior の計算
        MY_FLOAT_TYPE y_smoothness_indicators[2];
        y_smoothness_indicators[0]
            = (26.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][0]
             - 52.0 * z_interpolation_values_of_psi_substance_density[ix][1] * z_interpolation_values_of_psi_substance_density[ix][0]
             - 76.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 25.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][2]
             + 64.0 * z_interpolation_values_of_psi_substance_density[ix][1] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 13.0 * z_interpolation_values_of_psi_substance_density[ix][0] * z_interpolation_values_of_psi_substance_density[ix][0]) / 12.0;
        y_smoothness_indicators[1]
            = (26.0 * z_interpolation_values_of_psi_substance_density[ix][3] * z_interpolation_values_of_psi_substance_density[ix][1]
             - 52.0 * z_interpolation_values_of_psi_substance_density[ix][3] * z_interpolation_values_of_psi_substance_density[ix][2]
             - 76.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 25.0 * z_interpolation_values_of_psi_substance_density[ix][1] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 64.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][2]
             + 13.0 * z_interpolation_values_of_psi_substance_density[ix][3] * z_interpolation_values_of_psi_substance_density[ix][3]) / 12.0;

        //weight を作るための1時的な変数
        MY_FLOAT_TYPE y_alpha[2];
        for (int i = 0; i < 2; ++i) {
            MY_FLOAT_TYPE eps = 0.000001;
            y_alpha[i] = y_ideal_weights[i] / ((eps + y_smoothness_indicators[i]) * (eps + y_smoothness_indicators[i]));
        }
        //weights の計算
        MY_FLOAT_TYPE y_weights[2];
        for (int i = 0; i < 2; ++i) {
            y_weights[i] = y_alpha[i] / (y_alpha[0] + y_alpha[1]);
        }
        yz_interpolation_values_of_psi_substance_density[ix] = y_weights[0] * y_candidate_interpolants[0] + y_weights[1] * y_candidate_interpolants[1];
    }

    ////ここからx方向の補間
    // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
    MY_FLOAT_TYPE x_distance_from_face_to_position[4];
    for (int ix = 0; ix < 4; ++ix) {
        // (advected_index_x - 2 + ix + 0.5) の 0.5 は考えているfaceがy軸に垂直なものを考えているため、面の端から中心までのオフセット分の距離
        x_distance_from_face_to_position[ix] = position[0] - (advected_index_x - 1 + ix + 0.5) * all_grid._cell_length;
    }

    // candidate interpolants の計算
    MY_FLOAT_TYPE x_candidate_interpolants[2];
    x_candidate_interpolants[0]
        = yz_interpolation_values_of_psi_substance_density[1]
        + ((yz_interpolation_values_of_psi_substance_density[2] - yz_interpolation_values_of_psi_substance_density[0]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
        + ((yz_interpolation_values_of_psi_substance_density[2] - 2.0 * yz_interpolation_values_of_psi_substance_density[1] + yz_interpolation_values_of_psi_substance_density[0]) * (x_distance_from_face_to_position[1]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);
    x_candidate_interpolants[1]
        = yz_interpolation_values_of_psi_substance_density[1]
        + ((-yz_interpolation_values_of_psi_substance_density[3] + 4.0 * yz_interpolation_values_of_psi_substance_density[2] - 3.0 * yz_interpolation_values_of_psi_substance_density[1]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
        + ((yz_interpolation_values_of_psi_substance_density[3] - 2.0 * yz_interpolation_values_of_psi_substance_density[2] + yz_interpolation_values_of_psi_substance_density[1]) * (x_distance_from_face_to_position[1]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);

    //ideal weights の計算
    MY_FLOAT_TYPE x_ideal_weights[2];
    x_ideal_weights[0]
        = -x_distance_from_face_to_position[3] / (3.0 * all_grid._cell_length);
    x_ideal_weights[1]
        = x_distance_from_face_to_position[0] / (3.0 * all_grid._cell_length);

    // smoothness indicatior の計算
    MY_FLOAT_TYPE x_smoothness_indicators[2];
    x_smoothness_indicators[0]
        = (26.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[0]
         - 52.0 * yz_interpolation_values_of_psi_substance_density[1] * yz_interpolation_values_of_psi_substance_density[0]
         - 76.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[1]
         + 25.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[2]
         + 64.0 * yz_interpolation_values_of_psi_substance_density[1] * yz_interpolation_values_of_psi_substance_density[1]
         + 13.0 * yz_interpolation_values_of_psi_substance_density[0] * yz_interpolation_values_of_psi_substance_density[0]) / 12.0;
    x_smoothness_indicators[1]
        = (26.0 * yz_interpolation_values_of_psi_substance_density[3] * yz_interpolation_values_of_psi_substance_density[1]
         - 52.0 * yz_interpolation_values_of_psi_substance_density[3] * yz_interpolation_values_of_psi_substance_density[2]
         - 76.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[1]
         + 25.0 * yz_interpolation_values_of_psi_substance_density[1] * yz_interpolation_values_of_psi_substance_density[1]
         + 64.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[2]
         + 13.0 * yz_interpolation_values_of_psi_substance_density[3] * yz_interpolation_values_of_psi_substance_density[3]) / 12.0;

    //weight を作るための1時的な変数
    MY_FLOAT_TYPE x_alpha[2];
    for (int i = 0; i < 2; ++i) {
        MY_FLOAT_TYPE eps = 0.000001;
        x_alpha[i] = x_ideal_weights[i] / ((eps + x_smoothness_indicators[i]) * (eps + x_smoothness_indicators[i]));
    }
    //weights の計算
    MY_FLOAT_TYPE x_weights[2];
    for (int i = 0; i < 2; ++i) {
        x_weights[i] = x_alpha[i] / (x_alpha[0] + x_alpha[1]);
    }

    return x_weights[0] * x_candidate_interpolants[0] + x_weights[1] * x_candidate_interpolants[1];
}

MY_FLOAT_TYPE WENO4_interpolation_3D_cell_face_x_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
    const Grid_3D &all_grid
){
    int advected_index_x = floor((position[0]) / all_grid._cell_length);
    int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_z = floor((position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length);

    /////linear interpolation の処理
    //z方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE z_interpolation_values_of_psi_substance_density[4][4];
    //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    for (int ix = 0; ix < 4; ix++) {
        for(int iy = 0; iy < 4; ++iy){
            //補間に使う離散値をセット
            MY_FLOAT_TYPE psi_val[4];
            for (int iz = 0; iz < 4; ++iz) {
                int index_x = advected_index_x - 1 + ix;
                int index_y = advected_index_y - 1 + iy;
                int index_z = advected_index_z - 1 + iz;
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
            // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
            MY_FLOAT_TYPE z_distance_from_face_to_position[4];
            for (int i = 0; i < 4; ++i) {
                z_distance_from_face_to_position[i] = position[2] - (advected_index_z - 1 + i + 0.5) * all_grid._cell_length;
            }

            // candidate interpolants の計算
            MY_FLOAT_TYPE candidate_interpolants[2];
            candidate_interpolants[0]
                = psi_val[1]
                + ((psi_val[2] - psi_val[0]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
                + ((psi_val[2] - 2.0 * psi_val[1] + psi_val[0]) * (z_distance_from_face_to_position[1]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);
            candidate_interpolants[1]
                = psi_val[1]
                + ((-psi_val[3] + 4.0 * psi_val[2] - 3.0 * psi_val[1]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
                + ((psi_val[3] - 2.0 * psi_val[2] + psi_val[1]) * (z_distance_from_face_to_position[1]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);

            //ideal weights の計算
            MY_FLOAT_TYPE ideal_weights[2];
            ideal_weights[0]
                = -z_distance_from_face_to_position[3] / (3.0 * all_grid._cell_length);
            ideal_weights[1]
                = z_distance_from_face_to_position[0] / (3.0 * all_grid._cell_length);

            // smoothness indicatior の計算
            MY_FLOAT_TYPE smoothness_indicators[2];
            smoothness_indicators[0]
                = (26.0 * psi_val[2] * psi_val[0] - 52.0 * psi_val[1] * psi_val[0] - 76.0 * psi_val[2] * psi_val[1] + 25.0 * psi_val[2] * psi_val[2] + 64.0 * psi_val[1] * psi_val[1] + 13.0 * psi_val[0] * psi_val[0]) / 12.0;
            smoothness_indicators[1]
                = (26.0 * psi_val[3] * psi_val[1] - 52.0 * psi_val[3] * psi_val[2] - 76.0 * psi_val[2] * psi_val[1] + 25.0 * psi_val[1] * psi_val[1] + 64.0 * psi_val[2] * psi_val[2] + 13.0 * psi_val[3] * psi_val[3]) / 12.0;

            //weight を作るための1時的な変数
            MY_FLOAT_TYPE alpha[2];
            for (int i = 0; i < 2; ++i) {
                MY_FLOAT_TYPE eps = 0.000001;
                alpha[i] = ideal_weights[i] / ((eps + smoothness_indicators[i]) * (eps + smoothness_indicators[i]));
            }
            //weights の計算
            MY_FLOAT_TYPE weights[2];
            for (int i = 0; i < 2; ++i) {
                weights[i] = alpha[i] / (alpha[0] + alpha[1]);
            }
            z_interpolation_values_of_psi_substance_density[ix][iy] = weights[0] * candidate_interpolants[0] + weights[1] * candidate_interpolants[1];
        }
    }

    ////ここからy方向の補間
    //yz方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE yz_interpolation_values_of_psi_substance_density[4];
    for(int ix = 0; ix < 4; ++ix){
        // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
        MY_FLOAT_TYPE y_distance_from_face_to_position[4];
        for (int iy = 0; iy < 4; ++iy) {
            // (advected_index_x - 2 + ix + 0.5) の 0.5 は考えているfaceがy軸に垂直なものを考えているため、面の端から中心までのオフセット分の距離
            y_distance_from_face_to_position[iy] = position[1] - (advected_index_y - 1 + iy + 0.5) * all_grid._cell_length;
        }
        // candidate interpolants の計算
        MY_FLOAT_TYPE y_candidate_interpolants[2];
        y_candidate_interpolants[0]
            = z_interpolation_values_of_psi_substance_density[ix][1]
            + ((z_interpolation_values_of_psi_substance_density[ix][2] - z_interpolation_values_of_psi_substance_density[ix][0]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
            + ((z_interpolation_values_of_psi_substance_density[ix][2] - 2.0 * z_interpolation_values_of_psi_substance_density[ix][1] + z_interpolation_values_of_psi_substance_density[ix][0]) * (y_distance_from_face_to_position[1]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);
        y_candidate_interpolants[1]
            = z_interpolation_values_of_psi_substance_density[ix][1]
            + ((-z_interpolation_values_of_psi_substance_density[ix][3] + 4.0 * z_interpolation_values_of_psi_substance_density[ix][2] - 3.0 * z_interpolation_values_of_psi_substance_density[ix][1]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
            + ((z_interpolation_values_of_psi_substance_density[ix][3] - 2.0 * z_interpolation_values_of_psi_substance_density[ix][2] + z_interpolation_values_of_psi_substance_density[ix][1]) * (y_distance_from_face_to_position[1]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);

        //ideal weights の計算
        MY_FLOAT_TYPE y_ideal_weights[2];
        y_ideal_weights[0]
            = -y_distance_from_face_to_position[3] / (3.0 * all_grid._cell_length);
        y_ideal_weights[1]
            = y_distance_from_face_to_position[0] / (3.0 * all_grid._cell_length);

        // smoothness indicatior の計算
        MY_FLOAT_TYPE y_smoothness_indicators[2];
        y_smoothness_indicators[0]
            = (26.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][0]
             - 52.0 * z_interpolation_values_of_psi_substance_density[ix][1] * z_interpolation_values_of_psi_substance_density[ix][0]
             - 76.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 25.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][2]
             + 64.0 * z_interpolation_values_of_psi_substance_density[ix][1] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 13.0 * z_interpolation_values_of_psi_substance_density[ix][0] * z_interpolation_values_of_psi_substance_density[ix][0]) / 12.0;
        y_smoothness_indicators[1]
            = (26.0 * z_interpolation_values_of_psi_substance_density[ix][3] * z_interpolation_values_of_psi_substance_density[ix][1]
             - 52.0 * z_interpolation_values_of_psi_substance_density[ix][3] * z_interpolation_values_of_psi_substance_density[ix][2]
             - 76.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 25.0 * z_interpolation_values_of_psi_substance_density[ix][1] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 64.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][2]
             + 13.0 * z_interpolation_values_of_psi_substance_density[ix][3] * z_interpolation_values_of_psi_substance_density[ix][3]) / 12.0;

        //weight を作るための1時的な変数
        MY_FLOAT_TYPE y_alpha[2];
        for (int i = 0; i < 2; ++i) {
            MY_FLOAT_TYPE eps = 0.000001;
            y_alpha[i] = y_ideal_weights[i] / ((eps + y_smoothness_indicators[i]) * (eps + y_smoothness_indicators[i]));
        }
        //weights の計算
        MY_FLOAT_TYPE y_weights[2];
        for (int i = 0; i < 2; ++i) {
            y_weights[i] = y_alpha[i] / (y_alpha[0] + y_alpha[1]);
        }
        yz_interpolation_values_of_psi_substance_density[ix] = y_weights[0] * y_candidate_interpolants[0] + y_weights[1] * y_candidate_interpolants[1];
    }

    ////ここからx方向の補間
    // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
    MY_FLOAT_TYPE x_distance_from_face_to_position[4];
    for (int ix = 0; ix < 4; ++ix) {
        // (advected_index_x - 2 + ix + 0.5) の 0.5 は考えているfaceがy軸に垂直なものを考えているため、面の端から中心までのオフセット分の距離
        x_distance_from_face_to_position[ix] = position[0] - (advected_index_x - 1 + ix) * all_grid._cell_length;
    }

    // candidate interpolants の計算
    MY_FLOAT_TYPE x_candidate_interpolants[2];
    x_candidate_interpolants[0]
        = yz_interpolation_values_of_psi_substance_density[1]
        + ((yz_interpolation_values_of_psi_substance_density[2] - yz_interpolation_values_of_psi_substance_density[0]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
        + ((yz_interpolation_values_of_psi_substance_density[2] - 2.0 * yz_interpolation_values_of_psi_substance_density[1] + yz_interpolation_values_of_psi_substance_density[0]) * (x_distance_from_face_to_position[1]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);
    x_candidate_interpolants[1]
        = yz_interpolation_values_of_psi_substance_density[1]
        + ((-yz_interpolation_values_of_psi_substance_density[3] + 4.0 * yz_interpolation_values_of_psi_substance_density[2] - 3.0 * yz_interpolation_values_of_psi_substance_density[1]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
        + ((yz_interpolation_values_of_psi_substance_density[3] - 2.0 * yz_interpolation_values_of_psi_substance_density[2] + yz_interpolation_values_of_psi_substance_density[1]) * (x_distance_from_face_to_position[1]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);

    //ideal weights の計算
    MY_FLOAT_TYPE x_ideal_weights[2];
    x_ideal_weights[0]
        = -x_distance_from_face_to_position[3] / (3.0 * all_grid._cell_length);
    x_ideal_weights[1]
        = x_distance_from_face_to_position[0] / (3.0 * all_grid._cell_length);

    // smoothness indicatior の計算
    MY_FLOAT_TYPE x_smoothness_indicators[2];
    x_smoothness_indicators[0]
        = (26.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[0]
         - 52.0 * yz_interpolation_values_of_psi_substance_density[1] * yz_interpolation_values_of_psi_substance_density[0]
         - 76.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[1]
         + 25.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[2]
         + 64.0 * yz_interpolation_values_of_psi_substance_density[1] * yz_interpolation_values_of_psi_substance_density[1]
         + 13.0 * yz_interpolation_values_of_psi_substance_density[0] * yz_interpolation_values_of_psi_substance_density[0]) / 12.0;
    x_smoothness_indicators[1]
        = (26.0 * yz_interpolation_values_of_psi_substance_density[3] * yz_interpolation_values_of_psi_substance_density[1]
         - 52.0 * yz_interpolation_values_of_psi_substance_density[3] * yz_interpolation_values_of_psi_substance_density[2]
         - 76.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[1]
         + 25.0 * yz_interpolation_values_of_psi_substance_density[1] * yz_interpolation_values_of_psi_substance_density[1]
         + 64.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[2]
         + 13.0 * yz_interpolation_values_of_psi_substance_density[3] * yz_interpolation_values_of_psi_substance_density[3]) / 12.0;

    //weight を作るための1時的な変数
    MY_FLOAT_TYPE x_alpha[2];
    for (int i = 0; i < 2; ++i) {
        MY_FLOAT_TYPE eps = 0.000001;
        x_alpha[i] = x_ideal_weights[i] / ((eps + x_smoothness_indicators[i]) * (eps + x_smoothness_indicators[i]));
    }
    //weights の計算
    MY_FLOAT_TYPE x_weights[2];
    for (int i = 0; i < 2; ++i) {
        x_weights[i] = x_alpha[i] / (x_alpha[0] + x_alpha[1]);
    }
    return x_weights[0] * x_candidate_interpolants[0] + x_weights[1] * x_candidate_interpolants[1];
}

MY_FLOAT_TYPE WENO4_interpolation_3D_cell_face_y_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_y_values,
    const Grid_3D &all_grid
){
    int advected_index_x = floor((position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_y = floor((position[1]) / all_grid._cell_length);
    int advected_index_z = floor((position[2] - 0.5 * all_grid._cell_length) / all_grid._cell_length);

    /////linear interpolation の処理
    //z方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE z_interpolation_values_of_psi_substance_density[4][4];
    //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    for (int ix = 0; ix < 4; ix++) {
        for(int iy = 0; iy < 4; ++iy){
            //補間に使う離散値をセット
            MY_FLOAT_TYPE psi_val[4];
            for (int iz = 0; iz < 4; ++iz) {
                int index_x = advected_index_x - 1 + ix;
                int index_y = advected_index_y - 1 + iy;
                int index_z = advected_index_z - 1 + iz;
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
            // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
            MY_FLOAT_TYPE z_distance_from_face_to_position[4];
            for (int i = 0; i < 4; ++i) {
                z_distance_from_face_to_position[i] = position[2] - (advected_index_z - 1 + i + 0.5) * all_grid._cell_length;
            }

            // candidate interpolants の計算
            MY_FLOAT_TYPE candidate_interpolants[2];
            candidate_interpolants[0]
                = psi_val[1]
                + ((psi_val[2] - psi_val[0]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
                + ((psi_val[2] - 2.0 * psi_val[1] + psi_val[0]) * (z_distance_from_face_to_position[1]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);
            candidate_interpolants[1]
                = psi_val[1]
                + ((-psi_val[3] + 4.0 * psi_val[2] - 3.0 * psi_val[1]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
                + ((psi_val[3] - 2.0 * psi_val[2] + psi_val[1]) * (z_distance_from_face_to_position[1]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);

            //ideal weights の計算
            MY_FLOAT_TYPE ideal_weights[2];
            ideal_weights[0]
                = -z_distance_from_face_to_position[3] / (3.0 * all_grid._cell_length);
            ideal_weights[1]
                = z_distance_from_face_to_position[0] / (3.0 * all_grid._cell_length);

            // smoothness indicatior の計算
            MY_FLOAT_TYPE smoothness_indicators[2];
            smoothness_indicators[0]
                = (26.0 * psi_val[2] * psi_val[0] - 52.0 * psi_val[1] * psi_val[0] - 76.0 * psi_val[2] * psi_val[1] + 25.0 * psi_val[2] * psi_val[2] + 64.0 * psi_val[1] * psi_val[1] + 13.0 * psi_val[0] * psi_val[0]) / 12.0;
            smoothness_indicators[1]
                = (26.0 * psi_val[3] * psi_val[1] - 52.0 * psi_val[3] * psi_val[2] - 76.0 * psi_val[2] * psi_val[1] + 25.0 * psi_val[1] * psi_val[1] + 64.0 * psi_val[2] * psi_val[2] + 13.0 * psi_val[3] * psi_val[3]) / 12.0;

            //weight を作るための1時的な変数
            MY_FLOAT_TYPE alpha[2];
            for (int i = 0; i < 2; ++i) {
                MY_FLOAT_TYPE eps = 0.000001;
                alpha[i] = ideal_weights[i] / ((eps + smoothness_indicators[i]) * (eps + smoothness_indicators[i]));
            }
            //weights の計算
            MY_FLOAT_TYPE weights[2];
            for (int i = 0; i < 2; ++i) {
                weights[i] = alpha[i] / (alpha[0] + alpha[1]);
            }
            z_interpolation_values_of_psi_substance_density[ix][iy] = weights[0] * candidate_interpolants[0] + weights[1] * candidate_interpolants[1];
        }
    }

    ////ここからy方向の補間
    //yz方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE yz_interpolation_values_of_psi_substance_density[4];
    for(int ix = 0; ix < 4; ++ix){
        // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
        MY_FLOAT_TYPE y_distance_from_face_to_position[4];
        for (int iy = 0; iy < 4; ++iy) {
            // (advected_index_x - 2 + ix + 0.5) の 0.5 は考えているfaceがy軸に垂直なものを考えているため、面の端から中心までのオフセット分の距離
            y_distance_from_face_to_position[iy] = position[1] - (advected_index_y - 1 + iy) * all_grid._cell_length;
        }
        // candidate interpolants の計算
        MY_FLOAT_TYPE y_candidate_interpolants[2];
        y_candidate_interpolants[0]
            = z_interpolation_values_of_psi_substance_density[ix][1]
            + ((z_interpolation_values_of_psi_substance_density[ix][2] - z_interpolation_values_of_psi_substance_density[ix][0]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
            + ((z_interpolation_values_of_psi_substance_density[ix][2] - 2.0 * z_interpolation_values_of_psi_substance_density[ix][1] + z_interpolation_values_of_psi_substance_density[ix][0]) * (y_distance_from_face_to_position[1]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);
        y_candidate_interpolants[1]
            = z_interpolation_values_of_psi_substance_density[ix][1]
            + ((-z_interpolation_values_of_psi_substance_density[ix][3] + 4.0 * z_interpolation_values_of_psi_substance_density[ix][2] - 3.0 * z_interpolation_values_of_psi_substance_density[ix][1]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
            + ((z_interpolation_values_of_psi_substance_density[ix][3] - 2.0 * z_interpolation_values_of_psi_substance_density[ix][2] + z_interpolation_values_of_psi_substance_density[ix][1]) * (y_distance_from_face_to_position[1]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);

        //ideal weights の計算
        MY_FLOAT_TYPE y_ideal_weights[2];
        y_ideal_weights[0]
            = -y_distance_from_face_to_position[3] / (3.0 * all_grid._cell_length);
        y_ideal_weights[1]
            = y_distance_from_face_to_position[0] / (3.0 * all_grid._cell_length);

        // smoothness indicatior の計算
        MY_FLOAT_TYPE y_smoothness_indicators[2];
        y_smoothness_indicators[0]
            = (26.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][0]
             - 52.0 * z_interpolation_values_of_psi_substance_density[ix][1] * z_interpolation_values_of_psi_substance_density[ix][0]
             - 76.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 25.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][2]
             + 64.0 * z_interpolation_values_of_psi_substance_density[ix][1] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 13.0 * z_interpolation_values_of_psi_substance_density[ix][0] * z_interpolation_values_of_psi_substance_density[ix][0]) / 12.0;
        y_smoothness_indicators[1]
            = (26.0 * z_interpolation_values_of_psi_substance_density[ix][3] * z_interpolation_values_of_psi_substance_density[ix][1]
             - 52.0 * z_interpolation_values_of_psi_substance_density[ix][3] * z_interpolation_values_of_psi_substance_density[ix][2]
             - 76.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 25.0 * z_interpolation_values_of_psi_substance_density[ix][1] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 64.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][2]
             + 13.0 * z_interpolation_values_of_psi_substance_density[ix][3] * z_interpolation_values_of_psi_substance_density[ix][3]) / 12.0;

        //weight を作るための1時的な変数
        MY_FLOAT_TYPE y_alpha[2];
        for (int i = 0; i < 2; ++i) {
            MY_FLOAT_TYPE eps = 0.000001;
            y_alpha[i] = y_ideal_weights[i] / ((eps + y_smoothness_indicators[i]) * (eps + y_smoothness_indicators[i]));
        }
        //weights の計算
        MY_FLOAT_TYPE y_weights[2];
        for (int i = 0; i < 2; ++i) {
            y_weights[i] = y_alpha[i] / (y_alpha[0] + y_alpha[1]);
        }
        yz_interpolation_values_of_psi_substance_density[ix] = y_weights[0] * y_candidate_interpolants[0] + y_weights[1] * y_candidate_interpolants[1];
    }

    ////ここからx方向の補間
    // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
    MY_FLOAT_TYPE x_distance_from_face_to_position[4];
    for (int ix = 0; ix < 4; ++ix) {
        // (advected_index_x - 2 + ix + 0.5) の 0.5 は考えているfaceがy軸に垂直なものを考えているため、面の端から中心までのオフセット分の距離
        x_distance_from_face_to_position[ix] = position[0] - (advected_index_x - 1 + ix + 0.5) * all_grid._cell_length;
    }

    // candidate interpolants の計算
    MY_FLOAT_TYPE x_candidate_interpolants[2];
    x_candidate_interpolants[0]
        = yz_interpolation_values_of_psi_substance_density[1]
        + ((yz_interpolation_values_of_psi_substance_density[2] - yz_interpolation_values_of_psi_substance_density[0]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
        + ((yz_interpolation_values_of_psi_substance_density[2] - 2.0 * yz_interpolation_values_of_psi_substance_density[1] + yz_interpolation_values_of_psi_substance_density[0]) * (x_distance_from_face_to_position[1]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);
    x_candidate_interpolants[1]
        = yz_interpolation_values_of_psi_substance_density[1]
        + ((-yz_interpolation_values_of_psi_substance_density[3] + 4.0 * yz_interpolation_values_of_psi_substance_density[2] - 3.0 * yz_interpolation_values_of_psi_substance_density[1]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
        + ((yz_interpolation_values_of_psi_substance_density[3] - 2.0 * yz_interpolation_values_of_psi_substance_density[2] + yz_interpolation_values_of_psi_substance_density[1]) * (x_distance_from_face_to_position[1]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);

    //ideal weights の計算
    MY_FLOAT_TYPE x_ideal_weights[2];
    x_ideal_weights[0]
        = -x_distance_from_face_to_position[3] / (3.0 * all_grid._cell_length);
    x_ideal_weights[1]
        = x_distance_from_face_to_position[0] / (3.0 * all_grid._cell_length);

    // smoothness indicatior の計算
    MY_FLOAT_TYPE x_smoothness_indicators[2];
    x_smoothness_indicators[0]
        = (26.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[0]
         - 52.0 * yz_interpolation_values_of_psi_substance_density[1] * yz_interpolation_values_of_psi_substance_density[0]
         - 76.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[1]
         + 25.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[2]
         + 64.0 * yz_interpolation_values_of_psi_substance_density[1] * yz_interpolation_values_of_psi_substance_density[1]
         + 13.0 * yz_interpolation_values_of_psi_substance_density[0] * yz_interpolation_values_of_psi_substance_density[0]) / 12.0;
    x_smoothness_indicators[1]
        = (26.0 * yz_interpolation_values_of_psi_substance_density[3] * yz_interpolation_values_of_psi_substance_density[1]
         - 52.0 * yz_interpolation_values_of_psi_substance_density[3] * yz_interpolation_values_of_psi_substance_density[2]
         - 76.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[1]
         + 25.0 * yz_interpolation_values_of_psi_substance_density[1] * yz_interpolation_values_of_psi_substance_density[1]
         + 64.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[2]
         + 13.0 * yz_interpolation_values_of_psi_substance_density[3] * yz_interpolation_values_of_psi_substance_density[3]) / 12.0;

    //weight を作るための1時的な変数
    MY_FLOAT_TYPE x_alpha[2];
    for (int i = 0; i < 2; ++i) {
        MY_FLOAT_TYPE eps = 0.000001;
        x_alpha[i] = x_ideal_weights[i] / ((eps + x_smoothness_indicators[i]) * (eps + x_smoothness_indicators[i]));
    }
    //weights の計算
    MY_FLOAT_TYPE x_weights[2];
    for (int i = 0; i < 2; ++i) {
        x_weights[i] = x_alpha[i] / (x_alpha[0] + x_alpha[1]);
    }
    return x_weights[0] * x_candidate_interpolants[0] + x_weights[1] * x_candidate_interpolants[1];
}

MY_FLOAT_TYPE WENO4_interpolation_3D_cell_face_z_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_z_values,
    const Grid_3D &all_grid
){
    int advected_index_x = floor((position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_z = floor((position[2]) / all_grid._cell_length);

    /////linear interpolation の処理
    //z方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE z_interpolation_values_of_psi_substance_density[4][4];
    //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    for (int ix = 0; ix < 4; ix++) {
        for(int iy = 0; iy < 4; ++iy){
            //補間に使う離散値をセット
            MY_FLOAT_TYPE psi_val[4];
            for (int iz = 0; iz < 4; ++iz) {
                int index_x = advected_index_x - 1 + ix;
                int index_y = advected_index_y - 1 + iy;
                int index_z = advected_index_z - 1 + iz;
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
            // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
            MY_FLOAT_TYPE z_distance_from_face_to_position[4];
            for (int i = 0; i < 4; ++i) {
                z_distance_from_face_to_position[i] = position[2] - (advected_index_z - 1 + i) * all_grid._cell_length;
            }

            // candidate interpolants の計算
            MY_FLOAT_TYPE candidate_interpolants[2];
            candidate_interpolants[0]
                = psi_val[1]
                + ((psi_val[2] - psi_val[0]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
                + ((psi_val[2] - 2.0 * psi_val[1] + psi_val[0]) * (z_distance_from_face_to_position[1]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);
            candidate_interpolants[1]
                = psi_val[1]
                + ((-psi_val[3] + 4.0 * psi_val[2] - 3.0 * psi_val[1]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
                + ((psi_val[3] - 2.0 * psi_val[2] + psi_val[1]) * (z_distance_from_face_to_position[1]) * (z_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);

            //ideal weights の計算
            MY_FLOAT_TYPE ideal_weights[2];
            ideal_weights[0]
                = -z_distance_from_face_to_position[3] / (3.0 * all_grid._cell_length);
            ideal_weights[1]
                = z_distance_from_face_to_position[0] / (3.0 * all_grid._cell_length);

            // smoothness indicatior の計算
            MY_FLOAT_TYPE smoothness_indicators[2];
            smoothness_indicators[0]
                = (26.0 * psi_val[2] * psi_val[0] - 52.0 * psi_val[1] * psi_val[0] - 76.0 * psi_val[2] * psi_val[1] + 25.0 * psi_val[2] * psi_val[2] + 64.0 * psi_val[1] * psi_val[1] + 13.0 * psi_val[0] * psi_val[0]) / 12.0;
            smoothness_indicators[1]
                = (26.0 * psi_val[3] * psi_val[1] - 52.0 * psi_val[3] * psi_val[2] - 76.0 * psi_val[2] * psi_val[1] + 25.0 * psi_val[1] * psi_val[1] + 64.0 * psi_val[2] * psi_val[2] + 13.0 * psi_val[3] * psi_val[3]) / 12.0;

            //weight を作るための1時的な変数
            MY_FLOAT_TYPE alpha[2];
            for (int i = 0; i < 2; ++i) {
                MY_FLOAT_TYPE eps = 0.000001;
                alpha[i] = ideal_weights[i] / ((eps + smoothness_indicators[i]) * (eps + smoothness_indicators[i]));
            }
            //weights の計算
            MY_FLOAT_TYPE weights[2];
            for (int i = 0; i < 2; ++i) {
                weights[i] = alpha[i] / (alpha[0] + alpha[1]);
            }
            z_interpolation_values_of_psi_substance_density[ix][iy] = weights[0] * candidate_interpolants[0] + weights[1] * candidate_interpolants[1];
        }
    }

    ////ここからy方向の補間
    //yz方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE yz_interpolation_values_of_psi_substance_density[4];
    for(int ix = 0; ix < 4; ++ix){
        // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
        MY_FLOAT_TYPE y_distance_from_face_to_position[4];
        for (int iy = 0; iy < 4; ++iy) {
            // (advected_index_x - 2 + ix + 0.5) の 0.5 は考えているfaceがy軸に垂直なものを考えているため、面の端から中心までのオフセット分の距離
            y_distance_from_face_to_position[iy] = position[1] - (advected_index_y - 1 + iy + 0.5) * all_grid._cell_length;
        }
        // candidate interpolants の計算
        MY_FLOAT_TYPE y_candidate_interpolants[2];
        y_candidate_interpolants[0]
            = z_interpolation_values_of_psi_substance_density[ix][1]
            + ((z_interpolation_values_of_psi_substance_density[ix][2] - z_interpolation_values_of_psi_substance_density[ix][0]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
            + ((z_interpolation_values_of_psi_substance_density[ix][2] - 2.0 * z_interpolation_values_of_psi_substance_density[ix][1] + z_interpolation_values_of_psi_substance_density[ix][0]) * (y_distance_from_face_to_position[1]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);
        y_candidate_interpolants[1]
            = z_interpolation_values_of_psi_substance_density[ix][1]
            + ((-z_interpolation_values_of_psi_substance_density[ix][3] + 4.0 * z_interpolation_values_of_psi_substance_density[ix][2] - 3.0 * z_interpolation_values_of_psi_substance_density[ix][1]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
            + ((z_interpolation_values_of_psi_substance_density[ix][3] - 2.0 * z_interpolation_values_of_psi_substance_density[ix][2] + z_interpolation_values_of_psi_substance_density[ix][1]) * (y_distance_from_face_to_position[1]) * (y_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);

        //ideal weights の計算
        MY_FLOAT_TYPE y_ideal_weights[2];
        y_ideal_weights[0]
            = -y_distance_from_face_to_position[3] / (3.0 * all_grid._cell_length);
        y_ideal_weights[1]
            = y_distance_from_face_to_position[0] / (3.0 * all_grid._cell_length);

        // smoothness indicatior の計算
        MY_FLOAT_TYPE y_smoothness_indicators[2];
        y_smoothness_indicators[0]
            = (26.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][0]
             - 52.0 * z_interpolation_values_of_psi_substance_density[ix][1] * z_interpolation_values_of_psi_substance_density[ix][0]
             - 76.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 25.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][2]
             + 64.0 * z_interpolation_values_of_psi_substance_density[ix][1] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 13.0 * z_interpolation_values_of_psi_substance_density[ix][0] * z_interpolation_values_of_psi_substance_density[ix][0]) / 12.0;
        y_smoothness_indicators[1]
            = (26.0 * z_interpolation_values_of_psi_substance_density[ix][3] * z_interpolation_values_of_psi_substance_density[ix][1]
             - 52.0 * z_interpolation_values_of_psi_substance_density[ix][3] * z_interpolation_values_of_psi_substance_density[ix][2]
             - 76.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 25.0 * z_interpolation_values_of_psi_substance_density[ix][1] * z_interpolation_values_of_psi_substance_density[ix][1]
             + 64.0 * z_interpolation_values_of_psi_substance_density[ix][2] * z_interpolation_values_of_psi_substance_density[ix][2]
             + 13.0 * z_interpolation_values_of_psi_substance_density[ix][3] * z_interpolation_values_of_psi_substance_density[ix][3]) / 12.0;

        //weight を作るための1時的な変数
        MY_FLOAT_TYPE y_alpha[2];
        for (int i = 0; i < 2; ++i) {
            MY_FLOAT_TYPE eps = 0.000001;
            y_alpha[i] = y_ideal_weights[i] / ((eps + y_smoothness_indicators[i]) * (eps + y_smoothness_indicators[i]));
        }
        //weights の計算
        MY_FLOAT_TYPE y_weights[2];
        for (int i = 0; i < 2; ++i) {
            y_weights[i] = y_alpha[i] / (y_alpha[0] + y_alpha[1]);
        }
        yz_interpolation_values_of_psi_substance_density[ix] = y_weights[0] * y_candidate_interpolants[0] + y_weights[1] * y_candidate_interpolants[1];
    }

    ////ここからx方向の補間
    // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
    MY_FLOAT_TYPE x_distance_from_face_to_position[4];
    for (int ix = 0; ix < 4; ++ix) {
        // (advected_index_x - 2 + ix + 0.5) の 0.5 は考えているfaceがy軸に垂直なものを考えているため、面の端から中心までのオフセット分の距離
        x_distance_from_face_to_position[ix] = position[0] - (advected_index_x - 1 + ix + 0.5) * all_grid._cell_length;
    }

    // candidate interpolants の計算
    MY_FLOAT_TYPE x_candidate_interpolants[2];
    x_candidate_interpolants[0]
        = yz_interpolation_values_of_psi_substance_density[1]
        + ((yz_interpolation_values_of_psi_substance_density[2] - yz_interpolation_values_of_psi_substance_density[0]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
        + ((yz_interpolation_values_of_psi_substance_density[2] - 2.0 * yz_interpolation_values_of_psi_substance_density[1] + yz_interpolation_values_of_psi_substance_density[0]) * (x_distance_from_face_to_position[1]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);
    x_candidate_interpolants[1]
        = yz_interpolation_values_of_psi_substance_density[1]
        + ((-yz_interpolation_values_of_psi_substance_density[3] + 4.0 * yz_interpolation_values_of_psi_substance_density[2] - 3.0 * yz_interpolation_values_of_psi_substance_density[1]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length)
        + ((yz_interpolation_values_of_psi_substance_density[3] - 2.0 * yz_interpolation_values_of_psi_substance_density[2] + yz_interpolation_values_of_psi_substance_density[1]) * (x_distance_from_face_to_position[1]) * (x_distance_from_face_to_position[1])) / (2.0 * all_grid._cell_length * all_grid._cell_length);

    //ideal weights の計算
    MY_FLOAT_TYPE x_ideal_weights[2];
    x_ideal_weights[0]
        = -x_distance_from_face_to_position[3] / (3.0 * all_grid._cell_length);
    x_ideal_weights[1]
        = x_distance_from_face_to_position[0] / (3.0 * all_grid._cell_length);

    // smoothness indicatior の計算
    MY_FLOAT_TYPE x_smoothness_indicators[2];
    x_smoothness_indicators[0]
        = (26.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[0]
         - 52.0 * yz_interpolation_values_of_psi_substance_density[1] * yz_interpolation_values_of_psi_substance_density[0]
         - 76.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[1]
         + 25.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[2]
         + 64.0 * yz_interpolation_values_of_psi_substance_density[1] * yz_interpolation_values_of_psi_substance_density[1]
         + 13.0 * yz_interpolation_values_of_psi_substance_density[0] * yz_interpolation_values_of_psi_substance_density[0]) / 12.0;
    x_smoothness_indicators[1]
        = (26.0 * yz_interpolation_values_of_psi_substance_density[3] * yz_interpolation_values_of_psi_substance_density[1]
         - 52.0 * yz_interpolation_values_of_psi_substance_density[3] * yz_interpolation_values_of_psi_substance_density[2]
         - 76.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[1]
         + 25.0 * yz_interpolation_values_of_psi_substance_density[1] * yz_interpolation_values_of_psi_substance_density[1]
         + 64.0 * yz_interpolation_values_of_psi_substance_density[2] * yz_interpolation_values_of_psi_substance_density[2]
         + 13.0 * yz_interpolation_values_of_psi_substance_density[3] * yz_interpolation_values_of_psi_substance_density[3]) / 12.0;

    //weight を作るための1時的な変数
    MY_FLOAT_TYPE x_alpha[2];
    for (int i = 0; i < 2; ++i) {
        MY_FLOAT_TYPE eps = 0.000001;
        x_alpha[i] = x_ideal_weights[i] / ((eps + x_smoothness_indicators[i]) * (eps + x_smoothness_indicators[i]));
    }
    //weights の計算
    MY_FLOAT_TYPE x_weights[2];
    for (int i = 0; i < 2; ++i) {
        x_weights[i] = x_alpha[i] / (x_alpha[0] + x_alpha[1]);
    }
    return x_weights[0] * x_candidate_interpolants[0] + x_weights[1] * x_candidate_interpolants[1];
}

}//namespace smoke_simulation
