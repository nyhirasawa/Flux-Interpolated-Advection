#include "WENO6_interpolation_2d.h"

#include <vector>
#include "utils.h"

namespace smoke_simulation{
MY_FLOAT_TYPE WENO6_interpolation_2D_cell_center_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_center_values,
    const Grid &all_grid
){
    int advected_index_x = floor((position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);

    //WENO6 interpolation の処理
    //y方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE y_interpolation_values_of_psi_substance_density[6];
    //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    for (int ix = 0; ix < 6; ix++) {
        MY_FLOAT_TYPE psi_y_val[6];
        for (int iy = 0; iy < 6; ++iy) {
            int index_x = advected_index_x - 2 + ix;
            int index_y = advected_index_y - 2 + iy;
            if (index_x < 0) {
                index_x = 0;
            }
            if (index_x > all_grid.Grid_num_x - 1) {
                int exceed = index_x - (all_grid.Grid_num_x - 1);
                //境界の値がずっと外側まで続く場合
                index_x = all_grid.Grid_num_x - 1;
                //境界を境に鏡のように値が反射する場合
                //index_x = all_grid.Grid_num_x - 1 - exceed + 1;
            }
            if (index_y < 0) {
                index_y = 0;
            }
            if (index_y > all_grid.Grid_num_y - 1) {
                int exceed = index_y - (all_grid.Grid_num_y - 1);
                //境界の値がずっと外側まで続く場合
                index_y = all_grid.Grid_num_y - 1;
                //境界を境に鏡のように値が反射する場合
                //index_y = all_grid.Grid_num_y - 1 - exceed + 1;
            }
            psi_y_val[iy] = cell_center_values[get_voxel_center_index(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)];
        }

        // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
        MY_FLOAT_TYPE y_distance_from_face_to_position[6];
        for (int i = 0; i < 6; ++i) {
            // (advected_index_y - 2 + i + 0.5) の 0.5 は面の端から中心までのオフセット分の距離
            y_distance_from_face_to_position[i] = position[1] - (advected_index_y - 2 + i + 0.5) * all_grid._cell_length;
        }

        // candidate interpolants の計算
        MY_FLOAT_TYPE candidate_interpolants[3];
        for (int i = 0; i < 3; ++i) {
            candidate_interpolants[i]
                = psi_y_val[i]
                + (psi_y_val[i + 1] - psi_y_val[i]) * y_distance_from_face_to_position[i] / (all_grid._cell_length)
                + (psi_y_val[i + 2] - 2 * psi_y_val[i + 1] + psi_y_val[i]) * y_distance_from_face_to_position[i] * y_distance_from_face_to_position[i + 1]
                / (2 * all_grid._cell_length * all_grid._cell_length)
                + (psi_y_val[i + 3] - 3 * psi_y_val[i + 2] + 3 * psi_y_val[i + 1] - psi_y_val[i])
                * y_distance_from_face_to_position[i] * y_distance_from_face_to_position[i + 1] * y_distance_from_face_to_position[i + 2]
                / (6 * all_grid._cell_length * all_grid._cell_length * all_grid._cell_length);
        }

        //ideal weights の計算
        MY_FLOAT_TYPE ideal_weights[3];
        ideal_weights[0]
            = y_distance_from_face_to_position[4] * y_distance_from_face_to_position[5]
            / (20 * all_grid._cell_length * all_grid._cell_length);
        ideal_weights[1]
            = -y_distance_from_face_to_position[5] * y_distance_from_face_to_position[0]
            / (10 * all_grid._cell_length * all_grid._cell_length);
        ideal_weights[2]
            = y_distance_from_face_to_position[0] * y_distance_from_face_to_position[1]
            / (20 * all_grid._cell_length * all_grid._cell_length);

        // smoothness indicatior の計算
        MY_FLOAT_TYPE smoothness_indicators[3];
        smoothness_indicators[0]
            = (-3579 * psi_y_val[3] * psi_y_val[2] + 2634 * psi_y_val[3] * psi_y_val[1] - 683 * psi_y_val[3] * psi_y_val[0]
                - 6927 * psi_y_val[2] * psi_y_val[1] + 1854 * psi_y_val[2] * psi_y_val[0] - 1659 * psi_y_val[1] * psi_y_val[0]
                + 814 * psi_y_val[3] * psi_y_val[3] + 4326 * psi_y_val[2] * psi_y_val[2] + 2976 * psi_y_val[1] * psi_y_val[1]
                + 244 * psi_y_val[0] * psi_y_val[0]) / 180.0;
        smoothness_indicators[1]
            = (-3777 * psi_y_val[3] * psi_y_val[2] + 1074 * psi_y_val[3] * psi_y_val[1] - 1269 * psi_y_val[2] * psi_y_val[1]
                + 1986 * psi_y_val[3] * psi_y_val[3] + 1986 * psi_y_val[2] * psi_y_val[2] + 244 * psi_y_val[1] * psi_y_val[1]
                + 244 * psi_y_val[4] * psi_y_val[4] - 1269 * psi_y_val[4] * psi_y_val[3] + 1074 * psi_y_val[4] * psi_y_val[2]
                - 293 * psi_y_val[4] * psi_y_val[1]) / 180.0;
        smoothness_indicators[2]
            = (-3579 * psi_y_val[3] * psi_y_val[2] + 4326 * psi_y_val[3] * psi_y_val[3] + 814 * psi_y_val[2] * psi_y_val[2]
                + 2976 * psi_y_val[4] * psi_y_val[4] + 244 * psi_y_val[5] * psi_y_val[5] - 683 * psi_y_val[5] * psi_y_val[2]
                - 6927 * psi_y_val[4] * psi_y_val[3] + 2634 * psi_y_val[4] * psi_y_val[2] - 1659 * psi_y_val[5] * psi_y_val[4]
                + 1854 * psi_y_val[5] * psi_y_val[3]) / 180.0;

        //weight を作るための1時的な変数
        MY_FLOAT_TYPE alpha[3];
        for (int i = 0; i < 3; ++i) {
            MY_FLOAT_TYPE eps = 0.000001;
            alpha[i] = ideal_weights[i] / ((eps + smoothness_indicators[i]) * (eps + smoothness_indicators[i]));
        }
        //weights の計算
        MY_FLOAT_TYPE weights[3];
        for (int i = 0; i < 3; ++i) {
            weights[i] = alpha[i] / (alpha[0] + alpha[1] + alpha[2]);
        }
        //return weights[0] * candidate_interpolants[0] + weights[1] * candidate_interpolants[1] + weights[2] * candidate_interpolants[2];
        y_interpolation_values_of_psi_substance_density[ix] = weights[0] * candidate_interpolants[0] + weights[1] * candidate_interpolants[1] + weights[2] * candidate_interpolants[2];
    }

    ////ここからx方向の補間
    // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
    MY_FLOAT_TYPE x_distance_from_face_to_position[6];
    for (int ix = 0; ix < 6; ++ix) {
        // (advected_index_x - 2 + ix + 0.5) の 0.5 は考えているfaceがy軸に垂直なものを考えているため、面の端から中心までのオフセット分の距離
        x_distance_from_face_to_position[ix] = position[0] - (advected_index_x - 2 + ix + 0.5) * all_grid._cell_length;
    }

    // candidate interpolants の計算
    MY_FLOAT_TYPE x_candidate_interpolants[3];
    for (int i = 0; i < 3; ++i) {
        x_candidate_interpolants[i]
            = y_interpolation_values_of_psi_substance_density[i]
            + (y_interpolation_values_of_psi_substance_density[i + 1] - y_interpolation_values_of_psi_substance_density[i]) * x_distance_from_face_to_position[i] / (all_grid._cell_length)
            + (y_interpolation_values_of_psi_substance_density[i + 2] - 2 * y_interpolation_values_of_psi_substance_density[i + 1] + y_interpolation_values_of_psi_substance_density[i]) * x_distance_from_face_to_position[i] * x_distance_from_face_to_position[i + 1]
            / (2 * all_grid._cell_length * all_grid._cell_length)
            + (y_interpolation_values_of_psi_substance_density[i + 3] - 3 * y_interpolation_values_of_psi_substance_density[i + 2] + 3 * y_interpolation_values_of_psi_substance_density[i + 1] - y_interpolation_values_of_psi_substance_density[i])
            * x_distance_from_face_to_position[i] * x_distance_from_face_to_position[i + 1] * x_distance_from_face_to_position[i + 2]
            / (6 * all_grid._cell_length * all_grid._cell_length * all_grid._cell_length);
    }

    //ideal weights の計算
    MY_FLOAT_TYPE x_ideal_weights[3];
    x_ideal_weights[0]
        = x_distance_from_face_to_position[4] * x_distance_from_face_to_position[5]
        / (20 * all_grid._cell_length * all_grid._cell_length);
    x_ideal_weights[1]
        = -x_distance_from_face_to_position[5] * x_distance_from_face_to_position[0]
        / (10 * all_grid._cell_length * all_grid._cell_length);
    x_ideal_weights[2]
        = x_distance_from_face_to_position[0] * x_distance_from_face_to_position[1]
        / (20 * all_grid._cell_length * all_grid._cell_length);

    // smoothness indicatior の計算
    MY_FLOAT_TYPE x_smoothness_indicators[3];
    x_smoothness_indicators[0]
        = (-3579 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[2] + 2634 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[1] - 683 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[0]
            - 6927 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[1] + 1854 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[0] - 1659 * y_interpolation_values_of_psi_substance_density[1] * y_interpolation_values_of_psi_substance_density[0]
            + 814 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[3] + 4326 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[2] + 2976 * y_interpolation_values_of_psi_substance_density[1] * y_interpolation_values_of_psi_substance_density[1]
            + 244 * y_interpolation_values_of_psi_substance_density[0] * y_interpolation_values_of_psi_substance_density[0]) / 180.0;
    x_smoothness_indicators[1]
        = (-3777 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[2] + 1074 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[1] - 1269 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[1]
            + 1986 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[3] + 1986 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[2] + 244 * y_interpolation_values_of_psi_substance_density[1] * y_interpolation_values_of_psi_substance_density[1]
            + 244 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[4] - 1269 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[3] + 1074 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[2]
            - 293 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[1]) / 180.0;
    x_smoothness_indicators[2]
        = (-3579 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[2] + 4326 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[3] + 814 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[2]
            + 2976 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[4] + 244 * y_interpolation_values_of_psi_substance_density[5] * y_interpolation_values_of_psi_substance_density[5] - 683 * y_interpolation_values_of_psi_substance_density[5] * y_interpolation_values_of_psi_substance_density[2]
            - 6927 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[3] + 2634 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[2] - 1659 * y_interpolation_values_of_psi_substance_density[5] * y_interpolation_values_of_psi_substance_density[4]
            + 1854 * y_interpolation_values_of_psi_substance_density[5] * y_interpolation_values_of_psi_substance_density[3]) / 180.0;

    //weight を作るための1時的な変数
    MY_FLOAT_TYPE x_alpha[3];
    for (int i = 0; i < 3; ++i) {
        MY_FLOAT_TYPE eps = 0.000001;
        x_alpha[i] = x_ideal_weights[i] / ((eps + x_smoothness_indicators[i]) * (eps + x_smoothness_indicators[i]));
    }
    //weights の計算
    MY_FLOAT_TYPE x_weights[3];
    for (int i = 0; i < 3; ++i) {
        x_weights[i] = x_alpha[i] / (x_alpha[0] + x_alpha[1] + x_alpha[2]);
    }
    return x_weights[0] * x_candidate_interpolants[0] + x_weights[1] * x_candidate_interpolants[1] + x_weights[2] * x_candidate_interpolants[2];
}

MY_FLOAT_TYPE WENO6_interpolation_2D_cell_face_x_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
    const Grid &all_grid
){
    int advected_index_x = floor((position[0]) / all_grid._cell_length);
    int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);

    /////WENO6 interpolation の処理
    //y方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE y_interpolation_values_of_psi_substance_density[6];
    //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    for (int ix = 0; ix < 6; ix++) {
        MY_FLOAT_TYPE psi_y_val[6];
        for (int iy = 0; iy < 6; ++iy) {
            int index_x = advected_index_x - 2 + ix;
            int index_y = advected_index_y - 2 + iy;
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
            psi_y_val[iy] = cell_face_x_values[get_voxel_face_index_x(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)];
        }

        // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
        MY_FLOAT_TYPE y_distance_from_face_to_position[6];
        for (int i = 0; i < 6; ++i) {
            y_distance_from_face_to_position[i] = position[1] - (advected_index_y - 2 + i + 0.5) * all_grid._cell_length;
        }

        // candidate interpolants の計算
        MY_FLOAT_TYPE candidate_interpolants[3];
        for (int i = 0; i < 3; ++i) {
            candidate_interpolants[i]
                = psi_y_val[i]
                + (psi_y_val[i + 1] - psi_y_val[i]) * y_distance_from_face_to_position[i] / (all_grid._cell_length)
                + (psi_y_val[i + 2] - 2 * psi_y_val[i + 1] + psi_y_val[i]) * y_distance_from_face_to_position[i] * y_distance_from_face_to_position[i + 1]
                / (2 * all_grid._cell_length * all_grid._cell_length)
                + (psi_y_val[i + 3] - 3 * psi_y_val[i + 2] + 3 * psi_y_val[i + 1] - psi_y_val[i])
                * y_distance_from_face_to_position[i] * y_distance_from_face_to_position[i + 1] * y_distance_from_face_to_position[i + 2]
                / (6 * all_grid._cell_length * all_grid._cell_length * all_grid._cell_length);
        }

        //ideal weights の計算
        MY_FLOAT_TYPE ideal_weights[3];
        ideal_weights[0]
            = y_distance_from_face_to_position[4] * y_distance_from_face_to_position[5]
            / (20 * all_grid._cell_length * all_grid._cell_length);
        ideal_weights[1]
            = -y_distance_from_face_to_position[5] * y_distance_from_face_to_position[0]
            / (10 * all_grid._cell_length * all_grid._cell_length);
        ideal_weights[2]
            = y_distance_from_face_to_position[0] * y_distance_from_face_to_position[1]
            / (20 * all_grid._cell_length * all_grid._cell_length);

        // smoothness indicatior の計算
        MY_FLOAT_TYPE smoothness_indicators[3];
        smoothness_indicators[0]
            = (-3579 * psi_y_val[3] * psi_y_val[2] + 2634 * psi_y_val[3] * psi_y_val[1] - 683 * psi_y_val[3] * psi_y_val[0]
                - 6927 * psi_y_val[2] * psi_y_val[1] + 1854 * psi_y_val[2] * psi_y_val[0] - 1659 * psi_y_val[1] * psi_y_val[0]
                + 814 * psi_y_val[3] * psi_y_val[3] + 4326 * psi_y_val[2] * psi_y_val[2] + 2976 * psi_y_val[1] * psi_y_val[1]
                + 244 * psi_y_val[0] * psi_y_val[0]) / 180.0;
        smoothness_indicators[1]
            = (-3777 * psi_y_val[3] * psi_y_val[2] + 1074 * psi_y_val[3] * psi_y_val[1] - 1269 * psi_y_val[2] * psi_y_val[1]
                + 1986 * psi_y_val[3] * psi_y_val[3] + 1986 * psi_y_val[2] * psi_y_val[2] + 244 * psi_y_val[1] * psi_y_val[1]
                + 244 * psi_y_val[4] * psi_y_val[4] - 1269 * psi_y_val[4] * psi_y_val[3] + 1074 * psi_y_val[4] * psi_y_val[2]
                - 293 * psi_y_val[4] * psi_y_val[1]) / 180.0;
        smoothness_indicators[2]
            = (-3579 * psi_y_val[3] * psi_y_val[2] + 4326 * psi_y_val[3] * psi_y_val[3] + 814 * psi_y_val[2] * psi_y_val[2]
                + 2976 * psi_y_val[4] * psi_y_val[4] + 244 * psi_y_val[5] * psi_y_val[5] - 683 * psi_y_val[5] * psi_y_val[2]
                - 6927 * psi_y_val[4] * psi_y_val[3] + 2634 * psi_y_val[4] * psi_y_val[2] - 1659 * psi_y_val[5] * psi_y_val[4]
                + 1854 * psi_y_val[5] * psi_y_val[3]) / 180.0;

        //weight を作るための1時的な変数
        MY_FLOAT_TYPE alpha[3];
        for (int i = 0; i < 3; ++i) {
            MY_FLOAT_TYPE eps = 0.000001;
            alpha[i] = ideal_weights[i] / ((eps + smoothness_indicators[i]) * (eps + smoothness_indicators[i]));
        }
        //weights の計算
        MY_FLOAT_TYPE weights[3];
        for (int i = 0; i < 3; ++i) {
            weights[i] = alpha[i] / (alpha[0] + alpha[1] + alpha[2]);
        }
        //return weights[0] * candidate_interpolants[0] + weights[1] * candidate_interpolants[1] + weights[2] * candidate_interpolants[2];
        y_interpolation_values_of_psi_substance_density[ix] = weights[0] * candidate_interpolants[0] + weights[1] * candidate_interpolants[1] + weights[2] * candidate_interpolants[2];
    }

    ////ここからx方向の補間
    // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
    MY_FLOAT_TYPE x_distance_from_face_to_position[6];
    for (int ix = 0; ix < 6; ++ix) {
        // (advected_index_x - 2 + ix + 0.5) の 0.5 は考えているfaceがy軸に垂直なものを考えているため、面の端から中心までのオフセット分の距離
        x_distance_from_face_to_position[ix] = position[0] - (advected_index_x - 2 + ix) * all_grid._cell_length;
    }

    // candidate interpolants の計算
    MY_FLOAT_TYPE x_candidate_interpolants[3];
    for (int i = 0; i < 3; ++i) {
        x_candidate_interpolants[i]
            = y_interpolation_values_of_psi_substance_density[i]
            + (y_interpolation_values_of_psi_substance_density[i + 1] - y_interpolation_values_of_psi_substance_density[i]) * x_distance_from_face_to_position[i] / (all_grid._cell_length)
            + (y_interpolation_values_of_psi_substance_density[i + 2] - 2 * y_interpolation_values_of_psi_substance_density[i + 1] + y_interpolation_values_of_psi_substance_density[i]) * x_distance_from_face_to_position[i] * x_distance_from_face_to_position[i + 1]
            / (2 * all_grid._cell_length * all_grid._cell_length)
            + (y_interpolation_values_of_psi_substance_density[i + 3] - 3 * y_interpolation_values_of_psi_substance_density[i + 2] + 3 * y_interpolation_values_of_psi_substance_density[i + 1] - y_interpolation_values_of_psi_substance_density[i])
            * x_distance_from_face_to_position[i] * x_distance_from_face_to_position[i + 1] * x_distance_from_face_to_position[i + 2]
            / (6 * all_grid._cell_length * all_grid._cell_length * all_grid._cell_length);
    }

    //ideal weights の計算
    MY_FLOAT_TYPE x_ideal_weights[3];
    x_ideal_weights[0]
        = x_distance_from_face_to_position[4] * x_distance_from_face_to_position[5]
        / (20 * all_grid._cell_length * all_grid._cell_length);
    x_ideal_weights[1]
        = -x_distance_from_face_to_position[5] * x_distance_from_face_to_position[0]
        / (10 * all_grid._cell_length * all_grid._cell_length);
    x_ideal_weights[2]
        = x_distance_from_face_to_position[0] * x_distance_from_face_to_position[1]
        / (20 * all_grid._cell_length * all_grid._cell_length);

    // smoothness indicatior の計算
    MY_FLOAT_TYPE x_smoothness_indicators[3];
    x_smoothness_indicators[0]
        = (-3579 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[2] + 2634 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[1] - 683 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[0]
            - 6927 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[1] + 1854 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[0] - 1659 * y_interpolation_values_of_psi_substance_density[1] * y_interpolation_values_of_psi_substance_density[0]
            + 814 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[3] + 4326 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[2] + 2976 * y_interpolation_values_of_psi_substance_density[1] * y_interpolation_values_of_psi_substance_density[1]
            + 244 * y_interpolation_values_of_psi_substance_density[0] * y_interpolation_values_of_psi_substance_density[0]) / 180.0;
    x_smoothness_indicators[1]
        = (-3777 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[2] + 1074 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[1] - 1269 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[1]
            + 1986 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[3] + 1986 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[2] + 244 * y_interpolation_values_of_psi_substance_density[1] * y_interpolation_values_of_psi_substance_density[1]
            + 244 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[4] - 1269 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[3] + 1074 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[2]
            - 293 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[1]) / 180.0;
    x_smoothness_indicators[2]
        = (-3579 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[2] + 4326 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[3] + 814 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[2]
            + 2976 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[4] + 244 * y_interpolation_values_of_psi_substance_density[5] * y_interpolation_values_of_psi_substance_density[5] - 683 * y_interpolation_values_of_psi_substance_density[5] * y_interpolation_values_of_psi_substance_density[2]
            - 6927 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[3] + 2634 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[2] - 1659 * y_interpolation_values_of_psi_substance_density[5] * y_interpolation_values_of_psi_substance_density[4]
            + 1854 * y_interpolation_values_of_psi_substance_density[5] * y_interpolation_values_of_psi_substance_density[3]) / 180.0;

    //weight を作るための1時的な変数
    MY_FLOAT_TYPE x_alpha[3];
    for (int i = 0; i < 3; ++i) {
        MY_FLOAT_TYPE eps = 0.000001;
        x_alpha[i] = x_ideal_weights[i] / ((eps + x_smoothness_indicators[i]) * (eps + x_smoothness_indicators[i]));
    }
    //weights の計算
    MY_FLOAT_TYPE x_weights[3];
    for (int i = 0; i < 3; ++i) {
        x_weights[i] = x_alpha[i] / (x_alpha[0] + x_alpha[1] + x_alpha[2]);
    }
    return x_weights[0] * x_candidate_interpolants[0] + x_weights[1] * x_candidate_interpolants[1] + x_weights[2] * x_candidate_interpolants[2];
}

MY_FLOAT_TYPE WENO6_interpolation_2D_cell_face_y_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_y_values,
    const Grid &all_grid
){
    int advected_index_x = floor((position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_y = floor((position[1]) / all_grid._cell_length);
    //////WENO6 interpolation の処理
    //y方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE y_interpolation_values_of_psi_substance_density[6];
    //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    for (int ix = 0; ix < 6; ix++) {
        MY_FLOAT_TYPE psi_y_val[6];
        for (int iy = 0; iy < 6; ++iy) {
            int index_x = advected_index_x - 2 + ix;
            int index_y = advected_index_y - 2 + iy;
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
            psi_y_val[iy] = cell_face_y_values[get_voxel_face_index_y(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)];
        }

        // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
        MY_FLOAT_TYPE y_distance_from_face_to_position[6];
        for (int i = 0; i < 6; ++i) {
            y_distance_from_face_to_position[i] = position[1] - (advected_index_y - 2 + i) * all_grid._cell_length;
        }

        // candidate interpolants の計算
        MY_FLOAT_TYPE candidate_interpolants[3];
        for (int i = 0; i < 3; ++i) {
            candidate_interpolants[i]
                = psi_y_val[i]
                + (psi_y_val[i + 1] - psi_y_val[i]) * y_distance_from_face_to_position[i] / (all_grid._cell_length)
                + (psi_y_val[i + 2] - 2 * psi_y_val[i + 1] + psi_y_val[i]) * y_distance_from_face_to_position[i] * y_distance_from_face_to_position[i + 1]
                / (2 * all_grid._cell_length * all_grid._cell_length)
                + (psi_y_val[i + 3] - 3 * psi_y_val[i + 2] + 3 * psi_y_val[i + 1] - psi_y_val[i])
                * y_distance_from_face_to_position[i] * y_distance_from_face_to_position[i + 1] * y_distance_from_face_to_position[i + 2]
                / (6 * all_grid._cell_length * all_grid._cell_length * all_grid._cell_length);
        }

        //ideal weights の計算
        MY_FLOAT_TYPE ideal_weights[3];
        ideal_weights[0]
            = y_distance_from_face_to_position[4] * y_distance_from_face_to_position[5]
            / (20 * all_grid._cell_length * all_grid._cell_length);
        ideal_weights[1]
            = -y_distance_from_face_to_position[5] * y_distance_from_face_to_position[0]
            / (10 * all_grid._cell_length * all_grid._cell_length);
        ideal_weights[2]
            = y_distance_from_face_to_position[0] * y_distance_from_face_to_position[1]
            / (20 * all_grid._cell_length * all_grid._cell_length);

        // smoothness indicatior の計算
        MY_FLOAT_TYPE smoothness_indicators[3];
        smoothness_indicators[0]
            = (-3579 * psi_y_val[3] * psi_y_val[2] + 2634 * psi_y_val[3] * psi_y_val[1] - 683 * psi_y_val[3] * psi_y_val[0]
                - 6927 * psi_y_val[2] * psi_y_val[1] + 1854 * psi_y_val[2] * psi_y_val[0] - 1659 * psi_y_val[1] * psi_y_val[0]
                + 814 * psi_y_val[3] * psi_y_val[3] + 4326 * psi_y_val[2] * psi_y_val[2] + 2976 * psi_y_val[1] * psi_y_val[1]
                + 244 * psi_y_val[0] * psi_y_val[0]) / 180.0;
        smoothness_indicators[1]
            = (-3777 * psi_y_val[3] * psi_y_val[2] + 1074 * psi_y_val[3] * psi_y_val[1] - 1269 * psi_y_val[2] * psi_y_val[1]
                + 1986 * psi_y_val[3] * psi_y_val[3] + 1986 * psi_y_val[2] * psi_y_val[2] + 244 * psi_y_val[1] * psi_y_val[1]
                + 244 * psi_y_val[4] * psi_y_val[4] - 1269 * psi_y_val[4] * psi_y_val[3] + 1074 * psi_y_val[4] * psi_y_val[2]
                - 293 * psi_y_val[4] * psi_y_val[1]) / 180.0;
        smoothness_indicators[2]
            = (-3579 * psi_y_val[3] * psi_y_val[2] + 4326 * psi_y_val[3] * psi_y_val[3] + 814 * psi_y_val[2] * psi_y_val[2]
                + 2976 * psi_y_val[4] * psi_y_val[4] + 244 * psi_y_val[5] * psi_y_val[5] - 683 * psi_y_val[5] * psi_y_val[2]
                - 6927 * psi_y_val[4] * psi_y_val[3] + 2634 * psi_y_val[4] * psi_y_val[2] - 1659 * psi_y_val[5] * psi_y_val[4]
                + 1854 * psi_y_val[5] * psi_y_val[3]) / 180.0;

        //weight を作るための1時的な変数
        MY_FLOAT_TYPE alpha[3];
        for (int i = 0; i < 3; ++i) {
            MY_FLOAT_TYPE eps = 0.000001;
            alpha[i] = ideal_weights[i] / ((eps + smoothness_indicators[i]) * (eps + smoothness_indicators[i]));
        }
        //weights の計算
        MY_FLOAT_TYPE weights[3];
        for (int i = 0; i < 3; ++i) {
            weights[i] = alpha[i] / (alpha[0] + alpha[1] + alpha[2]);
        }
        //return weights[0] * candidate_interpolants[0] + weights[1] * candidate_interpolants[1] + weights[2] * candidate_interpolants[2];
        y_interpolation_values_of_psi_substance_density[ix] = weights[0] * candidate_interpolants[0] + weights[1] * candidate_interpolants[1] + weights[2] * candidate_interpolants[2];
    }

    ////ここからx方向の補間
    // interpolationに使うfaceからposition[1]までの距離 (face が position[1] より上にある場合は負値)
    MY_FLOAT_TYPE x_distance_from_face_to_position[6];
    for (int ix = 0; ix < 6; ++ix) {
        // (advected_index_x - 2 + ix + 0.5) の 0.5 は考えているfaceがy軸に垂直なものを考えているため、面の端から中心までのオフセット分の距離
        x_distance_from_face_to_position[ix] = position[0] - (advected_index_x - 2 + ix + 0.5) * all_grid._cell_length;
    }

    // candidate interpolants の計算
    MY_FLOAT_TYPE x_candidate_interpolants[3];
    for (int i = 0; i < 3; ++i) {
        x_candidate_interpolants[i]
            = y_interpolation_values_of_psi_substance_density[i]
            + (y_interpolation_values_of_psi_substance_density[i + 1] - y_interpolation_values_of_psi_substance_density[i]) * x_distance_from_face_to_position[i] / (all_grid._cell_length)
            + (y_interpolation_values_of_psi_substance_density[i + 2] - 2 * y_interpolation_values_of_psi_substance_density[i + 1] + y_interpolation_values_of_psi_substance_density[i]) * x_distance_from_face_to_position[i] * x_distance_from_face_to_position[i + 1]
            / (2 * all_grid._cell_length * all_grid._cell_length)
            + (y_interpolation_values_of_psi_substance_density[i + 3] - 3 * y_interpolation_values_of_psi_substance_density[i + 2] + 3 * y_interpolation_values_of_psi_substance_density[i + 1] - y_interpolation_values_of_psi_substance_density[i])
            * x_distance_from_face_to_position[i] * x_distance_from_face_to_position[i + 1] * x_distance_from_face_to_position[i + 2]
            / (6 * all_grid._cell_length * all_grid._cell_length * all_grid._cell_length);
    }

    //ideal weights の計算
    MY_FLOAT_TYPE x_ideal_weights[3];
    x_ideal_weights[0]
        = x_distance_from_face_to_position[4] * x_distance_from_face_to_position[5]
        / (20 * all_grid._cell_length * all_grid._cell_length);
    x_ideal_weights[1]
        = -x_distance_from_face_to_position[5] * x_distance_from_face_to_position[0]
        / (10 * all_grid._cell_length * all_grid._cell_length);
    x_ideal_weights[2]
        = x_distance_from_face_to_position[0] * x_distance_from_face_to_position[1]
        / (20 * all_grid._cell_length * all_grid._cell_length);

    // smoothness indicatior の計算
    MY_FLOAT_TYPE x_smoothness_indicators[3];
    x_smoothness_indicators[0]
        = (-3579 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[2] + 2634 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[1] - 683 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[0]
            - 6927 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[1] + 1854 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[0] - 1659 * y_interpolation_values_of_psi_substance_density[1] * y_interpolation_values_of_psi_substance_density[0]
            + 814 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[3] + 4326 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[2] + 2976 * y_interpolation_values_of_psi_substance_density[1] * y_interpolation_values_of_psi_substance_density[1]
            + 244 * y_interpolation_values_of_psi_substance_density[0] * y_interpolation_values_of_psi_substance_density[0]) / 180.0;
    x_smoothness_indicators[1]
        = (-3777 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[2] + 1074 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[1] - 1269 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[1]
            + 1986 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[3] + 1986 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[2] + 244 * y_interpolation_values_of_psi_substance_density[1] * y_interpolation_values_of_psi_substance_density[1]
            + 244 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[4] - 1269 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[3] + 1074 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[2]
            - 293 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[1]) / 180.0;
    x_smoothness_indicators[2]
        = (-3579 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[2] + 4326 * y_interpolation_values_of_psi_substance_density[3] * y_interpolation_values_of_psi_substance_density[3] + 814 * y_interpolation_values_of_psi_substance_density[2] * y_interpolation_values_of_psi_substance_density[2]
            + 2976 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[4] + 244 * y_interpolation_values_of_psi_substance_density[5] * y_interpolation_values_of_psi_substance_density[5] - 683 * y_interpolation_values_of_psi_substance_density[5] * y_interpolation_values_of_psi_substance_density[2]
            - 6927 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[3] + 2634 * y_interpolation_values_of_psi_substance_density[4] * y_interpolation_values_of_psi_substance_density[2] - 1659 * y_interpolation_values_of_psi_substance_density[5] * y_interpolation_values_of_psi_substance_density[4]
            + 1854 * y_interpolation_values_of_psi_substance_density[5] * y_interpolation_values_of_psi_substance_density[3]) / 180.0;

    //weight を作るための1時的な変数
    MY_FLOAT_TYPE x_alpha[3];
    for (int i = 0; i < 3; ++i) {
        MY_FLOAT_TYPE eps = 0.000001;
        x_alpha[i] = x_ideal_weights[i] / ((eps + x_smoothness_indicators[i]) * (eps + x_smoothness_indicators[i]));
    }
    //weights の計算
    MY_FLOAT_TYPE x_weights[3];
    for (int i = 0; i < 3; ++i) {
        x_weights[i] = x_alpha[i] / (x_alpha[0] + x_alpha[1] + x_alpha[2]);
    }
    return x_weights[0] * x_candidate_interpolants[0] + x_weights[1] * x_candidate_interpolants[1] + x_weights[2] * x_candidate_interpolants[2];
}

}//namespace smoke_simulation
