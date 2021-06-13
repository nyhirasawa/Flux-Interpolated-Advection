#include "linear_interpolation_2d.h"

#include <vector>
#include "utils.h"
#include "gauss_quadrature_points.h"

namespace smoke_simulation{
    MY_FLOAT_TYPE linear_interpolation_2D(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid &all_grid,
        const VEC3_TYPE origin_pos,
        const int grid_num_x,
        const int grid_num_y,
        const MY_FLOAT_TYPE cell_length
    ){
        int advected_index_x = floor((position[0] - origin_pos[0]) / cell_length);
        int advected_index_y = floor((position[1] - origin_pos[1]) / cell_length);

        /////linear interpolation の処理
        //y方向のinterpolationの結果を格納するための変数
        MY_FLOAT_TYPE y_interpolation_values_of_psi_substance_density[2];
        //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
        for (int ix = 0; ix < 2; ix++) {
            //補間に使う離散値をセット
            MY_FLOAT_TYPE psi_y_val[2];
            for (int iy = 0; iy < 2; ++iy) {
                int index_x = advected_index_x + ix;
                int index_y = advected_index_y + iy;
                if (index_x < 0) {
                    index_x = 0;
                }
                if (index_x > grid_num_x - 1) {
                    int exceed = index_x - (grid_num_x - 1);
                    //境界の値が外側までずっと続く場合
                    index_x = grid_num_x - 1;
                    //境界を境に鏡のように値が反射する場合
                    //index_x = all_grid.Grid_num_x - exceed + 1;
                }
                if (index_y < 0) {
                    index_y = 0;
                }
                if (index_y > grid_num_y - 1) {
                    int exceed = index_y - (grid_num_y - 1);
                    //境界の値が外側までずっと続く場合
                    index_y = grid_num_y - 1;
                    //境界を境に鏡のように値が反射する場合
                    //index_y = all_grid.Grid_num_y - 1 - exceed + 1;
                }
                psi_y_val[iy] = cell_center_values[get_voxel_center_index(index_x, index_y, grid_num_x, grid_num_y)];
            }
            //interpolation の係数
            MY_FLOAT_TYPE b0, b1;
            b0 = (position[1] - origin_pos[1]) / cell_length - advected_index_y;
            b1 = 1.0 - b0;
            y_interpolation_values_of_psi_substance_density[ix] = b1 * psi_y_val[0] + b0 * psi_y_val[1];
        }

        ////ここからx方向の補間
        //interpolation の係数
        MY_FLOAT_TYPE a0, a1;
        a0 = (position[0] - origin_pos[0]) / cell_length - advected_index_x;
        a1 = 1.0 - a0;
        //結果
        return a1 * y_interpolation_values_of_psi_substance_density[0] + a0 * y_interpolation_values_of_psi_substance_density[1];
    }

MY_FLOAT_TYPE linear_interpolation_2D_cell_center_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_center_values,
    const Grid &all_grid
){
    return linear_interpolation_2D(
        position,
        cell_center_values,
        all_grid,
        VEC3_TYPE(0.5 * all_grid._cell_length, 0.5 * all_grid._cell_length, 0.0),
        all_grid.Grid_num_x,
        all_grid.Grid_num_y,
        all_grid._cell_length
    );
/*
    int advected_index_x = floor((position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);

    /////linear interpolation の処理
    //y方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE y_interpolation_values_of_psi_substance_density[2];
    //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    for (int ix = 0; ix < 2; ix++) {
        //補間に使う離散値をセット
        MY_FLOAT_TYPE psi_y_val[2];
        for (int iy = 0; iy < 2; ++iy) {
            int index_x = advected_index_x + ix;
            int index_y = advected_index_y + iy;
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
            psi_y_val[iy] = cell_center_values[get_voxel_center_index(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)];
        }
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_y;
        b1 = 1.0 - b0;
        y_interpolation_values_of_psi_substance_density[ix] = b1 * psi_y_val[0] + b0 * psi_y_val[1];
    }

    ////ここからx方向の補間
    //interpolation の係数
    MY_FLOAT_TYPE a0, a1;
    a0 = (position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_x;
    a1 = 1.0 - a0;
    //結果
    return a1 * y_interpolation_values_of_psi_substance_density[0] + a0 * y_interpolation_values_of_psi_substance_density[1];
*/
}

MY_FLOAT_TYPE linear_interpolation_2D_cell_face_x_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
    const Grid &all_grid
) {
    return linear_interpolation_2D(
        position,
        cell_face_x_values,
        all_grid,
        VEC3_TYPE(0.0, 0.5 * all_grid._cell_length, 0.0),
        all_grid.Grid_num_x + 1,
        all_grid.Grid_num_y,
        all_grid._cell_length
    );
/*
    int advected_index_x = floor((position[0]) / all_grid._cell_length);
    int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);

    /////linear interpolation の処理
    //y方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE y_interpolation_values_of_psi_substance_density[2];
    //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    for (int ix = 0; ix < 2; ix++) {
        //補間に使う離散値をセット
        MY_FLOAT_TYPE psi_y_val[2];
        for (int iy = 0; iy < 2; ++iy) {
            int index_x = advected_index_x + ix;
            int index_y = advected_index_y + iy;
            if (index_x < 0) {
                index_x = 0;
            }
            if (index_x > all_grid.Grid_num_x) {
                int exceed = index_x - (all_grid.Grid_num_x);
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
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_y;
        b1 = 1.0 - b0;
        y_interpolation_values_of_psi_substance_density[ix] = b1 * psi_y_val[0] + b0 * psi_y_val[1];
    }

    ////ここからx方向の補間
    //interpolation の係数
    MY_FLOAT_TYPE a0, a1;
    a0 = (position[0]) / all_grid._cell_length - advected_index_x;
    a1 = 1.0 - a0;
    //結果
    return a1 * y_interpolation_values_of_psi_substance_density[0] + a0 * y_interpolation_values_of_psi_substance_density[1];
*/
}

MY_FLOAT_TYPE linear_interpolation_2D_cell_face_y_values(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_y_values,
    const Grid &all_grid
) {
    return linear_interpolation_2D(
        position,
        cell_face_y_values,
        all_grid,
        VEC3_TYPE(0.5 * all_grid._cell_length, 0.0, 0.0),
        all_grid.Grid_num_x,
        all_grid.Grid_num_y + 1,
        all_grid._cell_length
    );
/*
    int advected_index_x = floor((position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
    int advected_index_y = floor((position[1]) / all_grid._cell_length);

    /////linear interpolation の処理
    //y方向のinterpolationの結果を格納するための変数
    MY_FLOAT_TYPE y_interpolation_values_of_psi_substance_density[2];
    //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    for (int ix = 0; ix < 2; ix++) {
        //補間に使う離散値をセット
        MY_FLOAT_TYPE psi_y_val[2];
        for (int iy = 0; iy < 2; ++iy) {
            int index_x = advected_index_x + ix;
            int index_y = advected_index_y + iy;
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
                int exceed = index_y - (all_grid.Grid_num_y);
                //境界の値が外側までずっと続く場合
                index_y = all_grid.Grid_num_y;
                //境界を境に鏡のように値が反射する場合
                //index_y = all_grid.Grid_num_y - exceed + 1;
            }
            psi_y_val[iy] = cell_face_y_values[get_voxel_face_index_y(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)];
        }
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1]) / all_grid._cell_length - advected_index_y;
        b1 = 1.0 - b0;
        y_interpolation_values_of_psi_substance_density[ix] = b1 * psi_y_val[0] + b0 * psi_y_val[1];
    }

    ////ここからx方向の補間
    //interpolation の係数
    MY_FLOAT_TYPE a0, a1;
    a0 = (position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_x;
    a1 = 1.0 - a0;
    //結果
    return a1 * y_interpolation_values_of_psi_substance_density[0] + a0 * y_interpolation_values_of_psi_substance_density[1];
*/
}

MY_FLOAT_TYPE linear_interpolation_2D_psi_velocity_x(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_vertex_values,
    const Grid &all_grid
) {
    return linear_interpolation_2D(
        position,
        cell_vertex_values,
        all_grid,
        VEC3_TYPE(0.0, 0.0, 0.0),
        all_grid.Grid_num_x + 1,
        all_grid.Grid_num_y + 1,
        all_grid._cell_length
    );
}


    MY_FLOAT_TYPE linear_interpolation_2D_psi_velocity_y(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid &all_grid
    ){
        return linear_interpolation_2D(
            position,
            cell_center_values,
            all_grid,
            VEC3_TYPE(0.5 * all_grid._cell_length, -0.5 * all_grid._cell_length, 0.0),
            all_grid.Grid_num_x,
            all_grid.Grid_num_y + 2,
            all_grid._cell_length
        );
    }

    // psi_density を積分で計算するとき用
    MY_FLOAT_TYPE linear_interpolation_2D_cell_face_y_values_use_integral(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &psi_cell_face_y,
        const std::vector<MY_FLOAT_TYPE> &density_cell_center,
        const Grid &all_grid,
        const int num_gauss_quadrature_point_for_integrate_density
    ) {
        int advected_index_x = floor(position[0] / all_grid._cell_length);
        int advected_index_y = floor(position[1] / all_grid._cell_length);

        // positionが下端の面からどれだけ離れた位置にあるかを[0, 1]で表した量
        MY_FLOAT_TYPE over_under_face = (position[1] / all_grid._cell_length) - advected_index_y;

        MY_FLOAT_TYPE interpolated_psi_density = 0.0;

        //// position の下面のpsiをまず足す(psiからの寄与)
        interpolated_psi_density += linear_interpolation_2D_cell_face_y_values(
            VEC3_TYPE(position[0], advected_index_y * all_grid._cell_length, 0.0),
            psi_cell_face_y,
            all_grid
        );

        //// density からの寄与を計算する
        // positionがセルの中心より下にある場合
        if(over_under_face < 0.5){
            //求積点の位置のリスト
            std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                    advected_index_y * all_grid._cell_length,
            		position[1],
            		num_gauss_quadrature_point_for_integrate_density
                );
            //求積点の重みのリスト
            std::vector<MY_FLOAT_TYPE> quadrature_weight_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_weights_1D(
            		num_gauss_quadrature_point_for_integrate_density
                );

            for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                //求積点での密度の値
                MY_FLOAT_TYPE density_at_quadrature_point
                    = linear_interpolation_2D_cell_center_values(
                        VEC3_TYPE(
                            position[0],
                            quadrature_position_list_1D_density_integral[i_point][0],
                            0.0
                        ),
                        density_cell_center,
                        all_grid
                    );
                // 密度からの寄与を足す
                interpolated_psi_density
                    += density_at_quadrature_point
                    * over_under_face
                    * all_grid._cell_length
                    * quadrature_weight_list_1D_density_integral[i_point];
            }
        }
        else{
            ////セルの下半分の求積点からの寄与
            //求積点の位置のリスト
            std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                    advected_index_y * all_grid._cell_length,
                    (advected_index_y + 0.5) * all_grid._cell_length,
            		num_gauss_quadrature_point_for_integrate_density
                );
            //求積点の重みのリスト
            std::vector<MY_FLOAT_TYPE> quadrature_weight_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_weights_1D(
            		num_gauss_quadrature_point_for_integrate_density
                );

            for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                //求積点での密度の値
                MY_FLOAT_TYPE density_at_quadrature_point
                    = linear_interpolation_2D_cell_center_values(
                        VEC3_TYPE(
                            position[0],
                            quadrature_position_list_1D_density_integral[i_point][0],
                            0.0
                        ),
                        density_cell_center,
                        all_grid
                    );
                // 密度からの寄与を足す
                interpolated_psi_density
                    += density_at_quadrature_point
                    * 0.5
                    * all_grid._cell_length
                    * quadrature_weight_list_1D_density_integral[i_point];
            }
            ////セルの上半分の求積点からの寄与
            //求積点の位置のリスト
            quadrature_position_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                    (advected_index_y + 0.5) * all_grid._cell_length,
                    position[1],
            		num_gauss_quadrature_point_for_integrate_density
                );
            //求積点の重みのリスト
            quadrature_weight_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_weights_1D(
            		num_gauss_quadrature_point_for_integrate_density
                );

            for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                //求積点での密度の値
                MY_FLOAT_TYPE density_at_quadrature_point
                    = linear_interpolation_2D_cell_center_values(
                        VEC3_TYPE(
                            position[0],
                            quadrature_position_list_1D_density_integral[i_point][0],
                            0.0
                        ),
                        density_cell_center,
                        all_grid
                    );
                // 密度からの寄与を足す
                interpolated_psi_density
                    += density_at_quadrature_point
                    * (over_under_face - 0.5)
                    * all_grid._cell_length
                    * quadrature_weight_list_1D_density_integral[i_point];
            }
        }
        return interpolated_psi_density;
    }

    // psi_velocity_x を積分で計算するとき用
    MY_FLOAT_TYPE linear_interpolation_2D_psi_velocity_x_use_integral(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &psi_cell_face_y,
        const std::vector<MY_FLOAT_TYPE> &velocity_x,
        const Grid &all_grid,
        const int num_gauss_quadrature_point_for_integrate_density
    ){
        int advected_index_x = floor((position[0] + 0.5 * all_grid._cell_length) / all_grid._cell_length);
//        int advected_index_x = floor((position[0] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
//        int advected_index_x = floor((position[0]) / all_grid._cell_length);
//        int advected_index_y = floor((position[1] + 0.5 * all_grid._cell_length) / all_grid._cell_length);
//        int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        int advected_index_y = floor(position[1] / all_grid._cell_length);

        // positionが下端の面からどれだけ離れた位置にあるかを[0, 1]で表した量
//        MY_FLOAT_TYPE over_under_face = ((position[1] + 0.5 * all_grid._cell_length) / all_grid._cell_length) - advected_index_y;
//        MY_FLOAT_TYPE over_under_face = ((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length) - advected_index_y;
        MY_FLOAT_TYPE over_under_face = (position[1] / all_grid._cell_length) - advected_index_y;

        MY_FLOAT_TYPE interpolated_psi_density = 0.0;

        //// position の下面のpsiをまず足す(psiからの寄与)
//        interpolated_psi_density += linear_interpolation_2D_psi_velocity_x(
//            VEC3_TYPE(position[0], (advected_index_y + 0.5) * all_grid._cell_length, 0.0),
//            psi_cell_face_y,
//            all_grid
//        );
//        interpolated_psi_density += linear_interpolation_2D_psi_velocity_x(
//            VEC3_TYPE(position[0], (advected_index_y - 0.5) * all_grid._cell_length, 0.0),
//            psi_cell_face_y,
//            all_grid
//        );
        interpolated_psi_density += linear_interpolation_2D_psi_velocity_x(
            VEC3_TYPE(position[0], advected_index_y * all_grid._cell_length, 0.0),
            psi_cell_face_y,
            all_grid
        );

        //// density からの寄与を計算する
        // positionがセルの中心より下にある場合
        if(over_under_face < 0.5){
            //求積点の位置のリスト
//            std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral
//                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
//                    (advected_index_y + 0.5) * all_grid._cell_length,
//            		position[1],
//            		num_gauss_quadrature_point_for_integrate_density
//                );
//            std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral
//                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
//                    (advected_index_y - 0.5) * all_grid._cell_length,
//            		position[1],
//            		num_gauss_quadrature_point_for_integrate_density
//                );
            std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                    advected_index_y * all_grid._cell_length,
            		position[1],
            		num_gauss_quadrature_point_for_integrate_density
                );
            //求積点の重みのリスト
            std::vector<MY_FLOAT_TYPE> quadrature_weight_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_weights_1D(
            		num_gauss_quadrature_point_for_integrate_density
                );

            for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                //求積点での密度の値
                MY_FLOAT_TYPE density_at_quadrature_point
                    = linear_interpolation_2D_cell_face_x_values(
                        VEC3_TYPE(
                            position[0],
                            quadrature_position_list_1D_density_integral[i_point][0],
                            0.0
                        ),
                        velocity_x,
                        all_grid
                    );
                // 密度からの寄与を足す
                interpolated_psi_density
                    += density_at_quadrature_point
                    * over_under_face
                    * all_grid._cell_length
                    * quadrature_weight_list_1D_density_integral[i_point];
            }
        }
        else{
            ////セルの下半分の求積点からの寄与
            //求積点の位置のリスト
//            std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral
//                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
//                    (advected_index_y - 0.5) * all_grid._cell_length,
//                    advected_index_y * all_grid._cell_length,
//            		num_gauss_quadrature_point_for_integrate_density
//                );
//            std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral
//                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
//                    (advected_index_y + 0.5) * all_grid._cell_length,
//                    (advected_index_y + 1.0) * all_grid._cell_length,
//                	num_gauss_quadrature_point_for_integrate_density
//                );
            std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                    advected_index_y * all_grid._cell_length,
                    (advected_index_y + 0.5) * all_grid._cell_length,
            		num_gauss_quadrature_point_for_integrate_density
                );
            //求積点の重みのリスト
            std::vector<MY_FLOAT_TYPE> quadrature_weight_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_weights_1D(
            		num_gauss_quadrature_point_for_integrate_density
                );

            for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                //求積点での密度の値
                MY_FLOAT_TYPE density_at_quadrature_point
                    = linear_interpolation_2D_cell_face_x_values(
                        VEC3_TYPE(
                            position[0],
                            quadrature_position_list_1D_density_integral[i_point][0],
                            0.0
                        ),
                        velocity_x,
                        all_grid
                    );
                // 密度からの寄与を足す
                interpolated_psi_density
                    += density_at_quadrature_point
                    * 0.5
                    * all_grid._cell_length
                    * quadrature_weight_list_1D_density_integral[i_point];
            }
            ////セルの上半分の求積点からの寄与
            //求積点の位置のリスト
//            quadrature_position_list_1D_density_integral
//                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
//                    advected_index_y * all_grid._cell_length,
//                    position[1],
//            		num_gauss_quadrature_point_for_integrate_density
//                );
//            quadrature_position_list_1D_density_integral
//                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
//                    (advected_index_y + 1.0) * all_grid._cell_length,
//                    position[1],
//            		num_gauss_quadrature_point_for_integrate_density
//                );
            quadrature_position_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                    (advected_index_y + 0.5) * all_grid._cell_length,
                    position[1],
            		num_gauss_quadrature_point_for_integrate_density
                );
            //求積点の重みのリスト
            quadrature_weight_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_weights_1D(
            		num_gauss_quadrature_point_for_integrate_density
                );

            for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                //求積点での密度の値
                MY_FLOAT_TYPE density_at_quadrature_point
                    = linear_interpolation_2D_cell_face_x_values(
                        VEC3_TYPE(
                            position[0],
                            quadrature_position_list_1D_density_integral[i_point][0],
                            0.0
                        ),
                        velocity_x,
                        all_grid
                    );
                // 密度からの寄与を足す
                interpolated_psi_density
                    += density_at_quadrature_point
                    * (over_under_face - 0.5)
                    * all_grid._cell_length
                    * quadrature_weight_list_1D_density_integral[i_point];
            }
        }
        return interpolated_psi_density;
    }

    // psi_velocity_x を積分で計算するとき用
    MY_FLOAT_TYPE linear_interpolation_2D_psi_velocity_y_use_integral(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &psi_cell_face_y,
        const std::vector<MY_FLOAT_TYPE> &velocity_y,
        const Grid &all_grid,
        const int num_gauss_quadrature_point_for_integrate_density
    ){
        int advected_index_x = floor(position[0] / all_grid._cell_length);
        int advected_index_y = floor((position[1] + 0.5 * all_grid._cell_length) / all_grid._cell_length);

        // positionが下端の面からどれだけ離れた位置にあるかを[0, 1]で表した量
        MY_FLOAT_TYPE over_under_face = ((position[1] + 0.5 * all_grid._cell_length) / all_grid._cell_length) - advected_index_y;

        MY_FLOAT_TYPE interpolated_psi_density = 0.0;

        //// position の下面のpsiをまず足す(psiからの寄与)
        interpolated_psi_density += linear_interpolation_2D_psi_velocity_y(
            VEC3_TYPE(position[0], (advected_index_y - 0.5) * all_grid._cell_length, 0.0),
            psi_cell_face_y,
            all_grid
        );

        //// density からの寄与を計算する
        // positionがセルの中心より下にある場合
        if(over_under_face < 0.5){
            //求積点の位置のリスト
            std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                    (advected_index_y - 0.5) * all_grid._cell_length,
            		position[1],
            		num_gauss_quadrature_point_for_integrate_density
                );
            //求積点の重みのリスト
            std::vector<MY_FLOAT_TYPE> quadrature_weight_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_weights_1D(
            		num_gauss_quadrature_point_for_integrate_density
                );

            for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                //求積点での密度の値
                MY_FLOAT_TYPE density_at_quadrature_point
                    = linear_interpolation_2D_cell_face_y_values(
                        VEC3_TYPE(
                            position[0],
                            quadrature_position_list_1D_density_integral[i_point][0],
                            0.0
                        ),
                        velocity_y,
                        all_grid
                    );
                // 密度からの寄与を足す
                interpolated_psi_density
                    += density_at_quadrature_point
                    * over_under_face
                    * all_grid._cell_length
                    * quadrature_weight_list_1D_density_integral[i_point];
            }
        }
        else{
            ////セルの下半分の求積点からの寄与
            //求積点の位置のリスト
            std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                    (advected_index_y - 0.5) * all_grid._cell_length,
                    advected_index_y * all_grid._cell_length,
            		num_gauss_quadrature_point_for_integrate_density
                );
            //求積点の重みのリスト
            std::vector<MY_FLOAT_TYPE> quadrature_weight_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_weights_1D(
            		num_gauss_quadrature_point_for_integrate_density
                );

            for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                //求積点での密度の値
                MY_FLOAT_TYPE density_at_quadrature_point
                    = linear_interpolation_2D_cell_face_y_values(
                        VEC3_TYPE(
                            position[0],
                            quadrature_position_list_1D_density_integral[i_point][0],
                            0.0
                        ),
                        velocity_y,
                        all_grid
                    );
                // 密度からの寄与を足す
                interpolated_psi_density
                    += density_at_quadrature_point
                    * 0.5
                    * all_grid._cell_length
                    * quadrature_weight_list_1D_density_integral[i_point];
            }
            ////セルの上半分の求積点からの寄与
            //求積点の位置のリスト
            quadrature_position_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                    advected_index_y * all_grid._cell_length,
                    position[1],
            		num_gauss_quadrature_point_for_integrate_density
                );
            //求積点の重みのリスト
            quadrature_weight_list_1D_density_integral
                = gauss_quadrature_points_1D::get_quadtarure_weights_1D(
            		num_gauss_quadrature_point_for_integrate_density
                );

            for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                //求積点での密度の値
                MY_FLOAT_TYPE density_at_quadrature_point
                    = linear_interpolation_2D_cell_face_y_values(
                        VEC3_TYPE(
                            position[0],
                            quadrature_position_list_1D_density_integral[i_point][0],
                            0.0
                        ),
                        velocity_y,
                        all_grid
                    );
                // 密度からの寄与を足す
                interpolated_psi_density
                    += density_at_quadrature_point
                    * (over_under_face - 0.5)
                    * all_grid._cell_length
                    * quadrature_weight_list_1D_density_integral[i_point];
            }
        }
        return interpolated_psi_density;
    }

}//namespace smoke_simulation
