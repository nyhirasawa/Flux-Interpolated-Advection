#pragma once

#include "gauss_quadrature_points.h"
#include "grid_3d.h"

namespace smoke_simulation{
    // i_thread 番目のスレッドで計算する求積点のリスト quadrature_point_list_of_i_thread から 移流後の質量密度場substance_density_after_advectへの寄与を計算する
    void calc_density_in_cell_from_quadrature_point_list(
        const quadrature_point_vector &quadrature_point_list,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &substance_density_after_advect,
        const int i_thread,
        const std::string interpolation_method
    );
    // i_thread 番目のスレッドで計算する求積点のリスト quadrature_point_list_of_i_thread から 移流後のセル体積 cell_volume_after_advect への寄与を計算する
    void calc_cell_volumes_from_quadrature_point_list(
        const quadrature_point_vector &quadrature_point_list,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &cell_volume_after_advect,
        const int i_thread
    );
}
