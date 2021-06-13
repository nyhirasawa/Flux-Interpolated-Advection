#pragma once

#include "gauss_quadrature_points.h"
#include "grid_3d.h"

namespace smoke_simulation{
    // i_thread 番目のスレッドで計算する求積点のリスト quadrature_point_list_of_i_thread から 移流後の質量密度場substance_density_after_advectへの寄与を計算する
    // psiの計算を補間でなく積分によって行うパターン
	// quadrature_point_list_interpolate_density には積分するdensityからの寄与
	// quadrature_point_list_interpolate_psi にはpsiの値からの寄与
    void calc_density_in_cell_from_quadrature_point_list_for_integrate_density(
        const quadrature_point_vector &quadrature_point_list_interpolate_density,
        const quadrature_point_vector &quadrature_point_list_interpolate_psi,
        const Grid_3D &all_grid,
        std::vector<MY_FLOAT_TYPE> &substance_density_after_advect,
        const int i_thread,
        const std::string interpolation_method,
        const std::vector<MY_FLOAT_TYPE> &advected_values
    );
}
