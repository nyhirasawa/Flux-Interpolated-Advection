#pragma once

#include <vector>
#include <algorithm>
#include <functional>
#include <Eigen/Dense>
#include "cell_vertex_3D.h"
#include "cell_face_3D.h"
#include "grid_3d.h"
#include "polygon_3D.h"

#include "define_float_type.h"

namespace smoke_simulation{
    //参照渡しで結果を返すバージョン
    void split_face_and_add_quadrature_points_3D_for_integrate_density_ref(
        quadrature_point_vector &quadrature_point_list_interpolate_density,
        quadrature_point_vector &quadrature_point_list_interpolate_psi,
        quadrature_point_vector &quadrature_point_list_for_calc_cell_volume,
        const Grid_3D& all_grid,
        const cell_face_3D& face,
        const std::string split_method,
        const int i_thread,
//        const bool negative_weight,
        const int included_cell_index_x,
        const int included_cell_index_y,
        const int included_cell_index_z,
        const std::string integral_method,
        const int num_gauss_quadrature_points,
        const MY_FLOAT_TYPE global_weight_factor,
	    const int num_gauss_quadrature_point_for_integrate_density,
        const bool enable_cell_volume_correction,
	    const VEC3_TYPE origin_pos_of_grid,
		const MY_FLOAT_TYPE minimal_area
    );

void split_face_and_add_quadrature_points_3D_ref(
    quadrature_point_vector &quadrature_point_list,
    const Grid_3D& all_grid,
    const cell_face_3D& face,
    const std::string split_method,
    const int i_thread,
//    const bool negative_weight,
    const int included_cell_index_x,
    const int included_cell_index_y,
    const int included_cell_index_z,
    const std::string integral_method,
    const int num_gauss_quadrature_points,
    const MY_FLOAT_TYPE global_weight_factor,
    const MY_FLOAT_TYPE minimal_area
);
// 引数の polygons_in_xi_eta_space をセルの中心で切断する関数
std::vector<polygon_3D> split_polygons_in_xi_eta_space_by_cell_center_3D(
    const cell_face_3D &face,
    const std::vector<polygon_3D> &polygons_in_xi_eta_space,
    std::string split_direction,
    const MY_FLOAT_TYPE cell_length
);
// 引数の polygons_in_xi_eta_space をセルの中心で切断する関数(参照渡しで値を返すバージョン)
void split_polygons_in_xi_eta_space_by_cell_center_3D_ref(
    const cell_face_3D face,
    std::vector<polygon_3D> &polygons_in_xi_eta_space,
    std::string split_direction,
    const MY_FLOAT_TYPE cell_length
);
// 引数の polygons_in_xi_eta_space をセルの面で切断する関数(参照渡しで値を返すバージョン)
void split_polygons_in_xi_eta_space_by_cell_face_3D_ref(
    const cell_face_3D face,
    std::vector<polygon_3D> &polygons_in_xi_eta_space,
    std::string split_direction,
    const MY_FLOAT_TYPE cell_length
);
// ヘルパー関数. この関数の返り値のvector を再度この関数の引数 original_polygon_list_in_xi_eta_space に代入し,
// 異なるiso_value で切断することを繰り返して織木なるのcell face を細かく切断する
// (初回の切断には [0, 1] * [0, 1] の正方形をoriginal_polygon_list_in_xi_eta_space に入れればいい)
// negative density artifact を解消するために,  face に含まれる polygon のリスト全てを軸に垂直な平面(x or y or z軸, およびそれら)で切断する関数
// psi が x方向に(psi の積分の方向をxにとるならy方向に)区分線型関数であることを仮定したときの分割の方法
// split_direction は切断の方向を設定する。例えば(split_direction == "x" ならx軸に垂直な面で切断する)
// iso_value は切断する位置を決める。例えば split_mode == "x" で iso_value==1.9 なら x = 1.9 の平面でfaceを切断する
std::vector<polygon_3D> split_cell_face_by_axis_aligned_plane(
    const cell_face_3D face,
    const std::vector<polygon_3D> original_polygon_list_in_xi_eta_space,
    const MY_FLOAT_TYPE iso_value,
    std::string split_direction
);

// split_cell_face_by_axis_aligned_plane のヘルパー関数
// point_on_line_in_xi_eta_space_0 と point_on_line_in_xi_eta_space_1 を通る直線で polygon を切断し、切断後のpolygonのリストを返す
// OK
std::vector<polygon_3D> split_polygon_in_xi_eta_space_by_line(
    const polygon_3D polygonin_xi_eta_space,
    const VEC3_TYPE point_on_line_in_xi_eta_space_0,
    const VEC3_TYPE point_on_line_in_xi_eta_space_1
);

}//namespace smoke_simulation
