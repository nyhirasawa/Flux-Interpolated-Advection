#pragma once

#include <vector>
#include <Eigen/Dense>
#include "cell_vertex_3D.h"
#include "cell_face_3D.h"
#include "triangle_3D.h"
#include "gauss_quadrature_points.h"
#include "define_float_type.h"
#include "grid_3d.h"

namespace smoke_simulation {
class polygon_3D {
public:
	std::vector<VEC3_TYPE> _vertex_pos_list;
	// どのグリッドセルの面を分割したポリゴンかを表すインデックス
	Eigen::Vector3i _included_cell_index;

    // デフォルトコンストラクタ
    polygon_3D();
    //cell_face から polygon をつくるコンストラクタ
    polygon_3D(cell_face_3D square_face);

    VEC3_TYPE calc_center_position() const;
    void make_polygon_3D_from_cell_face_3D(cell_face_3D square_face);
	//このポリゴンを重心で複数の三角形に分割する。
	std::vector<triangle_3D> split_this_polygon_to_triangles_by_barycenter() const ;
	//ポリゴン上の求積点をquadrature_point_listに追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
	static void add_quadrature_point_list_by_triangle_split_in_xi_eta_space_helper(
		const VEC3_TYPE position,
		const VEC3_TYPE scaled_normal,
		quadrature_point_vector &quadrature_point_list,
		const int i_thread,
//		const bool negative_weight,
		const cell_face_3D& face,
		const int included_cell_index_x,
		const int included_cell_index_y,
		const int included_cell_index_z,
		const MY_FLOAT_TYPE global_weight_factor
	);
	void add_quadrature_point_list_by_triangle_split_in_xi_eta_space(
		quadrature_point_vector &quadrature_point_list,
		const int i_thread,
//	    const bool negative_weight,
		const cell_face_3D& face,
		const int included_cell_index_x,
		const int included_cell_index_y,
		const int included_cell_index_z,
	    const MY_FLOAT_TYPE global_weight_factor
	) const;
	//ポリゴン上の求積点をquadrature_point_listに追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
	void add_quadrature_point_list_by_recursive_triangle_split_in_xi_eta_space(
		quadrature_point_vector &quadrature_point_list,
		const int i_thread,
//	    const bool negative_weight,
		const cell_face_3D& face,
		const int included_cell_index_x,
		const int included_cell_index_y,
		const int included_cell_index_z,
	    const MY_FLOAT_TYPE global_weight_factor
	) const;
	//ポリゴン上の求積点をquadrature_point_listに追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
	void add_quadrature_point_list_by_recursive3_triangle_split_in_xi_eta_space(
		quadrature_point_vector &quadrature_point_list,
		const int i_thread,
//	    const bool negative_weight,
		const cell_face_3D& face,
		const int included_cell_index_x,
		const int included_cell_index_y,
		const int included_cell_index_z,
	    const MY_FLOAT_TYPE global_weight_factor
	) const;
	//ポリゴン上の求積点をquadrature_point_listに追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
	void add_quadrature_point_list_by_recursive4_triangle_split_in_xi_eta_space(
		quadrature_point_vector &quadrature_point_list,
		const int i_thread,
//	    const bool negative_weight,
		const cell_face_3D& face,
		const int included_cell_index_x,
		const int included_cell_index_y,
		const int included_cell_index_z,
	    const MY_FLOAT_TYPE global_weight_factor
	) const;

	// ポリゴン上の求積点をquadrature_point_listに追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
	// psiの計算を補間でなく積分によって行うパターン
	// quadrature_point_list_interpolate_density には積分するdensityからの寄与
	// quadrature_point_list_interpolate_psi にはpsiの値からの寄与
	static void add_quadrature_point_list_of_density_and_psi_by_triangle_split_in_xi_eta_space_helper(
		const VEC3_TYPE position,
		const VEC3_TYPE scaled_normal,
		quadrature_point_vector &quadrature_point_list_interpolate_density,
		quadrature_point_vector &quadrature_point_list_interpolate_psi,
		quadrature_point_vector &quadrature_point_list_for_calc_cell_volume,
		const int i_thread,
//		const bool negative_weight,
		const cell_face_3D& face,
		const int included_cell_index_x,
		const int included_cell_index_y,
		const int included_cell_index_z,
		const MY_FLOAT_TYPE global_weight_factor,
		const Grid_3D &all_grid,
		const int num_gauss_quadrature_point_for_integrate_density,
        const VEC3_TYPE origin_pos
	);
	void add_quadrature_point_list_of_density_and_psi_by_triangle_split_in_xi_eta_space(
		quadrature_point_vector &quadrature_point_list_interpolate_density,
		quadrature_point_vector &quadrature_point_list_interpolate_psi,
		quadrature_point_vector &quadrature_point_list_for_calc_cell_volume,
		const int i_thread,
//	    const bool negative_weight,
		const cell_face_3D& face,
		const int included_cell_index_x,
		const int included_cell_index_y,
		const int included_cell_index_z,
	    const MY_FLOAT_TYPE global_weight_factor,
	    const Grid_3D &all_grid,
	    const int num_gauss_quadrature_point_for_integrate_density,
	    const VEC3_TYPE origin_pos_of_grid
	) const;
	// ポリゴン上の求積点をquadrature_point_listに追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
	// psiの計算を補間でなく積分によって行うパターン
	// quadrature_point_list_interpolate_density には積分するdensityからの寄与
	// quadrature_point_list_interpolate_psi にはpsiの値からの寄与
	void add_quadrature_point_list_of_density_and_psi_by_recursive_triangle_split_in_xi_eta_space(
		quadrature_point_vector &quadrature_point_list_interpolate_density,
		quadrature_point_vector &quadrature_point_list_interpolate_psi,
		const int i_thread,
//	    const bool negative_weight,
		const cell_face_3D& face,
		const int included_cell_index_x,
		const int included_cell_index_y,
		const int included_cell_index_z,
	    const MY_FLOAT_TYPE global_weight_factor,
	    const Grid_3D &all_grid,
	    const int num_gauss_quadrature_point_for_integrate_density,
	    const VEC3_TYPE origin_pos_of_grid
	) const;
	// ポリゴン上の求積点をquadrature_point_listに追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
	// psiの計算を補間でなく積分によって行うパターン
	// quadrature_point_list_interpolate_density には積分するdensityからの寄与
	// quadrature_point_list_interpolate_psi にはpsiの値からの寄与
	void add_quadrature_point_list_of_density_and_psi_by_recursive3_triangle_split_in_xi_eta_space(
		quadrature_point_vector &quadrature_point_list_interpolate_density,
		quadrature_point_vector &quadrature_point_list_interpolate_psi,
		const int i_thread,
//	    const bool negative_weight,
		const cell_face_3D& face,
		const int included_cell_index_x,
		const int included_cell_index_y,
		const int included_cell_index_z,
	    const MY_FLOAT_TYPE global_weight_factor,
	    const Grid_3D &all_grid
	) const;
	// ポリゴン上の求積点をquadrature_point_listに追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
	// psiの計算を補間でなく積分によって行うパターン
	// quadrature_point_list_interpolate_density には積分するdensityからの寄与
	// quadrature_point_list_interpolate_psi にはpsiの値からの寄与
	void add_quadrature_point_list_of_density_and_psi_by_recursive4_triangle_split_in_xi_eta_space(
		quadrature_point_vector &quadrature_point_list_interpolate_density,
		quadrature_point_vector &quadrature_point_list_interpolate_psi,
		const int i_thread,
//	    const bool negative_weight,
		const cell_face_3D& face,
		const int included_cell_index_x,
		const int included_cell_index_y,
		const int included_cell_index_z,
	    const MY_FLOAT_TYPE global_weight_factor,
	    const Grid_3D &all_grid
	) const;


	//ポリゴン上の求積点をquadrature_point_listに追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
	void add_quadrature_point_list_by_triangle_split(
		quadrature_point_vector &quadrature_point_list,
		const int i_thread,
//	    const bool negative_weight,
	    const MY_FLOAT_TYPE global_weight_factor
	) const;
	//ポリゴン上の求積点をquadrature_point_list[i_thread]に追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
	void add_quadrature_point_list_by_recursive_triangle_split(
	    quadrature_point_vector &quadrature_point_list,
	    const int i_thread,
//	    const bool negative_weight,
	    const MY_FLOAT_TYPE global_weight_factor
	) const;
	//ポリゴン上の求積点をquadrature_point_list[i_thread]に追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
	void add_quadrature_point_list_by_recursive3_triangle_split(
	    quadrature_point_vector &quadrature_point_list,
	    const int i_thread,
//	    const bool negative_weight,
	    const MY_FLOAT_TYPE global_weight_factor
	) const;
	//ポリゴン上の求積点をquadrature_point_list[i_thread]に追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
	void add_quadrature_point_list_by_recursive4_triangle_split(
	    quadrature_point_vector &quadrature_point_list,
	    const int i_thread,
//	    const bool negative_weight,
	    const MY_FLOAT_TYPE global_weight_factor
	) const;
	//ポリゴン上の求積点をquadrature_point_list[i_thread]に追加する関数(ガウス求積)
	void add_quadrature_point_list_by_gauss_quadrature(
	    quadrature_point_vector &quadrature_point_list,
	    const int i_thread,
	    const int num_gauss_quadrature_point,
//	    const bool negative_weight,
	    const MY_FLOAT_TYPE global_weight_factor
	) const;
	std::vector<polygon_3D> split_this_polygon_to_triangle_and_square() const;
};
}//namespace smoke_simulation
