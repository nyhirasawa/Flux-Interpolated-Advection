#include "polygon_3D.h"

#include <iostream>
#include "gauss_quadrature_points.h"
#include "triangle_3D.h"
#include "define_float_type.h"
#include "grid_3d.h"

namespace smoke_simulation {
polygon_3D::polygon_3D(): _vertex_pos_list(), _included_cell_index(){
//    _vertex_pos_list.reserve(8);
}

//cell_face から polygon をつくるコンストラクタ
polygon_3D::polygon_3D(cell_face_3D square_face): _included_cell_index(){
    _vertex_pos_list.resize(4);
    for(int i_vert = 0; i_vert < 4; ++i_vert){
        _vertex_pos_list[i_vert] = square_face._vertex_list[i_vert]._vertex_pos;
    }
}

VEC3_TYPE polygon_3D::calc_center_position() const{
    //polygon を構成する頂点の数
    const int num_vertex = _vertex_pos_list.size();
    //結果を格納する変数
    VEC3_TYPE center_pos(0.0, 0.0, 0.0);
    for(int i_vert = 0; i_vert < num_vertex; ++i_vert){
        center_pos += _vertex_pos_list[i_vert];
    }
    return center_pos / MY_FLOAT_TYPE(num_vertex);
}

void polygon_3D::make_polygon_3D_from_cell_face_3D(cell_face_3D square_face){
    _vertex_pos_list.resize(4);
    for(int i_vert = 0; i_vert < 4; ++i_vert){
        _vertex_pos_list[i_vert] = square_face._vertex_list[i_vert]._vertex_pos;
    }
}

//このポリゴンを重心で複数の三角形に分割する。
std::vector<triangle_3D> polygon_3D::split_this_polygon_to_triangles_by_barycenter() const {
	const int num_vertex_in_polygon = _vertex_pos_list.size();
	//結果を格納する変数
	std::vector<triangle_3D> splitted_triangle_list;
	VEC3_TYPE center = calc_center_position();
	for(int i_vert = 0; i_vert < num_vertex_in_polygon; ++i_vert){
		triangle_3D added_triangle;
		added_triangle._vertex_list[0] = _vertex_pos_list[i_vert];
		added_triangle._vertex_list[1] = _vertex_pos_list[(i_vert + 1) % num_vertex_in_polygon];
		added_triangle._vertex_list[2] = cell_vertex_3D(center);
		splitted_triangle_list.push_back(added_triangle);
	}
	return splitted_triangle_list;
}

//xi-eta空間のポリゴン上の求積点をquadrature_point_list[i_thread]に追加する関数(xi-eta空間のポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
void polygon_3D::add_quadrature_point_list_by_triangle_split_in_xi_eta_space_helper(
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
) {
	quadrature_point_list._quadtarure_position[i_thread].push_back(position);
//	if(negative_weight){
//		quadrature_point_list._weighted_area_normal[i_thread].push_back(- global_weight_factor * scaled_normal);
//	}
//	else{
		quadrature_point_list._weighted_area_normal[i_thread].push_back(global_weight_factor * scaled_normal);
//	}
	quadrature_point_list._included_cell_index[i_thread].push_back(Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z));
}
void polygon_3D::add_quadrature_point_list_by_triangle_split_in_xi_eta_space(
    quadrature_point_vector &quadrature_point_list,
    const int i_thread,
//    const bool negative_weight,
    const cell_face_3D& face,
    const int included_cell_index_x,
    const int included_cell_index_y,
    const int included_cell_index_z,
    const MY_FLOAT_TYPE global_weight_factor
) const{
	//
	//ポリゴンを重心で三角形に分割
	std::vector<triangle_3D> splitted_triangles_of_this_polygon_in_xi_eta_space = split_this_polygon_to_triangles_by_barycenter();
	//
	const int num_triangle = splitted_triangles_of_this_polygon_in_xi_eta_space.size();
	for (int i_tri = 0; i_tri < num_triangle; ++i_tri) {
		// xi-eta 空間での三角形
		triangle_3D splitted_triangle_in_xi_eta_space = splitted_triangles_of_this_polygon_in_xi_eta_space[i_tri];
		// 実空間での三角形
		triangle_3D splitted_triangle_in_real_space;
		for(int i_vert = 0; i_vert < 3; ++i_vert){
			splitted_triangle_in_real_space._vertex_list[i_vert]
				= cell_vertex_3D(face.get_real_space_coordinate_from_xi_eta_space(
					splitted_triangle_in_xi_eta_space._vertex_list[i_vert]._vertex_pos[0],
					splitted_triangle_in_xi_eta_space._vertex_list[i_vert]._vertex_pos[1]));
		}
		//
		// 求積点の位置とスケールされた法線を求める
		const VEC3_TYPE position = splitted_triangle_in_real_space.calc_center_position();
		const VEC3_TYPE scaled_normal = splitted_triangle_in_real_space.calc_area_scaled_normal();
		//
		// 求積点を追加
		add_quadrature_point_list_by_triangle_split_in_xi_eta_space_helper(
			position,
			scaled_normal,
			quadrature_point_list,
			i_thread,
//			negative_weight,
			face,
			included_cell_index_x,
			included_cell_index_y,
			included_cell_index_z,
			global_weight_factor
		);
	}
}
//ポリゴン上の求積点をquadrature_point_listに追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
void polygon_3D::add_quadrature_point_list_by_recursive_triangle_split_in_xi_eta_space(
    quadrature_point_vector &quadrature_point_list,
    const int i_thread,
//    const bool negative_weight,
    const cell_face_3D& face,
    const int included_cell_index_x,
    const int included_cell_index_y,
    const int included_cell_index_z,
    const MY_FLOAT_TYPE global_weight_factor
) const{
    //ポリゴンを重心で三角形に分割
    std::vector<triangle_3D> splitted_triangles_of_this_polygon_in_xi_eta_space = split_this_polygon_to_triangles_by_barycenter();
    const int num_triangle = splitted_triangles_of_this_polygon_in_xi_eta_space.size();
    for (int i_tri = 0; i_tri < num_triangle; ++i_tri) {
        std::vector<triangle_3D> second_splitted_triangles = splitted_triangles_of_this_polygon_in_xi_eta_space[i_tri].split_this_polygon_to_triangles_by_barycenter();
        for(int i_tri_2= 0; i_tri_2 < second_splitted_triangles.size(); ++i_tri_2){
            // xi-eta 空間での三角形
            triangle_3D splitted_triangle_in_xi_eta_space = second_splitted_triangles[i_tri_2];
            // 実空間での三角形
            triangle_3D splitted_triangle_in_real_space;
            for(int i_vert = 0; i_vert < 3; ++i_vert){
                splitted_triangle_in_real_space._vertex_list[i_vert]
                    = cell_vertex_3D(face.get_real_space_coordinate_from_xi_eta_space(
                        splitted_triangle_in_xi_eta_space._vertex_list[i_vert]._vertex_pos[0],
                        splitted_triangle_in_xi_eta_space._vertex_list[i_vert]._vertex_pos[1]));
            }
            //求積点を追加
            quadrature_point_list._quadtarure_position[i_thread].push_back(splitted_triangle_in_real_space.calc_center_position());
//            if(negative_weight){
//                quadrature_point_list._weighted_area_normal[i_thread].push_back(- global_weight_factor * splitted_triangle_in_real_space.calc_area() * splitted_triangle_in_real_space.calc_normal());
//            }
//            else{
                quadrature_point_list._weighted_area_normal[i_thread].push_back( global_weight_factor * splitted_triangle_in_real_space.calc_area() * splitted_triangle_in_real_space.calc_normal());
//            }
            quadrature_point_list._included_cell_index[i_thread].push_back(Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z));
        }
    }
}

//ポリゴン上の求積点をquadrature_point_listに追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
void polygon_3D::add_quadrature_point_list_by_recursive3_triangle_split_in_xi_eta_space(
    quadrature_point_vector &quadrature_point_list,
    const int i_thread,
//    const bool negative_weight,
    const cell_face_3D& face,
    const int included_cell_index_x,
    const int included_cell_index_y,
    const int included_cell_index_z,
    const MY_FLOAT_TYPE global_weight_factor
) const {
    //ポリゴンを重心で三角形に分割
    std::vector<triangle_3D> splitted_triangles_of_this_polygon_in_xi_eta_space = split_this_polygon_to_triangles_by_barycenter();
    const int num_triangle = splitted_triangles_of_this_polygon_in_xi_eta_space.size();
    for (int i_tri = 0; i_tri < num_triangle; ++i_tri) {
        std::vector<triangle_3D> second_splitted_triangles = splitted_triangles_of_this_polygon_in_xi_eta_space[i_tri].split_this_polygon_to_triangles_by_barycenter();
        for(int i_tri_2= 0; i_tri_2 < second_splitted_triangles.size(); ++i_tri_2){
            std::vector<triangle_3D> third_splitted_triangles = second_splitted_triangles[i_tri_2].split_this_polygon_to_triangles_by_barycenter();
            for(int i_tri_3= 0; i_tri_3 < third_splitted_triangles.size(); ++i_tri_3){
                // xi-eta 空間での三角形
                triangle_3D splitted_triangle_in_xi_eta_space = third_splitted_triangles[i_tri_3];
                // 実空間での三角形
                triangle_3D splitted_triangle_in_real_space;
                for(int i_vert = 0; i_vert < 3; ++i_vert){
                    splitted_triangle_in_real_space._vertex_list[i_vert]
                        = cell_vertex_3D(face.get_real_space_coordinate_from_xi_eta_space(
                            splitted_triangle_in_xi_eta_space._vertex_list[i_vert]._vertex_pos[0],
                            splitted_triangle_in_xi_eta_space._vertex_list[i_vert]._vertex_pos[1]));
                }
                //求積点を追加
                quadrature_point_list._quadtarure_position[i_thread].push_back(splitted_triangle_in_real_space.calc_center_position());
//                if(negative_weight){
//                    quadrature_point_list._weighted_area_normal[i_thread].push_back(- global_weight_factor * splitted_triangle_in_real_space.calc_area() * splitted_triangle_in_real_space.calc_normal());
//                }
//                else{
                    quadrature_point_list._weighted_area_normal[i_thread].push_back( global_weight_factor * splitted_triangle_in_real_space.calc_area() * splitted_triangle_in_real_space.calc_normal());
//                }
                quadrature_point_list._included_cell_index[i_thread].push_back(Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z));
            }
        }
    }
}
//ポリゴン上の求積点をquadrature_point_listに追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
void polygon_3D::add_quadrature_point_list_by_recursive4_triangle_split_in_xi_eta_space(
    quadrature_point_vector &quadrature_point_list,
    const int i_thread,
//    const bool negative_weight,
    const cell_face_3D& face,
    const int included_cell_index_x,
    const int included_cell_index_y,
    const int included_cell_index_z,
    const MY_FLOAT_TYPE global_weight_factor
) const {
    //ポリゴンを重心で三角形に分割
    std::vector<triangle_3D> splitted_triangles_of_this_polygon_in_xi_eta_space = split_this_polygon_to_triangles_by_barycenter();
    const int num_triangle = splitted_triangles_of_this_polygon_in_xi_eta_space.size();
    for (int i_tri = 0; i_tri < num_triangle; ++i_tri) {
        std::vector<triangle_3D> second_splitted_triangles = splitted_triangles_of_this_polygon_in_xi_eta_space[i_tri].split_this_polygon_to_triangles_by_barycenter();
        for(int i_tri_2= 0; i_tri_2 < second_splitted_triangles.size(); ++i_tri_2){
            std::vector<triangle_3D> third_splitted_triangles = second_splitted_triangles[i_tri_2].split_this_polygon_to_triangles_by_barycenter();
            for(int i_tri_3= 0; i_tri_3 < third_splitted_triangles.size(); ++i_tri_3){
                std::vector<triangle_3D> fourth_splitted_triangles = third_splitted_triangles[i_tri_3].split_this_polygon_to_triangles_by_barycenter();
                for(int i_tri_4= 0; i_tri_4 < fourth_splitted_triangles.size(); ++i_tri_4){
                    // xi-eta 空間での三角形
                    triangle_3D splitted_triangle_in_xi_eta_space = fourth_splitted_triangles[i_tri_4];
                    // 実空間での三角形
                    triangle_3D splitted_triangle_in_real_space;
                    for(int i_vert = 0; i_vert < 3; ++i_vert){
                        splitted_triangle_in_real_space._vertex_list[i_vert]
                            = cell_vertex_3D(face.get_real_space_coordinate_from_xi_eta_space(
                                splitted_triangle_in_xi_eta_space._vertex_list[i_vert]._vertex_pos[0],
                                splitted_triangle_in_xi_eta_space._vertex_list[i_vert]._vertex_pos[1]));
                    }
                    //求積点を追加
                    quadrature_point_list._quadtarure_position[i_thread].push_back(splitted_triangle_in_real_space.calc_center_position());
//                    if(negative_weight){
//                        quadrature_point_list._weighted_area_normal[i_thread].push_back(- global_weight_factor * splitted_triangle_in_real_space.calc_area() * splitted_triangle_in_real_space.calc_normal());
//                    }
//                    else{
                        quadrature_point_list._weighted_area_normal[i_thread].push_back(global_weight_factor * splitted_triangle_in_real_space.calc_area() * splitted_triangle_in_real_space.calc_normal());
//                    }
                    quadrature_point_list._included_cell_index[i_thread].push_back(Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z));
                }
            }
        }
    }
}


// ポリゴン上の求積点をquadrature_point_listに追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
// psiの計算を補間でなく積分によって行うパターン
// quadrature_point_list_interpolate_density には積分するdensityからの寄与
// quadrature_point_list_interpolate_psi にはpsiの値からの寄与
void polygon_3D::add_quadrature_point_list_of_density_and_psi_by_triangle_split_in_xi_eta_space_helper(
		const VEC3_TYPE center_pos_of_triangle,
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
) {
	//三角形の中心が含まれるセルのインデックス
	int advected_index_y = floor((center_pos_of_triangle[1] - origin_pos[1]) / all_grid._cell_length);

	////////////////////////////////////////
	////densityについての求積点を追加
	////////////////////////////////////////
	// positionが下端の面からどれだけ離れた位置にあるかを[0, 1]で表した量
	MY_FLOAT_TYPE over_under_face = ((center_pos_of_triangle[1] - origin_pos[1]) / all_grid._cell_length) - advected_index_y;
//        if(over_under_face < 0.0 || over_under_face > 1.0){
//            std::cout<<"over_under_face: "<<over_under_face<<std::endl;
//        }
	// positionがセルの中心より下にある場合
	if(over_under_face < 0.5){
		std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral
			= gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
				advected_index_y * all_grid._cell_length + origin_pos[1],
				center_pos_of_triangle[1],
				num_gauss_quadrature_point_for_integrate_density
			);
		std::vector<MY_FLOAT_TYPE> quadrature_weight_list_1D_density_integral
			= gauss_quadrature_points_1D::get_quadtarure_weights_1D(
				num_gauss_quadrature_point_for_integrate_density
			);
/*
		if(abs((advected_index_y + over_under_face / 2.0) * all_grid._cell_length - quadrature_position_list_1D_density_integral[0][0])>0.0001){
			std::cout<<"(advected_index_y + over_under_face / 2.0) * all_grid._cell_length: "<<(advected_index_y + over_under_face / 2.0) * all_grid._cell_length<<std::endl;
			std::cout<<"quadrature_position_list_1D_density_integral[i_quad][0]: "<<quadrature_position_list_1D_density_integral[0][0]<<std::endl;
		}
*/
		for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
			quadrature_point_list_interpolate_density._quadtarure_position[i_thread].push_back(
				VEC3_TYPE(
					center_pos_of_triangle[0],
					quadrature_position_list_1D_density_integral[i_point][0],
//                        (advected_index_y + over_under_face / 2.0) * all_grid._cell_length,
					center_pos_of_triangle[2]
				)
			);
//			if(negative_weight){
//				quadrature_point_list_interpolate_density._weighted_area_normal[i_thread].push_back(
//					- global_weight_factor
//					* scaled_normal
//					* over_under_face
//					* all_grid._cell_length
//					* quadrature_weight_list_1D_density_integral[i_point]
//				);
//			}
//			else{
				quadrature_point_list_interpolate_density._weighted_area_normal[i_thread].push_back(
					global_weight_factor
					* scaled_normal
					* over_under_face
					* all_grid._cell_length
					* quadrature_weight_list_1D_density_integral[i_point]
				);
//			}
			quadrature_point_list_interpolate_density._included_cell_index[i_thread].push_back(
				Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z)
			);
		}
/*
		quadrature_point_list_interpolate_density._quadtarure_position[i_thread].push_back(
			VEC3_TYPE(
				center_pos_of_triangle[0],
				(advected_index_y + over_under_face / 2.0) * all_grid._cell_length,
				center_pos_of_triangle[2]
			)
		);
		if(negative_weight){
			quadrature_point_list_interpolate_density._weighted_area_normal[i_thread].push_back(
				- global_weight_factor
				* splitted_triangle_in_real_space.calc_area()
				* splitted_triangle_in_real_space.calc_normal()
				* over_under_face
				* all_grid._cell_length
			);
		}
		else{
			quadrature_point_list_interpolate_density._weighted_area_normal[i_thread].push_back(
				global_weight_factor
				* splitted_triangle_in_real_space.calc_area()
				* splitted_triangle_in_real_space.calc_normal()
				* over_under_face
				* all_grid._cell_length
			);
		}
		quadrature_point_list_interpolate_density._included_cell_index[i_thread].push_back(
			Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z)
		);
*/
	}
	// positionがセルの中心より上にある場合
	else{
		//求積点の位置
		std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral;
		//求積の重み
		std::vector<MY_FLOAT_TYPE> quadrature_weight_list_1D_density_integral;
		quadrature_weight_list_1D_density_integral = gauss_quadrature_points_1D::get_quadtarure_weights_1D(
				num_gauss_quadrature_point_for_integrate_density
			);
		//// positionがセルの中心より上にある場合はセル中心より下の1点と上の1点の
		//// 2つの求積点でセルでの積分ができる
		//セル中心より下の1点
		//求積点の位置を計算
		quadrature_position_list_1D_density_integral
			= gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
				advected_index_y * all_grid._cell_length + origin_pos[1],
				(advected_index_y + 0.5) * all_grid._cell_length + origin_pos[1],
				num_gauss_quadrature_point_for_integrate_density
			);
		//求積点を追加
		for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
			quadrature_point_list_interpolate_density._quadtarure_position[i_thread].push_back(
				VEC3_TYPE(
					center_pos_of_triangle[0],
					quadrature_position_list_1D_density_integral[i_point][0],
					center_pos_of_triangle[2]
				)
			);
//			if(negative_weight){
//				quadrature_point_list_interpolate_density._weighted_area_normal[i_thread].push_back(
//					- global_weight_factor
//					* scaled_normal
//					* 0.5
//					* all_grid._cell_length
//					* quadrature_weight_list_1D_density_integral[i_point]
//				);
//			}
//			else{
				quadrature_point_list_interpolate_density._weighted_area_normal[i_thread].push_back(
					global_weight_factor
					* scaled_normal
					* 0.5
					* all_grid._cell_length
					* quadrature_weight_list_1D_density_integral[i_point]
				);
//			}
			quadrature_point_list_interpolate_density._included_cell_index[i_thread].push_back(
				Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z)
			);
		}
		//セル中心より上の1点
		//求積点の位置を計算
		quadrature_position_list_1D_density_integral
			= gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
				(advected_index_y + 0.5) * all_grid._cell_length + origin_pos[1],
				(advected_index_y + over_under_face) * all_grid._cell_length + origin_pos[1],
				num_gauss_quadrature_point_for_integrate_density
			);
		//求積点を追加
		for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
			quadrature_point_list_interpolate_density._quadtarure_position[i_thread].push_back(
				VEC3_TYPE(
					center_pos_of_triangle[0],
					quadrature_position_list_1D_density_integral[i_point][0],
					center_pos_of_triangle[2]
				)
			);
//			if(negative_weight){
//				quadrature_point_list_interpolate_density._weighted_area_normal[i_thread].push_back(
//					- global_weight_factor
//					* scaled_normal
//					* (over_under_face - 0.5)
//					* all_grid._cell_length
//					* quadrature_weight_list_1D_density_integral[i_point]
//				);
//			}
//			else{
				quadrature_point_list_interpolate_density._weighted_area_normal[i_thread].push_back(
					global_weight_factor
					* scaled_normal
					* (over_under_face - 0.5)
					* all_grid._cell_length
					* quadrature_weight_list_1D_density_integral[i_point]
				);
//			}
			quadrature_point_list_interpolate_density._included_cell_index[i_thread].push_back(
				Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z)
			);
		}
	}

	////////////////////////////////////////
	////psiについての求積点を追加
	////////////////////////////////////////
//        quadrature_point_list_interpolate_psi._quadtarure_position[i_thread].push_back(VEC3_TYPE(advected_index_x * all_grid._cell_length, center_pos_of_triangle[1], advected_index_z * all_grid._cell_length));
//        quadrature_point_list_interpolate_psi._quadtarure_position[i_thread].push_back(VEC3_TYPE((advected_index_x + 0.5) * all_grid._cell_length, advected_index_y * all_grid._cell_length, (advected_index_z + 0.5) * all_grid._cell_length));
	quadrature_point_list_interpolate_psi._quadtarure_position[i_thread].push_back(
        VEC3_TYPE(center_pos_of_triangle[0], advected_index_y * all_grid._cell_length + origin_pos[1], center_pos_of_triangle[2])
    );
//	if(negative_weight){
//		quadrature_point_list_interpolate_psi._weighted_area_normal[i_thread].push_back(
//			- global_weight_factor
//			* scaled_normal
//		);
//	}
//	else{
		quadrature_point_list_interpolate_psi._weighted_area_normal[i_thread].push_back(
			global_weight_factor
			* scaled_normal
		);
//	}
	quadrature_point_list_interpolate_psi._included_cell_index[i_thread].push_back(
		Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z)
	);
    ////////////////////////////////////////
	////セル体積についての求積点を追加
	////////////////////////////////////////
    quadrature_point_list_for_calc_cell_volume._quadtarure_position[i_thread].push_back(center_pos_of_triangle);
    quadrature_point_list_for_calc_cell_volume._weighted_area_normal[i_thread].push_back( global_weight_factor * scaled_normal);
    quadrature_point_list_for_calc_cell_volume._included_cell_index[i_thread].push_back(
        Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z)
    );
}
void polygon_3D::add_quadrature_point_list_of_density_and_psi_by_triangle_split_in_xi_eta_space(
    quadrature_point_vector &quadrature_point_list_interpolate_density,
    quadrature_point_vector &quadrature_point_list_interpolate_psi,
    quadrature_point_vector &quadrature_point_list_for_calc_cell_volume,
    const int i_thread,
//    const bool negative_weight,
    const cell_face_3D& face,
    const int included_cell_index_x,
    const int included_cell_index_y,
    const int included_cell_index_z,
    const MY_FLOAT_TYPE global_weight_factor,
    const Grid_3D &all_grid,
    const int num_gauss_quadrature_point_for_integrate_density,
    const VEC3_TYPE origin_pos
) const{
	//ポリゴンを重心で三角形に分割
	std::vector<triangle_3D> splitted_triangles_of_this_polygon_in_xi_eta_space = split_this_polygon_to_triangles_by_barycenter();
	const int num_triangle = splitted_triangles_of_this_polygon_in_xi_eta_space.size();
	for (int i_tri = 0; i_tri < num_triangle; ++i_tri) {
		// xi-eta 空間での三角形
		triangle_3D splitted_triangle_in_xi_eta_space = splitted_triangles_of_this_polygon_in_xi_eta_space[i_tri];
		// 実空間での三角形
		triangle_3D splitted_triangle_in_real_space;
		for(int i_vert = 0; i_vert < 3; ++i_vert){
			splitted_triangle_in_real_space._vertex_list[i_vert]
				= cell_vertex_3D(face.get_real_space_coordinate_from_xi_eta_space(
					splitted_triangle_in_xi_eta_space._vertex_list[i_vert]._vertex_pos[0],
					splitted_triangle_in_xi_eta_space._vertex_list[i_vert]._vertex_pos[1]));
		}
		//
		// 三角形の中心の位置とスケールされた法線を求める
		const VEC3_TYPE center_pos_of_triangle = splitted_triangle_in_real_space.calc_center_position();
		const VEC3_TYPE scaled_normal = splitted_triangle_in_real_space.calc_area_scaled_normal();

		add_quadrature_point_list_of_density_and_psi_by_triangle_split_in_xi_eta_space_helper(
			center_pos_of_triangle,
			scaled_normal,
			quadrature_point_list_interpolate_density,
			quadrature_point_list_interpolate_psi,
            quadrature_point_list_for_calc_cell_volume,
			i_thread,
//			negative_weight,
			face,
			included_cell_index_x,
			included_cell_index_y,
			included_cell_index_z,
			global_weight_factor,
			all_grid,
			num_gauss_quadrature_point_for_integrate_density,
            origin_pos
		);
	}
}


// ポリゴン上の求積点をquadrature_point_listに追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
// psiの計算を補間でなく積分によって行うパターン
// quadrature_point_list_interpolate_density には積分するdensityからの寄与
// quadrature_point_list_interpolate_psi にはpsiの値からの寄与
void polygon_3D::add_quadrature_point_list_of_density_and_psi_by_recursive_triangle_split_in_xi_eta_space(
    quadrature_point_vector &quadrature_point_list_interpolate_density,
    quadrature_point_vector &quadrature_point_list_interpolate_psi,
    const int i_thread,
//    const bool negative_weight,
    const cell_face_3D& face,
    const int included_cell_index_x,
    const int included_cell_index_y,
    const int included_cell_index_z,
    const MY_FLOAT_TYPE global_weight_factor,
    const Grid_3D &all_grid,
    const int num_gauss_quadrature_point_for_integrate_density,
    const VEC3_TYPE origin_pos_of_grid
) const{
    //ポリゴンを重心で三角形に分割
    std::vector<triangle_3D> splitted_triangles_of_this_polygon_in_xi_eta_space = split_this_polygon_to_triangles_by_barycenter();
    const int num_triangle = splitted_triangles_of_this_polygon_in_xi_eta_space.size();
    for (int i_tri = 0; i_tri < num_triangle; ++i_tri) {
        std::vector<triangle_3D> second_splitted_triangles = splitted_triangles_of_this_polygon_in_xi_eta_space[i_tri].split_this_polygon_to_triangles_by_barycenter();
        for(int i_tri_2= 0; i_tri_2 < second_splitted_triangles.size(); ++i_tri_2){
            // xi-eta 空間での三角形
            triangle_3D splitted_triangle_in_xi_eta_space = second_splitted_triangles[i_tri_2];
            // 実空間での三角形
            triangle_3D splitted_triangle_in_real_space;
            for(int i_vert = 0; i_vert < 3; ++i_vert){
                splitted_triangle_in_real_space._vertex_list[i_vert]
                    = cell_vertex_3D(face.get_real_space_coordinate_from_xi_eta_space(
                        splitted_triangle_in_xi_eta_space._vertex_list[i_vert]._vertex_pos[0],
                        splitted_triangle_in_xi_eta_space._vertex_list[i_vert]._vertex_pos[1]));
            }
            //三角形の中心の位置
            VEC3_TYPE center_pos_of_triangle = splitted_triangle_in_real_space.calc_center_position();
            //三角形の中心が含まれるセルのインデックス
            int advected_index_y = floor((center_pos_of_triangle[1]) / all_grid._cell_length);

            ////////////////////////////////////////
            ////densityについての求積点を追加
            ////////////////////////////////////////
            // positionが下端の面からどれだけ離れた位置にあるかを[0, 1]で表した量
            MY_FLOAT_TYPE over_under_face = ((center_pos_of_triangle[1]) / all_grid._cell_length) - advected_index_y;
    //        if(over_under_face < 0.0 || over_under_face > 1.0){
    //            std::cout<<"over_under_face: "<<over_under_face<<std::endl;
    //        }
            // positionがセルの中心より下にある場合
            if(over_under_face < 0.5){
                std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral
                    = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                        advected_index_y * all_grid._cell_length,
                		center_pos_of_triangle[1],
                		num_gauss_quadrature_point_for_integrate_density
                    );
                std::vector<MY_FLOAT_TYPE> quadrature_weight_list_1D_density_integral
                    = gauss_quadrature_points_1D::get_quadtarure_weights_1D(
                		num_gauss_quadrature_point_for_integrate_density
                    );
    /*
                if(abs((advected_index_y + over_under_face / 2.0) * all_grid._cell_length - quadrature_position_list_1D_density_integral[0][0])>0.0001){
                    std::cout<<"(advected_index_y + over_under_face / 2.0) * all_grid._cell_length: "<<(advected_index_y + over_under_face / 2.0) * all_grid._cell_length<<std::endl;
                    std::cout<<"quadrature_position_list_1D_density_integral[i_quad][0]: "<<quadrature_position_list_1D_density_integral[0][0]<<std::endl;
                }
    */
                for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                    quadrature_point_list_interpolate_density._quadtarure_position[i_thread].push_back(
                        VEC3_TYPE(
                            center_pos_of_triangle[0],
                            quadrature_position_list_1D_density_integral[i_point][0],
    //                        (advected_index_y + over_under_face / 2.0) * all_grid._cell_length,
                            center_pos_of_triangle[2]
                        )
                    );
//                    if(negative_weight){
//                        quadrature_point_list_interpolate_density._weighted_area_normal[i_thread].push_back(
//                            - global_weight_factor
//                            * splitted_triangle_in_real_space.calc_area()
//                            * splitted_triangle_in_real_space.calc_normal()
//                            * over_under_face
//                            * all_grid._cell_length
//                            * quadrature_weight_list_1D_density_integral[i_point]
//                        );
//                    }
//                    else{
                        quadrature_point_list_interpolate_density._weighted_area_normal[i_thread].push_back(
                            global_weight_factor
                            * splitted_triangle_in_real_space.calc_area()
                            * splitted_triangle_in_real_space.calc_normal()
                            * over_under_face
                            * all_grid._cell_length
                            * quadrature_weight_list_1D_density_integral[i_point]
                        );
//                    }
                    quadrature_point_list_interpolate_density._included_cell_index[i_thread].push_back(
                        Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z)
                    );
                }
            }
            // positionがセルの中心より上にある場合
            else{
                //求積点の位置
                std::vector<VEC3_TYPE> quadrature_position_list_1D_density_integral;
                //求積の重み
                std::vector<MY_FLOAT_TYPE> quadrature_weight_list_1D_density_integral;
                quadrature_weight_list_1D_density_integral = gauss_quadrature_points_1D::get_quadtarure_weights_1D(
                		num_gauss_quadrature_point_for_integrate_density
                    );
                //// positionがセルの中心より上にある場合はセル中心より下の1点と上の1点の
                //// 2つの求積点でセルでの積分ができる
                //セル中心より下の1点
                //求積点の位置を計算
                quadrature_position_list_1D_density_integral
                    = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                        advected_index_y * all_grid._cell_length,
                		(advected_index_y + 0.5) * all_grid._cell_length,
                		num_gauss_quadrature_point_for_integrate_density
                    );
                //求積点を追加
                for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                    quadrature_point_list_interpolate_density._quadtarure_position[i_thread].push_back(
                        VEC3_TYPE(
                            center_pos_of_triangle[0],
                            quadrature_position_list_1D_density_integral[i_point][0],
                            center_pos_of_triangle[2]
                        )
                    );
//                    if(negative_weight){
//                        quadrature_point_list_interpolate_density._weighted_area_normal[i_thread].push_back(
//                            - global_weight_factor
//                            * splitted_triangle_in_real_space.calc_area()
//                            * splitted_triangle_in_real_space.calc_normal()
//                            * 0.5
//                            * all_grid._cell_length
//                            * quadrature_weight_list_1D_density_integral[i_point]
//                        );
//                    }
//                    else{
                        quadrature_point_list_interpolate_density._weighted_area_normal[i_thread].push_back(
                            global_weight_factor
                            * splitted_triangle_in_real_space.calc_area()
                            * splitted_triangle_in_real_space.calc_normal()
                            * 0.5
                            * all_grid._cell_length
                            * quadrature_weight_list_1D_density_integral[i_point]
                        );
//                    }
                    quadrature_point_list_interpolate_density._included_cell_index[i_thread].push_back(
                        Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z)
                    );
                }
                //セル中心より上の1点
                //求積点の位置を計算
                quadrature_position_list_1D_density_integral
                    = gauss_quadrature_points_1D::get_quadtarure_positions_on_1D_interval(
                        (advected_index_y + 0.5) * all_grid._cell_length,
                        (advected_index_y + over_under_face) * all_grid._cell_length,
                		num_gauss_quadrature_point_for_integrate_density
                    );
                //求積点を追加
                for(int i_point = 0; i_point < num_gauss_quadrature_point_for_integrate_density; ++i_point){
                    quadrature_point_list_interpolate_density._quadtarure_position[i_thread].push_back(
                        VEC3_TYPE(
                            center_pos_of_triangle[0],
                            quadrature_position_list_1D_density_integral[i_point][0],
                            center_pos_of_triangle[2]
                        )
                    );
//                    if(negative_weight){
//                        quadrature_point_list_interpolate_density._weighted_area_normal[i_thread].push_back(
//                            - global_weight_factor
//                            * splitted_triangle_in_real_space.calc_area()
//                            * splitted_triangle_in_real_space.calc_normal()
//                            * (over_under_face - 0.5)
//                            * all_grid._cell_length
//                            * quadrature_weight_list_1D_density_integral[i_point]
//                        );
//                    }
//                    else{
                        quadrature_point_list_interpolate_density._weighted_area_normal[i_thread].push_back(
                            global_weight_factor
                            * splitted_triangle_in_real_space.calc_area()
                            * splitted_triangle_in_real_space.calc_normal()
                            * (over_under_face - 0.5)
                            * all_grid._cell_length
                            * quadrature_weight_list_1D_density_integral[i_point]
                        );
//                    }
                    quadrature_point_list_interpolate_density._included_cell_index[i_thread].push_back(
                        Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z)
                    );
                }
            }
            ////////////////////////////////////////
            ////psiについての求積点を追加
            ////////////////////////////////////////
    //        quadrature_point_list_interpolate_psi._quadtarure_position[i_thread].push_back(VEC3_TYPE(advected_index_x * all_grid._cell_length, center_pos_of_triangle[1], advected_index_z * all_grid._cell_length));
    //        quadrature_point_list_interpolate_psi._quadtarure_position[i_thread].push_back(VEC3_TYPE((advected_index_x + 0.5) * all_grid._cell_length, advected_index_y * all_grid._cell_length, (advected_index_z + 0.5) * all_grid._cell_length));
            quadrature_point_list_interpolate_psi._quadtarure_position[i_thread].push_back(VEC3_TYPE(center_pos_of_triangle[0], advected_index_y * all_grid._cell_length, center_pos_of_triangle[2]));
//            if(negative_weight){
//                quadrature_point_list_interpolate_psi._weighted_area_normal[i_thread].push_back(
//                    - global_weight_factor
//                    * splitted_triangle_in_real_space.calc_area()
//                    * splitted_triangle_in_real_space.calc_normal()
//                );
//            }
//            else{
                quadrature_point_list_interpolate_psi._weighted_area_normal[i_thread].push_back(
                    global_weight_factor
                    * splitted_triangle_in_real_space.calc_area()
                    * splitted_triangle_in_real_space.calc_normal()
                );
//            }
            quadrature_point_list_interpolate_psi._included_cell_index[i_thread].push_back(
                Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z)
            );
		}
	}
}


//ポリゴン上の求積点をquadrature_point_list[i_thread]に追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
void polygon_3D::add_quadrature_point_list_by_triangle_split(
    quadrature_point_vector &quadrature_point_list,
    const int i_thread,
//    const bool negative_weight,
    const MY_FLOAT_TYPE global_weight_factor
) const{
    //ポリゴンを重心で三角形に分割
    std::vector<triangle_3D> splitted_triangles_of_this_polygon = split_this_polygon_to_triangles_by_barycenter();
    const int num_triangle = splitted_triangles_of_this_polygon.size();
    for (int i_tri = 0; i_tri < num_triangle; ++i_tri) {
        quadrature_point_list._quadtarure_position[i_thread].push_back(splitted_triangles_of_this_polygon[i_tri].calc_center_position());
//        if(negative_weight){
//            quadrature_point_list._normal[i_thread].push_back(splitted_triangles_of_this_polygon[i_tri].calc_normal());
//            quadrature_point_list._quadtarure_weight[i_thread].push_back(-global_weight_factor * splitted_triangles_of_this_polygon[i_tri].calc_area());
//            quadrature_point_list._weighted_area_normal[i_thread].push_back(- global_weight_factor * splitted_triangles_of_this_polygon[i_tri].calc_area() * splitted_triangles_of_this_polygon[i_tri].calc_normal());
//        }
//        else{
//            quadrature_point_list._normal[i_thread].push_back(splitted_triangles_of_this_polygon[i_tri].calc_normal());
//            quadrature_point_list._quadtarure_weight[i_thread].push_back(global_weight_factor * splitted_triangles_of_this_polygon[i_tri].calc_area());
            quadrature_point_list._weighted_area_normal[i_thread].push_back( global_weight_factor * splitted_triangles_of_this_polygon[i_tri].calc_area() * splitted_triangles_of_this_polygon[i_tri].calc_normal());
//        }
        quadrature_point_list._included_cell_index[i_thread].push_back(_included_cell_index);
    }
}

//ポリゴン上の求積点をquadrature_point_list[i_thread]に追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
void polygon_3D::add_quadrature_point_list_by_recursive_triangle_split(
    quadrature_point_vector &quadrature_point_list,
    const int i_thread,
//    const bool negative_weight,
    const MY_FLOAT_TYPE global_weight_factor
) const{
    //ポリゴンを重心で三角形に分割
    std::vector<triangle_3D> splitted_triangles_of_this_polygon = split_this_polygon_to_triangles_by_barycenter();
    const int num_triangle = splitted_triangles_of_this_polygon.size();
    for(int i_tri = 0; i_tri < num_triangle; ++i_tri){
        std::vector<triangle_3D> second_splitted_triangles = splitted_triangles_of_this_polygon[i_tri].split_this_polygon_to_triangles_by_barycenter();
        for(int i_tri_2= 0; i_tri_2 < second_splitted_triangles.size(); ++i_tri_2){
            quadrature_point_list._quadtarure_position[i_thread].push_back(second_splitted_triangles[i_tri_2].calc_center_position());
//            if(negative_weight){
//                quadrature_point_list._weighted_area_normal[i_thread].push_back(- global_weight_factor * second_splitted_triangles[i_tri_2].calc_area() * second_splitted_triangles[i_tri_2].calc_normal());
//            }
//            else{
                quadrature_point_list._weighted_area_normal[i_thread].push_back( global_weight_factor * second_splitted_triangles[i_tri_2].calc_area() * second_splitted_triangles[i_tri_2].calc_normal());
//            }
            quadrature_point_list._included_cell_index[i_thread].push_back(_included_cell_index);
        }
    }
}

//ポリゴン上の求積点をquadrature_point_list[i_thread]に追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
void polygon_3D::add_quadrature_point_list_by_recursive3_triangle_split(
    quadrature_point_vector &quadrature_point_list,
    const int i_thread,
//    const bool negative_weight,
    const MY_FLOAT_TYPE global_weight_factor
) const{
    //ポリゴンを重心で三角形に分割
    std::vector<triangle_3D> splitted_triangles_of_this_polygon = split_this_polygon_to_triangles_by_barycenter();
    const int num_triangle = splitted_triangles_of_this_polygon.size();
    for(int i_tri = 0; i_tri < num_triangle; ++i_tri){
        std::vector<triangle_3D> second_splitted_triangles = splitted_triangles_of_this_polygon[i_tri].split_this_polygon_to_triangles_by_barycenter();
        for(int i_tri_2= 0; i_tri_2 < second_splitted_triangles.size(); ++i_tri_2){
            std::vector<triangle_3D> third_splitted_triangles = second_splitted_triangles[i_tri_2].split_this_polygon_to_triangles_by_barycenter();
            for(int i_tri_3 = 0; i_tri_3 < third_splitted_triangles.size(); ++i_tri_3){
                quadrature_point_list._quadtarure_position[i_thread].push_back(third_splitted_triangles[i_tri_3].calc_center_position());
//                if(negative_weight){
//                    quadrature_point_list._weighted_area_normal[i_thread].push_back(- global_weight_factor * third_splitted_triangles[i_tri_3].calc_area() * third_splitted_triangles[i_tri_3].calc_normal());
//                }
//                else{
                    quadrature_point_list._weighted_area_normal[i_thread].push_back( global_weight_factor * third_splitted_triangles[i_tri_3].calc_area() * third_splitted_triangles[i_tri_3].calc_normal());
//                }
                quadrature_point_list._included_cell_index[i_thread].push_back(_included_cell_index);
            }
        }
    }
}

//ポリゴン上の求積点をquadrature_point_list[i_thread]に追加する関数(ポリゴンの中心で三角形分割し、それぞれの三角形の重心を求積点として使う)
void polygon_3D::add_quadrature_point_list_by_recursive4_triangle_split(
    quadrature_point_vector &quadrature_point_list,
    const int i_thread,
//    const bool negative_weight,
    const MY_FLOAT_TYPE global_weight_factor
) const{
    //ポリゴンを重心で三角形に分割
    std::vector<triangle_3D> splitted_triangles_of_this_polygon = split_this_polygon_to_triangles_by_barycenter();
    const int num_triangle = splitted_triangles_of_this_polygon.size();
    for(int i_tri = 0; i_tri < num_triangle; ++i_tri){
        std::vector<triangle_3D> second_splitted_triangles = splitted_triangles_of_this_polygon[i_tri].split_this_polygon_to_triangles_by_barycenter();
        for(int i_tri_2= 0; i_tri_2 < second_splitted_triangles.size(); ++i_tri_2){
            std::vector<triangle_3D> third_splitted_triangles = second_splitted_triangles[i_tri_2].split_this_polygon_to_triangles_by_barycenter();
            for(int i_tri_3 = 0; i_tri_3 < third_splitted_triangles.size(); ++i_tri_3){
                std::vector<triangle_3D> fourth_splitted_triangles = third_splitted_triangles[i_tri_3].split_this_polygon_to_triangles_by_barycenter();
                for(int i_tri_4 = 0; i_tri_4 < fourth_splitted_triangles.size(); ++i_tri_4){
                quadrature_point_list._quadtarure_position[i_thread].push_back(fourth_splitted_triangles[i_tri_4].calc_center_position());
//                    if(negative_weight){
//                        quadrature_point_list._weighted_area_normal[i_thread].push_back(- global_weight_factor * fourth_splitted_triangles[i_tri_4].calc_area() * fourth_splitted_triangles[i_tri_4].calc_normal());
//                    }
//                    else{
                        quadrature_point_list._weighted_area_normal[i_thread].push_back( global_weight_factor * fourth_splitted_triangles[i_tri_4].calc_area() * fourth_splitted_triangles[i_tri_4].calc_normal());
//                    }
                    quadrature_point_list._included_cell_index[i_thread].push_back(_included_cell_index);
                }
            }
        }
    }
}

//ポリゴン上の求積点をquadrature_point_list[i_thread]に追加する関数(ガウス求積)
void polygon_3D::add_quadrature_point_list_by_gauss_quadrature(
    quadrature_point_vector &quadrature_point_list,
    const int i_thread,
    const int num_gauss_quadrature_point,
//    const bool negative_weight,
    const MY_FLOAT_TYPE global_weight_factor
) const{
//    std::cout<<"aaa"<<std::endl;
    // [0, 1] * [0, 1] の正方形上のガウス求積点の位置のリスト
    std::vector<VEC3_TYPE> quadtarure_position_list_in_xi_eta_space;
    // [0, 1] * [0, 1] の正方形上のガウス求積点の重みのリスト
    std::vector<MY_FLOAT_TYPE> quadtarure_weight_list_in_xi_eta_space;
    if(num_gauss_quadrature_point == 1){
        quadtarure_position_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_positions_2D_1point;
        quadtarure_weight_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_weights_2D_1point;
    }
    else if(num_gauss_quadrature_point == 4){
        quadtarure_position_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_positions_2D_4point;
        quadtarure_weight_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_weights_2D_4point;
    }
    else if(num_gauss_quadrature_point == 9){
        quadtarure_position_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_positions_2D_9point;
        quadtarure_weight_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_weights_2D_9point;
    }
    else if(num_gauss_quadrature_point == 16){
        quadtarure_position_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_positions_2D_16point;
        quadtarure_weight_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_weights_2D_16point;
    }
    else if(num_gauss_quadrature_point == 25){
        quadtarure_position_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_positions_2D_25point;
        quadtarure_weight_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_weights_2D_25point;
    }
    else if(num_gauss_quadrature_point == 36){
        quadtarure_position_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_positions_2D_36point;
        quadtarure_weight_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_weights_2D_36point;
    }
    else if(num_gauss_quadrature_point == 49){
        quadtarure_position_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_positions_2D_49point;
        quadtarure_weight_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_weights_2D_49point;
    }
    else if(num_gauss_quadrature_point == 64){
        quadtarure_position_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_positions_2D_64point;
        quadtarure_weight_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_weights_2D_64point;
    }
    else if(num_gauss_quadrature_point == 81){
        quadtarure_position_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_positions_2D_81point;
        quadtarure_weight_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_weights_2D_81point;
    }
    else if(num_gauss_quadrature_point == 100){
        quadtarure_position_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_positions_2D_100point;
        quadtarure_weight_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_weights_2D_100point;
    }
    else if(num_gauss_quadrature_point == 121){
        quadtarure_position_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_positions_2D_121point;
        quadtarure_weight_list_in_xi_eta_space = gauss_quadrature_points_2D::quadtarure_weights_2D_121point;
    }
    else{
        std::cerr << "number of quadrature points is incorrect"<<std::endl;
    }
    //四角形の場合の処理
    if(_vertex_pos_list.size() == 4){
        for(int i_point = 0; i_point < num_gauss_quadrature_point; ++i_point){
            MY_FLOAT_TYPE xi = quadtarure_position_list_in_xi_eta_space[i_point][0];
            MY_FLOAT_TYPE eta = quadtarure_position_list_in_xi_eta_space[i_point][1];
            VEC3_TYPE quadtarure_position
                = (1.0 - xi) * (1.0 - eta) * _vertex_pos_list[0]
                +        xi  * (1.0 - eta) * _vertex_pos_list[1]
                +        xi  *        eta  * _vertex_pos_list[2]
                + (1.0 - xi) *        eta  * _vertex_pos_list[3];
            quadrature_point_list._quadtarure_position[i_thread].push_back(quadtarure_position);
//            if(negative_weight){
//            quadrature_point_list._normal[i_thread].push_back(splitted_triangles_of_this_polygon[i_tri].calc_normal());
//            quadrature_point_list._quadtarure_weight[i_thread].push_back(-global_weight_factor * splitted_triangles_of_this_polygon[i_tri].calc_area());
//                quadrature_point_list._weighted_area_normal[i_thread].push_back(
//                    -1 * global_weight_factor * quadtarure_weight_list_in_xi_eta_space[i_point]
//                    *((1.0 - xi) * (1.0 - eta) * (_vertex_pos_list[1] - _vertex_pos_list[0]).cross(_vertex_pos_list[3] - _vertex_pos_list[0])
//                    +        xi  * (1.0 - eta) * (_vertex_pos_list[1] - _vertex_pos_list[0]).cross(_vertex_pos_list[2] - _vertex_pos_list[1])
//                    +        xi  *        eta  * (_vertex_pos_list[2] - _vertex_pos_list[3]).cross(_vertex_pos_list[3] - _vertex_pos_list[0])
//                    + (1.0 - xi) *        eta  * (_vertex_pos_list[2] - _vertex_pos_list[3]).cross(_vertex_pos_list[2] - _vertex_pos_list[1]))
//                );
//            }
//            else{
//            quadrature_point_list._normal[i_thread].push_back(splitted_triangles_of_this_polygon[i_tri].calc_normal());
//            quadrature_point_list._quadtarure_weight[i_thread].push_back(global_weight_factor * splitted_triangles_of_this_polygon[i_tri].calc_area());
                quadrature_point_list._weighted_area_normal[i_thread].push_back(
                     global_weight_factor * quadtarure_weight_list_in_xi_eta_space[i_point]
                    *((1.0 - xi) * (1.0 - eta) * (_vertex_pos_list[1] - _vertex_pos_list[0]).cross(_vertex_pos_list[3] - _vertex_pos_list[0])
                    +        xi  * (1.0 - eta) * (_vertex_pos_list[1] - _vertex_pos_list[0]).cross(_vertex_pos_list[2] - _vertex_pos_list[1])
                    +        xi  *        eta  * (_vertex_pos_list[2] - _vertex_pos_list[3]).cross(_vertex_pos_list[3] - _vertex_pos_list[0])
                    + (1.0 - xi) *        eta  * (_vertex_pos_list[2] - _vertex_pos_list[3]).cross(_vertex_pos_list[2] - _vertex_pos_list[1]))
                );
//            }
            quadrature_point_list._included_cell_index[i_thread].push_back(_included_cell_index);
        }
    }
    //三角形の場合
    else if(_vertex_pos_list.size() == 3){
        triangle_3D tri(_vertex_pos_list[0], _vertex_pos_list[1], _vertex_pos_list[2]);
        MY_FLOAT_TYPE area = tri.calc_area();
        VEC3_TYPE normal = tri.calc_normal();
        for(int i_point = 0; i_point < num_gauss_quadrature_point; ++i_point){
            MY_FLOAT_TYPE xi = quadtarure_position_list_in_xi_eta_space[i_point][0];
            MY_FLOAT_TYPE eta = quadtarure_position_list_in_xi_eta_space[i_point][1];

            MAT3_TYPE transform_mat_from_xi_eta_to_real_space;
            transform_mat_from_xi_eta_to_real_space(0, 0) = _vertex_pos_list[1][0] - _vertex_pos_list[0][0]; transform_mat_from_xi_eta_to_real_space(0, 1) = _vertex_pos_list[2][0] - _vertex_pos_list[0][0]; transform_mat_from_xi_eta_to_real_space(0, 2) = 0.0;
            transform_mat_from_xi_eta_to_real_space(1, 0) = _vertex_pos_list[1][1] - _vertex_pos_list[0][1]; transform_mat_from_xi_eta_to_real_space(1, 1) = _vertex_pos_list[2][1] - _vertex_pos_list[0][1]; transform_mat_from_xi_eta_to_real_space(1, 2) = 0.0;
            transform_mat_from_xi_eta_to_real_space(2, 0) = _vertex_pos_list[1][2] - _vertex_pos_list[0][2]; transform_mat_from_xi_eta_to_real_space(2, 1) = _vertex_pos_list[2][2] - _vertex_pos_list[0][2]; transform_mat_from_xi_eta_to_real_space(2, 2) = 0.0;
            VEC3_TYPE position_in_reference_triangle;
            position_in_reference_triangle[0] = xi * (1.0 - eta);
            position_in_reference_triangle[1] = eta;
            position_in_reference_triangle[2] = 0.0;
            VEC3_TYPE position_in_real_space;
            position_in_real_space = transform_mat_from_xi_eta_to_real_space * position_in_reference_triangle + _vertex_pos_list[0];

            quadrature_point_list._quadtarure_position[i_thread].push_back(position_in_real_space);
//            if(negative_weight){
//            quadrature_point_list._normal[i_thread].push_back(splitted_triangles_of_this_polygon[i_tri].calc_normal());
//            quadrature_point_list._quadtarure_weight[i_thread].push_back(-global_weight_factor * splitted_triangles_of_this_polygon[i_tri].calc_area());
//                quadrature_point_list._weighted_area_normal[i_thread].push_back(
//                    - global_weight_factor * quadtarure_weight_list_in_xi_eta_space[i_point] * area * normal
//                );
//            }
//            else{
//            quadrature_point_list._normal[i_thread].push_back(splitted_triangles_of_this_polygon[i_tri].calc_normal());
//            quadrature_point_list._quadtarure_weight[i_thread].push_back(global_weight_factor * splitted_triangles_of_this_polygon[i_tri].calc_area());
                quadrature_point_list._weighted_area_normal[i_thread].push_back(
                    global_weight_factor * quadtarure_weight_list_in_xi_eta_space[i_point] * area * normal
                );
//            }
            quadrature_point_list._included_cell_index[i_thread].push_back(_included_cell_index);
        }
    }
    //5角形以上の場合は三角形と四角形に分割してそれぞれでガウス求積点を求める
    else {
        std::vector<polygon_3D> splitted_polygon_list = split_this_polygon_to_triangle_and_square();
        int num_splitted_polygon_list = splitted_polygon_list.size();
        for(int i_poly = 0; i_poly < num_splitted_polygon_list; ++i_poly){
            splitted_polygon_list[i_poly].add_quadrature_point_list_by_gauss_quadrature(
                quadrature_point_list,
                i_thread,
                num_gauss_quadrature_point,
//                negative_weight,
                global_weight_factor
            );
        }
//        return add_quadrature_point_list_by_triangle_split(quadrature_point_list, i_thread);
    }
}
std::vector<polygon_3D> polygon_3D::split_this_polygon_to_triangle_and_square() const{
    //頂点の数
    const int num_vertex = _vertex_pos_list.size();
    //結果を格納するstd::vector
    std::vector<polygon_3D> result;
    //三角形か四角形の場合は切断は行わずそのまま
    if(num_vertex < 5){
        polygon_3D tmp;
        tmp._vertex_pos_list = _vertex_pos_list;
        tmp._included_cell_index = _included_cell_index;
        result.push_back(tmp);
        return result;
    }
    //5角形以上の場合の処理
    int i_vert = 0;
    while(true){
        polygon_3D tmp;
        //最初の4角形を切り出す
        if(i_vert == 0){
            for(int i = 0; i<4; ++i){
                tmp._vertex_pos_list.push_back(_vertex_pos_list[i_vert + i]);
            }
            tmp._included_cell_index=_included_cell_index;
            result.push_back(tmp);
            i_vert += 3;
        }
        else{
            for(int i = 0; i<3; ++i){
                tmp._vertex_pos_list.push_back(_vertex_pos_list[i_vert + i]);
            }
            tmp._vertex_pos_list.push_back(_vertex_pos_list[0]);
            tmp._included_cell_index=_included_cell_index;
            result.push_back(tmp);
            i_vert += 2;
        }
        //これ以上4角形を切り出せない場合
        if(i_vert + 3 > num_vertex - 1){
            polygon_3D residual;
            for(int i = 0;; ++i){
                if(i_vert + i > num_vertex - 1){
                    break;
                }
                residual._vertex_pos_list.push_back(_vertex_pos_list[i_vert + i]);
            }
            residual._vertex_pos_list.push_back(_vertex_pos_list[0]);
//            residual._vertex_pos_list.push_back(_vertex_pos_list[i_vert]);
            residual._included_cell_index=_included_cell_index;
            result.push_back(residual);
            break;
        }
    }
    return result;
}
}//namespace smoke_simulation
