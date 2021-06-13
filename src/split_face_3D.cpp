#include "split_face_3D.h"
#include "facecutter.h"

#include <iostream>
#include <limits>
#include <algorithm>
#include <functional>

#include "grid_3d.h"
#include "define_float_type.h"

#define TMP_POLYGON_VERT_NUM	32
namespace smoke_simulation{

	const static bool use_optimized_split = true;
	//
	static void split_face_optimized (
		const cell_face_3D& face,
		const std::string split_method,
		MY_FLOAT_TYPE cell_length,
		std::function<void(const point2 *polygon, unsigned char num, const point2 &center, const float area)> func,
		const MY_FLOAT_TYPE minimal_area
	) {
		//
		point2 color [4];
		for( int i=0; i<4; ++i ) {
			if(split_method == "xz-CellFaceAligned"){
				color[i].x[0] = face._vertex_list[i]._vertex_pos[0] / cell_length - 0.5;
				color[i].x[1] = face._vertex_list[i]._vertex_pos[2] / cell_length - 0.5;
			} else if(split_method == "xz-AxisAligned"){
				color[i].x[0] = face._vertex_list[i]._vertex_pos[0] / cell_length;
				color[i].x[1] = face._vertex_list[i]._vertex_pos[2] / cell_length;
			}
			else if(split_method == "x-CellFaceAligned_z-AxisAligned"){
				color[i].x[0] = face._vertex_list[i]._vertex_pos[0] / cell_length - 0.5;
				color[i].x[1] = face._vertex_list[i]._vertex_pos[2] / cell_length;
			}
			else if(split_method == "x-AxisAligned_z-CellFaceAligned"){
				color[i].x[0] = face._vertex_list[i]._vertex_pos[0] / cell_length;
				color[i].x[1] = face._vertex_list[i]._vertex_pos[2] / cell_length - 0.5;
			}

		}
		//
		const vertex2 polygon[4] = {
			vertex2(point2(0.0,0.0),color[0]),
			vertex2(point2(1.0,0.0),color[1]),
			vertex2(point2(1.0,1.0),color[2]),
			vertex2(point2(0.0,1.0),color[3])
		};
		//
		cut_face(polygon,4,[&](const vertex2 *polygon, int num) {
			if( num >= 3 ) {
				assert( num <= TMP_POLYGON_VERT_NUM );
				point2 tmp[TMP_POLYGON_VERT_NUM];
				for( int i=0; i<num; ++i ) {
					tmp[i] = polygon[i].p;
				}
				//
//				const float minimal_area (1.0); // 1.0 は面の３角形分割の切断をしない
//				const float minimal_area (0.1); // 1.0 は面の３角形分割の切断をしない
				subdivide(tmp,num,minimal_area,func);
			}
		});
	}
	//
	static void get_3D_center_position_and_scaled_normal(
		const cell_face_3D& face,
		const point2 *polygon, unsigned char num, const point2 &center,
		VEC3_TYPE &position, VEC3_TYPE& scaled_normal, float area ) {
		//
		position = face.get_real_space_coordinate_from_xi_eta_space(center.x[0],center.x[1]);
		assert( num <= TMP_POLYGON_VERT_NUM );
		VEC3_TYPE tmp_positions[TMP_POLYGON_VERT_NUM];
		for( unsigned char n=0; n<num; ++n ) {
			tmp_positions[n] = face.get_real_space_coordinate_from_xi_eta_space(polygon[n].x[0],polygon[n].x[1]);
		}
		if( num == 3 ) {
			const VEC3_TYPE &p0 = tmp_positions[0];
			const VEC3_TYPE &p1 = tmp_positions[1];
			const VEC3_TYPE &p2 = tmp_positions[2];
			scaled_normal = 0.5 * (p1-p0).cross(p2-p0);
		} else if( num > 3 ) {
			scaled_normal = VEC3_TYPE(.0,.0,.0);
			for( unsigned char n=0; n<num; ++n ) {
				const VEC3_TYPE &p0 = position;
				const VEC3_TYPE &p1 = tmp_positions[n];
				const VEC3_TYPE &p2 = tmp_positions[(n+1) % num];
				scaled_normal += (p1-p0).cross(p2-p0);
			}
			scaled_normal *= 0.5;
		}
	}

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
    ){
    	//
		if( use_optimized_split ) {
			split_face_optimized(face,split_method,all_grid._cell_length,
			[&](const point2 *polygon, unsigned char num, const point2 &center, const float area) {
				//
				VEC3_TYPE position, scaled_normal;
				get_3D_center_position_and_scaled_normal(face,polygon,num,center,position,scaled_normal,area);
				polygon_3D::add_quadrature_point_list_of_density_and_psi_by_triangle_split_in_xi_eta_space_helper(
					position,
					scaled_normal,
					quadrature_point_list_interpolate_density,
					quadrature_point_list_interpolate_psi,
					quadrature_point_list_for_calc_cell_volume,
					i_thread,
//					negative_weight,
					face,
					included_cell_index_x,
					included_cell_index_y,
					included_cell_index_z,
					global_weight_factor,
					all_grid,
					num_gauss_quadrature_point_for_integrate_density,
					origin_pos_of_grid
				);
			},
			minimal_area);
		} else {
			//
	        // xi-eta 空間での切断されたポリゴンを格納する変数
	        std::vector<polygon_3D> splitted_faces_in_xi_eta_space;
	        //
	        // xi-eta 空間での face(つまり[0, 1] * [0, 1] の正方形)
	        polygon_3D face_in_xi_eta_space;
	        face_in_xi_eta_space._vertex_pos_list.push_back(VEC3_TYPE(0.0, 0.0, 0.0));
	        face_in_xi_eta_space._vertex_pos_list.push_back(VEC3_TYPE(1.0, 0.0, 0.0));
	        face_in_xi_eta_space._vertex_pos_list.push_back(VEC3_TYPE(1.0, 1.0, 0.0));
	        face_in_xi_eta_space._vertex_pos_list.push_back(VEC3_TYPE(0.0, 1.0, 0.0));

	        splitted_faces_in_xi_eta_space.push_back(face_in_xi_eta_space);
	        // セル中心でx,z軸に垂直な方向に切断
	        if(split_method == "xz-CellFaceAligned"){
	            // xに垂直な方向の切断
	            split_polygons_in_xi_eta_space_by_cell_center_3D_ref(face, splitted_faces_in_xi_eta_space, "x", all_grid._cell_length);
	            // zに垂直な方向の切断
	            split_polygons_in_xi_eta_space_by_cell_center_3D_ref(face, splitted_faces_in_xi_eta_space, "z", all_grid._cell_length);
	        }
	        // セル中心でx,z軸に垂直な方向に切断した後, セルの面でy軸に垂直な切断
	        else if(split_method == "xz-AxisAligned"){
	            // xに垂直な方向の切断
	            split_polygons_in_xi_eta_space_by_cell_face_3D_ref(face, splitted_faces_in_xi_eta_space, "x", all_grid._cell_length);
	            // zに垂直な方向の切断
	            split_polygons_in_xi_eta_space_by_cell_face_3D_ref(face, splitted_faces_in_xi_eta_space, "z", all_grid._cell_length);
	        }
	        // セル中心でx,z軸に垂直な方向に切断した後, セルの面でy軸に垂直な切断
	        else if(split_method == "x-AxisAligned_z-CellFaceAligned"){
	            // xに垂直な方向の切断
	            split_polygons_in_xi_eta_space_by_cell_face_3D_ref(face, splitted_faces_in_xi_eta_space, "x", all_grid._cell_length);
	            // zに垂直な方向の切断
	            split_polygons_in_xi_eta_space_by_cell_center_3D_ref(face, splitted_faces_in_xi_eta_space, "z", all_grid._cell_length);
	        }
	        // セル中心でx,z軸に垂直な方向に切断した後, セルの面でy軸に垂直な切断
	        else if(split_method == "x-CellFaceAligned_z-AxisAligned"){
	            // xに垂直な方向の切断
	            split_polygons_in_xi_eta_space_by_cell_center_3D_ref(face, splitted_faces_in_xi_eta_space, "x", all_grid._cell_length);
	            // zに垂直な方向の切断
	            split_polygons_in_xi_eta_space_by_cell_face_3D_ref(face, splitted_faces_in_xi_eta_space, "z", all_grid._cell_length);
	        }
	        //面を分割しない
	        else if(split_method == "no-split"){
	        }

	        // 各ポリゴンから求積点を計算し、リストに追加
	        int num_splitted_faces = splitted_faces_in_xi_eta_space.size();
	        for(size_t i_poly = 0; i_poly < num_splitted_faces; ++i_poly){
	            splitted_faces_in_xi_eta_space[i_poly]._included_cell_index = Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z);
	            if(integral_method == "triangle"){
	                //三角形分割で求積点を求める
	                splitted_faces_in_xi_eta_space[i_poly].add_quadrature_point_list_of_density_and_psi_by_triangle_split_in_xi_eta_space(
	                    quadrature_point_list_interpolate_density,
	                    quadrature_point_list_interpolate_psi,
	                    quadrature_point_list_for_calc_cell_volume,
	                    i_thread,
//	                    negative_weight,
	                    face,
	                    included_cell_index_x,
	                    included_cell_index_y,
	                    included_cell_index_z,
	                    global_weight_factor,
	                    all_grid,
	                    num_gauss_quadrature_point_for_integrate_density,
//	                    enable_cell_volume_correction,
	                    origin_pos_of_grid
	                );
	            }
	            else if(integral_method == "triangle-rec2"){
	                //三角形分割で求積点を求める
	                splitted_faces_in_xi_eta_space[i_poly].add_quadrature_point_list_of_density_and_psi_by_recursive_triangle_split_in_xi_eta_space(
	                    quadrature_point_list_interpolate_density,
	                    quadrature_point_list_interpolate_psi,
//	                    quadrature_point_list_for_calc_cell_volume,
	                    i_thread,
//	                    negative_weight,
	                    face,
	                    included_cell_index_x,
	                    included_cell_index_y,
	                    included_cell_index_z,
	                    global_weight_factor,
	                    all_grid,
	                    num_gauss_quadrature_point_for_integrate_density,
//	                    enable_cell_volume_correction,
	                    origin_pos_of_grid
	                );
	            }
	/*
	            else if(integral_method == "triangle-rec3"){
	                //三角形分割で求積点を求める
	                splitted_faces_in_xi_eta_space[i_poly].add_quadrature_point_list_of_density_and_psi_by_recursive3_triangle_split_in_xi_eta_space(
	                    quadrature_point_list_interpolate_density,
	                    quadrature_point_list_interpolate_psi,
	                    i_thread,
	                    negative_weight,
	                    face,
	                    included_cell_index_x,
	                    included_cell_index_y,
	                    included_cell_index_z,
	                    global_weight_factor,
	                    all_grid
	                );
	            }
	            else if(integral_method == "triangle-rec4"){
	                //三角形分割で求積点を求める
	                splitted_faces_in_xi_eta_space[i_poly].add_quadrature_point_list_of_density_and_psi_by_recursive4_triangle_split_in_xi_eta_space(
	                    quadrature_point_list_interpolate_density,
	                    quadrature_point_list_interpolate_psi,
	                    i_thread,
	                    negative_weight,
	                    face,
	                    included_cell_index_x,
	                    included_cell_index_y,
	                    included_cell_index_z,
	                    global_weight_factor,
	                    all_grid
	                );
	            }
	*/
	        }
		}
    }

    //参照渡しで結果を返すバージョン
    void split_face_and_add_quadrature_points_3D_ref(
        quadrature_point_vector &quadrature_point_list,
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
		const MY_FLOAT_TYPE minimal_area
    ){
    	//
		if( use_optimized_split ) {
			split_face_optimized(face,split_method,all_grid._cell_length,
			[&](const point2 *polygon, unsigned char num, const point2 &center, const float area) {
				//
				VEC3_TYPE position, scaled_normal;
				get_3D_center_position_and_scaled_normal(face,polygon,num,center,position,scaled_normal,area);
				polygon_3D::add_quadrature_point_list_by_triangle_split_in_xi_eta_space_helper(
					position,
					scaled_normal,
					quadrature_point_list,
					i_thread,
//					negative_weight,
					face,
					included_cell_index_x,
					included_cell_index_y,
					included_cell_index_z,
					global_weight_factor
				);
			},
			minimal_area);
		} else {
	        // xi-eta 空間での切断されたポリゴンを格納する変数
	        std::vector<polygon_3D> splitted_faces_in_xi_eta_space;
	        //
	        // xi-eta 空間での face(つまり[0, 1] * [0, 1] の正方形)
	        polygon_3D face_in_xi_eta_space;
	        face_in_xi_eta_space._vertex_pos_list.push_back(VEC3_TYPE(0.0, 0.0, 0.0));
	        face_in_xi_eta_space._vertex_pos_list.push_back(VEC3_TYPE(1.0, 0.0, 0.0));
	        face_in_xi_eta_space._vertex_pos_list.push_back(VEC3_TYPE(1.0, 1.0, 0.0));
	        face_in_xi_eta_space._vertex_pos_list.push_back(VEC3_TYPE(0.0, 1.0, 0.0));

	        splitted_faces_in_xi_eta_space.push_back(face_in_xi_eta_space);
	        // セル中心でx,z軸に垂直な方向に切断
	        if(split_method == "xz-CellFaceAligned"){
	            // xに垂直な方向の切断
	            split_polygons_in_xi_eta_space_by_cell_center_3D_ref(face, splitted_faces_in_xi_eta_space, "x", all_grid._cell_length);
	            // zに垂直な方向の切断
	            split_polygons_in_xi_eta_space_by_cell_center_3D_ref(face, splitted_faces_in_xi_eta_space, "z", all_grid._cell_length);
	        }
	        // セル中心でx,z軸に垂直な方向に切断した後, セルの面でy軸に垂直な切断
	        else if(split_method == "xz-AxisAligned"){
	            // xに垂直な方向の切断
	            split_polygons_in_xi_eta_space_by_cell_face_3D_ref(face, splitted_faces_in_xi_eta_space, "x", all_grid._cell_length);
	            // zに垂直な方向の切断
	            split_polygons_in_xi_eta_space_by_cell_face_3D_ref(face, splitted_faces_in_xi_eta_space, "z", all_grid._cell_length);
	        }
	        //面を分割しない
	        else if(split_method == "no-split"){
	        }

	        // 各ポリゴンから求積点を計算し、リストに追加
	        int num_splitted_faces = splitted_faces_in_xi_eta_space.size();
	        for(size_t i_poly = 0; i_poly < num_splitted_faces; ++i_poly){
	            splitted_faces_in_xi_eta_space[i_poly]._included_cell_index = Eigen::Vector3i(included_cell_index_x, included_cell_index_y, included_cell_index_z);
	            if(integral_method == "triangle"){
	                //三角形分割で求積点を求める
	                splitted_faces_in_xi_eta_space[i_poly].add_quadrature_point_list_by_triangle_split_in_xi_eta_space(
	                    quadrature_point_list,
	                    i_thread,
//	                    negative_weight,
	                    face,
	                    included_cell_index_x,
	                    included_cell_index_y,
	                    included_cell_index_z,
	                    global_weight_factor
	                );
	            }
	            else if(integral_method == "triangle-rec2"){
	                //三角形分割で求積点を求める
	                splitted_faces_in_xi_eta_space[i_poly].add_quadrature_point_list_by_recursive_triangle_split_in_xi_eta_space(
	                    quadrature_point_list,
	                    i_thread,
//	                    negative_weight,
	                    face,
	                    included_cell_index_x,
	                    included_cell_index_y,
	                    included_cell_index_z,
	                    global_weight_factor
	                );
	            }
	            else if(integral_method == "triangle-rec3"){
	                //三角形分割で求積点を求める
	                splitted_faces_in_xi_eta_space[i_poly].add_quadrature_point_list_by_recursive3_triangle_split_in_xi_eta_space(
	                    quadrature_point_list,
	                    i_thread,
//	                    negative_weight,
	                    face,
	                    included_cell_index_x,
	                    included_cell_index_y,
	                    included_cell_index_z,
	                    global_weight_factor
	                );
	            }
	            else if(integral_method == "triangle-rec4"){
	                //三角形分割で求積点を求める
	                splitted_faces_in_xi_eta_space[i_poly].add_quadrature_point_list_by_recursive4_triangle_split_in_xi_eta_space(
	                    quadrature_point_list,
	                    i_thread,
//	                    negative_weight,
	                    face,
	                    included_cell_index_x,
	                    included_cell_index_y,
	                    included_cell_index_z,
	                    global_weight_factor
	                );
	            }
	        }
		}
    }

// 引数の polygons_in_xi_eta_space をセルの中心で切断する関数
std::vector<polygon_3D> split_polygons_in_xi_eta_space_by_cell_center_3D(
    const cell_face_3D &face,
    const std::vector<polygon_3D> &polygons_in_xi_eta_space,
    std::string split_direction,
    const MY_FLOAT_TYPE cell_length
){
    //split_directionの成分にアクセスするためのインデックス
    int split_direction_index;
    if(split_direction == "x"){
        split_direction_index = 0;
    }
    else if(split_direction == "y"){
        split_direction_index = 1;
    }
    else if(split_direction == "z"){
        split_direction_index = 2;
    }
    //face のvertex position の split_direction成分のうち最大のものと最小のもの
    MY_FLOAT_TYPE max_xyz_value_on_face_vertex, min_xyz_value_on_face_vertex;
    // max_xyz_value_on_face_vertex の計算
    max_xyz_value_on_face_vertex = std::max(face._vertex_list[0]._vertex_pos[split_direction_index], face._vertex_list[1]._vertex_pos[split_direction_index]);
    max_xyz_value_on_face_vertex = std::max(max_xyz_value_on_face_vertex, face._vertex_list[2]._vertex_pos[split_direction_index]);
    max_xyz_value_on_face_vertex = std::max(max_xyz_value_on_face_vertex, face._vertex_list[3]._vertex_pos[split_direction_index]);
    // max_xyz_value_on_face_vertex の計算
    min_xyz_value_on_face_vertex = std::min(face._vertex_list[0]._vertex_pos[split_direction_index], face._vertex_list[1]._vertex_pos[split_direction_index]);
    min_xyz_value_on_face_vertex = std::min(min_xyz_value_on_face_vertex, face._vertex_list[2]._vertex_pos[split_direction_index]);
    min_xyz_value_on_face_vertex = std::min(min_xyz_value_on_face_vertex, face._vertex_list[3]._vertex_pos[split_direction_index]);
    int min_cell_center_index = floor((min_xyz_value_on_face_vertex - 0.5 * cell_length) / cell_length);
    int max_cell_center_index = floor((max_xyz_value_on_face_vertex - 0.5 * cell_length) / cell_length);
    //結果を格納する変数
    std::vector<polygon_3D> splitted_polygons = polygons_in_xi_eta_space;
    for(int i = 0;;++i){
        if(min_cell_center_index + i > max_cell_center_index){
            break;
        }
        splitted_polygons
            = split_cell_face_by_axis_aligned_plane(
                face,
                splitted_polygons,
                (min_cell_center_index + i + 0.5) * cell_length,
                split_direction);
    }
    return splitted_polygons;
}

// 引数の polygons_in_xi_eta_space をセルの中心で切断する関数(参照渡しで値を返すバージョン)
void split_polygons_in_xi_eta_space_by_cell_center_3D_ref(
    const cell_face_3D face,
    std::vector<polygon_3D> &polygons_in_xi_eta_space,
    std::string split_direction,
    const MY_FLOAT_TYPE cell_length
){
    //split_directionの成分にアクセスするためのインデックス
    int split_direction_index;
    if(split_direction == "x"){
        split_direction_index = 0;
    }
    else if(split_direction == "y"){
        split_direction_index = 1;
    }
    else if(split_direction == "z"){
        split_direction_index = 2;
    }
    //face のvertex position の split_direction成分のうち最大のものと最小のもの
    MY_FLOAT_TYPE max_xyz_value_on_face_vertex, min_xyz_value_on_face_vertex;
    // max_xyz_value_on_face_vertex の計算
    max_xyz_value_on_face_vertex = std::max(face._vertex_list[0]._vertex_pos[split_direction_index], face._vertex_list[1]._vertex_pos[split_direction_index]);
    max_xyz_value_on_face_vertex = std::max(max_xyz_value_on_face_vertex, face._vertex_list[2]._vertex_pos[split_direction_index]);
    max_xyz_value_on_face_vertex = std::max(max_xyz_value_on_face_vertex, face._vertex_list[3]._vertex_pos[split_direction_index]);
    // max_xyz_value_on_face_vertex の計算
    min_xyz_value_on_face_vertex = std::min(face._vertex_list[0]._vertex_pos[split_direction_index], face._vertex_list[1]._vertex_pos[split_direction_index]);
    min_xyz_value_on_face_vertex = std::min(min_xyz_value_on_face_vertex, face._vertex_list[2]._vertex_pos[split_direction_index]);
    min_xyz_value_on_face_vertex = std::min(min_xyz_value_on_face_vertex, face._vertex_list[3]._vertex_pos[split_direction_index]);
    int min_cell_center_index = floor((min_xyz_value_on_face_vertex - 0.5 * cell_length) / cell_length);
    int max_cell_center_index = floor((max_xyz_value_on_face_vertex - 0.5 * cell_length) / cell_length);
    //結果を格納する変数
    for(int i = 0;;++i){
        if(min_cell_center_index + i > max_cell_center_index){
            break;
        }
        polygons_in_xi_eta_space
            = split_cell_face_by_axis_aligned_plane(
                face,
                polygons_in_xi_eta_space,
                (min_cell_center_index + i + 0.5) * cell_length,
                split_direction);
    }
}

// 引数の polygons_in_xi_eta_space をセルの面で切断する関数(参照渡しで値を返すバージョン)
void split_polygons_in_xi_eta_space_by_cell_face_3D_ref(
    const cell_face_3D face,
    std::vector<polygon_3D> &polygons_in_xi_eta_space,
    std::string split_direction,
    const MY_FLOAT_TYPE cell_length
){
    //split_directionの成分にアクセスするためのインデックス
    int split_direction_index;
    if(split_direction == "x"){
        split_direction_index = 0;
    }
    else if(split_direction == "y"){
        split_direction_index = 1;
    }
    else if(split_direction == "z"){
        split_direction_index = 2;
    }
    //face のvertex position の split_direction成分のうち最大のものと最小のもの
    MY_FLOAT_TYPE max_xyz_value_on_face_vertex, min_xyz_value_on_face_vertex;
    // max_xyz_value_on_face_vertex の計算
    max_xyz_value_on_face_vertex = std::max(face._vertex_list[0]._vertex_pos[split_direction_index], face._vertex_list[1]._vertex_pos[split_direction_index]);
    max_xyz_value_on_face_vertex = std::max(max_xyz_value_on_face_vertex, face._vertex_list[2]._vertex_pos[split_direction_index]);
    max_xyz_value_on_face_vertex = std::max(max_xyz_value_on_face_vertex, face._vertex_list[3]._vertex_pos[split_direction_index]);
    // max_xyz_value_on_face_vertex の計算
    min_xyz_value_on_face_vertex = std::min(face._vertex_list[0]._vertex_pos[split_direction_index], face._vertex_list[1]._vertex_pos[split_direction_index]);
    min_xyz_value_on_face_vertex = std::min(min_xyz_value_on_face_vertex, face._vertex_list[2]._vertex_pos[split_direction_index]);
    min_xyz_value_on_face_vertex = std::min(min_xyz_value_on_face_vertex, face._vertex_list[3]._vertex_pos[split_direction_index]);
    int min_cell_center_index = floor((min_xyz_value_on_face_vertex) / cell_length);
    int max_cell_center_index = floor((max_xyz_value_on_face_vertex) / cell_length);
    //結果を格納する変数
    for(int i = 0;;++i){
        if(min_cell_center_index + i > max_cell_center_index){
            break;
        }
        polygons_in_xi_eta_space
            = split_cell_face_by_axis_aligned_plane(
                face,
                polygons_in_xi_eta_space,
                (min_cell_center_index + i) * cell_length,
                split_direction);
    }
}

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
) {
    //original polygon を切断した後の polygon のリスト(結果を格納するための変数)
    std::vector<polygon_3D> splitted_polygons;
    //////正規化されたパラメータ空間 \xi - \eta 空間で考える//////
    ////\xi - \eta 空間上で切断を行う直線を求める
    // face の 4つの頂点の位置のsplit_direction成分を格納する(例えばsplit_direction=="x"ならface の 4つの頂点の位置のx成分を格納する)
    MY_FLOAT_TYPE xyz_value_on_cell_vertex[4];
    for(int i_vert = 0; i_vert < 4; ++i_vert){
        if(split_direction == "x"){
            xyz_value_on_cell_vertex[i_vert] = face._vertex_list[i_vert]._vertex_pos[0];
        }
        else if(split_direction == "y"){
            xyz_value_on_cell_vertex[i_vert] = face._vertex_list[i_vert]._vertex_pos[1];
        }
        else if(split_direction == "z"){
            xyz_value_on_cell_vertex[i_vert] = face._vertex_list[i_vert]._vertex_pos[2];
        }
    }
    //各頂点の位置の split_direction 成分 (xyz_value_on_cell_vertexのこと) の値が iso_value より大きいか
    bool xyz_value_on_cell_vertex_is_above_iso_value[4];
    for(int i_vert = 0; i_vert < 4; ++i_vert){
        xyz_value_on_cell_vertex_is_above_iso_value[i_vert] = (xyz_value_on_cell_vertex[i_vert] > iso_value);
    }
    //直線の本数とパターンを判別する0から15までの数(https://en.wikipedia.org/wiki/Marching_squares のLook-up table contour lines の 0 から 15 番を参照)
    int split_case_number = 0;
    if(xyz_value_on_cell_vertex_is_above_iso_value[0]){
        split_case_number += 1;
    }
    if(xyz_value_on_cell_vertex_is_above_iso_value[1]){
        split_case_number += 2;
    }
    if(xyz_value_on_cell_vertex_is_above_iso_value[2]){
        split_case_number += 4;
    }
    if(xyz_value_on_cell_vertex_is_above_iso_value[3]){
        split_case_number += 8;
    }
    //polygon を全く切断しなくて良い場合
    if(split_case_number == 0 || split_case_number == 15){
        return original_polygon_list_in_xi_eta_space;
    }
    //polygon を切断する直線が2本の場合
    else if(split_case_number == 5 || split_case_number == 10){
        //\xi - \eta 空間上の直線を定義する2点を face の辺上でiso value をとる位置から求める
        VEC3_TYPE points_on_line_in_xi_eta_space[4];
        //points_on_line_in_xi_eta_space を設定するためのindex
        int i_point = 0;
        // face の4辺を回るループ
        for(int i_edge = 0; i_edge < 4; ++i_edge){
            ////直線と考えている face の辺の交点を\xi - \eta 空間上で求める
            //考えている辺を構成する頂点の位置を\xi - \eta 空間表示したもの
            VEC3_TYPE vertex_pos_of_iedge_in_xi_eta_space_0;
            VEC3_TYPE vertex_pos_of_iedge_in_xi_eta_space_1;
            if(i_edge == 0){
                vertex_pos_of_iedge_in_xi_eta_space_0[0] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_0[1] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_0[2] = 0.0;
                vertex_pos_of_iedge_in_xi_eta_space_1[0] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_1[1] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_1[2] = 0.0;
                ////線型補間により points_on_line_in_xi_eta_space を求める
                //線型補間の重み
                MY_FLOAT_TYPE c0, c1;
                c0 = (abs(iso_value - xyz_value_on_cell_vertex[1]))/(abs(xyz_value_on_cell_vertex[0] - xyz_value_on_cell_vertex[1]));
                c1 = 1.0 - c0;
                points_on_line_in_xi_eta_space[i_point] = c0 * vertex_pos_of_iedge_in_xi_eta_space_0 + c1 * vertex_pos_of_iedge_in_xi_eta_space_1;
                ++i_point;
            }
            else if(i_edge == 1){
                vertex_pos_of_iedge_in_xi_eta_space_0[0] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_0[1] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_0[2] = 0.0;
                vertex_pos_of_iedge_in_xi_eta_space_1[0] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_1[1] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_1[2] = 0.0;
                ////線型補間により points_on_line_in_xi_eta_space を求める
                //線型補間の重み
                MY_FLOAT_TYPE c0, c1;
                c0 = (abs(iso_value - xyz_value_on_cell_vertex[2]))/(abs(xyz_value_on_cell_vertex[1] - xyz_value_on_cell_vertex[2]));
                c1 = 1.0 - c0;
                points_on_line_in_xi_eta_space[i_point] = c0 * vertex_pos_of_iedge_in_xi_eta_space_0 + c1 * vertex_pos_of_iedge_in_xi_eta_space_1;
                ++i_point;
            }
            else if(i_edge == 2){
                vertex_pos_of_iedge_in_xi_eta_space_0[0] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_0[1] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_0[2] = 0.0;
                vertex_pos_of_iedge_in_xi_eta_space_1[0] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_1[1] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_1[2] = 0.0;
                ////線型補間により points_on_line_in_xi_eta_space を求める
                //線型補間の重み
                MY_FLOAT_TYPE c0, c1;
                c0 = (abs(iso_value - xyz_value_on_cell_vertex[3]))/(abs(xyz_value_on_cell_vertex[2] - xyz_value_on_cell_vertex[3]));
                c1 = 1.0 - c0;
                points_on_line_in_xi_eta_space[i_point] = c0 * vertex_pos_of_iedge_in_xi_eta_space_0 + c1 * vertex_pos_of_iedge_in_xi_eta_space_1;
                ++i_point;
            }
            else if(i_edge == 3){
                vertex_pos_of_iedge_in_xi_eta_space_0[0] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_0[1] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_0[2] = 0.0;
                vertex_pos_of_iedge_in_xi_eta_space_1[0] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_1[1] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_1[2] = 0.0;
                ////線型補間により points_on_line_in_xi_eta_space を求める
                //線型補間の重み
                MY_FLOAT_TYPE c0, c1;
                c0 = (abs(iso_value - xyz_value_on_cell_vertex[0]))/(abs(xyz_value_on_cell_vertex[3] - xyz_value_on_cell_vertex[0]));
                c1 = 1.0 - c0;
                points_on_line_in_xi_eta_space[i_point] = c0 * vertex_pos_of_iedge_in_xi_eta_space_0 + c1 * vertex_pos_of_iedge_in_xi_eta_space_1;
                ++i_point;
            }
        }
        // 切断を行う線の数が2つなのに切断を行う線とエッジの交差点が4つでない場合は正しくない
        if(i_point =! 4){
            std::cerr << "split line number is 2. But number of intersections of lines and edge is not 4."<<std::endl;
        }
        ////4つの, エッジ上の交差点がどのような組み合わせで線を作るかを求める
        //1つ目の線の上の点のリスト
        VEC3_TYPE points_on_line0_in_xi_eta_space[2];
        VEC3_TYPE points_on_line1_in_xi_eta_space[2];
        //cell face 中心での xyz_value
        const MY_FLOAT_TYPE xyz_value_at_cell_center
            = (xyz_value_on_cell_vertex[0] + xyz_value_on_cell_vertex[1] + xyz_value_on_cell_vertex[2] + xyz_value_on_cell_vertex[3]) / 4.0;
            //cell face 中心での xyz_value が iso_value より大きいかどうか
        const bool xyz_value_at_cell_center_is_above_isovalue = (xyz_value_at_cell_center > iso_value);
        if(split_case_number == 5 && !xyz_value_at_cell_center_is_above_isovalue){
            points_on_line0_in_xi_eta_space[0] = points_on_line_in_xi_eta_space[0];
            points_on_line0_in_xi_eta_space[1] = points_on_line_in_xi_eta_space[3];
            points_on_line1_in_xi_eta_space[0] = points_on_line_in_xi_eta_space[1];
            points_on_line1_in_xi_eta_space[1] = points_on_line_in_xi_eta_space[2];
        }
        else if(split_case_number == 5 && xyz_value_at_cell_center_is_above_isovalue){
            points_on_line0_in_xi_eta_space[0] = points_on_line_in_xi_eta_space[0];
            points_on_line0_in_xi_eta_space[1] = points_on_line_in_xi_eta_space[1];
            points_on_line1_in_xi_eta_space[0] = points_on_line_in_xi_eta_space[2];
            points_on_line1_in_xi_eta_space[1] = points_on_line_in_xi_eta_space[3];
        }
        else if(split_case_number == 10 && !xyz_value_at_cell_center_is_above_isovalue){
            points_on_line0_in_xi_eta_space[0] = points_on_line_in_xi_eta_space[0];
            points_on_line0_in_xi_eta_space[1] = points_on_line_in_xi_eta_space[1];
            points_on_line1_in_xi_eta_space[0] = points_on_line_in_xi_eta_space[2];
            points_on_line1_in_xi_eta_space[1] = points_on_line_in_xi_eta_space[3];
        }
        else if(split_case_number == 10 && xyz_value_at_cell_center_is_above_isovalue){
            points_on_line0_in_xi_eta_space[0] = points_on_line_in_xi_eta_space[0];
            points_on_line0_in_xi_eta_space[1] = points_on_line_in_xi_eta_space[3];
            points_on_line1_in_xi_eta_space[0] = points_on_line_in_xi_eta_space[1];
            points_on_line1_in_xi_eta_space[1] = points_on_line_in_xi_eta_space[2];
        }

        //// points_on_line_in_xi_eta_space の2点から定義される直線で original_polygon_list_in_xi_eta_space の中の各polygon を分割する
        const int polygon_num = original_polygon_list_in_xi_eta_space.size();
        //0番目の線で切断
//        std::cout<<"split by line 0"<<std::endl;
        std::vector<polygon_3D> splitted_polygons_tmp;
        for(int i_poly = 0; i_poly < polygon_num; ++i_poly){
            std::vector<polygon_3D> splitted_polygons_of_i_poly
                = split_polygon_in_xi_eta_space_by_line(
                    original_polygon_list_in_xi_eta_space[i_poly],
                    points_on_line0_in_xi_eta_space[0],
                    points_on_line0_in_xi_eta_space[1]
                );
            const int num_split = splitted_polygons_of_i_poly.size();
            for(int i_split = 0; i_split < num_split; ++i_split){
                splitted_polygons_tmp.push_back(splitted_polygons_of_i_poly[i_split]);
            }
        }
        //1番目の線で切断

//        std::cout<<"split by line 1"<<std::endl;
        for(int i_poly = 0; i_poly < splitted_polygons_tmp.size(); ++i_poly){
            std::vector<polygon_3D> splitted_polygons_of_i_poly
                = split_polygon_in_xi_eta_space_by_line(
                    splitted_polygons_tmp[i_poly],
                    points_on_line1_in_xi_eta_space[0],
                    points_on_line1_in_xi_eta_space[1]
                );
//            std::cout<<"splitted_polygons_of_i_poly.size(): "<<splitted_polygons_of_i_poly.size()<<std::endl;
            const int num_split = splitted_polygons_of_i_poly.size();
            for(int i_split = 0; i_split < num_split; ++i_split){
                splitted_polygons.push_back(splitted_polygons_of_i_poly[i_split]);
            }
        }
    }
    //polygon を切断する直線が1本の場合
    else{
        //\xi - \eta 空間上の直線を定義する2点を face の辺上でiso value をとる位置から求める
        VEC3_TYPE points_on_line_in_xi_eta_space[2];
        //points_on_line_in_xi_eta_space を設定するためのindex
        int i_point = 0;
        // face の4辺を回るループ
        for(int i_edge = 0; i_edge < 4; ++i_edge){
            //直線と考えている face の辺が交点を持たない場合は次の辺に移る
            if(xyz_value_on_cell_vertex_is_above_iso_value[i_edge] == xyz_value_on_cell_vertex_is_above_iso_value[(i_edge + 1) % 4]){
                continue;
            }
            ////直線と考えている face の辺の交点を\xi - \eta 空間上で求める
            //考えている辺を構成する頂点の位置を\xi - \eta 空間表示したもの
            VEC3_TYPE vertex_pos_of_iedge_in_xi_eta_space_0;
            VEC3_TYPE vertex_pos_of_iedge_in_xi_eta_space_1;
            if(i_edge == 0){
                vertex_pos_of_iedge_in_xi_eta_space_0[0] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_0[1] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_0[2] = 0.0;
                vertex_pos_of_iedge_in_xi_eta_space_1[0] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_1[1] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_1[2] = 0.0;
                ////線型補間により points_on_line_in_xi_eta_space を求める
                //線型補間の重み
                MY_FLOAT_TYPE c0, c1;
                c0 = (abs(iso_value - xyz_value_on_cell_vertex[1]))/(abs(xyz_value_on_cell_vertex[0] - xyz_value_on_cell_vertex[1]));
                c1 = 1.0 - c0;
                points_on_line_in_xi_eta_space[i_point] = c0 * vertex_pos_of_iedge_in_xi_eta_space_0 + c1 * vertex_pos_of_iedge_in_xi_eta_space_1;
                ++i_point;
            }
            else if(i_edge == 1){
                vertex_pos_of_iedge_in_xi_eta_space_0[0] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_0[1] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_0[2] = 0.0;
                vertex_pos_of_iedge_in_xi_eta_space_1[0] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_1[1] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_1[2] = 0.0;
                ////線型補間により points_on_line_in_xi_eta_space を求める
                //線型補間の重み
                MY_FLOAT_TYPE c0, c1;
                c0 = (abs(iso_value - xyz_value_on_cell_vertex[2]))/(abs(xyz_value_on_cell_vertex[1] - xyz_value_on_cell_vertex[2]));
                c1 = 1.0 - c0;
                points_on_line_in_xi_eta_space[i_point] = c0 * vertex_pos_of_iedge_in_xi_eta_space_0 + c1 * vertex_pos_of_iedge_in_xi_eta_space_1;
                ++i_point;
            }
            else if(i_edge == 2){
                vertex_pos_of_iedge_in_xi_eta_space_0[0] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_0[1] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_0[2] = 0.0;
                vertex_pos_of_iedge_in_xi_eta_space_1[0] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_1[1] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_1[2] = 0.0;
                ////線型補間により points_on_line_in_xi_eta_space を求める
                //線型補間の重み
                MY_FLOAT_TYPE c0, c1;
                c0 = (abs(iso_value - xyz_value_on_cell_vertex[3]))/(abs(xyz_value_on_cell_vertex[2] - xyz_value_on_cell_vertex[3]));
                c1 = 1.0 - c0;
                points_on_line_in_xi_eta_space[i_point] = c0 * vertex_pos_of_iedge_in_xi_eta_space_0 + c1 * vertex_pos_of_iedge_in_xi_eta_space_1;
                ++i_point;
            }
            else if(i_edge == 3){
                vertex_pos_of_iedge_in_xi_eta_space_0[0] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_0[1] = 1.0; vertex_pos_of_iedge_in_xi_eta_space_0[2] = 0.0;
                vertex_pos_of_iedge_in_xi_eta_space_1[0] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_1[1] = 0.0; vertex_pos_of_iedge_in_xi_eta_space_1[2] = 0.0;
                ////線型補間により points_on_line_in_xi_eta_space を求める
                //線型補間の重み
                MY_FLOAT_TYPE c0, c1;
                c0 = (abs(iso_value - xyz_value_on_cell_vertex[0]))/(abs(xyz_value_on_cell_vertex[3] - xyz_value_on_cell_vertex[0]));
                c1 = 1.0 - c0;
                points_on_line_in_xi_eta_space[i_point] = c0 * vertex_pos_of_iedge_in_xi_eta_space_0 + c1 * vertex_pos_of_iedge_in_xi_eta_space_1;
                ++i_point;
            }
        }
        //// points_on_line_in_xi_eta_space の2点から定義される直線で original_polygon_list_in_xi_eta_space の中の各polygon を分割する
        const int polygon_num = original_polygon_list_in_xi_eta_space.size();
        for(int i_poly = 0; i_poly < polygon_num; ++i_poly){
            std::vector<polygon_3D> splitted_polygons_of_i_poly
                = split_polygon_in_xi_eta_space_by_line(
                    original_polygon_list_in_xi_eta_space[i_poly],
                    points_on_line_in_xi_eta_space[0],
                    points_on_line_in_xi_eta_space[1]
                );
            const int num_split = splitted_polygons_of_i_poly.size();
            for(int i_split = 0; i_split < num_split; ++i_split){
                splitted_polygons.push_back(splitted_polygons_of_i_poly[i_split]);
            }
        }
    }
    return splitted_polygons;


}

// split_cell_face_by_axis_aligned_plane のヘルパー関数
// point_on_line_in_xi_eta_space_0 と point_on_line_in_xi_eta_space_1 を通る直線で polygon を切断し、切断後のpolygonのリストを返す
// OK
std::vector<polygon_3D> split_polygon_in_xi_eta_space_by_line(
    const polygon_3D polygonin_xi_eta_space,
    const VEC3_TYPE point_on_line_in_xi_eta_space_0,
    const VEC3_TYPE point_on_line_in_xi_eta_space_1
){
    //point_on_line_in_xi_eta_space_0 と point_on_line_in_xi_eta_space_1 を通る直線 と polygonのedge の交点の情報を持つための構造体
    struct intersection_of_line_and_edge{
        int vertex_index_of_polygon_edge_0;
        int vertex_index_of_polygon_edge_1;
        VEC3_TYPE intersection_position;
    };
    //polygon のエッジと直線の交差点の情報を格納していくリスト
    std::vector<intersection_of_line_and_edge> intersection_list;
    ////各エッジが考えている直線と交差するか判定する
    const int edge_num = polygonin_xi_eta_space._vertex_pos_list.size();
    for(int i_edge = 0 ;i_edge < edge_num; ++i_edge){
        ////微小量での割り算が発生しないように直線の方程式の作り方を選ぶ
        // \eta = slope * \xi + intersept の表式を使う場合
        if(abs(point_on_line_in_xi_eta_space_1[0] - point_on_line_in_xi_eta_space_0[0])
            > abs(point_on_line_in_xi_eta_space_1[1] - point_on_line_in_xi_eta_space_0[1])) {
                //直線の傾き
                const MY_FLOAT_TYPE slope
                    = (point_on_line_in_xi_eta_space_1[1] - point_on_line_in_xi_eta_space_0[1])
                    / (point_on_line_in_xi_eta_space_1[0] - point_on_line_in_xi_eta_space_0[0]);
                //直線の切片
                const MY_FLOAT_TYPE intercept = point_on_line_in_xi_eta_space_0[1] - slope * point_on_line_in_xi_eta_space_0[0];
                //polygon の i_edge 番目の edge の2つの頂点の位置
                VEC3_TYPE vertex_pos_of_i_edge_0 = polygonin_xi_eta_space._vertex_pos_list[i_edge];
                VEC3_TYPE vertex_pos_of_i_edge_1 = polygonin_xi_eta_space._vertex_pos_list[(i_edge + 1) % edge_num];
                //線型系を解いて交点を求める
                Eigen::Matrix2d coefficient;
                coefficient(0, 0) = vertex_pos_of_i_edge_1[0] - vertex_pos_of_i_edge_0[0]; coefficient(0, 1) = -1;
                coefficient(1, 0) = vertex_pos_of_i_edge_1[1] - vertex_pos_of_i_edge_0[1]; coefficient(1, 1) = -slope;
                Eigen::Vector2d rhs_vector;
                rhs_vector(0) = - vertex_pos_of_i_edge_0[0];
                rhs_vector(1) = intercept - vertex_pos_of_i_edge_0[1];
                Eigen::Vector2d t_xi = coefficient.inverse() * rhs_vector;
                //交差しない場合
                if(t_xi[0] < 0 || t_xi[0] > 1){
                    continue;
                }
                //交差する場合, 交差点の情報を intersection_list に格納していく
                else{
                    VEC3_TYPE intersection_position = vertex_pos_of_i_edge_0 + t_xi[0] * (vertex_pos_of_i_edge_1 - vertex_pos_of_i_edge_0);
                    // 交差点がinf や nan を持つ場合は除外
                    bool have_nan_inf_value = false;
                    for(int i_xyz = 0; i_xyz < 3; ++i_xyz){
                        if(std::isnan(intersection_position[i_xyz])){
                            have_nan_inf_value = true;
                        }
                        if(std::isinf(intersection_position[i_xyz])){
                            have_nan_inf_value = true;
                        }
                    }
                    if(!have_nan_inf_value){
                        intersection_of_line_and_edge intersection;
                        intersection.vertex_index_of_polygon_edge_0 = i_edge;
                        intersection.vertex_index_of_polygon_edge_1 = (i_edge + 1) % edge_num;
                        intersection.intersection_position = intersection_position;
                        intersection_list.push_back(intersection);
                    }
                }
        }
        // \xi = slope * \eta + intersept の表式を使う場合
        else {
            //直線の傾き
            const MY_FLOAT_TYPE slope
                = (point_on_line_in_xi_eta_space_1[0] - point_on_line_in_xi_eta_space_0[0])
                / (point_on_line_in_xi_eta_space_1[1] - point_on_line_in_xi_eta_space_0[1]);
            //直線の切片
            const MY_FLOAT_TYPE intercept = point_on_line_in_xi_eta_space_0[0] - slope * point_on_line_in_xi_eta_space_0[1];
            //polygon の i_edge 番目の edge の2つの頂点の位置
            VEC3_TYPE vertex_pos_of_i_edge_0 = polygonin_xi_eta_space._vertex_pos_list[i_edge];
            VEC3_TYPE vertex_pos_of_i_edge_1 = polygonin_xi_eta_space._vertex_pos_list[(i_edge + 1) % edge_num];
            //線型系を解いて交点を求める
            Eigen::Matrix2d coefficient;
            coefficient(0, 0) = vertex_pos_of_i_edge_1[0] - vertex_pos_of_i_edge_0[0]; coefficient(0, 1) = -slope;
            coefficient(1, 0) = vertex_pos_of_i_edge_1[1] - vertex_pos_of_i_edge_0[1]; coefficient(1, 1) = -1;
            Eigen::Vector2d rhs_vector;
            rhs_vector(0) = intercept - vertex_pos_of_i_edge_0[0];
            rhs_vector(1) = - vertex_pos_of_i_edge_0[1];
            Eigen::Vector2d t_eta = coefficient.inverse() * rhs_vector;
            //交差しない場合
            if (t_eta[0] < 0 || t_eta[0] > 1){
                continue;
            }
            //交差する場合, 交差点の情報を intersection_list に格納していく
            else {
                VEC3_TYPE intersection_position = vertex_pos_of_i_edge_0 + t_eta[0] * (vertex_pos_of_i_edge_1 - vertex_pos_of_i_edge_0);
                // 交差点がinf や nan を持つ場合は除外
                bool have_nan_inf_value = false;
                for(int i_xyz = 0; i_xyz < 3; ++i_xyz){
                    if(std::isnan(intersection_position[i_xyz])){
                        have_nan_inf_value = true;
                    }
                    if(std::isinf(intersection_position[i_xyz])){
                        have_nan_inf_value = true;
                    }
                }
                if(!have_nan_inf_value){
                    intersection_of_line_and_edge intersection;
                    intersection.vertex_index_of_polygon_edge_0 = i_edge;
                    intersection.vertex_index_of_polygon_edge_1 = (i_edge + 1) % edge_num;
                    intersection.intersection_position = intersection_position;
                    intersection_list.push_back(intersection);
                }
            }
        }
    }
    ////直線とポリゴンのエッジの交点のリスト intersection_list からポリゴン polygonin_xi_eta_space を切断する
    //結果を格納する変数
    std::vector<polygon_3D> splitted_polygons;
    //直線とポリゴンが全く交差しない場合は元のポリゴンを1つ格納したリストをそのまま返す
    if(intersection_list.size() == 0){
        splitted_polygons.push_back(polygonin_xi_eta_space);
        return splitted_polygons;
    }
    //交差点の数が1の場合(幾何学的に間違っている or 直線で切断する位置がポリゴンの端過ぎて数値誤差が出てる) や
    //交差点の数が2より大きい場合(polygon が凸ではない)は今回のシミュレーションでは起こってはいけないのでエラーとして別で処理する(一応もとのポリゴンを切断せずに返す)
    else if(intersection_list.size() == 1 || intersection_list.size() > 2){
//        std::cerr << "intersection point number is " <<intersection_list.size()<<" . not 2."<<std::endl;
        splitted_polygons.push_back(polygonin_xi_eta_space);
        return splitted_polygons;
    }
    //交差点の数が2で正しい場合
    else{
        ////切断後1つ目のポリゴン
        polygon_3D polygon_0;
        //// polygon_0 にvertex を追加していく
        // polygon_0 に追加する1つ目の vertex のindex
        int added_vertex_index = std::max(intersection_list[0].vertex_index_of_polygon_edge_0, intersection_list[0].vertex_index_of_polygon_edge_1);
        // 0 と num_vertex-1 を結ぶエッジの場合, 0のインデックスは num_vertex を表すのでadded_vertex_indexには0番目(つまりnum_vertex番目)にインデックスを追加する
        if((intersection_list[0].vertex_index_of_polygon_edge_0 == 0
            &&intersection_list[0].vertex_index_of_polygon_edge_1 == polygonin_xi_eta_space._vertex_pos_list.size() - 1)
        ||(intersection_list[0].vertex_index_of_polygon_edge_1 == 0
            &&intersection_list[0].vertex_index_of_polygon_edge_0 == polygonin_xi_eta_space._vertex_pos_list.size() - 1)){
            added_vertex_index = 0;
        }
        while(true){
            polygon_0._vertex_pos_list.push_back(polygonin_xi_eta_space._vertex_pos_list[added_vertex_index]);
            bool added_vertex_is_included_in_another_split_edge
                =( added_vertex_index == intersection_list[1].vertex_index_of_polygon_edge_0
                || added_vertex_index == intersection_list[1].vertex_index_of_polygon_edge_1);
            if(added_vertex_is_included_in_another_split_edge){
                break;
            }
            added_vertex_index = (added_vertex_index + 1) % edge_num;
        }
        //交差点を追加
        polygon_0._vertex_pos_list.push_back(intersection_list[1].intersection_position);
        polygon_0._vertex_pos_list.push_back(intersection_list[0].intersection_position);

        ////切断後2つ目のポリゴン
        polygon_3D polygon_1;
        //// polygon_1 にvertex を追加していく
        // polygon_1 に追加する1つ目の vertex のindex
        added_vertex_index = std::max(intersection_list[1].vertex_index_of_polygon_edge_0, intersection_list[1].vertex_index_of_polygon_edge_1);
        // 0 と num_vertex-1 を結ぶエッジの場合, 0のインデックスは num_vertex を表すのでadded_vertex_indexには0番目(つまりnum_vertex番目)にインデックスを追加する
        if((intersection_list[1].vertex_index_of_polygon_edge_0 == 0
            &&intersection_list[1].vertex_index_of_polygon_edge_1 == polygonin_xi_eta_space._vertex_pos_list.size() - 1)
        ||(intersection_list[1].vertex_index_of_polygon_edge_1 == 0
            &&intersection_list[1].vertex_index_of_polygon_edge_0 == polygonin_xi_eta_space._vertex_pos_list.size() - 1)){
            added_vertex_index = 0;
        }
        while(true){
            polygon_1._vertex_pos_list.push_back(polygonin_xi_eta_space._vertex_pos_list[added_vertex_index]);
            bool added_vertex_is_included_in_another_split_edge
                =( added_vertex_index == intersection_list[0].vertex_index_of_polygon_edge_0
                || added_vertex_index == intersection_list[0].vertex_index_of_polygon_edge_1);
            if(added_vertex_is_included_in_another_split_edge){
                break;
            }
            added_vertex_index = (added_vertex_index + 1) % edge_num;
        }
        //交差点を追加
        polygon_1._vertex_pos_list.push_back(intersection_list[0].intersection_position);
        polygon_1._vertex_pos_list.push_back(intersection_list[1].intersection_position);

        //polygon_0 と polygon_1 を結果を格納するための変数に追加
        splitted_polygons.push_back(polygon_0);
        splitted_polygons.push_back(polygon_1);
    }
    return splitted_polygons;
}

}//namespace smoke_simulation
