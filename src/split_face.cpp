#include "split_face.h"

namespace smoke_simulation {
	std::vector<cell_face> split_face(const Grid& all_grid, const cell_face& face, const std::string split_method){
		std::vector<cell_face> splitted_faces;
		// x軸, y軸のどちらにも沿った切断
		if(split_method == "xy-AxisAligned"){
			std::vector<cell_face> splitted_faces_x = all_grid.split_cell_face_by_axis(face, "x");
			for (int i_face_0 = 0; i_face_0 < splitted_faces_x.size(); ++i_face_0) {
				std::vector<cell_face> splitted_faces_xy_tmp = all_grid.split_cell_face_by_axis(splitted_faces_x[i_face_0], "y");
				for (int i_face_1 = 0; i_face_1 < splitted_faces_xy_tmp.size(); ++i_face_1) {
					splitted_faces.push_back(splitted_faces_xy_tmp[i_face_1]);
				}
			}
		}
		// セルの面でx軸に沿った切断
		else if(split_method == "x-AxisAligned"){
			splitted_faces = all_grid.split_cell_face_by_axis(face, "x");
		}
		// セルの面でy軸に沿った切断
		else if(split_method == "y-AxisAligned"){
			splitted_faces = all_grid.split_cell_face_by_axis(face, "y");
		}
		// セル中心でx軸に沿った切断(未実装)
//		else if(split_method == "x-CellFaceAligned"){
//			splitted_faces = all_grid.split_cell_face_by_cell_center(face, "x");
//		}
		// セル中心でy軸に沿った切断
		else if(split_method == "y-CellFaceAligned"){
			splitted_faces = all_grid.split_cell_face_by_cell_center(face, "y");
		}
		// セル中心でx軸, y軸のどちらにも沿った切断(未実装)
//		else if(split_method == "xy-CellFaceAligned"){
//			std::vector<cell_face> splitted_faces_x = all_grid.split_cell_face_by_cell_center(face, "x");
//			for (int i_face_0 = 0; i_face_0 < splitted_faces_x.size(); ++i_face_0) {
//				std::vector<cell_face> splitted_faces_xy_tmp = all_grid.split_cell_face_by_cell_center(splitted_faces_x[i_face_0], "y");
//				for (int i_face_1 = 0; i_face_1 < splitted_faces_xy_tmp.size(); ++i_face_1) {
//					splitted_faces.push_back(splitted_faces_xy_tmp[i_face_1]);
//				}
//			}
//		}
		// セルの面でx軸に沿った切断をした後にセル中心でy軸に沿った切断
		else if(split_method == "x-AxisAligned_y-CellFaceAligned"){
			std::vector<cell_face> splitted_faces_x = all_grid.split_cell_face_by_axis(face, "x");
			for (int i_face_0 = 0; i_face_0 < splitted_faces_x.size(); ++i_face_0) {
				std::vector<cell_face> splitted_faces_xy_tmp = all_grid.split_cell_face_by_cell_center(splitted_faces_x[i_face_0], "y");
				for (int i_face_1 = 0; i_face_1 < splitted_faces_xy_tmp.size(); ++i_face_1) {
					splitted_faces.push_back(splitted_faces_xy_tmp[i_face_1]);
				}
			}
		}
		//uniform に切断
		//splitted_faces = all_grid.split_cell_face_uniformly(face, physical_const::kSplit_num_uniformly);
		// セルの頂点の位置で切断
		// 考えるcellを構成するvertex のindex
		/*
		std::vector<int> vertex_index_1{ ix + dx[0], iy + dy[0] };
		std::vector<int> vertex_index_2{ ix + dx[1], iy + dy[1] };
		std::vector<int> vertex_index_3{ ix + dx[2], iy + dy[2] };
		std::vector<int> vertex_index_4{ ix + dx[3], iy + dy[3] };
		std::vector<cell_vertex> backtrace_cell_vertex_list = calc_backtrace_vertex_of_cell(all_grid, vertex_index_1, vertex_index_2, vertex_index_3, vertex_index_4, time_step_length);
		std::vector<MY_FLOAT_TYPE> x_coordinates_of_backtrace_cell_vertex;
		x_coordinates_of_backtrace_cell_vertex.push_back(backtrace_cell_vertex_list[0]._vertex_pos[0]);
		x_coordinates_of_backtrace_cell_vertex.push_back(backtrace_cell_vertex_list[1]._vertex_pos[0]);
		x_coordinates_of_backtrace_cell_vertex.push_back(backtrace_cell_vertex_list[2]._vertex_pos[0]);
		x_coordinates_of_backtrace_cell_vertex.push_back(backtrace_cell_vertex_list[3]._vertex_pos[0]);

		std::sort(x_coordinates_of_backtrace_cell_vertex.begin(), x_coordinates_of_backtrace_cell_vertex.end());
		std::vector<MY_FLOAT_TYPE> x_coordinates_of_split_positions;
		x_coordinates_of_split_positions.push_back(x_coordinates_of_backtrace_cell_vertex[1]);
		x_coordinates_of_split_positions.push_back(x_coordinates_of_backtrace_cell_vertex[2]);
		splitted_faces = all_grid.split_cell_face_by_axis_aligned_plane(face, x_coordinates_of_split_positions, "y");
		*/

		//面を分割しない
		if(split_method == "no-split"){
			splitted_faces.push_back(face);
		}
		return splitted_faces;
	}
}
