#include "grid.h"

#include <algorithm>

#include "utils.h"
#include "physical_const.h"
#include "define_float_type.h"
#include "calc_psi.h"
#include "linear_interpolation_2d.h"
#include "WENO6_interpolation_2d.h"
#include "linear_interpolation_1d.h"
#include "WENO6_interpolation_1d.h"

namespace smoke_simulation {

	Grid::Grid(int nx, int ny, MY_FLOAT_TYPE cell_length, const std::string interpolation_method)
	: Grid_num_x(nx),
	  Grid_num_y(ny),
	  min_pos_x(0.0),
	  min_pos_y(0.0),
	  _cell_length(cell_length),
	  _cell_volume(cell_length * cell_length)//,
//	  _interpolation_method(interpolation_method)
	  {
		//substance_densityのメモリ確保
		substance_density.resize(nx * ny);
		//substance_density のy方向の累積和のメモリ確保
		psi_substance_density_cell_face_x.resize((nx + 1) * ny);
		psi_substance_density_cell_face_y.resize(nx * (ny + 1));

		psi_velocity_cell_vertex_x.resize((nx + 1) * (ny + 1));
		psi_velocity_cell_center_y.resize(nx * (ny + 2));

		//velocityのメモリ確保
		velocity_in_voxel_face_x.resize((nx + 1) * ny);
		velocity_in_voxel_face_y.resize(nx * (ny + 1));

		//pressureのメモリ確保
		pressure.resize(nx * ny);

		//セル体積のメモリ確保
		cell_volume_cell_center.resize(nx * ny);
		cell_volume_cell_face_x.resize((nx + 1) * ny);
		cell_volume_cell_face_y.resize(nx * (ny + 1));
	}

	Grid::~Grid() {
	}

	//(vertex_index_x, vertex_index_y)番目の頂点での速度を、隣り合うface の値からinterpolationによって計算
	VEC3_TYPE Grid::calc_vertex_velocity(int vertex_index_x, int vertex_index_y) const {
		MY_FLOAT_TYPE vertex_velocity_x, vertex_velocity_y;
		if (vertex_index_x == 0) {
			vertex_velocity_y = velocity_in_voxel_face_y[get_voxel_face_index_y(vertex_index_x, vertex_index_y, Grid_num_x, Grid_num_y)];
		}
		else if (vertex_index_x == Grid_num_x) {
			vertex_velocity_y = velocity_in_voxel_face_y[get_voxel_face_index_y(vertex_index_x - 1, vertex_index_y, Grid_num_x, Grid_num_y)];
		}
		else {
			vertex_velocity_y
				= (velocity_in_voxel_face_y[get_voxel_face_index_y(vertex_index_x - 1, vertex_index_y, Grid_num_x, Grid_num_y)]
					+ velocity_in_voxel_face_y[get_voxel_face_index_y(vertex_index_x, vertex_index_y, Grid_num_x, Grid_num_y)])
				/ 2.0;
		}

		if (vertex_index_y == 0) {
			vertex_velocity_x = velocity_in_voxel_face_x[get_voxel_face_index_x(vertex_index_x, vertex_index_y, Grid_num_x, Grid_num_y)];
		}
		else if (vertex_index_y == Grid_num_y) {
			vertex_velocity_x = velocity_in_voxel_face_x[get_voxel_face_index_x(vertex_index_x, vertex_index_y - 1, Grid_num_x, Grid_num_y)];
		}
		else {
			vertex_velocity_x
				= (velocity_in_voxel_face_x[get_voxel_face_index_x(vertex_index_x, vertex_index_y, Grid_num_x, Grid_num_y)]
					+ velocity_in_voxel_face_x[get_voxel_face_index_x(vertex_index_x, vertex_index_y - 1, Grid_num_x, Grid_num_y)])
				/ 2.0;
		}
		return VEC3_TYPE(vertex_velocity_x, vertex_velocity_y, 0.0);
	}

	//position での psi(deisity) の値を周辺の cell face での値から interpolation して求める。
	VEC3_TYPE Grid::calc_psi_substance_density_by_interpolation_cell_face(
		VEC3_TYPE position,
		const std::string interpolation_method,
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density
	) const {
		//結果を格納する変数(それぞれpsi のx, y, z成分)
		MY_FLOAT_TYPE interpolated_psi_x, interpolated_psi_y, interpolated_psi_z;
		if(use_integral) {
			// psi のうち x 成分の計算
			interpolated_psi_x = 0.0;
			if (interpolation_method == "linear") {
				interpolated_psi_y =  linear_interpolation_2D_cell_face_y_values_use_integral(
					position,
					psi_substance_density_cell_face_y,
					substance_density,
					*this,
					num_gauss_quadrature_point_for_integrate_density
				);
/*
				interpolated_psi_y =  linear_interpolation_2D_cell_face_y_values(
					position,
					psi_substance_density_cell_face_y,
					*this
				);
*/
			}
			else if (interpolation_method == "WENO6") {
//				interpolated_psi_y = WENO6_interpolation_2D_cell_face_y_values_use_integral(
//				    position,
//				    psi_substance_density_cell_face_y,
//				    *this
//				);
			}
			else if(interpolation_method == "1Dy_linear"){
				interpolated_psi_y = linear_interpolation_1Dy_cell_face_y_values(
				    position,
				    psi_substance_density_cell_face_y,
				    *this
				);
			}
			else if(interpolation_method == "1Dy_WENO6"){
				interpolated_psi_y = WENO6_interpolation_1Dy_cell_face_y_values(
				    position,
				    psi_substance_density_cell_face_y,
				    *this
				);
			}
			// psi のうち z 成分の計算
			interpolated_psi_z = 0.0;
		}
		else{
			// psi のうち x 成分の計算
			interpolated_psi_x = 0.0;
			if (interpolation_method == "linear") {
				interpolated_psi_y =  linear_interpolation_2D_cell_face_y_values(
					position,
					psi_substance_density_cell_face_y,
					*this
				);
			}
			else if (interpolation_method == "WENO6") {
				interpolated_psi_y = WENO6_interpolation_2D_cell_face_y_values(
				    position,
				    psi_substance_density_cell_face_y,
				    *this
				);
			}
			else if(interpolation_method == "1Dy_linear"){
				interpolated_psi_y = linear_interpolation_1Dy_cell_face_y_values(
				    position,
				    psi_substance_density_cell_face_y,
				    *this
				);
			}
			else if(interpolation_method == "1Dy_WENO6"){
				interpolated_psi_y = WENO6_interpolation_1Dy_cell_face_y_values(
				    position,
				    psi_substance_density_cell_face_y,
				    *this
				);
			}
			// psi のうち z 成分の計算
			interpolated_psi_z = 0.0;
		}
		if (physical_const::kPsi_definition == "x") {
			return VEC3_TYPE(interpolated_psi_x, interpolated_psi_y, interpolated_psi_z);
		}
		else if (physical_const::kPsi_definition == "y") {
			return VEC3_TYPE(interpolated_psi_x, interpolated_psi_y, interpolated_psi_z);
		}
		else if (physical_const::kPsi_definition == "xy") {
			return VEC3_TYPE(interpolated_psi_x / 2.0, interpolated_psi_y / 2.0, interpolated_psi_z);
		}
	}

	//position での psi(deisity) の値を周辺の cell face での値から interpolation して求める。
	VEC3_TYPE Grid::calc_psi_velocity_x_by_interpolation(
		VEC3_TYPE position,
		const std::string interpolation_method,
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density

	) const {
		//結果を格納する変数(それぞれpsi のx, y, z成分)
		MY_FLOAT_TYPE interpolated_psi_x, interpolated_psi_y, interpolated_psi_z;
		if(use_integral) {
			// psi のうち x 成分の計算
			interpolated_psi_x = 0.0;
			if (interpolation_method == "linear") {
/*
				interpolated_psi_y =  linear_interpolation_2D_psi_velocity_x_use_integral(
					position,
					psi_velocity_cell_vertex_x,
					velocity_in_voxel_face_x,
					*this,
					num_gauss_quadrature_point_for_integrate_density
				);
*/
				interpolated_psi_y =  linear_interpolation_2D_psi_velocity_x(
					position,
					psi_velocity_cell_vertex_x,
					*this
				);
			}
			else if (interpolation_method == "WENO6") {
	//			interpolated_psi_y = WENO6_interpolation_2D_psi_velocity_x_use_integral(
	//				position,
	//				psi_velocity_cell_vertex_x,
	//				*this
	//			);
			}
			else if(interpolation_method == "1Dy_linear"){
				interpolated_psi_y = linear_interpolation_1Dy_psi_velocity_x(
					position,
					psi_velocity_cell_vertex_x,
					*this
				);
			}
			else if(interpolation_method == "1Dy_WENO6"){
	//			interpolated_psi_y = WENO6_interpolation_1Dy_psi_velocity_x(
	//				position,
	//				psi_velocity_cell_vertex_x,
	//				*this
	//			);
			}

			// psi のうち z 成分の計算
			interpolated_psi_z = 0.0;
		}
		else{
			// psi のうち x 成分の計算
			interpolated_psi_x = 0.0;

			if (interpolation_method == "linear") {
				interpolated_psi_y =  linear_interpolation_2D_psi_velocity_x(
					position,
					psi_velocity_cell_vertex_x,
					*this
				);
			}
			else if (interpolation_method == "WENO6") {
	//			interpolated_psi_y = WENO6_interpolation_2D_psi_velocity_x(
	//				position,
	//				psi_velocity_cell_vertex_x,
	//				*this
	//			);
			}
			else if(interpolation_method == "1Dy_linear"){
				interpolated_psi_y = linear_interpolation_1Dy_psi_velocity_x(
					position,
					psi_velocity_cell_vertex_x,
					*this
				);
			}
			else if(interpolation_method == "1Dy_WENO6"){
	//			interpolated_psi_y = WENO6_interpolation_1Dy_psi_velocity_x(
	//				position,
	//				psi_velocity_cell_vertex_x,
	//				*this
	//			);
			}

			// psi のうち z 成分の計算
			interpolated_psi_z = 0.0;
		}
		return VEC3_TYPE(interpolated_psi_x, interpolated_psi_y, interpolated_psi_z);
	}

	//position での psi(deisity) の値を周辺の cell face での値から interpolation して求める。
	VEC3_TYPE Grid::calc_psi_velocity_y_by_interpolation(
		VEC3_TYPE position,
		const std::string interpolation_method,
		const bool use_integral,
		const int num_gauss_quadrature_point_for_integrate_density
	) const {
		//結果を格納する変数(それぞれpsi のx, y, z成分)
		MY_FLOAT_TYPE interpolated_psi_x, interpolated_psi_y, interpolated_psi_z;
		// psi のうち x 成分の計算
		interpolated_psi_x = 0.0;
		if(use_integral) {
			if (interpolation_method == "linear") {
/*
				interpolated_psi_y =  linear_interpolation_2D_psi_velocity_y_use_integral(
					position,
					psi_velocity_cell_center_y,
					velocity_in_voxel_face_y,
					*this,
					num_gauss_quadrature_point_for_integrate_density
				);
*/
				interpolated_psi_y =  linear_interpolation_2D_psi_velocity_y(
					position,
					psi_velocity_cell_center_y,
					*this
				);
			}
			else if (interpolation_method == "WENO6") {
	//			interpolated_psi_y = WENO6_interpolation_2D_psi_velocity_x(
	//				position,
	//				psi_velocity_cell_center_y,
	//				*this
	//			);
			}
			else if(interpolation_method == "1Dy_linear"){
				interpolated_psi_y = linear_interpolation_1Dy_psi_velocity_y(
					position,
					psi_velocity_cell_center_y,
					*this
				);
			}
			else if(interpolation_method == "1Dy_WENO6"){
	//			interpolated_psi_y = WENO6_interpolation_1Dy_psi_velocity_x(
	//				position,
	//				psi_velocity_cell_center_y,
	//				*this
	//			);
			}
			// psi のうち z 成分の計算
			interpolated_psi_z = 0.0;
		}
		else{
			if (interpolation_method == "linear") {
				interpolated_psi_y =  linear_interpolation_2D_psi_velocity_y(
					position,
					psi_velocity_cell_center_y,
					*this
				);
			}
			else if (interpolation_method == "WENO6") {
	//			interpolated_psi_y = WENO6_interpolation_2D_psi_velocity_x(
	//				position,
	//				psi_velocity_cell_center_y,
	//				*this
	//			);
			}
			else if(interpolation_method == "1Dy_linear"){
				interpolated_psi_y = linear_interpolation_1Dy_psi_velocity_y(
					position,
					psi_velocity_cell_center_y,
					*this
				);
			}
			else if(interpolation_method == "1Dy_WENO6"){
	//			interpolated_psi_y = WENO6_interpolation_1Dy_psi_velocity_x(
	//				position,
	//				psi_velocity_cell_center_y,
	//				*this
	//			);
			}

			// psi のうち z 成分の計算
			interpolated_psi_z = 0.0;
		}
		return VEC3_TYPE(interpolated_psi_x, interpolated_psi_y, interpolated_psi_z);
	}

	//position での deisity の値を周辺の cell center での値から interpolation して求める。
	MY_FLOAT_TYPE Grid::calc_substance_density_by_interpolation(
		VEC3_TYPE position,
		const std::string interpolation_method) const {
		// linear interpolation を使う場合
		if (interpolation_method == "linear") {
			return linear_interpolation_2D_cell_center_values(
			    position,
			    substance_density,
			    *this
			);
		}
		// WENO6 interpolation を使う場合
		else if (interpolation_method == "WENO6") {
			return WENO6_interpolation_2D_cell_center_values(
			    position,
			    substance_density,
			    *this
			);
		}
		// 1d linear interpolation を使う場合
		else if (interpolation_method == "1Dy_linear") {
			return linear_interpolation_1Dy_cell_center_values(
			    position,
			    substance_density,
			    *this
			);
		}
		// 1d linear interpolation を使う場合
		else if (interpolation_method == "1Dy_WENO6") {
			return WENO6_interpolation_1Dy_cell_center_values(
			    position,
			    substance_density,
			    *this
			);
		}
	}

	// negative density artifact を解消するために, バックトレースしたセルの面を軸(x軸またはy軸, およびその両方)で切断する関数
	std::vector<cell_face> Grid::split_cell_face_by_axis(const cell_face face, std::string split_mode) const {
		if (split_mode == "x") {
			// 1つ目の vertex が含まれるセルのインデックスを格納する変数(0要素目がx軸, 1要素目がy軸方向のインデックス)
			std::vector<int> included_cell_index_of_face_vertex_1(2);
			// 2つ目の vertex が含まれるセルのインデックスを格納する変数(0要素目がx軸, 1要素目がy軸方向のインデックス)
			std::vector<int> included_cell_index_of_face_vertex_2(2);

			included_cell_index_of_face_vertex_1[0] = floor((face._vertex_list[0]._vertex_pos[0] - min_pos_x) / _cell_length);
			included_cell_index_of_face_vertex_1[1] = floor((face._vertex_list[0]._vertex_pos[1] - min_pos_y) / _cell_length);
			included_cell_index_of_face_vertex_2[0] = floor((face._vertex_list[1]._vertex_pos[0] - min_pos_x) / _cell_length);
			included_cell_index_of_face_vertex_2[1] = floor((face._vertex_list[1]._vertex_pos[1] - min_pos_y) / _cell_length);

			//結果を格納する変数
			std::vector<cell_face> splitted_faces_list;
			//2つのvertex が含まれるセルのインデックスのx成分が異なる場合、
			//それらを結ぶ線はグリッドの面が作る線のうち, y軸に平行な線のいくつかをまたがっている
			if (included_cell_index_of_face_vertex_1[1] != included_cell_index_of_face_vertex_2[1]) {
				//2つの vertex の位置のうちx座標が小さいものと大きいもの
				VEC3_TYPE vertex_pos_xmin, vertex_pos_xmax;
				//2つの vertex を含むセルのインデックスのうちx座標が小さいものと大きいもの
				std::vector<int> vertex_index_xmin(2), vertex_index_xmax(2);
				if (included_cell_index_of_face_vertex_1[1] < included_cell_index_of_face_vertex_2[1]) {
					vertex_pos_xmin = face._vertex_list[0]._vertex_pos;
					vertex_pos_xmax = face._vertex_list[1]._vertex_pos;
					vertex_index_xmin = included_cell_index_of_face_vertex_1;
					vertex_index_xmax = included_cell_index_of_face_vertex_2;
				}
				else {
					vertex_pos_xmin = face._vertex_list[1]._vertex_pos;
					vertex_pos_xmax = face._vertex_list[0]._vertex_pos;
					vertex_index_xmin = included_cell_index_of_face_vertex_2;
					vertex_index_xmax = included_cell_index_of_face_vertex_1;
				}

				//// 2つの vertex を結ぶ線の方程式をつくる
				//// 以下の傾きと切片を用いて x = slope * y + intercept のように表せる
				// 2つの vertex を結ぶ線の傾きを求める
				MY_FLOAT_TYPE slope = (vertex_pos_xmax[0] - vertex_pos_xmin[0]) / (vertex_pos_xmax[1] - vertex_pos_xmin[1]);
				// 2つの vertex を結ぶ線の切片を求める
				MY_FLOAT_TYPE intercept = vertex_pos_xmin[0] - slope * vertex_pos_xmin[1];
				////引数の face を y軸に平行な grid face で切断していく
				//切断後の面を作る頂点を集めたリスト
				std::vector<cell_vertex> vertex_list;
				//引数の面の始点を追加する
				vertex_list.push_back(cell_vertex(vertex_pos_xmin));
				for (int i_index_y = vertex_index_xmin[1] + 1;; ++i_index_y) {
					if (i_index_y > vertex_index_xmax[1]) {
						break;
					}
					//切断する面の x座標と y座標
					MY_FLOAT_TYPE face_pos_y = min_pos_y + i_index_y * _cell_length;
					MY_FLOAT_TYPE face_pos_x = slope * face_pos_y + intercept;
					//2頂点を結んだ直線とy軸に平行な grid face の交点の位置
					VEC3_TYPE intersection_pos(face_pos_x, face_pos_y, 0.0);
					//引数の面の始点を追加する
					vertex_list.push_back(cell_vertex(intersection_pos));

					//// 切断した点(2頂点を結んだ直線とy軸に平行な grid face の交点) が2頂点を結んだ直線の始点と面を作る場合
					//if (i_index_x == vertex_index_xmin[0] + 1) {
					//}
					//// 切断した点(2頂点を結んだ直線とy軸に平行な grid face の交点) が2頂点を結んだ直線の終点と面を作る場合
					//else if (i_index_x == vertex_index_xmax[0]) {
					//}
					//// 切断した点(2頂点を結んだ直線とy軸に平行な grid face の交点) 同士が面を作る場合
					//else {
					//}
				}
				//引数の面の終点を追加する
				vertex_list.push_back(cell_vertex(vertex_pos_xmax));
				//切断後の面を追加していく
				for (int i_vert = 0; i_vert < vertex_list.size() - 1; ++i_vert) {
					cell_face tmp;
					tmp._vertex_list[0] = vertex_list[i_vert];
					tmp._vertex_list[1] = vertex_list[i_vert + 1];
					splitted_faces_list.push_back(tmp);
				}
				//cell face を構成する頂点を 大きい座標 -> 小さい座標 の順番で格納すべき場合
				if (included_cell_index_of_face_vertex_2[1] < included_cell_index_of_face_vertex_1[1]) {
					for (auto itr_face = splitted_faces_list.begin(); itr_face != splitted_faces_list.end(); ++itr_face) {
						//std::cout << "flipped" << std::endl;
						itr_face->invert_vertex();
					}
				}
			}
			//2つのvertex が含まれるセルのインデックスのx成分が同じ場合は、
			//それらを結ぶ線はグリッドの面が作る線のうち, y軸に平行な線をまたがらないので元の面を返すだけでいい
			else {
				splitted_faces_list.push_back(face);
			}
			return splitted_faces_list;
		}
		else if (split_mode == "y") {
			// 1つ目の vertex が含まれるセルのインデックスを格納する変数(0要素目がx軸, 1要素目がy軸方向のインデックス)
			std::vector<int> included_cell_index_of_face_vertex_1(2);
			// 2つ目の vertex が含まれるセルのインデックスを格納する変数(0要素目がx軸, 1要素目がy軸方向のインデックス)
			std::vector<int> included_cell_index_of_face_vertex_2(2);

			included_cell_index_of_face_vertex_1[0] = floor((face._vertex_list[0]._vertex_pos[0] - min_pos_x) / _cell_length);
			included_cell_index_of_face_vertex_1[1] = floor((face._vertex_list[0]._vertex_pos[1] - min_pos_y) / _cell_length);
			included_cell_index_of_face_vertex_2[0] = floor((face._vertex_list[1]._vertex_pos[0] - min_pos_x) / _cell_length);
			included_cell_index_of_face_vertex_2[1] = floor((face._vertex_list[1]._vertex_pos[1] - min_pos_y) / _cell_length);

			//結果を格納する変数
			std::vector<cell_face> splitted_faces_list;
			//2つのvertex が含まれるセルのインデックスのx成分が異なる場合、
			//それらを結ぶ線はグリッドの面が作る線のうち, y軸に平行な線のいくつかをまたがっている
			if (included_cell_index_of_face_vertex_1[0] != included_cell_index_of_face_vertex_2[0]) {
				//2つの vertex の位置のうちx座標が小さいものと大きいもの
				VEC3_TYPE vertex_pos_xmin, vertex_pos_xmax;
				//2つの vertex を含むセルのインデックスのうちx座標が小さいものと大きいもの
				std::vector<int> vertex_index_xmin(2), vertex_index_xmax(2);
				if (included_cell_index_of_face_vertex_1[0] < included_cell_index_of_face_vertex_2[0]) {
					vertex_pos_xmin = face._vertex_list[0]._vertex_pos;
					vertex_pos_xmax = face._vertex_list[1]._vertex_pos;
					vertex_index_xmin = included_cell_index_of_face_vertex_1;
					vertex_index_xmax = included_cell_index_of_face_vertex_2;
				}
				else {
					vertex_pos_xmin = face._vertex_list[1]._vertex_pos;
					vertex_pos_xmax = face._vertex_list[0]._vertex_pos;
					vertex_index_xmin = included_cell_index_of_face_vertex_2;
					vertex_index_xmax = included_cell_index_of_face_vertex_1;
				}

				//// 2つの vertex を結ぶ線の方程式をつくる
				//// 以下の傾きと切片を用いて y = slope * x + intercept のように表せる
				// 2つの vertex を結ぶ線の傾きを求める
				MY_FLOAT_TYPE slope = (vertex_pos_xmax[1] - vertex_pos_xmin[1]) / (vertex_pos_xmax[0] - vertex_pos_xmin[0]);
				// 2つの vertex を結ぶ線の切片を求める
				MY_FLOAT_TYPE intercept = vertex_pos_xmin[1] - slope * vertex_pos_xmin[0];
				////引数の face を y軸に平行な grid face で切断していく
				//切断後の面を作る頂点を集めたリスト
				std::vector<cell_vertex> vertex_list;
				//引数の面の始点を追加する
				vertex_list.push_back(cell_vertex(vertex_pos_xmin));
				for (int i_index_x = vertex_index_xmin[0] + 1;; ++i_index_x) {
					if (i_index_x > vertex_index_xmax[0]) {
						break;
					}
					//切断する面の x座標と y座標
					MY_FLOAT_TYPE face_pos_x = min_pos_x + i_index_x * _cell_length;
					MY_FLOAT_TYPE face_pos_y = slope * face_pos_x + intercept;
					//2頂点を結んだ直線とy軸に平行な grid face の交点の位置
					VEC3_TYPE intersection_pos(face_pos_x, face_pos_y, 0.0);
					//引数の面の始点を追加する
					vertex_list.push_back(cell_vertex(intersection_pos));

					//// 切断した点(2頂点を結んだ直線とy軸に平行な grid face の交点) が2頂点を結んだ直線の始点と面を作る場合
					//if (i_index_x == vertex_index_xmin[0] + 1) {
					//}
					//// 切断した点(2頂点を結んだ直線とy軸に平行な grid face の交点) が2頂点を結んだ直線の終点と面を作る場合
					//else if (i_index_x == vertex_index_xmax[0]) {
					//}
					//// 切断した点(2頂点を結んだ直線とy軸に平行な grid face の交点) 同士が面を作る場合
					//else {
					//}
				}
				//引数の面の終点を追加する
				vertex_list.push_back(cell_vertex(vertex_pos_xmax));
				//切断後の面を追加していく
				for (int i_vert = 0; i_vert < vertex_list.size() - 1; ++i_vert) {
					cell_face tmp;
					tmp._vertex_list[0] = vertex_list[i_vert];
					tmp._vertex_list[1] = vertex_list[i_vert + 1];
					splitted_faces_list.push_back(tmp);
				}
				//cell face を構成する頂点を 大きい座標 -> 小さい座標 の順番で格納すべき場合
				if (included_cell_index_of_face_vertex_2[0] < included_cell_index_of_face_vertex_1[0]) {
					for (auto itr_face = splitted_faces_list.begin(); itr_face != splitted_faces_list.end(); ++itr_face) {
						//std::cout << "flipped" << std::endl;
						itr_face->invert_vertex();
					}
				}
			}
			//2つのvertex が含まれるセルのインデックスのx成分が同じ場合は、
			//それらを結ぶ線はグリッドの面が作る線のうち, y軸に平行な線をまたがらないので元の面を返すだけでいい
			else {
				splitted_faces_list.push_back(face);
			}
			return splitted_faces_list;
		}

	}

	// negative density artifact を解消するために, バックトレースしたセルの面を軸(x軸またはy軸, およびその両方)で切断する関数
	// psi が x方向に(psi の積分の方向をxにとるならy方向に)区分線型関数であることを仮定したときの分割の方法
	// split_mode == "x" は未実装
	std::vector<cell_face> Grid::split_cell_face_by_cell_center(const cell_face face, std::string split_mode) const {
		if (split_mode == "x") {
			// 1つ目の vertex が含まれるセルのインデックスを格納する変数(0要素目がx軸, 1要素目がy軸方向のインデックス)
			std::vector<int> included_cell_index_of_face_vertex_1(2);
			// 2つ目の vertex が含まれるセルのインデックスを格納する変数(0要素目がx軸, 1要素目がy軸方向のインデックス)
			std::vector<int> included_cell_index_of_face_vertex_2(2);

			included_cell_index_of_face_vertex_1[0] = floor((face._vertex_list[0]._vertex_pos[0] - min_pos_x) / _cell_length);
			included_cell_index_of_face_vertex_1[1] = floor((face._vertex_list[0]._vertex_pos[1] - min_pos_y) / _cell_length - 0.5);
			included_cell_index_of_face_vertex_2[0] = floor((face._vertex_list[1]._vertex_pos[0] - min_pos_x) / _cell_length);
			included_cell_index_of_face_vertex_2[1] = floor((face._vertex_list[1]._vertex_pos[1] - min_pos_y) / _cell_length - 0.5);

			//結果を格納する変数
			std::vector<cell_face> splitted_faces_list;
			//2つのvertex が含まれるセルのインデックスのx成分が異なる場合、
			//それらを結ぶ線はグリッドの面が作る線のうち, y軸に平行な線のいくつかをまたがっている
			if (included_cell_index_of_face_vertex_1[1] != included_cell_index_of_face_vertex_2[1]) {
				//2つの vertex の位置のうちx座標が小さいものと大きいもの
				VEC3_TYPE vertex_pos_xmin, vertex_pos_xmax;
				//2つの vertex を含むセルのインデックスのうちx座標が小さいものと大きいもの
				std::vector<int> vertex_index_xmin(2), vertex_index_xmax(2);
				if (included_cell_index_of_face_vertex_1[1] < included_cell_index_of_face_vertex_2[1]) {
					vertex_pos_xmin = face._vertex_list[0]._vertex_pos;
					vertex_pos_xmax = face._vertex_list[1]._vertex_pos;
					vertex_index_xmin = included_cell_index_of_face_vertex_1;
					vertex_index_xmax = included_cell_index_of_face_vertex_2;
				}
				else {
					vertex_pos_xmin = face._vertex_list[1]._vertex_pos;
					vertex_pos_xmax = face._vertex_list[0]._vertex_pos;
					vertex_index_xmin = included_cell_index_of_face_vertex_2;
					vertex_index_xmax = included_cell_index_of_face_vertex_1;
				}

				//// 2つの vertex を結ぶ線の方程式をつくる
				//// 以下の傾きと切片を用いて x = slope * y + intercept のように表せる
				// 2つの vertex を結ぶ線の傾きを求める
				MY_FLOAT_TYPE slope = (vertex_pos_xmax[0] - vertex_pos_xmin[0]) / (vertex_pos_xmax[1] - vertex_pos_xmin[1]);
				// 2つの vertex を結ぶ線の切片を求める
				MY_FLOAT_TYPE intercept = vertex_pos_xmin[0] - slope * vertex_pos_xmin[1];
				////引数の face を y軸に平行な grid face で切断していく
				//切断後の面を作る頂点を集めたリスト
				std::vector<cell_vertex> vertex_list;
				//引数の面の始点を追加する
				vertex_list.push_back(cell_vertex(vertex_pos_xmin));
				for (int i_index_y = vertex_index_xmin[1] + 1;; ++i_index_y) {
					if (i_index_y > vertex_index_xmax[1]) {
						break;
					}
					//切断する面の x座標と y座標
					MY_FLOAT_TYPE face_pos_y = min_pos_y + (i_index_y + 0.5) * _cell_length;
					MY_FLOAT_TYPE face_pos_x = slope * face_pos_y + intercept;
					//2頂点を結んだ直線とy軸に平行な grid face の交点の位置
					VEC3_TYPE intersection_pos(face_pos_x, face_pos_y, 0.0);
					//引数の面の始点を追加する
					vertex_list.push_back(cell_vertex(intersection_pos));

					//// 切断した点(2頂点を結んだ直線とy軸に平行な grid face の交点) が2頂点を結んだ直線の始点と面を作る場合
					//if (i_index_x == vertex_index_xmin[0] + 1) {
					//}
					//// 切断した点(2頂点を結んだ直線とy軸に平行な grid face の交点) が2頂点を結んだ直線の終点と面を作る場合
					//else if (i_index_x == vertex_index_xmax[0]) {
					//}
					//// 切断した点(2頂点を結んだ直線とy軸に平行な grid face の交点) 同士が面を作る場合
					//else {
					//}
				}
				//引数の面の終点を追加する
				vertex_list.push_back(cell_vertex(vertex_pos_xmax));
				//切断後の面を追加していく
				for (int i_vert = 0; i_vert < vertex_list.size() - 1; ++i_vert) {
					cell_face tmp;
					tmp._vertex_list[0] = vertex_list[i_vert];
					tmp._vertex_list[1] = vertex_list[i_vert + 1];
					splitted_faces_list.push_back(tmp);
				}
				//cell face を構成する頂点を 大きい座標 -> 小さい座標 の順番で格納すべき場合
				if (included_cell_index_of_face_vertex_2[1] < included_cell_index_of_face_vertex_1[1]) {
					for (auto itr_face = splitted_faces_list.begin(); itr_face != splitted_faces_list.end(); ++itr_face) {
						//std::cout << "flipped" << std::endl;
						itr_face->invert_vertex();
					}
				}
			}
			//2つのvertex が含まれるセルのインデックスのx成分が同じ場合は、
			//それらを結ぶ線はグリッドの面が作る線のうち, y軸に平行な線をまたがらないので元の面を返すだけでいい
			else {
				splitted_faces_list.push_back(face);
			}
			return splitted_faces_list;
		}
		else if (split_mode == "y") {
			// 1つ目の vertex が含まれるセルのインデックスを格納する変数(0要素目がx軸, 1要素目がy軸方向のインデックス)
			std::vector<int> included_cell_index_of_face_vertex_1(2);
			// 2つ目の vertex が含まれるセルのインデックスを格納する変数(0要素目がx軸, 1要素目がy軸方向のインデックス)
			std::vector<int> included_cell_index_of_face_vertex_2(2);
//			std::cout<<"face._vertex_list[0]._vertex_pos: "<<face._vertex_list[0]._vertex_pos<<std::endl;
//			std::cout<<"face._vertex_list[1]._vertex_pos: "<<face._vertex_list[1]._vertex_pos<<std::endl;
			included_cell_index_of_face_vertex_1[0] = floor((face._vertex_list[0]._vertex_pos[0] - min_pos_x) / _cell_length - 0.5);
			included_cell_index_of_face_vertex_1[1] = floor((face._vertex_list[0]._vertex_pos[1] - min_pos_y) / _cell_length);
			included_cell_index_of_face_vertex_2[0] = floor((face._vertex_list[1]._vertex_pos[0] - min_pos_x) / _cell_length - 0.5);
			included_cell_index_of_face_vertex_2[1] = floor((face._vertex_list[1]._vertex_pos[1] - min_pos_y) / _cell_length);

			//結果を格納する変数
			std::vector<cell_face> splitted_faces_list;
			//2つのvertex が含まれるセルのインデックスのx成分が異なる場合、
			//それらを結ぶ線はグリッドの面が作る線のうち, y軸に平行な線のいくつかをまたがっている
			if (included_cell_index_of_face_vertex_1[0] != included_cell_index_of_face_vertex_2[0]) {
				//2つの vertex の位置のうちx座標が小さいものと大きいもの
				VEC3_TYPE vertex_pos_xmin, vertex_pos_xmax;
				//2つの vertex を含むセルのインデックスのうちx座標が小さいものと大きいもの
				std::vector<int> vertex_index_xmin(2), vertex_index_xmax(2);
				if (included_cell_index_of_face_vertex_1[0] < included_cell_index_of_face_vertex_2[0]) {
					vertex_pos_xmin = face._vertex_list[0]._vertex_pos;
					vertex_pos_xmax = face._vertex_list[1]._vertex_pos;
					vertex_index_xmin = included_cell_index_of_face_vertex_1;
					vertex_index_xmax = included_cell_index_of_face_vertex_2;
				}
				else {
					vertex_pos_xmin = face._vertex_list[1]._vertex_pos;
					vertex_pos_xmax = face._vertex_list[0]._vertex_pos;
					vertex_index_xmin = included_cell_index_of_face_vertex_2;
					vertex_index_xmax = included_cell_index_of_face_vertex_1;
				}

				//// 2つの vertex を結ぶ線の方程式をつくる
				//// 以下の傾きと切片を用いて y = slope * x + intercept のように表せる
				// 2つの vertex を結ぶ線の傾きを求める
				MY_FLOAT_TYPE slope = (vertex_pos_xmax[1] - vertex_pos_xmin[1]) / (vertex_pos_xmax[0] - vertex_pos_xmin[0]);
				// 2つの vertex を結ぶ線の切片を求める
				MY_FLOAT_TYPE intercept = vertex_pos_xmin[1] - slope * vertex_pos_xmin[0];
				////引数の face を y軸に平行な grid face で切断していく
				//切断後の面を作る頂点を集めたリスト
				std::vector<cell_vertex> vertex_list;
				//引数の面の始点を追加する
				vertex_list.push_back(cell_vertex(vertex_pos_xmin));
				for (int i_index_x = vertex_index_xmin[0] + 1;; ++i_index_x) {
					if (i_index_x > vertex_index_xmax[0]) {
						break;
					}
					//切断する面の x座標と y座標
					MY_FLOAT_TYPE face_pos_x = min_pos_x + (i_index_x + 0.5) * _cell_length;
					MY_FLOAT_TYPE face_pos_y = slope * face_pos_x + intercept;
					//2頂点を結んだ直線とy軸に平行な grid face の交点の位置
					VEC3_TYPE intersection_pos(face_pos_x, face_pos_y, 0.0);
					//引数の面の始点を追加する
					vertex_list.push_back(cell_vertex(intersection_pos));
				}
				//引数の面の終点を追加する
				vertex_list.push_back(cell_vertex(vertex_pos_xmax));
				//切断後の面を追加していく
				for (int i_vert = 0; i_vert < vertex_list.size() - 1; ++i_vert) {
					cell_face tmp;
					tmp._vertex_list[0] = vertex_list[i_vert];
					tmp._vertex_list[1] = vertex_list[i_vert + 1];
					splitted_faces_list.push_back(tmp);
				}
				//cell face を構成する頂点を 大きい座標 -> 小さい座標 の順番で格納すべき場合
				if (included_cell_index_of_face_vertex_2[0] < included_cell_index_of_face_vertex_1[0]) {
					for (auto itr_face = splitted_faces_list.begin(); itr_face != splitted_faces_list.end(); ++itr_face) {
						itr_face->invert_vertex();
					}
				}
			}
			//2つのvertex が含まれるセルのインデックスのx成分が同じ場合は、
			//それらを結ぶ線はグリッドの面が作る線のうち, y軸に平行な線をまたがらないので元の面を返すだけでいい
			else {
				splitted_faces_list.push_back(face);
			}
			return splitted_faces_list;
		}

	}


	// バックトレースしたセルの面を等間隔で切断する関数
	std::vector<cell_face> Grid::split_cell_face_uniformly(cell_face face, const int split_face_num) const {
		VEC3_TYPE vertex_pos_begin = face._vertex_list[0]._vertex_pos;
		VEC3_TYPE vertex_pos_end = face._vertex_list[1]._vertex_pos;
		VEC3_TYPE edge_direction = (vertex_pos_end - vertex_pos_begin).normalized();
		//切断前のfaceの面積
		const MY_FLOAT_TYPE edge_length = (vertex_pos_end - vertex_pos_begin).norm();
		//切断後の1つのfaceの面積
		const MY_FLOAT_TYPE splitted_edge_length = edge_length / split_face_num;

		//結果を格納する変数
		std::vector<cell_face> splitted_faces_list;
		for (int i_face = 0; i_face < split_face_num; ++i_face) {
			VEC3_TYPE splitted_vertex_pos_begin = vertex_pos_begin + i_face * splitted_edge_length * edge_direction;
			VEC3_TYPE splitted_vertex_pos_end = vertex_pos_begin + (i_face + 1) * splitted_edge_length * edge_direction;
			//追加する面
			cell_face tmp;
			tmp._vertex_list[0] = cell_vertex(splitted_vertex_pos_begin);
			tmp._vertex_list[1] = cell_vertex(splitted_vertex_pos_end);
			splitted_faces_list.push_back(tmp);
		}
		return splitted_faces_list;
	}

	// バックトレースしたセルの面を軸(x軸またはy軸, およびその両方)に平行な平面で切断する関数
	// split_mode == "x" の時は x軸に平行な, split_mode == "y" の時は y軸に平行な平面で切断する.
	// split_positions は切断する平面の位置を表す。例えばsplit_mode == "y"でsplit_positions={1.4, 1.8, 2.5}の時はx=1.4, 1.8, 2.5 の三つの平面で切断する。
	// split_mode == "x" は未実装
	std::vector<cell_face> Grid::split_cell_face_by_axis_aligned_plane(const cell_face face, std::vector<MY_FLOAT_TYPE> split_positions, std::string split_mode) const {
		if (split_mode == "x") {
			// 1つ目の vertex が含まれるセルのインデックスを格納する変数(0要素目がx軸, 1要素目がy軸方向のインデックス)
			std::vector<int> included_cell_index_of_face_vertex_1(2);
			// 2つ目の vertex が含まれるセルのインデックスを格納する変数(0要素目がx軸, 1要素目がy軸方向のインデックス)
			std::vector<int> included_cell_index_of_face_vertex_2(2);

			included_cell_index_of_face_vertex_1[0] = floor((face._vertex_list[0]._vertex_pos[0] - min_pos_x) / _cell_length);
			included_cell_index_of_face_vertex_1[1] = floor((face._vertex_list[0]._vertex_pos[1] - min_pos_y) / _cell_length);
			included_cell_index_of_face_vertex_2[0] = floor((face._vertex_list[1]._vertex_pos[0] - min_pos_x) / _cell_length);
			included_cell_index_of_face_vertex_2[1] = floor((face._vertex_list[1]._vertex_pos[1] - min_pos_y) / _cell_length);

			//結果を格納する変数
			std::vector<cell_face> splitted_faces_list;
			//2つのvertex が含まれるセルのインデックスのx成分が異なる場合、
			//それらを結ぶ線はグリッドの面が作る線のうち, y軸に平行な線のいくつかをまたがっている
			if (included_cell_index_of_face_vertex_1[1] != included_cell_index_of_face_vertex_2[1]) {
				//2つの vertex の位置のうちx座標が小さいものと大きいもの
				VEC3_TYPE vertex_pos_xmin, vertex_pos_xmax;
				//2つの vertex を含むセルのインデックスのうちx座標が小さいものと大きいもの
				std::vector<int> vertex_index_xmin(2), vertex_index_xmax(2);
				if (included_cell_index_of_face_vertex_1[1] < included_cell_index_of_face_vertex_2[1]) {
					vertex_pos_xmin = face._vertex_list[0]._vertex_pos;
					vertex_pos_xmax = face._vertex_list[1]._vertex_pos;
					vertex_index_xmin = included_cell_index_of_face_vertex_1;
					vertex_index_xmax = included_cell_index_of_face_vertex_2;
				}
				else {
					vertex_pos_xmin = face._vertex_list[1]._vertex_pos;
					vertex_pos_xmax = face._vertex_list[0]._vertex_pos;
					vertex_index_xmin = included_cell_index_of_face_vertex_2;
					vertex_index_xmax = included_cell_index_of_face_vertex_1;
				}

				//// 2つの vertex を結ぶ線の方程式をつくる
				//// 以下の傾きと切片を用いて x = slope * y + intercept のように表せる
				// 2つの vertex を結ぶ線の傾きを求める
				MY_FLOAT_TYPE slope = (vertex_pos_xmax[0] - vertex_pos_xmin[0]) / (vertex_pos_xmax[1] - vertex_pos_xmin[1]);
				// 2つの vertex を結ぶ線の切片を求める
				MY_FLOAT_TYPE intercept = vertex_pos_xmin[0] - slope * vertex_pos_xmin[1];
				////引数の face を y軸に平行な grid face で切断していく
				//切断後の面を作る頂点を集めたリスト
				std::vector<cell_vertex> vertex_list;
				//引数の面の始点を追加する
				vertex_list.push_back(cell_vertex(vertex_pos_xmin));
				for (int i_index_y = vertex_index_xmin[1] + 1;; ++i_index_y) {
					if (i_index_y > vertex_index_xmax[1]) {
						break;
					}
					//切断する面の x座標と y座標
					MY_FLOAT_TYPE face_pos_y = min_pos_y + i_index_y * _cell_length;
					MY_FLOAT_TYPE face_pos_x = slope * face_pos_y + intercept;
					//2頂点を結んだ直線とy軸に平行な grid face の交点の位置
					VEC3_TYPE intersection_pos(face_pos_x, face_pos_y, 0.0);
					//引数の面の始点を追加する
					vertex_list.push_back(cell_vertex(intersection_pos));
				}
				//引数の面の終点を追加する
				vertex_list.push_back(cell_vertex(vertex_pos_xmax));
				//切断後の面を追加していく
				for (int i_vert = 0; i_vert < vertex_list.size() - 1; ++i_vert) {
					cell_face tmp;
					tmp._vertex_list[0] = vertex_list[i_vert];
					tmp._vertex_list[1] = vertex_list[i_vert + 1];
					splitted_faces_list.push_back(tmp);
				}
				//cell face を構成する頂点を 大きい座標 -> 小さい座標 の順番で格納すべき場合
				if (included_cell_index_of_face_vertex_2[1] < included_cell_index_of_face_vertex_1[1]) {
					for (auto itr_face = splitted_faces_list.begin(); itr_face != splitted_faces_list.end(); ++itr_face) {
						//std::cout << "flipped" << std::endl;
						itr_face->invert_vertex();
					}
				}
			}
			//2つのvertex が含まれるセルのインデックスのx成分が同じ場合は、
			//それらを結ぶ線はグリッドの面が作る線のうち, y軸に平行な線をまたがらないので元の面を返すだけでいい
			else {
				splitted_faces_list.push_back(face);
			}
			return splitted_faces_list;
		}
		else if (split_mode == "y") {
			//切断する面の位置を昇順にソート
			std::sort(split_positions.begin(), split_positions.end());

			//2つの vertex の位置のうちx座標が小さいものと大きいもの
			VEC3_TYPE vertex_pos_xmin, vertex_pos_xmax;
			if (face._vertex_list[0]._vertex_pos[0] < face._vertex_list[1]._vertex_pos[0]) {
				vertex_pos_xmin = face._vertex_list[0]._vertex_pos;
				vertex_pos_xmax = face._vertex_list[1]._vertex_pos;
			}
			else {
				vertex_pos_xmin = face._vertex_list[1]._vertex_pos;
				vertex_pos_xmax = face._vertex_list[0]._vertex_pos;
			}
			//// 2つの vertex を結ぶ線の方程式をつくる
			//// 以下の傾きと切片を用いて y = slope * x + intercept のように表せる
			// 2つの vertex を結ぶ線の傾きを求める
			MY_FLOAT_TYPE slope = (vertex_pos_xmax[1] - vertex_pos_xmin[1]) / (vertex_pos_xmax[0] - vertex_pos_xmin[0]);
			// 2つの vertex を結ぶ線の切片を求める
			MY_FLOAT_TYPE intercept = vertex_pos_xmin[1] - slope * vertex_pos_xmin[0];

			//切断後の面を作る頂点を集めたリスト
			std::vector<cell_vertex> vertex_list;
			//引数の面の始点を追加する
			vertex_list.push_back(cell_vertex(vertex_pos_xmin));
			//分割する平面の数
			const int split_plane_num = split_positions.size();
			//各平面を走査するループ
			for (int i_plane = 0; i_plane < split_plane_num; ++i_plane) {
				// 切断する面の位置が二つの頂点の間の位置にあるとき、考えている面は切断される
				if (vertex_pos_xmin[0] < split_positions[i_plane] && split_positions[i_plane] < vertex_pos_xmax[0]) {
					//切断する面の x座標と y座標
					MY_FLOAT_TYPE face_pos_x = split_positions[i_plane];
					MY_FLOAT_TYPE face_pos_y = slope * face_pos_x + intercept;
					//2頂点を結んだ直線とy軸に平行な grid face の交点の位置
					VEC3_TYPE intersection_pos(face_pos_x, face_pos_y, 0.0);
					//切断する面と引数の面の交差点を追加する
					vertex_list.push_back(cell_vertex(intersection_pos));
				}
			}
			//引数の面の終点を追加する
			vertex_list.push_back(cell_vertex(vertex_pos_xmax));

			//結果を格納する変数
			std::vector<cell_face> splitted_faces_list;
			//切断後の面を追加していく
			for (int i_vert = 0; i_vert < vertex_list.size() - 1; ++i_vert) {
				cell_face tmp;
				tmp._vertex_list[0] = vertex_list[i_vert];
				tmp._vertex_list[1] = vertex_list[i_vert + 1];
				splitted_faces_list.push_back(tmp);
			}
			//cell face を構成する頂点を 大きい座標 -> 小さい座標 の順番で格納すべき場合
			if (face._vertex_list[1]._vertex_pos[0] < face._vertex_list[0]._vertex_pos[0]) {
				for (auto itr_face = splitted_faces_list.begin(); itr_face != splitted_faces_list.end(); ++itr_face) {
					//std::cout << "flipped" << std::endl;
					itr_face->invert_vertex();
				}
			}
			return splitted_faces_list;
		}

	}

	//position での速度をinterpolationによって計算
	VEC3_TYPE Grid::calc_interpolated_velocity(
		VEC3_TYPE position,
		const std::string interpolation_method) const {
		VEC3_TYPE interpolated_velocity;
		if (interpolation_method == "linear") {
			interpolated_velocity[0] = linear_interpolation_2D_cell_face_x_values(
			    position,
			    velocity_in_voxel_face_x,
			    *this
			);
			interpolated_velocity[1] = linear_interpolation_2D_cell_face_y_values(
				   position,
				   velocity_in_voxel_face_y,
				   *this
			   );
			interpolated_velocity[2] = 0.0;
		}
		else if (interpolation_method == "WENO6") {
			interpolated_velocity[0] = WENO6_interpolation_2D_cell_face_x_values(
			    position,
			    velocity_in_voxel_face_x,
			    *this
			);
			interpolated_velocity[1] = WENO6_interpolation_2D_cell_face_y_values(
			    position,
			    velocity_in_voxel_face_y,
			    *this
			);
			interpolated_velocity[2] = 0.0;
		}
		else if (interpolation_method == "1Dy_linear") {
			interpolated_velocity[0] = linear_interpolation_1Dy_cell_face_x_values(
			    position,
			    velocity_in_voxel_face_x,
			    *this
			);
			interpolated_velocity[1] = linear_interpolation_1Dy_cell_face_y_values(
				   position,
				   velocity_in_voxel_face_y,
				   *this
			   );
			interpolated_velocity[2] = 0.0;
		}
		else if (interpolation_method == "1Dy_WENO6") {
			interpolated_velocity[0] = WENO6_interpolation_1Dy_cell_face_x_values(
			    position,
			    velocity_in_voxel_face_x,
			    *this
			);
			interpolated_velocity[1] = WENO6_interpolation_1Dy_cell_face_y_values(
				   position,
				   velocity_in_voxel_face_y,
				   *this
			   );
			interpolated_velocity[2] = 0.0;
		}
		return interpolated_velocity;
	}

	//position での速度をinterpolationによって計算
	MY_FLOAT_TYPE Grid::calc_interpolated_velocity_x(
		VEC3_TYPE position,
		const std::string interpolation_method
	) const {
		MY_FLOAT_TYPE interpolated_velocity_x;
		if (interpolation_method == "linear") {
			interpolated_velocity_x = linear_interpolation_2D_cell_face_x_values(
			    position,
			    velocity_in_voxel_face_x,
			    *this
			);
		}
		else if (interpolation_method == "WENO6") {
			interpolated_velocity_x = WENO6_interpolation_2D_cell_face_x_values(
			    position,
			    velocity_in_voxel_face_x,
			    *this
			);
		}
		else if (interpolation_method == "1Dy_linear") {
			interpolated_velocity_x = linear_interpolation_1Dy_cell_face_x_values(
			    position,
			    velocity_in_voxel_face_x,
			    *this
			);
		}
		else if (interpolation_method == "1Dy_WENO6") {
			interpolated_velocity_x = WENO6_interpolation_1Dy_cell_face_x_values(
			    position,
			    velocity_in_voxel_face_x,
			    *this
			);
		}
		return interpolated_velocity_x;
	}
	//position での速度をinterpolationによって計算
	MY_FLOAT_TYPE Grid::calc_interpolated_velocity_y(
		VEC3_TYPE position,
		const std::string interpolation_method
	) const {
		MY_FLOAT_TYPE interpolated_velocity_y;
		if (interpolation_method == "linear") {
			interpolated_velocity_y = linear_interpolation_2D_cell_face_y_values(
			    position,
			    velocity_in_voxel_face_y,
			    *this
			);
		}
		else if (interpolation_method == "WENO6") {
			interpolated_velocity_y = WENO6_interpolation_2D_cell_face_y_values(
			    position,
			    velocity_in_voxel_face_y,
			    *this
			);
		}
		else if (interpolation_method == "1Dy_linear") {
			interpolated_velocity_y = linear_interpolation_1Dy_cell_face_y_values(
			    position,
			    velocity_in_voxel_face_y,
			    *this
			);
		}
		else if (interpolation_method == "1Dy_WENO6") {
			interpolated_velocity_y = WENO6_interpolation_1Dy_cell_face_y_values(
			    position,
			    velocity_in_voxel_face_y,
			    *this
			);
		}
		return interpolated_velocity_y;
	}


	// cell center で定義される量cell_center_valuesのpositionでの補間を tri-linear interpolation で計算する
    MY_FLOAT_TYPE Grid::interpolate_cell_center_defined_values_trilinear(
        const VEC3_TYPE &position,
        const std::vector<MY_FLOAT_TYPE>& cell_center_values) const {
        int advected_index_x = floor((position[0] - 0.5 * _cell_length) / _cell_length);
        int advected_index_y = floor((position[1] - 0.5 * _cell_length) / _cell_length);
        /////linear interpolation の処理
        //y方向のinterpolationの結果を格納するための変数
        MY_FLOAT_TYPE y_interpolation_values_of_substance_density[2];
        //x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
        for (int ix = 0; ix < 2; ix++) {
            //補間に使う離散値をセット
            MY_FLOAT_TYPE psi_val[2];
            for (int iy = 0; iy < 2; ++iy) {
                int index_x = advected_index_x + ix;
                int index_y = advected_index_y + iy;
                // グリッドの外側を参照しようとしたときの処理(x方向)
                if (index_x < 0) {
                    index_x = 0;
                }
                if (index_x > Grid_num_x - 1) {
                    int exceed = index_x - (Grid_num_x - 1);
                    //境界の値が外側までずっと続く場合
                    index_x = Grid_num_x - 1;
                	//境界を境に鏡のように値が反射する場合
                    //index_x = Grid_num_x - exceed + 1;
                }
                // グリッドの外側を参照しようとしたときの処理(y方向)
                if (index_y < 0) {
                    index_y = 0;
                }
                if (index_y > Grid_num_y - 1) {
                    int exceed = index_y - (Grid_num_y - 1);
                    //境界の値が外側までずっと続く場合
                    index_y = Grid_num_y - 1;
                    //境界を境に鏡のように値が反射する場合
                    //index_y = Grid_num_y - 1 - exceed + 1;
                }
                psi_val[iy] = cell_center_values[get_voxel_center_index(index_x, index_y, Grid_num_x, Grid_num_y)];
            }
            //interpolation の係数
            MY_FLOAT_TYPE b0, b1;
            b0 = (position[1] - 0.5 * _cell_length) / _cell_length - advected_index_y;
            b1 = 1.0 - b0;
            y_interpolation_values_of_substance_density[ix] = b1 * psi_val[0] + b0 * psi_val[1];
        }

        ////ここからx方向の補間
        //interpolation の係数
        MY_FLOAT_TYPE a0, a1;
        a0 = (position[0] - 0.5 * _cell_length) / _cell_length - advected_index_x;
        a1 = 1.0 - a0;
        //結果
        return a1 * y_interpolation_values_of_substance_density[0] + a0 * y_interpolation_values_of_substance_density[1];
    }

}// namespace smoke_simulation
