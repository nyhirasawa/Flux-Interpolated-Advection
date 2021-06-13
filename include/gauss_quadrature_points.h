#pragma once

#include <vector>
#include <thread>
#include <Eigen/Dense>

#include "cell_face.h"
#include "define_float_type.h"

namespace smoke_simulation {
	class quadrature_point_vector{
	public:
		std::vector<std::vector<VEC3_TYPE>> _quadtarure_position;
//		std::vector<std::vector<VEC3_TYPE>> _normal;
//		std::vector<std::vector<MY_FLOAT_TYPE>> _quadtarure_weight;
		//求積点の重み * 面積 * 法線
		std::vector<std::vector<VEC3_TYPE>> _weighted_area_normal;
		std::vector<std::vector<Eigen::Vector3i>> _included_cell_index;

		quadrature_point_vector();
		//quadrature_point(VEC3_TYPE position, VEC3_TYPE normal, MY_FLOAT_TYPE weight, Eigen::Vector3i included_cell_index);
	};

	class gauss_quadrature_points_1D {
		// 1次元ガウス求積の評価点の位置と重みを保持しておくための定数
		// 定数名 quadtarure_positions から始まる定数は1次元ガウス求積の評価点の位置を保持する.
		// 評価点の位置は1次元空間中の位置なので, VEC3_TYPE型の第0成分のみが非0で、第1, 2成分は0になっている.
		// 定数名 quadtarure_weights から始まる定数は1次元ガウス求積の評価点の重みを保持する.
		// 1次元ガウス求積の評価点の位置は区間[-1, 1]の積分を評価するために必要な評価点の位置を保持しておく.
		// 1次元ガウス求積の評価点の重みは和が1になるように正規化したうえで保持しておく.
		static const std::vector<VEC3_TYPE> quadtarure_positions_1D_1point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_1D_3point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_1D_5point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_1D_7point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_1D_9point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_1D_11point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_1D_13point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_1D_25point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_1D_49point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_1D_1point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_1D_3point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_1D_5point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_1D_7point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_1D_9point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_1D_11point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_1D_13point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_1D_25point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_1D_49point;

	public:
		gauss_quadrature_points_1D();
		static const std::vector<VEC3_TYPE> get_quadtarure_positions_1D(const int num_quadrature_points);
		static const std::vector<MY_FLOAT_TYPE> get_quadtarure_weights_1D(const int num_quadrature_points);
		static const std::vector<VEC3_TYPE> calc_quadtarure_positions_on_cell_face(const cell_face face, const int num_quadrature_points);
		static const std::vector<VEC3_TYPE> get_quadtarure_positions_on_1D_interval(const MY_FLOAT_TYPE begin_of_interval, const MY_FLOAT_TYPE end_of_interval, const int num_quadrature_points);
	};

	class gauss_quadrature_points_2D {
	public:
		// 2次元ガウス求積の評価点の位置と重みを保持しておくための定数
		// 定数名 quadtarure_positions から始まる定数は2次元ガウス求積の評価点の位置を保持する.
		// 評価点の位置は2次元空間中の位置なので, VEC3_TYPE型の第0, 1成分のみが非0で、第2成分は0になっている.
		// 定数名 quadtarure_weights から始まる定数は2次元ガウス求積の評価点の重みを保持する.
		// 2次元ガウス求積の評価点の位置は区間[0, 1] * [0, 1]の積分を評価するために必要な評価点の位置を保持しておく.
		// 2次元ガウス求積の評価点の重みは和が1になるように正規化したうえで保持しておく.
		static const std::vector<VEC3_TYPE> quadtarure_positions_2D_1point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_2D_4point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_2D_9point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_2D_16point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_2D_25point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_2D_36point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_2D_49point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_2D_64point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_2D_81point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_2D_100point;
		static const std::vector<VEC3_TYPE> quadtarure_positions_2D_121point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_2D_1point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_2D_4point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_2D_9point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_2D_16point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_2D_25point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_2D_36point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_2D_49point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_2D_64point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_2D_81point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_2D_100point;
		static const std::vector<MY_FLOAT_TYPE> quadtarure_weights_2D_121point;
	};

}// namespace smoke_simulation
