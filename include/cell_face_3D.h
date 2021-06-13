#pragma once

#include <vector>
#include <Eigen/Dense>
#include "cell_vertex_3D.h"

#include "define_float_type.h"

class cell_face_3D {
public:
	std::vector<cell_vertex_3D> _vertex_list;
//	bool is_boundary_face;

	cell_face_3D();
	VEC3_TYPE calc_approximate_normal() const;
    //2つの三角形に分割して四辺形の面積を計算(一般に3次元空間の四辺形は平面にならないからこれは近似的な面積)
	MY_FLOAT_TYPE calc_approximate_face_area() const;
	//face の中心を計算
	VEC3_TYPE calc_face_center() const;
	//\xi - \eta パラメータ空間から実空間への座標変換
	VEC3_TYPE get_real_space_coordinate_from_xi_eta_space(const MY_FLOAT_TYPE xi, const MY_FLOAT_TYPE eta) const;
};
