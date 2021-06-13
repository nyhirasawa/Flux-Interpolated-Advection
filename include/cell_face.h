#pragma once

#include <vector>
#include <Eigen/Dense>
#include "cell_vertex.h"

#include "define_float_type.h"

class cell_face {
public:
	std::vector<cell_vertex> _vertex_list;
	bool is_boundary_face;

	cell_face();
	VEC3_TYPE calc_normal() const;
	//cell face の面積を計算する(2次元の場合は線の長さになる)
	MY_FLOAT_TYPE calc_face_area() const;
	//face の中心を計算
	VEC3_TYPE calc_face_center() const;
	//面を構成する2点をスワップする
	void invert_vertex();
};
