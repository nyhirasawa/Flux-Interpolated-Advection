#pragma once

#include <vector>
#include <Eigen/Dense>
#include "cell_vertex_3D.h"
#include "cell_face_3D.h"

#include "define_float_type.h"

class triangle_3D {
public:
	std::vector<cell_vertex_3D> _vertex_list;

    // デフォルトコンストラクタ
    triangle_3D();
	//3つの頂点座標から三角形を構成
	triangle_3D(VEC3_TYPE vertex_pos_0, VEC3_TYPE vertex_pos_1, VEC3_TYPE vertex_pos_2);


    VEC3_TYPE calc_center_position() const;
	//ヘロンの公式を使って面積を計算する
	MY_FLOAT_TYPE calc_area() const;
	VEC3_TYPE calc_normal() const;
	VEC3_TYPE calc_area_scaled_normal() const;
	//このポリゴンを重心で複数の三角形に分割する。
	std::vector<triangle_3D> split_this_polygon_to_triangles_by_barycenter() const;
};
