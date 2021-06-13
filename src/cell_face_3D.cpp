#include "cell_face_3D.h"

#include "define_float_type.h"

cell_face_3D::cell_face_3D()/*: is_boundary_face()*/ {
	_vertex_list.resize(4);
}

VEC3_TYPE cell_face_3D::calc_approximate_normal() const {
    VEC3_TYPE edge_direction_0 = _vertex_list[1]._vertex_pos - _vertex_list[0]._vertex_pos;
    VEC3_TYPE edge_direction_1 = _vertex_list[2]._vertex_pos - _vertex_list[1]._vertex_pos;
    VEC3_TYPE edge_direction_2 = _vertex_list[3]._vertex_pos - _vertex_list[2]._vertex_pos;
    VEC3_TYPE edge_direction_3 = _vertex_list[0]._vertex_pos - _vertex_list[3]._vertex_pos;
    VEC3_TYPE normal_0 = edge_direction_0.cross(edge_direction_1).normalized();
    VEC3_TYPE normal_1 = edge_direction_1.cross(edge_direction_2).normalized();
    VEC3_TYPE normal_2 = edge_direction_2.cross(edge_direction_3).normalized();
    VEC3_TYPE normal_3 = edge_direction_3.cross(edge_direction_0).normalized();
	return (normal_0 + normal_1 + normal_2 +normal_3).normalized();
}

//2つの三角形に分割して四辺形の面積を計算(一般に3次元空間の四辺形は平面にならないからこれは近似的な面積)
MY_FLOAT_TYPE cell_face_3D::calc_approximate_face_area() const {
    MY_FLOAT_TYPE edge_length_0, edge_length_1, edge_length_2, edge_length_3;
    edge_length_0 = (_vertex_list[0]._vertex_pos - _vertex_list[1]._vertex_pos).norm();
    edge_length_1 = (_vertex_list[1]._vertex_pos - _vertex_list[2]._vertex_pos).norm();
    edge_length_2 = (_vertex_list[2]._vertex_pos - _vertex_list[3]._vertex_pos).norm();
    edge_length_3 = (_vertex_list[3]._vertex_pos - _vertex_list[0]._vertex_pos).norm();
    MY_FLOAT_TYPE diagonal_length_0, diagonal_length_1;
    diagonal_length_0 = (_vertex_list[1]._vertex_pos - _vertex_list[3]._vertex_pos).norm();
    diagonal_length_1 = (_vertex_list[0]._vertex_pos - _vertex_list[2]._vertex_pos).norm();

    MY_FLOAT_TYPE s0, s1, s2, s3;
    s0 = (edge_length_0 + edge_length_3 + diagonal_length_0) / 2.0;
    s1 = (edge_length_1 + edge_length_2 + diagonal_length_0) / 2.0;
    s2 = (edge_length_0 + edge_length_1 + diagonal_length_1) / 2.0;
    s3 = (edge_length_2 + edge_length_3 + diagonal_length_1) / 2.0;
    MY_FLOAT_TYPE triangle_area_0, triangle_area_1, triangle_area_2, triangle_area_3;
    triangle_area_0 = sqrt(s0 * (s0 - edge_length_0) * (s0 - edge_length_3) * (s0 - diagonal_length_0));
    triangle_area_1 = sqrt(s1 * (s1 - edge_length_1) * (s1 - edge_length_2) * (s1 - diagonal_length_0));
    triangle_area_2 = sqrt(s2 * (s2 - edge_length_0) * (s2 - edge_length_1) * (s2 - diagonal_length_1));
    triangle_area_3 = sqrt(s3 * (s3 - edge_length_2) * (s3 - edge_length_3) * (s3 - diagonal_length_1));
    return((triangle_area_0 + triangle_area_1) + (triangle_area_2 + triangle_area_3)) / 2.0;
}

//face の中心を計算
VEC3_TYPE cell_face_3D::calc_face_center() const {
	return (_vertex_list[0]._vertex_pos
          + _vertex_list[1]._vertex_pos
          + _vertex_list[2]._vertex_pos
          + _vertex_list[3]._vertex_pos) / 4.0;
}

//\xi - \eta パラメータ空間から実空間への座標変換
VEC3_TYPE cell_face_3D::get_real_space_coordinate_from_xi_eta_space(const MY_FLOAT_TYPE xi, const MY_FLOAT_TYPE eta) const{
	return (1.0 - xi) * (1.0 - eta) * _vertex_list[0]._vertex_pos
			+ xi * (1.0 - eta) * _vertex_list[1]._vertex_pos
			+ xi * eta * _vertex_list[2]._vertex_pos
			+ (1.0 - xi) * eta * _vertex_list[3]._vertex_pos;
}
