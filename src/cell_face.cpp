#include "cell_face.h"

#include "define_float_type.h"

cell_face::cell_face(): is_boundary_face() {
	_vertex_list.resize(2);
}

VEC3_TYPE cell_face::calc_normal() const {
	VEC3_TYPE edge_direction = _vertex_list[1]._vertex_pos - _vertex_list[0]._vertex_pos;
	VEC3_TYPE normal = edge_direction.cross(VEC3_TYPE(0.0, 0.0, 1.0)).normalized();
	return normal;
}

MY_FLOAT_TYPE cell_face::calc_face_area() const {
	VEC3_TYPE edge = _vertex_list[1]._vertex_pos - _vertex_list[0]._vertex_pos;
	return edge.norm();
}

//face の中心を計算
VEC3_TYPE cell_face::calc_face_center() const {
	return (_vertex_list[0]._vertex_pos + _vertex_list[1]._vertex_pos) / 2.0;
}

//面を構成する2点をスワップする
void cell_face::invert_vertex() {
	cell_vertex tmp = _vertex_list[0];
	_vertex_list[0] = _vertex_list[1];
	_vertex_list[1] = tmp;
}
