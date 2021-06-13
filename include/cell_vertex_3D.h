#pragma once

#include <Eigen/Dense>
#include <vector>

#include "define_float_type.h"

class cell_vertex_3D {
public:
	Eigen::Vector3i _vertex_index;
	VEC3_TYPE _vertex_pos;

	cell_vertex_3D();
	cell_vertex_3D(VEC3_TYPE vertex_pos);
	cell_vertex_3D(Eigen::Vector3i vertex_index, const MY_FLOAT_TYPE cell_length);
	cell_vertex_3D(Eigen::Vector3i vertex_index, VEC3_TYPE vertex_pos);
};
