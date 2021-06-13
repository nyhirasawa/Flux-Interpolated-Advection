#pragma once

#include <Eigen/Dense>
#include <vector>

#include "define_float_type.h"

class cell_vertex {
public:
	std::vector<int> _vertex_index;
	VEC3_TYPE _vertex_pos;

	cell_vertex();
	cell_vertex(VEC3_TYPE vertex_pos);
	cell_vertex(std::vector<int> vertex_index, const MY_FLOAT_TYPE cell_length);
	cell_vertex(std::vector<int> vertex_index, VEC3_TYPE vertex_pos);
};
