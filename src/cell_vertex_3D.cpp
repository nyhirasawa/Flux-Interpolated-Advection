#include "cell_vertex_3D.h"

#include <vector>
#include "physical_const.h"
#include "define_float_type.h"

cell_vertex_3D::cell_vertex_3D() : _vertex_pos() { _vertex_index.resize(3); }
cell_vertex_3D::cell_vertex_3D(VEC3_TYPE vertex_pos) : _vertex_pos(vertex_pos) { _vertex_index.resize(3); }
cell_vertex_3D::cell_vertex_3D(Eigen::Vector3i vertex_index, const MY_FLOAT_TYPE cell_length) : _vertex_index(vertex_index) {
	_vertex_pos[0] = vertex_index[0] * cell_length;
	_vertex_pos[1] = vertex_index[1] * cell_length;
	_vertex_pos[2] = vertex_index[2] * cell_length;
}
cell_vertex_3D::cell_vertex_3D(Eigen::Vector3i vertex_index, VEC3_TYPE vertex_pos) : _vertex_index(vertex_index), _vertex_pos(vertex_pos) {}
