#include "triangle_3D.h"

#include "define_float_type.h"

triangle_3D::triangle_3D(){
    _vertex_list.resize(3);
}

triangle_3D::triangle_3D(VEC3_TYPE vertex_pos_0, VEC3_TYPE vertex_pos_1, VEC3_TYPE vertex_pos_2){
    _vertex_list.resize(3);
    _vertex_list[1]._vertex_pos = vertex_pos_1;
    _vertex_list[0]._vertex_pos = vertex_pos_0;
    _vertex_list[2]._vertex_pos = vertex_pos_2;
}


VEC3_TYPE triangle_3D::calc_center_position() const{
    return (_vertex_list[0]._vertex_pos + _vertex_list[1]._vertex_pos + _vertex_list[2]._vertex_pos)/3.0;
}

//ヘロンの公式を使って面積を計算する
MY_FLOAT_TYPE triangle_3D::calc_area() const{
	/*
    MY_FLOAT_TYPE edge_length_a = (_vertex_list[1]._vertex_pos - _vertex_list[0]._vertex_pos).norm();
    MY_FLOAT_TYPE edge_length_b = (_vertex_list[2]._vertex_pos - _vertex_list[1]._vertex_pos).norm();
    MY_FLOAT_TYPE edge_length_c = (_vertex_list[0]._vertex_pos - _vertex_list[2]._vertex_pos).norm();
    MY_FLOAT_TYPE s = (edge_length_a + edge_length_b + edge_length_c) / 2.0;

    MY_FLOAT_TYPE area = sqrt(s * (s - edge_length_a) * (s - edge_length_b) * (s - edge_length_c));
    //三角形の面積が0になるようなつぶれた三角形にヘロンの公式を使うと nan や inf が発生してしまうみたいなので、nan や inf が発生した場合は0を返すことにする
    if(std::isfinite(area)){
        return area;
    }
    else{
        return 0.0;
    }
	*/
    VEC3_TYPE edge_0 = _vertex_list[1]._vertex_pos - _vertex_list[0]._vertex_pos;
    VEC3_TYPE edge_1 = _vertex_list[2]._vertex_pos - _vertex_list[0]._vertex_pos;
    return 0.5 * edge_0.cross(edge_1).norm();
}

VEC3_TYPE triangle_3D::calc_normal() const{
    VEC3_TYPE edge_0 = _vertex_list[1]._vertex_pos - _vertex_list[0]._vertex_pos;
    VEC3_TYPE edge_1 = _vertex_list[2]._vertex_pos - _vertex_list[0]._vertex_pos;
    return edge_0.cross(edge_1).normalized();
}

VEC3_TYPE triangle_3D::calc_area_scaled_normal() const {
    VEC3_TYPE edge_0 = _vertex_list[1]._vertex_pos - _vertex_list[0]._vertex_pos;
    VEC3_TYPE edge_1 = _vertex_list[2]._vertex_pos - _vertex_list[0]._vertex_pos;
    return 0.5 * edge_0.cross(edge_1);
}


//このポリゴンを重心で複数の三角形に分割する。
std::vector<triangle_3D> triangle_3D::split_this_polygon_to_triangles_by_barycenter() const {
    const int num_vertex_in_polygon = _vertex_list.size();
    //結果を格納する変数
    std::vector<triangle_3D> splitted_triangle_list;
    for(int i_vert = 0; i_vert < num_vertex_in_polygon; ++i_vert){
        triangle_3D added_triangle;
        added_triangle._vertex_list[0] = _vertex_list[i_vert];
        added_triangle._vertex_list[1] = _vertex_list[(i_vert + 1) % num_vertex_in_polygon];
        added_triangle._vertex_list[2] = cell_vertex_3D(calc_center_position());
        splitted_triangle_list.push_back(added_triangle);
    }
    return splitted_triangle_list;
}
