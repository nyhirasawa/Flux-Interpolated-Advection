#ifndef SMOKE_SIMULATION_GRID_3D_H
#define SMOKE_SIMULATION_GRID_3D_H

#include <vector>
#include <string>
#include <Eigen/Dense>

//#include "polygon_3D.h"

#include "define_float_type.h"

namespace smoke_simulation{
    //系の物理量が乗るグリッドの定義
    class Grid_3D{
    public:
        //position での psi_density の x 成分の値を周辺の cell face での値から linear interpolation して求める。
        //(calc_interpolated_velocity のヘルプ関数)
//        MY_FLOAT_TYPE calc_psi_density_x_by_1Dy_linear_integral_cell_face_x(VEC3_TYPE position) const;
        //position での psi_density の y 成分の値を周辺の cell face での値から linear interpolation して求める。
        //(calc_interpolated_velocity のヘルプ関数)
        MY_FLOAT_TYPE calc_psi_density_y_by_1Dy_linear_integral_cell_face_y(VEC3_TYPE position) const;
        //position での psi_density の z 成分の値を周辺の cell face での値から linear interpolation して求める。
        //(calc_interpolated_velocity のヘルプ関数)
//        MY_FLOAT_TYPE calc_psi_density_z_by_1Dy_linear_integral_cell_face_z(VEC3_TYPE position) const;
        //position での psi_density の x 成分の値を周辺の cell face での値から linear interpolation して求める。
        //(calc_interpolated_velocity のヘルプ関数)
//        MY_FLOAT_TYPE calc_psi_density_x_by_linear_integral_cell_face_x(VEC3_TYPE position) const;
        //position での psi_density の y 成分の値を周辺の cell face での値から linear interpolation して求める。
        //(calc_interpolated_velocity のヘルプ関数)
        MY_FLOAT_TYPE calc_psi_density_y_by_linear_integral_cell_face_y(VEC3_TYPE position) const;
        //position での psi_density の z 成分の値を周辺の cell face での値から linear interpolation して求める。
        //(calc_interpolated_velocity のヘルプ関数)
//        MY_FLOAT_TYPE calc_psi_density_z_by_linear_integral_cell_face_z(VEC3_TYPE position) const;

        const int Grid_num_x, Grid_num_y, Grid_num_z;
        const MY_FLOAT_TYPE min_pos_x, min_pos_y, min_pos_z;
        const MY_FLOAT_TYPE _cell_length, _cell_face_area, _cell_volume;
//        const std::string _interpolation_method;
//        std::string _interpolation_method;
        std::vector<MY_FLOAT_TYPE> velocity_in_voxel_face_x, velocity_in_voxel_face_y, velocity_in_voxel_face_z;
        std::vector<MY_FLOAT_TYPE> pressure;
        std::vector<MY_FLOAT_TYPE> substance_density;

        std::vector<MY_FLOAT_TYPE> psi_substance_density_cell_face_x, psi_substance_density_cell_face_y, psi_substance_density_cell_face_z;
        std::vector<MY_FLOAT_TYPE> cell_volume_cell_center;

        std::vector<MY_FLOAT_TYPE> psi_velocity_x_cell_face;
        std::vector<MY_FLOAT_TYPE> psi_velocity_y_cell_face;
        std::vector<MY_FLOAT_TYPE> psi_velocity_z_cell_face;

        //元時刻までで, 壁から逃げていった運動量の合計(初期時刻から元時刻までの間に逃げていった運動量の合計を全て足し合わせたもの)
        MY_FLOAT_TYPE total_loss_of_momentum_x;
        MY_FLOAT_TYPE total_loss_of_momentum_y;
        MY_FLOAT_TYPE total_loss_of_momentum_z;

//        std::vector<std::vector<MY_FLOAT_TYPE>> external_force_field;

        Grid_3D(int nx, int ny, int nz, MY_FLOAT_TYPE cell_length);
        ~Grid_3D();
        //全cellのcell faceでのpsi(deisity)を計算する
//        void calc_psi_substance_density_cell_face_3D();
        //position での速度をinterpolationによって計算
        VEC3_TYPE calc_interpolated_velocity(VEC3_TYPE position, const std::string interpolation_method) const;
        //position での psi_density をinterpolationによって計算
        VEC3_TYPE calc_interpolated_psi_density(VEC3_TYPE position, const std::string interpolation_method) const;
        //position での psi_density をinterpolationによって計算
        VEC3_TYPE calc_interpolated_psi_velocity(const int dim, VEC3_TYPE position, const std::string interpolation_method) const;
        //position での psi_density をinterpolationによって計算
        VEC3_TYPE calc_interpolated_psi_velocity_x(VEC3_TYPE position, const std::string interpolation_method) const;
        //position での psi_density をinterpolationによって計算
        VEC3_TYPE calc_interpolated_psi_velocity_y(VEC3_TYPE position, const std::string interpolation_method) const;
        //position での psi_density をinterpolationによって計算
        VEC3_TYPE calc_interpolated_psi_velocity_z(VEC3_TYPE position, const std::string interpolation_method) const;
        //position での deisity の値を周辺の cell center での値から interpolation して求める。
        MY_FLOAT_TYPE calc_substance_density_by_interpolation(const VEC3_TYPE &position, const std::string interpolation_method) const;
        // cell center で定義される量cell_center_valuesのpositionでの補間を tri-linear interpolation で計算する
//        MY_FLOAT_TYPE interpolate_cell_center_defined_values_trilinear(const VEC3_TYPE &position, const std::vector<MY_FLOAT_TYPE>& cell_center_values) const;
        // cell center で定義される量cell_center_valuesのpositionでの補間を y方向の1次元linear interpolation で計算する
        MY_FLOAT_TYPE interpolate_cell_center_defined_values_y_direction_1d_linear(const VEC3_TYPE &position, const std::vector<MY_FLOAT_TYPE>& cell_center_values) const;
        // セルの頂点の速度を隣接する面での速度の平均によって計算する
        VEC3_TYPE calc_cell_vertex_velocity(const Eigen::Vector3i cell_vertex_index) const;

        // negative density artifact を解消するために, バックトレースしたセルの面を軸に垂直な平面(x or y or z軸, およびそれら)で切断する関数
        // psi が x方向に(psi の積分の方向をxにとるならy方向に)区分線型関数であることを仮定したときの分割の方法
        // split_mode は切断の方向を設定する。例えば(split_mode == "x" ならx軸に垂直な面で切断する)
        // iso_value は切断する位置を決める。例えば split_mode == "x" で iso_value==1.9 なら x = 1.9 の平面でfaceを切断する
//        std::vector<polygon_3D> split_cell_face_by_axis_aligned_plane(const cell_face_3D face, const MY_FLOAT_TYPE iso_value, std::string split_mode) const;
    };
}
#endif//SMOKE_SIMULATION_GRID_3D_H
