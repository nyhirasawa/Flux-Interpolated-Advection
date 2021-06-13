#ifndef SMOKE_SIMULATION_GRID_H
#define SMOKE_SIMULATION_GRID_H

#include <vector>
#include <Eigen/Dense>

#include "cell_face.h"

#include "define_float_type.h"

namespace smoke_simulation {
	//セルの境界条件
	enum BoundaryCondition {
		FLUID = 0,
		WALL = 1,
	};

	//系の物理量が乗るグリッドの定義
	class Grid {
	public:
		const int Grid_num_x, Grid_num_y;
		const MY_FLOAT_TYPE min_pos_x, min_pos_y;
		const MY_FLOAT_TYPE _cell_length, _cell_volume;
//		const std::string _interpolation_method;
		std::vector<MY_FLOAT_TYPE> velocity_in_voxel_face_x, velocity_in_voxel_face_y;
		std::vector<MY_FLOAT_TYPE> pressure;
		std::vector<MY_FLOAT_TYPE> substance_density;
		std::vector<MY_FLOAT_TYPE> psi_substance_density_cell_face_x, psi_substance_density_cell_face_y;
		std::vector<MY_FLOAT_TYPE> psi_velocity_cell_vertex_x;
		std::vector<MY_FLOAT_TYPE> psi_velocity_cell_center_y;
		std::vector<MY_FLOAT_TYPE> cell_volume_cell_center;
		std::vector<MY_FLOAT_TYPE> cell_volume_cell_face_x;
		std::vector<MY_FLOAT_TYPE> cell_volume_cell_face_y;

		Grid(int nx, int ny, MY_FLOAT_TYPE cell_length, const std::string interpolation_method);
		~Grid();
		//全cellのcell faceでのpsi(deisity)を計算する
//		void calc_psi_substance_density_cell_face();
		//(vertex_index_x, vertex_index_y)番目の頂点での速度を、隣り合うface の値からinterpolationによって計算
		VEC3_TYPE calc_vertex_velocity(int vertex_index_x, int vertex_index_y) const;
		//position での psi(deisity) の値を周辺の cell face での値から interpolation して求める。
		VEC3_TYPE calc_psi_substance_density_by_interpolation_cell_face(
			VEC3_TYPE position,
			const std::string interpolation_method,
			const bool use_integral,
			const int num_gauss_quadrature_point_for_integrate_density
		) const;

		//position での psi(velocity_x) の値を周辺の離散値から interpolation して求める。
		VEC3_TYPE calc_psi_velocity_x_by_interpolation(
			VEC3_TYPE position,
			const std::string interpolation_method,
			const bool use_integral,
			const int num_gauss_quadrature_point_for_integrate_density
		) const;
		//position での psi(velocity_y) の値を周辺の離散値から interpolation して求める。
		VEC3_TYPE calc_psi_velocity_y_by_interpolation(
			VEC3_TYPE position,
			const std::string interpolation_method,
			const bool use_integral,
			const int num_gauss_quadrature_point_for_integrate_density
		) const;

		//position での psi(deisity) の値を周辺の cell center での値から interpolation して求める。
		//VEC3_TYPE calc_psi_substance_density_by_interpolation(VEC3_TYPE position, bool& is_out_fo_range);
		//position での deisity の値を周辺の cell center での値から interpolation して求める。
		MY_FLOAT_TYPE calc_substance_density_by_interpolation(
			VEC3_TYPE position,
			const std::string interpolation_method) const;

		// negative density artifact を解消するために, バックトレースしたセルの面をx軸(またはy軸, およびその両方)で切断する関数
		std::vector<cell_face> split_cell_face_by_axis(cell_face face, std::string split_mode) const;
		// negative density artifact を解消するために, バックトレースしたセルの面を軸(x軸またはy軸, およびその両方)で切断する関数
		// psi が x方向に(psi の積分の方向をxにとるならy方向に)区分線型関数であることを仮定したときの分割の方法
		// split_mode == "x" は未実装
		std::vector<cell_face> split_cell_face_by_cell_center(const cell_face face, std::string split_mode) const;
		// バックトレースしたセルの面を等間隔で split_face_num 個に切断する関数
		std::vector<cell_face> split_cell_face_uniformly(cell_face face, const int split_face_num) const;
		// バックトレースしたセルの面を軸(x軸またはy軸, およびその両方)に平行な平面で切断する関数
		// split_mode == "x" の時は x軸に平行な, split_mode == "y" の時は y軸に平行な平面で切断する.
		// split_positions は切断する平面の位置を表す。例えばsplit_mode == "y"でsplit_positions={1.4, 1.8, 2.5}の時はx=1.4, 1.8, 2.5 の三つの平面で切断する。
		// split_mode == "x" は未実装
		std::vector<cell_face> split_cell_face_by_axis_aligned_plane(const cell_face face, std::vector<MY_FLOAT_TYPE> split_positions, std::string split_mode) const;

		//position での速度をinterpolationによって計算
		VEC3_TYPE calc_interpolated_velocity(
			VEC3_TYPE position,
			const std::string interpolation_method
		) const;
		//position での速度をinterpolationによって計算
		MY_FLOAT_TYPE calc_interpolated_velocity_x(
			VEC3_TYPE position,
			const std::string interpolation_method
		) const;
		//position での速度をinterpolationによって計算
		MY_FLOAT_TYPE calc_interpolated_velocity_y(
			VEC3_TYPE position,
			const std::string interpolation_method
		) const;

		// cell center で定義される量cell_center_valuesのpositionでの補間を tri-linear interpolation で計算する
	    MY_FLOAT_TYPE interpolate_cell_center_defined_values_trilinear(const VEC3_TYPE &position, const std::vector<MY_FLOAT_TYPE>& cell_center_values) const;
	};
}// namespace smoke_simulation
#endif//SMOKE_SIMULATION_GRID_H
