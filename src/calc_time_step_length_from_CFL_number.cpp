#include "calc_time_step_length_from_CFL_number.h"

#include "grid.h"
#include "grid_3d.h"
#include "grid_1d.h"
#include "physical_const.h"
#include "utils.h"
#include "define_float_type.h"

namespace smoke_simulation {
	MY_FLOAT_TYPE calc_time_step_length_from_CFL_number_1D(
		const Grid_1D& all_grid,
		const MY_FLOAT_TYPE now_time,
		const int i_movie_frame,
		bool& write_frame_to_movie_file,
		const int CFL_number,
		const MY_FLOAT_TYPE movie_fps){
		//速度の絶対値の最大値を求める
		MY_FLOAT_TYPE max_absolute_velocity = 0.0;
		for (int iy = 0; iy < all_grid.Grid_num_y + 1; ++iy) {
			if (abs(all_grid.velocity_in_voxel_face_y[iy]) > max_absolute_velocity) {
				max_absolute_velocity = abs(all_grid.velocity_in_voxel_face_y[iy]);
			}
		}
		//CFL number から計算したtime step length
		MY_FLOAT_TYPE tiem_step_length_from_CFL_number = CFL_number * all_grid._cell_length / max_absolute_velocity;
		//次の動画のフレームの時刻
		MY_FLOAT_TYPE next_movie_frame_time = (i_movie_frame + 1) / movie_fps;
		// time step length として tiem_step_length_from_CFL_number を採用したときの次時刻が、
		// 次の動画のフレームの時刻 を超える場合は, 次時刻が次の動画のフレームの時刻になるように調節する
		if (now_time + tiem_step_length_from_CFL_number > next_movie_frame_time) {
			write_frame_to_movie_file = true;
			return (next_movie_frame_time - now_time);
		}
		else {
			write_frame_to_movie_file = false;
			return tiem_step_length_from_CFL_number;
		}
	}

	MY_FLOAT_TYPE calc_time_step_length_from_CFL_number(
		const Grid& all_grid,
		const MY_FLOAT_TYPE now_time,
		const int i_movie_frame,
		bool& write_frame_to_movie_file,
		const int CFL_number,
		const MY_FLOAT_TYPE movie_fps
	) {
		//速度の絶対値の最大値を求める
		MY_FLOAT_TYPE max_absolute_velocity = 0.0;
		for (int ix = 0; ix < all_grid.Grid_num_x + 1; ++ix) {
			for (int iy = 0; iy < all_grid.Grid_num_y; ++iy) {
				if (abs(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]) > max_absolute_velocity) {
					max_absolute_velocity = abs(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]);
				}
			}
		}
		for (int ix = 0; ix < all_grid.Grid_num_x; ++ix) {
			for (int iy = 0; iy < all_grid.Grid_num_y + 1; ++iy) {
				if (abs(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]) > max_absolute_velocity) {
					max_absolute_velocity = abs(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]);
				}
			}
		}
		//CFL number から計算したtime step length
		MY_FLOAT_TYPE tiem_step_length_from_CFL_number = CFL_number * all_grid._cell_length / max_absolute_velocity;
		//次の動画のフレームの時刻
		MY_FLOAT_TYPE next_movie_frame_time = (i_movie_frame + 1) / movie_fps;
		// time step length として tiem_step_length_from_CFL_number を採用したときの次時刻が、
		// 次の動画のフレームの時刻 を超える場合は, 次時刻が次の動画のフレームの時刻になるように調節する
		if (now_time + tiem_step_length_from_CFL_number > next_movie_frame_time) {
			write_frame_to_movie_file = true;
			return (next_movie_frame_time - now_time);
		}
		else {
			write_frame_to_movie_file = false;
			return tiem_step_length_from_CFL_number;
		}
	}

	MY_FLOAT_TYPE calc_time_step_length_from_CFL_number_3D(
		const Grid_3D& all_grid,
		const MY_FLOAT_TYPE now_time,
		const int i_movie_frame,
		bool& write_frame_to_movie_file,
		const int CFL_number,
		const MY_FLOAT_TYPE movie_fps
	) {
		//速度の絶対値の最大値を求める
		MY_FLOAT_TYPE max_absolute_velocity = 0.0;
		for (int ix = 0; ix < all_grid.Grid_num_x + 1; ++ix) {
			for (int iy = 0; iy < all_grid.Grid_num_y; ++iy) {
				for (int iz = 0; iz < all_grid.Grid_num_z; ++iz) {
					if (abs(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]) > max_absolute_velocity) {
						max_absolute_velocity = abs(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]);
					}
				}
			}
		}
		for (int ix = 0; ix < all_grid.Grid_num_x ; ++ix) {
			for (int iy = 0; iy < all_grid.Grid_num_y + 1; ++iy) {
				for (int iz = 0; iz < all_grid.Grid_num_z; ++iz) {
					if (abs(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]) > max_absolute_velocity) {
						max_absolute_velocity = abs(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]);
					}
				}
			}
		}
		for (int ix = 0; ix < all_grid.Grid_num_x; ++ix) {
			for (int iy = 0; iy < all_grid.Grid_num_y; ++iy) {
				for (int iz = 0; iz < all_grid.Grid_num_z + 1; ++iz) {
					if (abs(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]) > max_absolute_velocity) {
						max_absolute_velocity = abs(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]);
					}
				}
			}
		}
		//CFL number から計算したtime step length
		MY_FLOAT_TYPE tiem_step_length_from_CFL_number = CFL_number * all_grid._cell_length / max_absolute_velocity;
		//次の動画のフレームの時刻
		MY_FLOAT_TYPE next_movie_frame_time = (i_movie_frame + 1) / movie_fps;
		// time step length として tiem_step_length_from_CFL_number を採用したときの次時刻が、
		// 次の動画のフレームの時刻 を超える場合は, 次時刻が次の動画のフレームの時刻になるように調節する
		if (now_time + tiem_step_length_from_CFL_number > next_movie_frame_time) {
			write_frame_to_movie_file = true;
			return (next_movie_frame_time - now_time);
		}
		else {
			write_frame_to_movie_file = false;
			return tiem_step_length_from_CFL_number;
		}
	}

	MY_FLOAT_TYPE calc_time_step_length_from_fixed_time_step_length(
		const MY_FLOAT_TYPE now_time,
		const int i_movie_frame,
		bool& write_frame_to_movie_file,
		const MY_FLOAT_TYPE movie_fps
	) {
		// 固定の time step length
		MY_FLOAT_TYPE fixed_tiem_step_length = physical_const::kDt_fixed;
		//次の動画のフレームの時刻
		MY_FLOAT_TYPE next_movie_frame_time = (i_movie_frame + 1) / movie_fps;
		// time step length として fixed_tiem_step_length を採用したときの次時刻が、
		// 次の動画のフレームの時刻 を超える場合は, 次時刻が次の動画のフレームの時刻になるように調節する
		if (now_time + fixed_tiem_step_length > next_movie_frame_time) {
			write_frame_to_movie_file = true;
			return (next_movie_frame_time - now_time);
		}
		else {
			write_frame_to_movie_file = false;
			return fixed_tiem_step_length;
		}
	}

}//namespace smoke_simulation
