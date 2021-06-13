#include "grid.h"
#include "grid_3d.h"
#include "initialize_grid.h"
#include "physical_const.h"
#include "utils.h"
#include "define_float_type.h"

#include <math.h>

namespace smoke_simulation {
	void initialize_grid(
        Grid_1D& all_grid
    ){
		//速度の初期化
		//velocity_in_voxel_face_xの初期化
		for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
			all_grid.velocity_in_voxel_face_y[iy] = 1.0;
		}
		//密度の初期化
		// substance_density の初期化
		for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
			MY_FLOAT_TYPE cell_center_pos = (0.5 + iy) * all_grid._cell_length;
			if(cell_center_pos > 0.1 && cell_center_pos < 0.2){
				all_grid.substance_density[iy] = 1.0;
			}
			else {
				all_grid.substance_density[iy] = 0.0;
			}
		}
	}

	//グリッドを初期化する関数
	void initialize_grid(
		Grid& all_grid,
		const bool fix_velocity,
		const bool initialize_velocity,
		const bool initialize_density
	) {
		//四角形に密度, 速度場を配置
		if (physical_const::kInitial_density_shape == "square") {
			const int width_of_density = (all_grid.Grid_num_x / 9)/2;
//			const int half_index = 32;
			const int half_index = all_grid.Grid_num_x / 2;
			//速度の初期化
			if(initialize_velocity){
				//velocity_in_voxel_face_xの初期化
				for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
					for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
						if (   ix >= half_index - width_of_density
							&& ix < half_index + width_of_density+1
							&& iy >= half_index - width_of_density
							&& iy < half_index + width_of_density+1) {
							all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 0.0;
						}
						else {
							all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 0.0;
						}
					}
				}
				//velocity_in_voxel_face_yの初期化
				for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
					for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
						if (   ix >= half_index - width_of_density
							&& ix < half_index + width_of_density+1
							&& iy >= half_index - width_of_density
							&& iy < half_index + width_of_density+1) {
							all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 10.0;
						}
						else {
							all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 0.0;
						}
					}
				}
			}
			//密度の初期化
			if(initialize_density){
				// substance_density の初期化
				for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
					for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
						if (   ix >= half_index - width_of_density
							&& ix < half_index + width_of_density+1
							&& iy >= half_index - width_of_density
							&& iy < half_index + width_of_density+1) {
							all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 1.0;
						}
						else {
							all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 0.0;
						}
					}
				}
			}
		}
		//円形に密度, 速度場を配置
		else if (physical_const::kInitial_density_shape == "circle") {
			//半径
			const MY_FLOAT_TYPE radius = 0.1;
			// center pos of grid
			VEC3_TYPE center_pos_of_grid(
				all_grid.min_pos_x + (all_grid.Grid_num_x / 2) * all_grid._cell_length,
				all_grid.min_pos_y + (all_grid.Grid_num_y / 2) * all_grid._cell_length,
				0.0);
			//速度の初期化
			if(initialize_velocity){
				//velocity_in_voxel_face_xの初期化
				for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
					for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
						VEC3_TYPE face_pos_ixiy(
							all_grid.min_pos_x +  ix        * all_grid._cell_length,
							all_grid.min_pos_y + (iy + 0.5) * all_grid._cell_length,
							0.0);

						if ((face_pos_ixiy - center_pos_of_grid).norm() < radius) {
							all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 10.0;
						}
						else {
							all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 0.0;
						}
					}
				}
				//velocity_in_voxel_face_yの初期化
				for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
					for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
						VEC3_TYPE face_pos_ixiy(
							all_grid.min_pos_x + (ix + 0.5) * all_grid._cell_length,
							all_grid.min_pos_y + iy * all_grid._cell_length,
							0.0);

						if ((face_pos_ixiy - center_pos_of_grid).norm() < radius) {
							all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 0.0;
						}
						else {
							all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 0.0;
						}
					}
				}
			}
			//密度の初期化
			if(initialize_density){
				// substance_density の初期化
				for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
					for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
						VEC3_TYPE face_pos_ixiy(
							all_grid.min_pos_x + ix * all_grid._cell_length,
							all_grid.min_pos_y + iy * all_grid._cell_length,
							0.0);
						if ((face_pos_ixiy - center_pos_of_grid).norm() < radius) {
							all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 1.0;
						}
						else {
							all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 0.0;
						}
					}
				}
			}
		}
		else if(physical_const::kInitial_density_shape == "Tylor-Green"){
			//速度の初期化
			if(initialize_velocity){
				//velocity_in_voxel_face_xの初期化
				for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
					for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
						//考えているセルの面の位置座標
						const MY_FLOAT_TYPE position_x = ix * all_grid._cell_length;
						const MY_FLOAT_TYPE position_y = (iy + 0.5) * all_grid._cell_length;
//						all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//							= cos(2.0 * smoke_simulation::physical_const::kPI * position_x) * sin(2.0 * smoke_simulation::physical_const::kPI * position_y);
						all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
							= sin(2.0 * smoke_simulation::physical_const::kPI * position_x) * cos(2.0 * smoke_simulation::physical_const::kPI * position_y);
					}
				}
				//velocity_in_voxel_face_yの初期化
				for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
					for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
						//考えているセルの面の位置座標
						const MY_FLOAT_TYPE position_x = (ix + 0.5) * all_grid._cell_length;
						const MY_FLOAT_TYPE position_y = iy * all_grid._cell_length;
						all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
							= -cos(2.0 * smoke_simulation::physical_const::kPI * position_x) * sin(2.0 * smoke_simulation::physical_const::kPI * position_y);
//						all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//							= -sin(2.0 * smoke_simulation::physical_const::kPI * position_x) * cos(2.0 * smoke_simulation::physical_const::kPI * position_y);
					}
				}
			}
			//密度の初期化
			if(initialize_density){
				for(int ix = 0; ix < all_grid.Grid_num_x; ++ix){
					for(int iy = 0; iy < all_grid.Grid_num_y; ++iy){
						all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 0.5;
					}
				}
			}
		}
		else if (physical_const::kInitial_density_shape == "vortex_sheet_2d") {
			////速度場の初期化
			//velocity_in_voxel_face_xの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					VEC3_TYPE p(ix * all_grid._cell_length, (iy + 0.5) * all_grid._cell_length, 0.0);
					MY_FLOAT_TYPE distance_from_center = (p - VEC3_TYPE(0.5, 0.5, 0.0)).norm();
					if( distance_from_center < 0.25) {
						all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
							= 10.0 * distance_from_center * (-p[1] + 0.5);
					}
					else{
						all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
							= 0.0;
					}
				}
			}
			//velocity_in_voxel_face_yの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
					VEC3_TYPE p((ix + 0.5) * all_grid._cell_length, iy * all_grid._cell_length, 0.0);
					MY_FLOAT_TYPE distance_from_center = (p - VEC3_TYPE(0.5, 0.5, 0.0)).norm();
					if( distance_from_center < 0.25) {
						all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
							= 10.0 * distance_from_center * (p[0] - 0.5);
					}
					else{
						all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
							= 0.0;
					}
				}
			}
			////密度場の初期化
			// substance_density の初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					VEC3_TYPE p = VEC3_TYPE(ix + 0.5, iy + 0.5, 0.0) * all_grid._cell_length;
					if( p[0] < 0.5 && p[1] < 0.5 ) {
						all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 0.0;
					}
					if( p[0] < 0.5 && p[1] > 0.5 ) {
						all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 1.0;
					}
					if( p[0] > 0.5 && p[1] < 0.5 ) {
						all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 1.0;
					}
					if( p[0] > 0.5 && p[1] > 0.5 ) {
						all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 0.0;
					}
				}
			}
		}
		//固定した速度場を使う場合
		if(fix_velocity && physical_const::kInitial_density_shape != "Tylor-Green") {
			if(initialize_velocity){
				//velocity_in_voxel_face_xの初期化
				for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
					for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
						all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 0.0;
					}
				}

				//velocity_in_voxel_face_yの初期化
				for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
					for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
						all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = 0.0;
					}
				}
			}
		}
		// セル体積の初期化
		for(int ix = 0; ix < all_grid.Grid_num_x; ++ix){
			for(int iy = 0; iy < all_grid.Grid_num_y; ++iy){
				all_grid.cell_volume_cell_center[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)] = all_grid._cell_volume;
			}
		}
	}

	//グリッドを初期化する関数
	void initialize_Grid_3D(Grid_3D& all_grid, const bool fix_velocity) {
		//四角形に密度, 速度場を配置
		if (physical_const::kInitial_density_shape == "square") {
		    const int width_of_density = (all_grid.Grid_num_x / 9)/2;
//			const int half_index = 32;
			const int half_index = all_grid.Grid_num_x / 2;

			//velocity_in_voxel_face_xの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						if (   ix >= half_index - width_of_density
							&& ix < half_index + width_of_density + 1
							&& iy >= half_index - width_of_density
							&& iy < half_index + width_of_density + 1
							&& iz >= half_index - width_of_density
							&& iz < half_index + width_of_density + 1) {
							all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = -10.0;
						}
						else {
							all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
					}
				}
			}
			//velocity_in_voxel_face_yの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
					}
				}
			}
			//velocity_in_voxel_face_zの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z + 1; iz++) {
						all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
					}
				}
			}
/*
			//velocity_in_voxel_face_xの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
					}
				}
			}
			//velocity_in_voxel_face_yの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						if (   ix >= half_index - width_of_density
							&& ix < half_index + width_of_density + 1
							&& iy >= half_index - width_of_density
							&& iy < half_index + width_of_density + 1
							&& iz >= half_index - width_of_density
							&& iz < half_index + width_of_density + 1) {
							all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 10.0;
						}
						else {
							all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
					}
				}
			}
			//velocity_in_voxel_face_zの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z + 1; iz++) {
						all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
					}
				}
			}
*/
/*
			//velocity_in_voxel_face_xの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						if (   ix >= half_index - width_of_density
							&& ix < half_index + width_of_density + 1
							&& iy >= half_index - width_of_density
							&& iy < half_index + width_of_density + 1
							&& iz >= half_index - width_of_density
							&& iz < half_index + width_of_density + 1) {
							all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
						else {
							all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
					}
				}
			}
			//velocity_in_voxel_face_yの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						if (   ix >= half_index - width_of_density
							&& ix < half_index + width_of_density + 1
							&& iy >= half_index - width_of_density
							&& iy < half_index + width_of_density + 1
							&& iz >= half_index - width_of_density
							&& iz < half_index + width_of_density + 1) {
							all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 10.0;
						}
						else {
							all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
					}
				}
			}
			//velocity_in_voxel_face_zの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z + 1; iz++) {
						if (   ix >= half_index - width_of_density
							&& ix < half_index + width_of_density + 1
							&& iy >= half_index - width_of_density
							&& iy < half_index + width_of_density + 1
							&& iz >= half_index - width_of_density
							&& iz < half_index + width_of_density + 1) {
							all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
						else {
							all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
					}
				}
			}
*/
			// substance_density の初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
//						if (   ix >= half_index - width_of_density-1
//							&& ix < half_index + width_of_density
						if (   ix >= half_index - width_of_density
							&& ix < half_index + width_of_density + 1
							&& iy >= half_index - width_of_density
							&& iy < half_index + width_of_density + 1
							&& iz >= half_index - width_of_density
							&& iz < half_index + width_of_density + 1) {
							all_grid.substance_density[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 1.0;
						}
						else {
							all_grid.substance_density[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
					}
				}
			}
		}
		//円形に密度, 速度場を配置
		else if (physical_const::kInitial_density_shape == "circle") {
			//半径
			const MY_FLOAT_TYPE radius = 0.1;
			// center pos of grid
			VEC3_TYPE center_pos_of_grid(
				all_grid.min_pos_x + ((all_grid.Grid_num_x / 2)) * all_grid._cell_length,
				all_grid.min_pos_y + ((all_grid.Grid_num_y / 2)) * all_grid._cell_length,
				all_grid.min_pos_z + ((all_grid.Grid_num_z / 2)) * all_grid._cell_length);

			//velocity_in_voxel_face_xの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						VEC3_TYPE face_pos_ixiyiz(
							all_grid.min_pos_x +  ix        * all_grid._cell_length,
							all_grid.min_pos_y + (iy + 0.5) * all_grid._cell_length,
							all_grid.min_pos_z + (iz + 0.5) * all_grid._cell_length);
						if ((face_pos_ixiyiz - center_pos_of_grid).norm() < radius) {
							all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
						else {
							all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
					}
				}
			}
			//velocity_in_voxel_face_yの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						VEC3_TYPE face_pos_ixiyiz(
						all_grid.min_pos_x + (ix + 0.5) * all_grid._cell_length,
						all_grid.min_pos_y + iy * all_grid._cell_length,
						all_grid.min_pos_z + (iz + 0.5) * all_grid._cell_length);
						if ((face_pos_ixiyiz - center_pos_of_grid).norm() < radius) {
							all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 10.0;
						}
						else {
							all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
					}
				}
			}
			//velocity_in_voxel_face_zの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z + 1; iz++) {
						VEC3_TYPE face_pos_ixiyiz(
						all_grid.min_pos_x + (ix + 0.5) * all_grid._cell_length,
						all_grid.min_pos_y + (iy + 0.5) * all_grid._cell_length,
						all_grid.min_pos_z + iz * all_grid._cell_length);
						if ((face_pos_ixiyiz - center_pos_of_grid).norm() < radius) {
							all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
						else {
							all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
					}
				}
			}
			// substance_density の初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						VEC3_TYPE face_pos_ixiyiz(
							all_grid.min_pos_x + ix * all_grid._cell_length,
							all_grid.min_pos_y + iy * all_grid._cell_length,
							all_grid.min_pos_z + iz * all_grid._cell_length);

						if ((face_pos_ixiyiz - center_pos_of_grid).norm() < radius) {
							all_grid.substance_density[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 1.0;
						}
						else {
							all_grid.substance_density[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
					}
				}
			}
		}
		else if(physical_const::kInitial_density_shape == "Tylor-Green"){
			//速度の初期化
			//velocity_in_voxel_face_xの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						//考えているセルの面の位置座標
						const MY_FLOAT_TYPE position_x = ix * all_grid._cell_length;
						const MY_FLOAT_TYPE position_y = (iy + 0.5) * all_grid._cell_length;
						const MY_FLOAT_TYPE position_z = (iz + 0.5) * all_grid._cell_length;
						all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
							= sin(2.0 * smoke_simulation::physical_const::kPI * position_x)
							* cos(2.0 * smoke_simulation::physical_const::kPI * position_y)
							* cos(2.0 * smoke_simulation::physical_const::kPI * position_z);
					}
				}
			}
			//velocity_in_voxel_face_yの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						//考えているセルの面の位置座標
						const MY_FLOAT_TYPE position_x = (ix + 0.5) * all_grid._cell_length;
						const MY_FLOAT_TYPE position_y = iy * all_grid._cell_length;
						const MY_FLOAT_TYPE position_z = (iz + 0.5) * all_grid._cell_length;
						all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
							= -cos(2.0 * smoke_simulation::physical_const::kPI * position_x)
							*  sin(2.0 * smoke_simulation::physical_const::kPI * position_y)
							*  cos(2.0 * smoke_simulation::physical_const::kPI * position_z);
					}
				}
			}
			//velocity_in_voxel_face_yの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z + 1; iz++) {
						//考えているセルの面の位置座標
						const MY_FLOAT_TYPE position_x = (ix + 0.5) * all_grid._cell_length;
						const MY_FLOAT_TYPE position_y = (iy + 0.5) * all_grid._cell_length;
						const MY_FLOAT_TYPE position_z = iz * all_grid._cell_length;
						all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
							= 0.0 * cos(2.0 * smoke_simulation::physical_const::kPI * position_x)
							*  cos(2.0 * smoke_simulation::physical_const::kPI * position_y)
							*  sin(2.0 * smoke_simulation::physical_const::kPI * position_z);
					}
				}
			}
			//密度の初期化
			for(int ix = 0; ix < all_grid.Grid_num_x; ++ix){
				for(int iy = 0; iy < all_grid.Grid_num_y; ++iy){
					for(int iz = 0; iz < all_grid.Grid_num_z; ++iz){
						all_grid.substance_density[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.01;
					}
				}
			}
		}

		//固定した速度場を使う場合
		if (fix_velocity && physical_const::kInitial_density_shape != "Tylor-Green") {
/*
			//velocity_in_voxel_face_xの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						if (ix > all_grid.Grid_num_x / 2) {
							all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 1.0;
						}
						else {
							all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 1.0;
						}
					}
				}
			}
			//velocity_in_voxel_face_yの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						if (ix < all_grid.Grid_num_x / 2) {
							all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 1.0;
						}
						else {
							all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 1.0;
						}
					}
				}
			}
			//velocity_in_voxel_face_zの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z + 1; iz++) {
						if (ix < all_grid.Grid_num_x / 2) {
							all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
						else {
							all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
						}
					}
				}
			}
*/
			//velocity_in_voxel_face_xの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
					}
				}
			}
			//velocity_in_voxel_face_yの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
						all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
					}
				}
			}
			//velocity_in_voxel_face_zの初期化
			for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
				for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
					for (int iz = 0; iz < all_grid.Grid_num_z + 1; iz++) {
						all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = 0.0;
					}
				}
			}
		}
		// セル体積の初期化
		for(int ix = 0; ix < all_grid.Grid_num_x; ++ix){
			for(int iy = 0; iy < all_grid.Grid_num_y; ++iy){
				for(int iz = 0; iz < all_grid.Grid_num_z; ++iz){
					all_grid.cell_volume_cell_center[get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] = all_grid._cell_volume;
				}
			}
		}
	}
}//namespace smoke_simulation
