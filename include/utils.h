#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
//#include <opencv2/opencv.hpp>
#include <Eigen/Dense>

#include "physical_const.h"

#include "define_float_type.h"

namespace smoke_simulation {
	inline int get_voxel_face_index_x(int ix, int iy, const int grid_num_x, const int grid_num_y) {
/*
		if (ix < 0) {
			std::cout << "範囲外参照x0" << std::endl;
		}
		if (ix > grid_num_x) {
			std::cout << "範囲外参照x1" << std::endl;
		}
		if (iy < 0) {
			std::cout << "範囲外参照x2" << std::endl;
		}
		if (iy > grid_num_y - 1) {
			std::cout << "範囲外参照x3" << std::endl;
		}
*/
		return grid_num_y * ix + iy;
	}
	inline int get_voxel_face_index_y(int ix, int iy, const int grid_num_x, const int grid_num_y) {
/*
		if (ix < 0) {
			std::cout << "範囲外参照y0" << std::endl;
		}
		if (ix > grid_num_x - 1) {
			std::cout << "範囲外参照y1" << std::endl;
		}
		if (iy < 0) {
			std::cout << "範囲外参照y2" << std::endl;
		}
		if (iy > grid_num_y) {
			std::cout << "範囲外参照y3" << std::endl;
		}
*/
		return (grid_num_y + 1) * ix + iy;
	}
	inline int get_voxel_center_index(int ix, int iy, const int grid_num_x, const int grid_num_y) {
/*
		if (ix < 0) {
			std::cout << "範囲外参照c0" << std::endl;
		}
		if (ix > grid_num_x - 1) {
			std::cout << "範囲外参照c1" << std::endl;
		}
		if (iy < 0) {
			std::cout << "範囲外参照c2" << std::endl;
		}
		if (iy > grid_num_y - 1) {
			std::cout << "範囲外参照c3" << std::endl;
		}
*/
		return grid_num_y * ix + iy;
	}
	inline int get_voxel_face_index_x_3D(int ix, int iy, int iz, const int grid_num_x, const int grid_num_y, const int grid_num_z){
/*
	    if(ix<0){
//			std::cout<<"ix: "<<ix<<", iy: "<<iy<<", iz: "<<iz<<std::endl;
	        std::cout<<"範囲外参照x0"<<std::endl;
	    }
	    if(ix>grid_num_x){
//			std::cout<<"ix: "<<ix<<", iy: "<<iy<<", iz: "<<iz<<std::endl;
	        std::cout<<"範囲外参照x1"<<std::endl;
	    }
	    if(iy<0){
//			std::cout<<"ix: "<<ix<<", iy: "<<iy<<", iz: "<<iz<<std::endl;
	        std::cout<<"範囲外参照x2"<<std::endl;
	    }
	    if(iy>grid_num_y-1){
//			std::cout<<"ix: "<<ix<<", iy: "<<iy<<", iz: "<<iz<<std::endl;
	        std::cout<<"範囲外参照x3"<<std::endl;
	    }
	    if(iz<0){
//			std::cout<<"ix: "<<ix<<", iy: "<<iy<<", iz: "<<iz<<std::endl;
	        std::cout<<"範囲外参照x4"<<std::endl;
	    }
	    if(iz>grid_num_z-1){
//			std::cout<<"ix: "<<ix<<", iy: "<<iy<<", iz: "<<iz<<std::endl;
	        std::cout<<"範囲外参照x5"<<std::endl;
	    }
*/
	    return grid_num_z*(grid_num_y*ix+iy)+iz;
	}

	inline int get_voxel_face_index_y_3D(int ix, int iy, int iz, const int grid_num_x, const int grid_num_y, const int grid_num_z){
/*
	    if(ix<0){
	        std::cout<<"範囲外参照y0"<<std::endl;
	    }
	    if(ix>grid_num_x-1){
	        std::cout<<"範囲外参照y1"<<std::endl;
	    }
	    if(iy<0){
	        std::cout<<"範囲外参照y2"<<std::endl;
	    }
	    if(iy>grid_num_y){
	        std::cout<<"範囲外参照y3"<<std::endl;
	    }
	    if(iz<0){
	        std::cout<<"範囲外参照y4"<<std::endl;
	    }
	    if(iz>grid_num_z-1){
	        std::cout<<"範囲外参照y4"<<std::endl;
	    }
*/
	    return grid_num_z*((grid_num_y+1)*ix+iy)+iz;
	}

	inline int get_voxel_face_index_z_3D(int ix, int iy, int iz, const int grid_num_x, const int grid_num_y, const int grid_num_z){
/*
	    if(ix<0){
	        std::cout<<"範囲外参照z0"<<std::endl;
	    }
	    if(ix>grid_num_x-1){
	        std::cout<<"範囲外参照z1"<<std::endl;
	    }
	    if(iy<0){
	        std::cout<<"範囲外参照z2"<<std::endl;
	    }
	    if(iy>grid_num_y-1){
	        std::cout<<"範囲外参照z3"<<std::endl;
	    }
	    if(iz<0){
	        std::cout<<"範囲外参照z4"<<std::endl;
	    }
	    if(iz>grid_num_z){
	        std::cout<<"範囲外参照z5"<<std::endl;
	    }
*/
	    return (grid_num_z+1)*(grid_num_y*ix+iy)+iz;
	}

	inline int get_voxel_center_index_3D(int ix, int iy, int iz, const int grid_num_x, const int grid_num_y, const int grid_num_z){
/*
	    if(ix<0){
	        std::cout<<"範囲外参照c0"<<std::endl;
	    }
	    if(ix>grid_num_x-1){
	        std::cout<<"範囲外参照c1"<<std::endl;
	    }
	    if(iy<0){
	        std::cout<<"範囲外参照c2"<<std::endl;
	    }
	    if(iy>grid_num_y-1){
	        std::cout<<"範囲外参照c3"<<std::endl;
	    }
	    if(iz<0){
	        std::cout<<"範囲外参照c4"<<std::endl;
	    }
	    if(iz>grid_num_z-1){
	        std::cout<<"範囲外参照c5"<<std::endl;
	    }
*/
	    return grid_num_z*(grid_num_y*ix+iy)+iz;
	}

}//namespace smoke_simulation

#endif //INITIALIZE_GRID_H
