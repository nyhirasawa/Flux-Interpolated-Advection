#include "linear_interpolation_1d.h"

#include "utils.h"

namespace smoke_simulation{
    //////////////////////////////////////////////////
    ////////// 1次元の補間
    //////////////////////////////////////////////////
    MY_FLOAT_TYPE linear_interpolation_1D(
        MY_FLOAT_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_1D &all_grid,
        const MY_FLOAT_TYPE origin_pos,
        const int grid_num_y,
        const MY_FLOAT_TYPE cell_length
    ){
        int advected_index_y = floor((position - origin_pos) / cell_length);

        //補間に使う離散値をセット
        MY_FLOAT_TYPE psi_val[2];

        for(int iy = 0; iy < 2; ++ iy ){
            int index_y = advected_index_y + iy;
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_y < 0) {
                index_y = 0;
            }
            if (index_y > grid_num_y - 1) {
                int exceed = index_y - (grid_num_y - 1);
                //境界の値が外側までずっと続く場合
                index_y = grid_num_y - 1;
                //境界を境に鏡のように値が反射する場合
                //index_y = grid_num_y - 1 - exceed + 1;
            }
            psi_val[iy] = cell_center_values[index_y];
        }

        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position - origin_pos) / cell_length - advected_index_y;
        b1 = 1.0 - b0;
        return b1 * psi_val[0] + b0 * psi_val[1];
    }

    MY_FLOAT_TYPE linear_interpolation_1D_cell_center_values(
        MY_FLOAT_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_1D &all_grid
    ){
        return linear_interpolation_1D(
            position,
            cell_center_values,
            all_grid,
            0.5 * all_grid._cell_length,
            all_grid.Grid_num_y,
            all_grid._cell_length
        );
    }

    //////////////////////////////////////////////////
    ////////// 2次元の補間
    //////////////////////////////////////////////////
    MY_FLOAT_TYPE linear_interpolation_1Dy(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid &all_grid,
        const VEC3_TYPE origin_pos,
        const int grid_num_x,
        const int grid_num_y,
        const MY_FLOAT_TYPE cell_length
    ){
        int advected_index_x = floor((position[0] - (origin_pos[0] - 0.5 * cell_length)) / cell_length);
        int advected_index_y = floor((position[1] - origin_pos[1]) / cell_length);

        //補間に使う離散値をセット
        MY_FLOAT_TYPE psi_val[2];
        for (int iy = 0; iy < 2; ++iy) {
            int index_x = advected_index_x;
            int index_y = advected_index_y + iy;
            // グリッドの外側を参照しようとしたときの処理(x方向)
            if (index_x < 0) {
                index_x = 0;
            }
            if (index_x > grid_num_x - 1) {
                int exceed = index_x - (grid_num_x - 1);
                //境界の値が外側までずっと続く場合
                index_x = grid_num_x - 1;
                //境界を境に鏡のように値が反射する場合
                //index_x = grid_num_x - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_y < 0) {
                index_y = 0;
            }
            if (index_y > grid_num_y - 1) {
                int exceed = index_y - (grid_num_y - 1);
                //境界の値が外側までずっと続く場合
                index_y = grid_num_y - 1;
                //境界を境に鏡のように値が反射する場合
                //index_y = grid_num_y - 1 - exceed + 1;
            }
            psi_val[iy] = cell_center_values[get_voxel_center_index(index_x, index_y, grid_num_x, grid_num_y)];
        }
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1] - origin_pos[1]) / cell_length - advected_index_y;
        b1 = 1.0 - b0;
        return b1 * psi_val[0] + b0 * psi_val[1];
    }

    MY_FLOAT_TYPE linear_interpolation_1Dy_cell_center_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid &all_grid
    ){
        return linear_interpolation_1Dy(
            position,
            cell_center_values,
            all_grid,
            VEC3_TYPE(0.5 * all_grid._cell_length, 0.5 * all_grid._cell_length, 0.0),
            all_grid.Grid_num_x,
            all_grid.Grid_num_y,
            all_grid._cell_length
        );
/*
        int advected_index_x = floor((position[0]) / all_grid._cell_length);
        int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);

        //補間に使う離散値をセット
        MY_FLOAT_TYPE psi_val[2];
        for (int iy = 0; iy < 2; ++iy) {
            int index_x = advected_index_x;
            int index_y = advected_index_y + iy;
            // グリッドの外側を参照しようとしたときの処理(x方向)
            if (index_x < 0) {
                index_x = 0;
            }
            if (index_x > all_grid.Grid_num_x - 1) {
                int exceed = index_x - (all_grid.Grid_num_x - 1);
                //境界の値が外側までずっと続く場合
                index_x = all_grid.Grid_num_x - 1;
                //境界を境に鏡のように値が反射する場合
                //index_x = all_grid.Grid_num_x - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_y < 0) {
                index_y = 0;
            }
            if (index_y > all_grid.Grid_num_y - 1) {
                int exceed = index_y - (all_grid.Grid_num_y - 1);
                //境界の値が外側までずっと続く場合
                index_y = all_grid.Grid_num_y - 1;
                //境界を境に鏡のように値が反射する場合
                //index_y = all_grid.Grid_num_y - 1 - exceed + 1;
            }
            psi_val[iy] = cell_center_values[get_voxel_center_index(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)];
        }
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_y;
        b1 = 1.0 - b0;
        return b1 * psi_val[0] + b0 * psi_val[1];
*/
    }

    MY_FLOAT_TYPE linear_interpolation_1Dy_cell_face_x_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid &all_grid
    ){
        return linear_interpolation_1Dy(
            position,
            cell_face_x_values,
            all_grid,
            VEC3_TYPE(0.0, 0.5 * all_grid._cell_length, 0.0),
            all_grid.Grid_num_x + 1,
            all_grid.Grid_num_y,
            all_grid._cell_length
        );

/*
        int advected_index_x = floor((position[0]) / all_grid._cell_length);
        int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);

        //補間に使う離散値をセット
        MY_FLOAT_TYPE psi_val[2];
        for (int iy = 0; iy < 2; ++iy) {
        	int index_x = advected_index_x;
            int index_y = advected_index_y + iy;
            // グリッドの外側を参照しようとしたときの処理(x方向)
        	if (index_x < 0) {
        		index_x = 0;
        	}
        	if (index_x > all_grid.Grid_num_x) {
        		int exceed = index_x - all_grid.Grid_num_x;
        		//境界の値が外側までずっと続く場合
        		index_x = all_grid.Grid_num_x;
        		//境界を境に鏡のように値が反射する場合
        		//index_x = all_grid.Grid_num_x - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_y < 0) {
            	index_y = 0;
            }
            if (index_y > all_grid.Grid_num_y - 1) {
            	int exceed = index_y - (all_grid.Grid_num_y - 1);
            	//境界の値が外側までずっと続く場合
            	index_y = all_grid.Grid_num_y - 1;
            	//境界を境に鏡のように値が反射する場合
            	//index_y = all_grid.Grid_num_y - 1 - exceed + 1;
            }
            psi_val[iy] = cell_face_x_values[get_voxel_face_index_x(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)];
    	}
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_y;
        b1 = 1.0 - b0;
        return b1 * psi_val[0] + b0 * psi_val[1];
*/
    }

    MY_FLOAT_TYPE linear_interpolation_1Dy_cell_face_y_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_y_values,
        const Grid &all_grid
    ){
        return linear_interpolation_1Dy(
            position,
            cell_face_y_values,
            all_grid,
            VEC3_TYPE(0.5 * all_grid._cell_length, 0.0, 0.0),
            all_grid.Grid_num_x,
            all_grid.Grid_num_y + 1,
            all_grid._cell_length
        );
/*
        int advected_index_x = floor((position[0]) / all_grid._cell_length);
        int advected_index_y = floor((position[1]) / all_grid._cell_length);

    	/////linear interpolation の処理
        //補間に使う離散値をセット
        MY_FLOAT_TYPE psi_val[2];
        for (int iy = 0; iy < 2; ++iy) {
        	int index_x = advected_index_x;
            int index_y = advected_index_y + iy;
            // グリッドの外側を参照しようとしたときの処理(x方向)
            if (index_x < 0) {
            	index_x = 0;
            }
            if (index_x > all_grid.Grid_num_x - 1) {
            	int exceed = index_x - (all_grid.Grid_num_x - 1);
            	//境界の値が外側までずっと続く場合
            	index_x = all_grid.Grid_num_x - 1;
            	//境界を境に鏡のように値が反射する場合
            	//index_x = all_grid.Grid_num_x - 1 - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_y < 0) {
            	index_y = 0;
            }
            if (index_y > all_grid.Grid_num_y) {
            	int exceed = index_y - all_grid.Grid_num_y;
            	//境界の値が外側までずっと続く場合
            	index_y = all_grid.Grid_num_y;
            	//境界を境に鏡のように値が反射する場合
            	//index_y = all_grid.Grid_num_y - exceed + 1;
            }
            psi_val[iy] = cell_face_y_values[get_voxel_face_index_y(index_x, index_y, all_grid.Grid_num_x, all_grid.Grid_num_y)];
    	}

        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1]) / all_grid._cell_length - advected_index_y;
        b1 = 1.0 - b0;
        return b1 * psi_val[0] + b0 * psi_val[1];
*/
    }

    MY_FLOAT_TYPE linear_interpolation_1Dy_psi_velocity_x(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid &all_grid
    ){
        return linear_interpolation_1Dy(
            position,
            cell_face_x_values,
            all_grid,
            VEC3_TYPE(0.0, 0.0, 0.0),
            all_grid.Grid_num_x + 1,
            all_grid.Grid_num_y + 1,
            all_grid._cell_length
        );
    }

    MY_FLOAT_TYPE linear_interpolation_1Dy_psi_velocity_y(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid &all_grid
    ){
        return linear_interpolation_1Dy(
            position,
            cell_face_x_values,
            all_grid,
            VEC3_TYPE(0.5 * all_grid._cell_length, -0.5 * all_grid._cell_length, 0.0),
            all_grid.Grid_num_x,
            all_grid.Grid_num_y + 2,
            all_grid._cell_length
        );
    }

    //////////////////////////////////////////////////
    ////////// 3次元の補間
    //////////////////////////////////////////////////
    MY_FLOAT_TYPE linear_interpolation_1Dy(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_3D &all_grid,
        const VEC3_TYPE origin_pos,
        const int grid_num_x,
        const int grid_num_y,
        const int grid_num_z,
        const MY_FLOAT_TYPE cell_length
    ) {
        int advected_index_x = floor((position[0] - (origin_pos[0] - 0.5 * cell_length)) / cell_length);
        int advected_index_y = floor((position[1] - origin_pos[1]) / cell_length);
        int advected_index_z = floor((position[2] - (origin_pos[2] - 0.5 * cell_length)) / cell_length);

        //補間に使う離散値をセット
        MY_FLOAT_TYPE psi_val[2];
        for (int iy = 0; iy < 2; ++iy) {
            int index_x = advected_index_x;
            int index_y = advected_index_y + iy;
            int index_z = advected_index_z;
            // グリッドの外側を参照しようとしたときの処理(x方向)
            if (index_x < 0) {
                index_x = 0;
            }
            if (index_x > grid_num_x - 1) {
                int exceed = index_x - (grid_num_x - 1);
                //境界の値が外側までずっと続く場合
                index_x = grid_num_x - 1;
                //境界を境に鏡のように値が反射する場合
                //index_x = grid_num_x - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_y < 0) {
                index_y = 0;
            }
            if (index_y > grid_num_y - 1) {
                int exceed = index_y - (grid_num_y - 1);
                //境界の値が外側までずっと続く場合
                index_y = grid_num_y - 1;
                //境界を境に鏡のように値が反射する場合
                //index_y = grid_num_y - 1 - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_z < 0) {
                index_z = 0;
            }
            if (index_z > grid_num_z - 1) {
                int exceed = index_z - (grid_num_z - 1);
                //境界の値が外側までずっと続く場合
                index_z = grid_num_z - 1;
                //境界を境に鏡のように値が反射する場合
                //index_z = grid_num_z - 1 - exceed + 1;
            }
            psi_val[iy] = cell_center_values[get_voxel_center_index_3D(index_x, index_y, index_z, grid_num_x, grid_num_y, grid_num_z)];
        }
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1] - origin_pos[1]) / cell_length - advected_index_y;
        b1 = 1.0 - b0;
        return b1 * psi_val[0] + b0 * psi_val[1];
    }

    MY_FLOAT_TYPE linear_interpolation_1Dy_cell_center_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_center_values,
        const Grid_3D &all_grid
    ){
        linear_interpolation_1Dy(
            position,
            cell_center_values,
            all_grid,
            VEC3_TYPE(0.5 * all_grid._cell_length, 0.5 * all_grid._cell_length, 0.5 * all_grid._cell_length),
            all_grid.Grid_num_x,
            all_grid.Grid_num_y,
            all_grid.Grid_num_z,
            all_grid._cell_length
        );
/*
        int advected_index_x = floor((position[0]) / all_grid._cell_length);
        int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        int advected_index_z = floor((position[2]) / all_grid._cell_length);

        //補間に使う離散値をセット
        MY_FLOAT_TYPE psi_val[2];
        for (int iy = 0; iy < 2; ++iy) {
            int index_x = advected_index_x;
            int index_y = advected_index_y + iy;
            int index_z = advected_index_z;
            // グリッドの外側を参照しようとしたときの処理(x方向)
            if (index_x < 0) {
                index_x = 0;
            }
            if (index_x > all_grid.Grid_num_x - 1) {
                int exceed = index_x - (all_grid.Grid_num_x - 1);
                //境界の値が外側までずっと続く場合
                index_x = all_grid.Grid_num_x - 1;
                //境界を境に鏡のように値が反射する場合
                //index_x = all_grid.Grid_num_x - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_y < 0) {
                index_y = 0;
            }
            if (index_y > all_grid.Grid_num_y - 1) {
                int exceed = index_y - (all_grid.Grid_num_y - 1);
                //境界の値が外側までずっと続く場合
                index_y = all_grid.Grid_num_y - 1;
                //境界を境に鏡のように値が反射する場合
                //index_y = all_grid.Grid_num_y - 1 - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_z < 0) {
                index_z = 0;
            }
            if (index_z > all_grid.Grid_num_z - 1) {
                int exceed = index_z - (all_grid.Grid_num_z - 1);
                //境界の値が外側までずっと続く場合
                index_z = all_grid.Grid_num_z - 1;
                //境界を境に鏡のように値が反射する場合
                //index_z = all_grid.Grid_num_z - 1 - exceed + 1;
            }
            psi_val[iy] = cell_center_values[get_voxel_center_index_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
        }
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_y;
        b1 = 1.0 - b0;
        return b1 * psi_val[0] + b0 * psi_val[1];
*/
    }

    MY_FLOAT_TYPE linear_interpolation_1Dy_cell_face_x_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
        const Grid_3D &all_grid
    ){
        return linear_interpolation_1Dy(
            position,
            cell_face_x_values,
            all_grid,
            VEC3_TYPE(0.0, 0.5 * all_grid._cell_length, 0.5 * all_grid._cell_length),
            all_grid.Grid_num_x + 1,
            all_grid.Grid_num_y,
            all_grid.Grid_num_z,
            all_grid._cell_length
        );
/*
        int advected_index_x = floor((position[0]) / all_grid._cell_length);
        int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        int advected_index_z = floor((position[2]) / all_grid._cell_length);

        //補間に使う離散値をセット
        MY_FLOAT_TYPE psi_val[2];
        for (int iy = 0; iy < 2; ++iy) {
        	int index_x = advected_index_x;
            int index_y = advected_index_y + iy;
            int index_z = advected_index_z;
            // グリッドの外側を参照しようとしたときの処理(x方向)
        	if (index_x < 0) {
        		index_x = 0;
        	}
        	if (index_x > all_grid.Grid_num_x) {
        		int exceed = index_x - all_grid.Grid_num_x;
        		//境界の値が外側までずっと続く場合
        		index_x = all_grid.Grid_num_x;
        		//境界を境に鏡のように値が反射する場合
        		//index_x = all_grid.Grid_num_x - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_y < 0) {
            	index_y = 0;
            }
            if (index_y > all_grid.Grid_num_y - 1) {
            	int exceed = index_y - (all_grid.Grid_num_y - 1);
            	//境界の値が外側までずっと続く場合
            	index_y = all_grid.Grid_num_y - 1;
            	//境界を境に鏡のように値が反射する場合
            	//index_y = all_grid.Grid_num_y - 1 - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_z < 0) {
            	index_z = 0;
            }
            if (index_z > all_grid.Grid_num_z - 1) {
            	int exceed = index_z - (all_grid.Grid_num_z - 1);
            	//境界の値が外側までずっと続く場合
            	index_z = all_grid.Grid_num_z - 1;
            	//境界を境に鏡のように値が反射する場合
            	//index_z = all_grid.Grid_num_z - 1 - exceed + 1;
            }
            psi_val[iy] = cell_face_x_values[get_voxel_face_index_x_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
    	}
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_y;
        b1 = 1.0 - b0;
        return b1 * psi_val[0] + b0 * psi_val[1];
*/
    }

    MY_FLOAT_TYPE linear_interpolation_1Dy_cell_face_y_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_y_values,
        const Grid_3D &all_grid
    ){
        return linear_interpolation_1Dy(
            position,
            cell_face_y_values,
            all_grid,
            VEC3_TYPE(0.5 * all_grid._cell_length, 0.0, 0.5 * all_grid._cell_length),
            all_grid.Grid_num_x,
            all_grid.Grid_num_y + 1,
            all_grid.Grid_num_z,
            all_grid._cell_length
        );
/*
        int advected_index_x = floor((position[0]) / all_grid._cell_length);
        int advected_index_y = floor((position[1]) / all_grid._cell_length);
        int advected_index_z = floor((position[2]) / all_grid._cell_length);

    	/////linear interpolation の処理
        //補間に使う離散値をセット
        MY_FLOAT_TYPE psi_val[2];
        for (int iy = 0; iy < 2; ++iy) {
        	int index_x = advected_index_x;
            int index_y = advected_index_y + iy;
            int index_z = advected_index_z;
            // グリッドの外側を参照しようとしたときの処理(x方向)
            if (index_x < 0) {
            	index_x = 0;
            }
            if (index_x > all_grid.Grid_num_x - 1) {
            	int exceed = index_x - (all_grid.Grid_num_x - 1);
            	//境界の値が外側までずっと続く場合
            	index_x = all_grid.Grid_num_x - 1;
            	//境界を境に鏡のように値が反射する場合
            	//index_x = all_grid.Grid_num_x - 1 - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_y < 0) {
            	index_y = 0;
            }
            if (index_y > all_grid.Grid_num_y) {
            	int exceed = index_y - all_grid.Grid_num_y;
            	//境界の値が外側までずっと続く場合
            	index_y = all_grid.Grid_num_y;
            	//境界を境に鏡のように値が反射する場合
            	//index_y = all_grid.Grid_num_y - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_z < 0) {
            	index_z = 0;
            }
            if (index_z > all_grid.Grid_num_z - 1) {
            	int exceed = index_z - (all_grid.Grid_num_z - 1);
            	//境界の値が外側までずっと続く場合
            	index_z = all_grid.Grid_num_z - 1;
            	//境界を境に鏡のように値が反射する場合
            	//index_z = all_grid.Grid_num_z - 1 - exceed + 1;
            }
            psi_val[iy] = cell_face_y_values[get_voxel_face_index_y_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
    	}

        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1]) / all_grid._cell_length - advected_index_y;
        b1 = 1.0 - b0;
        return b1 * psi_val[0] + b0 * psi_val[1];
*/
    }

    MY_FLOAT_TYPE linear_interpolation_1Dy_cell_face_z_values(
        VEC3_TYPE position,
        const std::vector<MY_FLOAT_TYPE> &cell_face_z_values,
        const Grid_3D &all_grid
    ){
        return linear_interpolation_1Dy(
            position,
            cell_face_z_values,
            all_grid,
            VEC3_TYPE(0.5 * all_grid._cell_length, 0.5 * all_grid._cell_length, 0.0),
            all_grid.Grid_num_x,
            all_grid.Grid_num_y,
            all_grid.Grid_num_z + 1,
            all_grid._cell_length
        );
/*
        int advected_index_x = floor((position[0]) / all_grid._cell_length);
        int advected_index_y = floor((position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length);
        int advected_index_z = floor((position[2]) / all_grid._cell_length);
        //補間に使う離散値をセット
        MY_FLOAT_TYPE psi_val[2];
        for (int iy = 0; iy < 2; ++iy) {
        	int index_x = advected_index_x;
            int index_y = advected_index_y + iy;
            int index_z = advected_index_z;
            // グリッドの外側を参照しようとしたときの処理(x方向)
            if (index_x < 0) {
            	index_x = 0;
            }
            if (index_x > all_grid.Grid_num_x - 1) {
            	int exceed = index_x - (all_grid.Grid_num_x - 1);
            	//境界の値が外側までずっと続く場合
            	index_x = all_grid.Grid_num_x - 1;
            	//境界を境に鏡のように値が反射する場合
            	//index_x = all_grid.Grid_num_x - 1 - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_y < 0) {
            	index_y = 0;
            }
            if (index_y > all_grid.Grid_num_y - 1) {
            	int exceed = index_y - (all_grid.Grid_num_y - 1);
            	//境界の値が外側までずっと続く場合
            	index_y = all_grid.Grid_num_y - 1;
            	//境界を境に鏡のように値が反射する場合
            	//index_y = all_grid.Grid_num_y - 1 - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_z < 0) {
            	index_z = 0;
            }
            if (index_z > all_grid.Grid_num_z) {
            	int exceed = index_z - all_grid.Grid_num_z;
            	//境界の値が外側までずっと続く場合
            	index_z = all_grid.Grid_num_z;
            	//境界を境に鏡のように値が反射する場合
            	//index_z = all_grid.Grid_num_z - exceed + 1;
            }
            psi_val[iy] = cell_face_z_values[get_voxel_face_index_z_3D(index_x, index_y, index_z, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
        }
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1] - 0.5 * all_grid._cell_length) / all_grid._cell_length - advected_index_y;
        b1 = 1.0 - b0;
        return b1 * psi_val[0] + b0 * psi_val[1];
*/
}
MY_FLOAT_TYPE linear_interpolation_1Dy_psi_velocity_x(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
    const Grid_3D &all_grid
) {
    return linear_interpolation_1Dy(
        position,
        cell_face_x_values,
        all_grid,
        VEC3_TYPE(0.0, 0.0, 0.5 * all_grid._cell_length),
        all_grid.Grid_num_x + 1,
        all_grid.Grid_num_y + 1,
        all_grid.Grid_num_z,
        all_grid._cell_length
    );
}
MY_FLOAT_TYPE linear_interpolation_1Dy_psi_velocity_y(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
    const Grid_3D &all_grid
) {
    return linear_interpolation_1Dy(
        position,
        cell_face_x_values,
        all_grid,
        VEC3_TYPE(0.5 * all_grid._cell_length, -0.5 * all_grid._cell_length, 0.5 * all_grid._cell_length),
        all_grid.Grid_num_x,
        all_grid.Grid_num_y + 2,
        all_grid.Grid_num_z,
        all_grid._cell_length
    );
}
MY_FLOAT_TYPE linear_interpolation_1Dy_psi_velocity_z(
    VEC3_TYPE position,
    const std::vector<MY_FLOAT_TYPE> &cell_face_x_values,
    const Grid_3D &all_grid
) {
    return linear_interpolation_1Dy(
        position,
        cell_face_x_values,
        all_grid,
        VEC3_TYPE(0.5 * all_grid._cell_length, 0.0, 0.0),
        all_grid.Grid_num_x,
        all_grid.Grid_num_y + 1,
        all_grid.Grid_num_z + 1,
        all_grid._cell_length
    );
}


}//namespace smoke_simulation
