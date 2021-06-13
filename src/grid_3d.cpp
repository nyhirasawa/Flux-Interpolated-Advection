#include "grid_3d.h"

#include "utils.h"
#include "define_float_type.h"
#include "calc_psi.h"
#include "linear_interpolation_1d.h"
#include "linear_interpolation_3d.h"
#include "WENO4_interpolation_3d.h"
#include "WENO6_interpolation_1d.h"
#include "WENO6_interpolation_3d.h"

namespace smoke_simulation{
    Grid_3D::Grid_3D(int nx, int ny, int nz, MY_FLOAT_TYPE cell_length)
    : Grid_num_x(nx),
      Grid_num_y(ny),
      Grid_num_z(nz),
      min_pos_x(0.0),
      min_pos_y(0.0),
      min_pos_z(0.0),
      _cell_length(cell_length),
      _cell_face_area(cell_length*cell_length),
      _cell_volume(cell_length*cell_length*cell_length),
      total_loss_of_momentum_x(0.0),
      total_loss_of_momentum_y(0.0),
      total_loss_of_momentum_z(0.0){
        //substance_densityのメモリ確保
        substance_density.resize(nx*ny*nz);

        //velocityのメモリ確保
        velocity_in_voxel_face_x.resize((nx+1) * ny * nz);
        velocity_in_voxel_face_y.resize(nx * (ny+1) * nz);
        velocity_in_voxel_face_z.resize(nx * ny * (nz+1));

        psi_substance_density_cell_face_x.resize((nx + 1) * ny * nz);
        psi_substance_density_cell_face_y.resize(nx * (ny + 1) * nz);
        psi_substance_density_cell_face_z.resize(nx * ny * (nz + 1));

        psi_velocity_x_cell_face.resize((nx + 1) * (ny + 1) * nz);
        psi_velocity_y_cell_face.resize(nx * (ny + 2) * nz);
        psi_velocity_z_cell_face.resize(nx * (ny + 1) * (nz + 1));

        //pressureのメモリ確保
        pressure.resize(nx*ny*nz);

        //セル体積のメモリ確保
        cell_volume_cell_center.resize(nx * ny * nz);
    }

    Grid_3D::~Grid_3D(){
    }

    // position での psi_density の y 成分の値を周辺の cell face での値から linear interpolation して求める。
    // (calc_interpolated_velocity のヘルプ関数)
    MY_FLOAT_TYPE Grid_3D::calc_psi_density_y_by_1Dy_linear_integral_cell_face_y(VEC3_TYPE position) const{
        int advected_index_x = floor((position[0]) / _cell_length);
        int advected_index_y = floor((position[1]) / _cell_length);
        int advected_index_z = floor((position[2]) / _cell_length);

    	///// 積分の処理
        ////////////////////////////////////////
        ////psiからの寄与
        ////////////////////////////////////////
        // 下端の面でのpsiの値
//        MY_FLOAT_TYPE psi_under
//            = psi_substance_density_cell_face_y[
//                get_voxel_face_index_y_3D(
//                    advected_index_x, advected_index_y, advected_index_z,
//                    Grid_num_x, Grid_num_y, Grid_num_z)
//                ];
/*
        VEC3_TYPE psi_under_interpolated
            = calc_interpolated_psi_density(VEC3_TYPE((advected_index_x + 0.5) * _cell_length, (advected_index_y) * _cell_length, (advected_index_z + 0.5) * _cell_length), "linear");
        MY_FLOAT_TYPE psi_under = psi_under_interpolated[1];
*/


        ////////////////////////////////////////
        ////densityからの寄与
        ////////////////////////////////////////

        //// position を含むセルでの積分を計算
        MY_FLOAT_TYPE psi_under=0.0;
        //リーマン和で積分を計算
/*
        const int num_split = 100;
        MY_FLOAT_TYPE h = (position[1] - advected_index_y * _cell_length) / num_split;
        for(int i = 0; i < num_split; ++i){
            psi_under
                += linear_interpolation_1Dy_cell_center_values(
                    VEC3_TYPE(
                        position[0],
                        advected_index_y * _cell_length + (i+0.5)*h,
                        position[2]
                    ),
                    substance_density,
                    *this
                ) * h;
        }
*/
        ////積分を厳密に計算
        // positionが下端の面からどれだけ離れた位置にあるかを[0, 1]で表した量
        MY_FLOAT_TYPE over_under_face = ((position[1]) / _cell_length) - advected_index_y;
        // positionがセルの中心より下にある場合
        if(over_under_face < 0.5){
            // positionがセルの中心より下にある場合は1つの求積点でセルでの積分ができる
            psi_under
                += linear_interpolation_1Dy_cell_center_values(
                    VEC3_TYPE(
                        position[0],
                        (advected_index_y + over_under_face / 2.0) * _cell_length,
                        position[2]
                    ),
                    substance_density,
                    *this
                )
                * over_under_face
                * _cell_length;
//            psi_under
//                += linear_interpolation_3D_cell_center_values(
//                    VEC3_TYPE(
//                        position[0],
//                        (advected_index_y + over_under_face / 2.0) * _cell_length,
//                        position[2]
//                    ),
//                    substance_density,
//                    *this
//                )
//                * over_under_face
//                * _cell_length;
        }
        // positionがセルの中心より上にある場合
        else{
            //// positionがセルの中心より上にある場合はセル中心より下の1点と上の1点の
            //// 2つの求積点でセルでの積分ができる
            //セル中心より下の1点
            psi_under
                += linear_interpolation_1Dy_cell_center_values(
                    VEC3_TYPE(
                        position[0],
                        (advected_index_y + (0.5 / 2.0)) * _cell_length,
                        position[2]
                    ),
                    substance_density,
                    *this
                )
                * 0.5
                * _cell_length;
//            psi_under
//                += linear_interpolation_3D_cell_center_values(
//                    VEC3_TYPE(
//                        position[0],
//                        (advected_index_y + (0.5 / 2.0)) * _cell_length,
//                        position[2]
//                    ),
//                    substance_density,
//                    *this
//                )
//                * 0.5
//                * _cell_length;
            //セル中心より上の1点
            psi_under
                += linear_interpolation_1Dy_cell_center_values(
                    VEC3_TYPE(
                        position[0],
                        (advected_index_y + ((over_under_face - 0.5)) / 2.0 + 0.5) * _cell_length,
                        position[2]
                    ),
                    substance_density,
                    *this
                )
                * (over_under_face - 0.5)
                * _cell_length;
//            psi_under
//                += linear_interpolation_3D_cell_center_values(
//                    VEC3_TYPE(
//                        position[0],
//                        (advected_index_y + ((over_under_face - 0.5)) / 2.0 + 0.5) * _cell_length,
//                        position[2]
//                    ),
//                    substance_density,
//                    *this
//                )
//                * (over_under_face - 0.5)
//                * _cell_length;
        }
        return psi_under;
    }
    //position での psi_density の y 成分の値を周辺の cell face での値から linear interpolation して求める。
    //(calc_interpolated_velocity のヘルプ関数)
    MY_FLOAT_TYPE Grid_3D::calc_psi_density_y_by_linear_integral_cell_face_y(VEC3_TYPE position) const{
//        int advected_index_x = floor((position[0]) / _cell_length);
        int advected_index_y = floor((position[1]) / _cell_length);
//        int advected_index_z = floor((position[2]) / _cell_length);

    	///// 積分の処理
        ////////////////////////////////////////
        ////psiからの寄与
        ////////////////////////////////////////
        // 下端の面でのpsiの値
//        MY_FLOAT_TYPE psi_under
//            = psi_substance_density_cell_face_y[
//                get_voxel_face_index_y_3D(
//                    advected_index_x, advected_index_y, advected_index_z,
//                    Grid_num_x, Grid_num_y, Grid_num_z)
//                ];
        VEC3_TYPE psi_under_interpolated
            = calc_interpolated_psi_density(
                VEC3_TYPE(
                    position[0],
                    (advected_index_y) * _cell_length,
                    position[2]
                ),
                "linear");
        MY_FLOAT_TYPE psi_under = psi_under_interpolated[1];



        ////////////////////////////////////////
        ////densityからの寄与
        ////////////////////////////////////////

        //// position を含むセルでの積分を計算
//        MY_FLOAT_TYPE psi_under=0.0;
        //リーマン和で積分を計算
/*
        const int num_split = 100;
        MY_FLOAT_TYPE h = (position[1] - advected_index_y * _cell_length) / num_split;
        for(int i = 0; i < num_split; ++i){
            psi_under
                += linear_interpolation_1Dy_cell_center_values(
                    VEC3_TYPE(
                        position[0],
                        advected_index_y * _cell_length + (i+0.5)*h,
                        position[2]
                    ),
                    substance_density,
                    *this
                ) * h;
        }
*/
        ////積分を厳密に計算
        // positionが下端の面からどれだけ離れた位置にあるかを[0, 1]で表した量
        MY_FLOAT_TYPE over_under_face = ((position[1]) / _cell_length) - advected_index_y;
        // positionがセルの中心より下にある場合
        if(over_under_face < 0.5){
            // positionがセルの中心より下にある場合は1つの求積点でセルでの積分ができる
            psi_under
                += linear_interpolation_3D_cell_center_values(
                    VEC3_TYPE(
                        position[0],
                        (advected_index_y + over_under_face / 2.0) * _cell_length,
                        position[2]
                    ),
                    substance_density,
                    *this
                )
                * over_under_face
                * _cell_length;
        }
        // positionがセルの中心より上にある場合
        else{
            //// positionがセルの中心より上にある場合はセル中心より下の1点と上の1点の
            //// 2つの求積点でセルでの積分ができる
            //セル中心より下の1点
            psi_under
                += linear_interpolation_3D_cell_center_values(
                    VEC3_TYPE(
                        position[0],
                        (advected_index_y + (0.5 / 2.0)) * _cell_length,
                        position[2]
                    ),
                    substance_density,
                    *this
                )
                * 0.5
                * _cell_length;
            //セル中心より上の1点
            psi_under
                += linear_interpolation_3D_cell_center_values(
                    VEC3_TYPE(
                        position[0],
                        (advected_index_y + ((over_under_face - 0.5)) / 2.0 + 0.5) * _cell_length,
                        position[2]
                    ),
                    substance_density,
                    *this
                )
                * (over_under_face - 0.5)
                * _cell_length;
        }
        return psi_under;

/*
        int advected_index_x = floor((position[0] - 0.5 * _cell_length) / _cell_length);
        int advected_index_y = floor((position[1]) / _cell_length);
        int advected_index_z = floor((position[2] - 0.5 * _cell_length) / _cell_length);

    	/////linear interpolation の処理
    	//y方向の積分の結果を格納するための変数
    	MY_FLOAT_TYPE y_interpolation_values_of_psi_substance_density[2][2];
    	//x軸方向に走査するループ。ループ内の処理ではy軸方向に関するinterpolationを行う。
    	for (int ix = 0; ix < 2; ix++) {
            for(int iz = 0; iz < 2; ++iz){
            	int index_x = advected_index_x + ix;
                int index_z = advected_index_z + iz;
                // グリッドの外側を参照しようとしたときの処理(x方向)
            	if (index_x < 0) {
            		index_x = 0;
            	}
            	if (index_x > Grid_num_x - 1) {
            		int exceed = index_x - (Grid_num_x - 1);
            		//境界の値が外側までずっと続く場合
            		index_x = Grid_num_x - 1;
            		//境界を境に鏡のように値が反射する場合
            		//index_x = Grid_num_x - 1 - exceed + 1;
            	}
                // グリッドの外側を参照しようとしたときの処理(y方向)
                if (index_z < 0) {
            		index_z = 0;
            	}
            	if (index_z > Grid_num_z - 1) {
            		int exceed = index_z - (Grid_num_z - 1);
            		//境界の値が外側までずっと続く場合
            		index_z = Grid_num_z - 1;
            		//境界を境に鏡のように値が反射する場合
            		//index_z = Grid_num_z - 1 - exceed + 1;
            	}
    			y_interpolation_values_of_psi_substance_density[ix][iz]
                    = calc_psi_density_y_by_1Dy_linear_integral_cell_face_y(
                        VEC3_TYPE(
                            (index_x + 0.5) * _cell_length,
                            position[1],
                            (index_z + 0.5) * _cell_length
                        )
                    );
            }
    	}

    	////ここからz方向の補間
        //yz方向のinterpolationの結果を格納するための変数
        MY_FLOAT_TYPE yz_interpolation_values_of_psi_substance_density[2];
        for(int ix = 0; ix < 2; ++ix){
            //interpolation の係数
            MY_FLOAT_TYPE c0, c1;
            c0 = (position[2] - 0.5 * _cell_length) / _cell_length - advected_index_z;
            c1 = 1.0 - c0;
    		yz_interpolation_values_of_psi_substance_density[ix]
                = c1 * y_interpolation_values_of_psi_substance_density[ix][0] + c0 * y_interpolation_values_of_psi_substance_density[ix][1];
        }

        ////ここからx方向の補間
    	//interpolation の係数
    	MY_FLOAT_TYPE a0, a1;
    	a0 = (position[0] - 0.5 * _cell_length) / _cell_length - advected_index_x;
    	a1 = 1.0 - a0;
    	//結果
    	return a1 * yz_interpolation_values_of_psi_substance_density[0] + a0 * yz_interpolation_values_of_psi_substance_density[1];
*/
    }

    //position での psi_density をinterpolationによって計算
    VEC3_TYPE Grid_3D::calc_interpolated_psi_density(VEC3_TYPE position, const std::string interpolation_method) const {
        VEC3_TYPE interpolated_psi_density;
        if (interpolation_method == "linear") {
            interpolated_psi_density[0]
                = linear_interpolation_3D_cell_face_x_values(position, psi_substance_density_cell_face_x, *this);
            interpolated_psi_density[1]
                = linear_interpolation_3D_cell_face_y_values(position, psi_substance_density_cell_face_y, *this);
            interpolated_psi_density[2]
                = linear_interpolation_3D_cell_face_z_values(position, psi_substance_density_cell_face_z, *this);
        }
        else if (interpolation_method == "WENO4") {
            interpolated_psi_density[0]
                = WENO4_interpolation_3D_cell_face_x_values(position, psi_substance_density_cell_face_x, *this);
            interpolated_psi_density[1]
                = WENO4_interpolation_3D_cell_face_y_values(position, psi_substance_density_cell_face_y, *this);
            interpolated_psi_density[2]
                = WENO4_interpolation_3D_cell_face_z_values(position, psi_substance_density_cell_face_z, *this);
        }
        else if (interpolation_method == "WENO6") {
            interpolated_psi_density[0]
                = WENO6_interpolation_3D_cell_face_x_values(position, psi_substance_density_cell_face_x, *this);
            interpolated_psi_density[1]
                = WENO6_interpolation_3D_cell_face_y_values(position, psi_substance_density_cell_face_y, *this);
            interpolated_psi_density[2]
                = WENO6_interpolation_3D_cell_face_z_values(position, psi_substance_density_cell_face_z, *this);
        }
        else if (interpolation_method == "1Dy_linear") {
            if(physical_const::kPsi_definition == "y"){
                interpolated_psi_density[0] = 0.0;
            }
            else{
                interpolated_psi_density[0]
                    = linear_interpolation_1Dy_cell_face_x_values(position, psi_substance_density_cell_face_x, *this);
            }

            interpolated_psi_density[1]
                = linear_interpolation_1Dy_cell_face_y_values(position, psi_substance_density_cell_face_y, *this);

            if(physical_const::kPsi_definition == "y"){
                interpolated_psi_density[2] = 0.0;
            }
            else{
                interpolated_psi_density[2]
                    = linear_interpolation_1Dy_cell_face_z_values(position, psi_substance_density_cell_face_z, *this);
            }
        }
        else if (interpolation_method == "1Dy_WENO6") {
            interpolated_psi_density[0]
                = WENO6_interpolation_1Dy_cell_face_x_values(position, psi_substance_density_cell_face_x, *this);
            interpolated_psi_density[1]
                = WENO6_interpolation_1Dy_cell_face_y_values(position, psi_substance_density_cell_face_y, *this);
            interpolated_psi_density[2]
                = WENO6_interpolation_1Dy_cell_face_z_values(position, psi_substance_density_cell_face_z, *this);
        }
        else if(interpolation_method == "1Dy_linear_integral"){
            if(physical_const::kPsi_definition == "y"){
                interpolated_psi_density[0] = 0.0;
            }
            else{
//                interpolated_psi_density[0] = calc_psi_density_x_by_1Dy_linear_integral_cell_face_x(position);
            }

            interpolated_psi_density[1] = calc_psi_density_y_by_1Dy_linear_integral_cell_face_y(position);

            if(physical_const::kPsi_definition == "y"){
                interpolated_psi_density[2] = 0.0;
            }
            else{
//                interpolated_psi_density[2] = calc_psi_density_z_by_1Dy_linear_integral_cell_face_z(position);
            }
        }
        else if(interpolation_method == "linear_integral"){
            if(physical_const::kPsi_definition == "y"){
                interpolated_psi_density[0] = 0.0;
            }
            else{
//                interpolated_psi_density[0] = calc_psi_density_x_by_linear_integral_cell_face_x(position);
            }

            interpolated_psi_density[1] = calc_psi_density_y_by_linear_integral_cell_face_y(position);

            if(physical_const::kPsi_definition == "y"){
                interpolated_psi_density[2] = 0.0;
            }
            else{
//                interpolated_psi_density[2] = calc_psi_density_z_by_linear_integral_cell_face_z(position);
            }
        }
        return interpolated_psi_density;
    }

    //position での psi_density をinterpolationによって計算
    VEC3_TYPE Grid_3D::calc_interpolated_psi_velocity(const int dim, VEC3_TYPE position, const std::string interpolation_method) const {
        if(dim == 0) {
            return calc_interpolated_psi_velocity_x(position, interpolation_method);
        }
        else if(dim == 1) {
            return calc_interpolated_psi_velocity_y(position, interpolation_method);
        }
        else if(dim == 2) {
            return calc_interpolated_psi_velocity_z(position, interpolation_method);
        }
    }

    //position での psi_velocity_x を interpolation によって計算
    VEC3_TYPE Grid_3D::calc_interpolated_psi_velocity_x(VEC3_TYPE position, const std::string interpolation_method) const {
        VEC3_TYPE interpolated_psi_velocity_x;
        if (interpolation_method == "linear") {
            interpolated_psi_velocity_x[0] = 0.0;
            interpolated_psi_velocity_x[1] = linear_interpolation_3D_psi_velocity_x(position, psi_velocity_x_cell_face, *this);
            interpolated_psi_velocity_x[2] = 0.0;
        }
        else if (interpolation_method == "WENO6") {
//            interpolated_psi_velocity_x[0] = 0.0;
//            interpolated_psi_velocity_x[1] = WENO6_interpolation_3D_psi_velocity_x(position, psi_velocity_x_cell_face, *this);
//            interpolated_psi_velocity_x[2] = 0.0;
        }
        else if (interpolation_method == "1Dy_linear") {
            interpolated_psi_velocity_x[0] = 0.0;
            interpolated_psi_velocity_x[1] = linear_interpolation_1Dy_psi_velocity_x(position, psi_velocity_x_cell_face, *this);
            interpolated_psi_velocity_x[2] = 0.0;
        }
        else if (interpolation_method == "1Dy_WENO6") {
//            interpolated_psi_velocity_x[0] = 0.0;
//            interpolated_psi_velocity_x[1] = WENO6_interpolation_1Dy_psi_velocity_x(position, psi_velocity_x_cell_face, *this);
//            interpolated_psi_velocity_x[2] = 0.0;
        }
        return interpolated_psi_velocity_x;
    }

    //position での psi_velocity_x を interpolation によって計算
    VEC3_TYPE Grid_3D::calc_interpolated_psi_velocity_y(VEC3_TYPE position, const std::string interpolation_method) const {
        VEC3_TYPE interpolated_psi_velocity_y;
        if (interpolation_method == "linear") {
            interpolated_psi_velocity_y[0] = 0.0;
            interpolated_psi_velocity_y[1] = linear_interpolation_3D_psi_velocity_y(position, psi_velocity_y_cell_face, *this);
            interpolated_psi_velocity_y[2] = 0.0;
        }
        else if (interpolation_method == "WENO6") {
//            interpolated_psi_velocity_y[0] = 0.0;
//            interpolated_psi_velocity_y[1] = WENO6_interpolation_3D_psi_velocity_y(position, psi_velocity_y_cell_face, *this);
//            interpolated_psi_velocity_y[2] = 0.0;
        }
        else if (interpolation_method == "1Dy_linear") {
//            interpolated_psi_velocity_y[0] = 0.0;
//            interpolated_psi_velocity_y[1] = linear_interpolation_1Dy_psi_velocity_y(position, psi_velocity_y_cell_face, *this);
//            interpolated_psi_velocity_y[2] = 0.0;
        }
        else if (interpolation_method == "1Dy_WENO6") {
//            interpolated_psi_velocity_y[0] = 0.0;
//            interpolated_psi_velocity_y[1] = WENO6_interpolation_1Dy_psi_velocity_y(position, psi_velocity_y_cell_face, *this);
//            interpolated_psi_velocity_y[2] = 0.0;
        }
        return interpolated_psi_velocity_y;
    }

    //position での psi_velocity_x を interpolation によって計算
    VEC3_TYPE Grid_3D::calc_interpolated_psi_velocity_z(VEC3_TYPE position, const std::string interpolation_method) const {
        VEC3_TYPE interpolated_psi_velocity_z;
        if (interpolation_method == "linear") {
            interpolated_psi_velocity_z[0] = 0.0;
            interpolated_psi_velocity_z[1] = linear_interpolation_3D_psi_velocity_z(position, psi_velocity_z_cell_face, *this);
            interpolated_psi_velocity_z[2] = 0.0;
        }
        else if (interpolation_method == "WENO6") {
//            interpolated_psi_velocity_z[0] = 0.0;
//            interpolated_psi_velocity_z[1] = WENO6_interpolation_3D_psi_velocity_z(position, psi_velocity_z_cell_face, *this);
//            interpolated_psi_velocity_z[2] = 0.0;
        }
        else if (interpolation_method == "1Dy_linear") {
//            interpolated_psi_velocity_z[0] = 0.0;
//            interpolated_psi_velocity_z[1] = linear_interpolation_1Dy_psi_velocity_z(position, psi_velocity_z_cell_face, *this);
//            interpolated_psi_velocity_z[2] = 0.0;
        }
        else if (interpolation_method == "1Dy_WENO6") {
//            interpolated_psi_velocity_z[0] = 0.0;
//            interpolated_psi_velocity_z[1] = WENO6_interpolation_1Dy_psi_velocity_z(position, psi_velocity_z_cell_face, *this);
//            interpolated_psi_velocity_z[2] = 0.0;
        }
        return interpolated_psi_velocity_z;
    }

    //position での速度をinterpolationによって計算
    VEC3_TYPE Grid_3D::calc_interpolated_velocity(VEC3_TYPE position, const std::string interpolation_method) const {
        VEC3_TYPE interpolated_velocity;
        if (interpolation_method == "linear") {
            interpolated_velocity[0] = linear_interpolation_3D_cell_face_x_values(position, velocity_in_voxel_face_x, *this);
            interpolated_velocity[1] = linear_interpolation_3D_cell_face_y_values(position, velocity_in_voxel_face_y, *this);
            interpolated_velocity[2] = linear_interpolation_3D_cell_face_z_values(position, velocity_in_voxel_face_z, *this);
        }
        else if (interpolation_method == "WENO4") {
            interpolated_velocity[0] = WENO4_interpolation_3D_cell_face_x_values(position, velocity_in_voxel_face_x, *this);
            interpolated_velocity[1] = WENO4_interpolation_3D_cell_face_y_values(position, velocity_in_voxel_face_y, *this);
            interpolated_velocity[2] = WENO4_interpolation_3D_cell_face_z_values(position, velocity_in_voxel_face_z, *this);
        }
        else if (interpolation_method == "WENO6" || interpolation_method == "WENO6-optimized") {
            interpolated_velocity[0] = WENO6_interpolation_3D_cell_face_x_values(position, velocity_in_voxel_face_x, *this);
            interpolated_velocity[1] = WENO6_interpolation_3D_cell_face_y_values(position, velocity_in_voxel_face_y, *this);
            interpolated_velocity[2] = WENO6_interpolation_3D_cell_face_z_values(position, velocity_in_voxel_face_z, *this);
        }
        else if (interpolation_method == "1Dy_linear"){
            interpolated_velocity[0] = linear_interpolation_1Dy_cell_face_x_values(position, velocity_in_voxel_face_x, *this);
            interpolated_velocity[1] = linear_interpolation_1Dy_cell_face_y_values(position, velocity_in_voxel_face_y, *this);
            interpolated_velocity[2] = linear_interpolation_1Dy_cell_face_z_values(position, velocity_in_voxel_face_z, *this);
        }
        else if (interpolation_method == "1Dy_WENO6"){
            interpolated_velocity[0] = WENO6_interpolation_1Dy_cell_face_x_values(position, velocity_in_voxel_face_x, *this);
            interpolated_velocity[1] = WENO6_interpolation_1Dy_cell_face_y_values(position, velocity_in_voxel_face_y, *this);
            interpolated_velocity[2] = WENO6_interpolation_1Dy_cell_face_z_values(position, velocity_in_voxel_face_z, *this);
        }
        return interpolated_velocity;
    }

    //position での deisity の値を周辺の cell center での値から interpolation して求める。
    MY_FLOAT_TYPE Grid_3D::calc_substance_density_by_interpolation(const VEC3_TYPE &position, const std::string interpolation_method) const {
		//linear interpolation を使う場合
		if (interpolation_method == "linear") {
            return linear_interpolation_3D_cell_center_values(
                position,
                substance_density,
                *this
            );
//			return calc_substance_density_by_linear_interpolation(position);
		}
        //WENO4 interpolation を使う場合
		else if (interpolation_method == "WENO4") {
            return WENO4_interpolation_3D_cell_center_values(
                position,
                substance_density,
                *this
            );
		}
		//linear interpolation を使う場合
		else if (interpolation_method == "WENO6") {
            return WENO6_interpolation_3D_cell_center_values(
                position,
                substance_density,
                *this
            );
		}
    }

    // cell center で定義される量cell_center_valuesのpositionでの補間を y方向の1次元linear interpolation で計算する
    MY_FLOAT_TYPE Grid_3D::interpolate_cell_center_defined_values_y_direction_1d_linear(
        const VEC3_TYPE &position,
        const std::vector<MY_FLOAT_TYPE>& cell_center_values) const{
        int advected_index_x = floor((position[0]) / _cell_length);
        int advected_index_y = floor((position[1] - 0.5 * _cell_length) / _cell_length);
        int advected_index_z = floor((position[2]) / _cell_length);

        /////linear interpolation の処理
        //補間に使う離散値をセット
        MY_FLOAT_TYPE val[2];
        for (int iy = 0; iy < 2; ++iy) {
            int index_x = advected_index_x;
            int index_y = advected_index_y + iy;
            int index_z = advected_index_z;
            // グリッドの外側を参照しようとしたときの処理(x方向)
            if (index_x < 0) {
                index_x = 0;
            }
            if (index_x > Grid_num_x - 1) {
                int exceed = index_x - (Grid_num_x - 1);
                //境界の値が外側までずっと続く場合
                index_x = Grid_num_x - 1;
                //境界を境に鏡のように値が反射する場合
                //index_x = Grid_num_x - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_y < 0) {
                index_y = 0;
            }
            if (index_y > Grid_num_y - 1) {
                int exceed = index_y - (Grid_num_y - 1);
                //境界の値が外側までずっと続く場合
                index_y = Grid_num_y - 1;
                //境界を境に鏡のように値が反射する場合
                //index_y = Grid_num_y - 1 - exceed + 1;
            }
            // グリッドの外側を参照しようとしたときの処理(y方向)
            if (index_z < 0) {
                index_z = 0;
            }
            if (index_z > Grid_num_z - 1) {
                int exceed = index_z - (Grid_num_z - 1);
                //境界の値が外側までずっと続く場合
                index_z = Grid_num_z - 1;
                //境界を境に鏡のように値が反射する場合
                //index_z = Grid_num_z - 1 - exceed + 1;
            }
            val[iy] = cell_center_values[get_voxel_center_index_3D(index_x, index_y, index_z, Grid_num_x, Grid_num_y, Grid_num_z)];
        }
        //interpolation の係数
        MY_FLOAT_TYPE b0, b1;
        b0 = (position[1] - 0.5 * _cell_length) / _cell_length - advected_index_y;
        b1 = 1.0 - b0;
        return b1 * val[0] + b0 * val[1];
    }

    // セルの頂点の速度を隣接する面での速度の平均によって計算する
    VEC3_TYPE Grid_3D::calc_cell_vertex_velocity(const Eigen::Vector3i cell_vertex_index) const{
        const int ix = cell_vertex_index[0];
        const int iy = cell_vertex_index[1];
        const int iz = cell_vertex_index[2];
        // 速度のx成分の計算
        MY_FLOAT_TYPE velocity_x;
        if(iy == 0 && iz == 0){
            velocity_x
                = velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, Grid_num_x, Grid_num_y, Grid_num_z)];
        }
        else if(iy == 0 && iz == Grid_num_z){
            velocity_x
                = velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy    , iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)];
        }
        else if(iy == Grid_num_y && iz == 0){
            velocity_x
                = velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy - 1, iz    , Grid_num_x, Grid_num_y, Grid_num_z)];
        }
        else if(iy == Grid_num_y && iz == Grid_num_z){
            velocity_x
                = velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy - 1, iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)];
        }
        else if(iy == 0){
            velocity_x
                = (velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy    , iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy    , iz    , Grid_num_x, Grid_num_y, Grid_num_z)]) / 2.0;
        }
        else if(iy == Grid_num_y){
            velocity_x
                = (velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy - 1, iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy - 1, iz    , Grid_num_x, Grid_num_y, Grid_num_z)]) / 2.0;
        }
        else if(iz == 0){
            velocity_x
                = (velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy - 1, iz    , Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy    , iz    , Grid_num_x, Grid_num_y, Grid_num_z)]) / 2.0;
        }
        else if(iz == Grid_num_z){
            velocity_x
                = (velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy - 1, iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy    , iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)]) / 2.0;
        }
        else{
            velocity_x
                = (velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy - 1, iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy    , iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy - 1, iz    , Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy    , iz    , Grid_num_x, Grid_num_y, Grid_num_z)]) / 4.0;
        }
        // 速度のy成分の計算
        MY_FLOAT_TYPE velocity_y;
        if(ix == 0 && iz == 0){
            velocity_y
                = velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix    , iy, iz    , Grid_num_x, Grid_num_y, Grid_num_z)];
        }
        else if(ix == 0 && iz ==Grid_num_z){
            velocity_y
                = velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix    , iy, iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)];
        }
        else if(ix == Grid_num_x && iz == 0){
            velocity_y
                = velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix - 1, iy, iz    , Grid_num_x, Grid_num_y, Grid_num_z)];
        }
        else if(ix == Grid_num_x && iz == Grid_num_z){
            velocity_y
                = velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix - 1, iy, iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)];
        }
        else if(ix == 0){
            velocity_y
                = (velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix    , iy, iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix    , iy, iz    , Grid_num_x, Grid_num_y, Grid_num_z)]) / 2.0;
        }
        else if(ix == Grid_num_x){
            velocity_y
                = (velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix - 1, iy, iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix - 1, iy, iz    , Grid_num_x, Grid_num_y, Grid_num_z)]) / 2.0;
        }
        else if(iz == 0){
            velocity_y
                = (velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix - 1, iy, iz    , Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix    , iy, iz    , Grid_num_x, Grid_num_y, Grid_num_z)]) / 2.0;
        }
        else if(iz == Grid_num_z){
            velocity_y
                = (velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix - 1, iy, iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix    , iy, iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)]) / 2.0;
        }
        else{
            velocity_y
                = (velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix - 1, iy, iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix    , iy, iz - 1, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix - 1, iy, iz    , Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix    , iy, iz    , Grid_num_x, Grid_num_y, Grid_num_z)]) / 4.0;
        }
        // 速度のx成分の計算
        MY_FLOAT_TYPE velocity_z;
        if(ix == 0 && iy == 0){
            velocity_z
                = velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix    , iy    , iz, Grid_num_x, Grid_num_y, Grid_num_z)];
        }
        else if(ix == 0 && iy == Grid_num_y){
            velocity_z
                = velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix    , iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)];
        }
        else if(ix == Grid_num_x && iy == 0){
            velocity_z
                = velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix - 1, iy    , iz, Grid_num_x, Grid_num_y, Grid_num_z)];
        }
        else if(ix == Grid_num_x && iy == Grid_num_y){
            velocity_z
                = velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix - 1, iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)];
        }
        else if(ix == 0){
            velocity_z
                = (velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix    , iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix    , iy    , iz, Grid_num_x, Grid_num_y, Grid_num_z)]) / 2.0;
        }
        else if(ix == Grid_num_x){
            velocity_z
                = (velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix - 1, iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix - 1, iy    , iz, Grid_num_x, Grid_num_y, Grid_num_z)]) / 2.0;
        }
        else if(iy == 0){
            velocity_z
                = (velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix - 1, iy    , iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix    , iy    , iz, Grid_num_x, Grid_num_y, Grid_num_z)]) / 2.0;
        }
        else if(iy == Grid_num_y){
            velocity_z
                = (velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix - 1, iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix    , iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)]) / 2.0;
        }
        else{
            velocity_z
                = (velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix - 1, iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix - 1, iy    , iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix    , iy - 1, iz, Grid_num_x, Grid_num_y, Grid_num_z)]
                +  velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix    , iy    , iz, Grid_num_x, Grid_num_y, Grid_num_z)]) / 4.0;
        }

        return VEC3_TYPE(velocity_x, velocity_y, velocity_z);
    }

}
