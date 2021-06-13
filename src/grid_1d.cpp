#include "grid_1d.h"

#include "linear_interpolation_1d.h"
#include <iostream>

namespace smoke_simulation {

    Grid_1D::Grid_1D(
        int ny,
        MY_FLOAT_TYPE cell_length):
        Grid_num_y(ny),
        min_pos_y(0.0),
        _cell_length(cell_length)
        {
            velocity_in_voxel_face_y.resize(ny + 1);
            substance_density.resize(ny);
            psi_substance_density_cell_face_y.resize(ny + 1);
            cell_volume_cell_center.resize(ny);
    }
    Grid_1D::~Grid_1D(){

    }


    MY_FLOAT_TYPE Grid_1D::calc_interpolated_velocity(MY_FLOAT_TYPE pos) {
        return 0.1;
    }

    MY_FLOAT_TYPE Grid_1D::calc_interpolated_psi(MY_FLOAT_TYPE pos){
        return linear_interpolation_1D(
            pos,
            psi_substance_density_cell_face_y,
            *this,
            0.0,
            Grid_num_y,
            _cell_length
        );
/*
        MY_FLOAT_TYPE interpolated_psi = 0.0;

        // 下面のpsiの値からの寄与
        int advected_index = floor(pos / _cell_length);
//        std::cout << "psi_substance_density_cell_face_y[advected_index] begin" <<std::endl;
        if(advected_index < 0){
            advected_index = 0;
        }
        if(advected_index > Grid_num_y){
            advected_index = Grid_num_y;
        }
        interpolated_psi += psi_substance_density_cell_face_y[advected_index];
//        std::cout << "psi_substance_density_cell_face_y[advected_index] end" <<std::endl;

        MY_FLOAT_TYPE over_under_face = (pos / _cell_length) - advected_index;
        //posがセルの半分より下にあるとき
        if(over_under_face < 0.5 ){
            MY_FLOAT_TYPE density_at_quadrature_point = linear_interpolation_1D_cell_center_values(
                (pos + advected_index * _cell_length) / 2.0,
                substance_density,
                *this
            );
            interpolated_psi += density_at_quadrature_point * over_under_face * _cell_length;
        }
        //posがセルの半分より上にあるとき
        else{
            MY_FLOAT_TYPE density_at_quadrature_point = linear_interpolation_1D_cell_center_values(
                (advected_index * _cell_length + (advected_index + 0.5) * _cell_length) / 2.0,
                substance_density,
                *this
            );
            interpolated_psi += density_at_quadrature_point * 0.5 * _cell_length;

            density_at_quadrature_point = linear_interpolation_1D_cell_center_values(
                (pos + (advected_index + 0.5) * _cell_length) / 2.0,
                substance_density,
                *this
            );
            interpolated_psi += density_at_quadrature_point * (over_under_face - 0.5) * _cell_length;
        }
        return interpolated_psi;
*/
    }

} //namespace smoke_simulation
