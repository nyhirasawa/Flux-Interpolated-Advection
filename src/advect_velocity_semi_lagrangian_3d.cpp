#include "advect_velocity_semi_lagrangian_3d.h"

#include "utils.h"
#include "linear_interpolation_3d.h"

namespace smoke_simulation{
    //advect項の計算
    //速度場を時間 -dt だけバックトレースしてadvect項を計算する
    void advect_velocity_semi_lagrangian_x_3D(
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const std::string interpolation_method,
        std::vector<MY_FLOAT_TYPE> &velocity_after_advect_x
    ) {
        //velocityのx成分を計算
        for(int ix=0;ix<all_grid.Grid_num_x+1;ix++){
            for(int iy=0;iy<all_grid.Grid_num_y;iy++){
                for(int iz=0;iz<all_grid.Grid_num_z;iz++){
                    MY_FLOAT_TYPE velocity_y ,velocity_z;
                    //考えてるx面での速度のy成分の計算
                    if(ix<=0){
                        velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    }
                    else if(ix>=all_grid.Grid_num_x){
                        velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix-1,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix-1,iy+1,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    }
                    else{
                        velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix-1,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix-1,iy+1,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/4.0;
                    }
                    //考えてるx面での速度z成分の計算
                    if(ix<=0){
                        velocity_z=(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz+1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    }
                    else if(ix>=all_grid.Grid_num_x){
                        velocity_z=(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix-1,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix-1,iy,iz+1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    }
                    else{
                        velocity_z=(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix-1,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix-1,iy,iz+1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz+1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/4.0;
                    }

                    //バックトレース先の位置
                    MY_FLOAT_TYPE advected_pos_x=(ix) * all_grid._cell_length - (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] * time_step_length);
                    MY_FLOAT_TYPE advected_pos_y=(iy + 0.5) * all_grid._cell_length - (velocity_y * time_step_length);
                    MY_FLOAT_TYPE advected_pos_z=(iz + 0.5) * all_grid._cell_length - (velocity_z * time_step_length);
                    bool tmp;
                    velocity_after_advect_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        = linear_interpolation_3D_cell_face_x_values(
                            VEC3_TYPE(advected_pos_x, advected_pos_y, advected_pos_z),
                            all_grid.velocity_in_voxel_face_x,
                            all_grid
                        );
                }
            }
        }
    }

    //advect項の計算
    //速度場を時間 -dt だけバックトレースしてadvect項を計算する
    void advect_velocity_semi_lagrangian_y_3D(
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const std::string interpolation_method,
        std::vector<MY_FLOAT_TYPE> &velocity_after_advect_y
    ) {
        //velocityのy成分を計算
        for(int ix=0;ix<all_grid.Grid_num_x;ix++){
            for(int iy=0;iy<all_grid.Grid_num_y+1;iy++){
                for(int iz=0;iz<all_grid.Grid_num_z;iz++){
                    MY_FLOAT_TYPE velocity_x, velocity_z;
                    //考えてるy面での速度のx成分の計算
                    if(iy<=0){
                        velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    }
                    else if(iy>=all_grid.Grid_num_y){
                        velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy-1,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy-1,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    }
                    else{
                        velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy-1,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy-1,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/4.0;
                    }
                    //考えてるy面での速度のz成分の計算
                    if(iy<=0){
                        velocity_z=(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz+1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    }
                    else if(iy>=all_grid.Grid_num_y){
                        velocity_z=(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy-1,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy-1,iz+1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    }
                    else{
                        velocity_z=(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy-1,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy-1,iz+1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz+1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/4.0;
                    }

                    //バックトレース先の位置
                    MY_FLOAT_TYPE advected_pos_x=(ix + 0.5) * all_grid._cell_length - (velocity_x * time_step_length);
                    MY_FLOAT_TYPE advected_pos_y=(iy) * all_grid._cell_length - (all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] * time_step_length);
                    MY_FLOAT_TYPE advected_pos_z=(iz + 0.5) * all_grid._cell_length - (velocity_z * time_step_length);

                    bool tmp;
                    velocity_after_advect_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        = linear_interpolation_3D_cell_face_y_values(
                            VEC3_TYPE(advected_pos_x, advected_pos_y, advected_pos_z),
                            all_grid.velocity_in_voxel_face_y,
                            all_grid
                        );
                }
            }
        }
    }

    //advect項の計算
    //速度場を時間 -dt だけバックトレースしてadvect項を計算する
    void advect_velocity_semi_lagrangian_z_3D(
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const std::string interpolation_method,
        std::vector<MY_FLOAT_TYPE> &velocity_after_advect_z
    ) {
        //velocityのz成分を計算
        for(int ix=0;ix<all_grid.Grid_num_x;ix++){
            for(int iy=0;iy<all_grid.Grid_num_y;iy++){
                for(int iz=0;iz<all_grid.Grid_num_z+1;iz++){
                    MY_FLOAT_TYPE velocity_x, velocity_y;
                    //考えてるz面での速度のx成分の計算
                    if(iz<=0){
                        velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    }
                    else if(iz>=all_grid.Grid_num_z){
                        velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz-1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz-1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    }
                    else{
                        velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz-1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz-1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/4.0;
                    }
                    //考えてるz面での速度のy成分の計算
                    if(iz<=0){
                        velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    }
                    else if(iz>=all_grid.Grid_num_z){
                        velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz-1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz-1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/2.0;
                    }
                    else{
                        velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz-1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz-1, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                                   +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)])/4.0;
                    }

                    //バックトレース先の位置
                    MY_FLOAT_TYPE advected_pos_x=(ix + 0.5) * all_grid._cell_length - (velocity_x * time_step_length);
                    MY_FLOAT_TYPE advected_pos_y=(iy + 0.5) * all_grid._cell_length - (velocity_y * time_step_length);
                    MY_FLOAT_TYPE advected_pos_z=(iz) * all_grid._cell_length - (all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)] * time_step_length);

                    bool tmp;
                    velocity_after_advect_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        = linear_interpolation_3D_cell_face_z_values(
                            VEC3_TYPE(advected_pos_x, advected_pos_y, advected_pos_z),
                            all_grid.velocity_in_voxel_face_z,
                            all_grid
                        );
                }
            }
        }
    }

    //advect項の計算 (semi-Lagrangian)
    void advect_velocity_semi_lagrangian_3D(
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const std::string interpolation_method
    ) {
        //移流後の速度場
        std::vector<MY_FLOAT_TYPE> velocity_after_advect_x((all_grid.Grid_num_x + 1) * all_grid.Grid_num_y * all_grid.Grid_num_z);
        std::vector<MY_FLOAT_TYPE> velocity_after_advect_y(all_grid.Grid_num_x * (all_grid.Grid_num_y + 1) * all_grid.Grid_num_z);
        std::vector<MY_FLOAT_TYPE> velocity_after_advect_z(all_grid.Grid_num_x * all_grid.Grid_num_y * (all_grid.Grid_num_z + 1));
        //x成分の移流
        advect_velocity_semi_lagrangian_x_3D(
            all_grid,
            time_step_length,
            interpolation_method,
            velocity_after_advect_x
        );
        //y成分の移流
        advect_velocity_semi_lagrangian_y_3D(
            all_grid,
            time_step_length,
            interpolation_method,
            velocity_after_advect_y
        );
        //z成分の移流
        advect_velocity_semi_lagrangian_z_3D(
            all_grid,
            time_step_length,
            interpolation_method,
            velocity_after_advect_z
        );
        //計算結果をコピー
        for (int ix = 0; ix < all_grid.Grid_num_x + 1; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
                    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        = velocity_after_advect_x[get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                }
            }
        }
        for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y + 1; iy++) {
                for (int iz = 0; iz < all_grid.Grid_num_z; iz++) {
                    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        = velocity_after_advect_y[get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                }
            }
        }
        for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                for (int iz = 0; iz < all_grid.Grid_num_z + 1; iz++) {
                    all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
                        = velocity_after_advect_z[get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)];
                }
            }
        }
    }
} // namespace smoke_simulation
