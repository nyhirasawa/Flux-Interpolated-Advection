#ifndef SMOKE_SIMULATION_GRID_1D_H
#define SMOKE_SIMULATION_GRID_1D_H

#include <vector>
#include "define_float_type.h"

namespace smoke_simulation {
    //系の物理量が乗るグリッドの定義
    class Grid_1D {
    public:
        const int Grid_num_y;
        const MY_FLOAT_TYPE min_pos_y;
        const MY_FLOAT_TYPE _cell_length;
//		const std::string _interpolation_method;
        std::vector<MY_FLOAT_TYPE> velocity_in_voxel_face_y;
        std::vector<MY_FLOAT_TYPE> substance_density;
        std::vector<MY_FLOAT_TYPE> psi_substance_density_cell_face_y;
        std::vector<MY_FLOAT_TYPE> cell_volume_cell_center;

        Grid_1D(int ny, MY_FLOAT_TYPE cell_length);
        ~Grid_1D();

        MY_FLOAT_TYPE calc_interpolated_velocity(MY_FLOAT_TYPE pos);
        MY_FLOAT_TYPE calc_interpolated_psi(MY_FLOAT_TYPE pos);
    };

} //namespace smoke_simulation

#endif // SMOKE_SIMULATION_GRID_1D_H
