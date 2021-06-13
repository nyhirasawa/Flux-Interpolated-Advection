#pragma once

#include "grid_3d.h"

namespace smoke_simulation{
    //Poisson eq. を解くことによる圧力の計算
    void calc_pressure_3D(Grid_3D& all_grid);
    //pressure gradient termの計算
    void calc_pressure_gradient_term_3D(Grid_3D& all_grid);
} // namespace smoke_simulation
