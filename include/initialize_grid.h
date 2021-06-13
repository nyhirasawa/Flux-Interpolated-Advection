#ifndef INITIALIZE_GRID_H
#define INITIALIZE_GRID_H

#include "grid_1d.h"
#include "grid.h"
#include "grid_3d.h"

#include "define_float_type.h"

namespace smoke_simulation{
    //グリッドを初期化する関数
    void initialize_grid(
        Grid_1D& all_grid
    );
    void initialize_grid(
        Grid& all_grid,
        const bool fix_velocity,
        const bool initialize_velocity,
        const bool initialize_density
    );
    //グリッドを初期化する関数
    void initialize_Grid_3D(Grid_3D& all_grid, const bool fix_velocity);
}//namespace smoke_simulation

#endif //INITIALIZE_GRID_H
