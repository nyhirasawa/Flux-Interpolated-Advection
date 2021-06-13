#ifndef WRITE_SUBSTANCE_DENSITY_DATA_H
#define WRITE_SUBSTANCE_DENSITY_DATA_H

#include "grid_3d.h"

#include "define_float_type.h"

namespace smoke_simulation{

void write_substance_density_data(Grid_3D& all_grid, int file_number, const std::string density_data_path);

}//namespace smoke_simulation
#endif//WRITE_SUBSTANCE_DENSITY_DATA_H
