#include "write_substance_density_data.h"
#include <fstream>
#include <string>
#include <sstream>
#include "physical_const.h"
#include "utils.h"
#include "define_float_type.h"

namespace smoke_simulation{

void write_substance_density_data(Grid_3D& all_grid, int file_number, const std::string density_data_path){
    std::ostringstream filename;
    filename<<density_data_path<<"/substance_density_"<<file_number<<".dat"<<std::flush;
//    std::string filename = "test.txt";
    std::ofstream writing_file;
    writing_file.open(filename.str(), std::ios::out);
    writing_file<<all_grid.Grid_num_x<<" "
                <<all_grid.Grid_num_y<<" "
                <<all_grid.Grid_num_z<<std::endl;
    for(int ix=0;ix<all_grid.Grid_num_x;ix++){
        for(int iy=0;iy<all_grid.Grid_num_y;iy++){
            for(int iz=0;iz<all_grid.Grid_num_z;iz++){
                writing_file<<all_grid.substance_density[smoke_simulation::get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]<<std::endl;
            }
        }
    }
}

}
