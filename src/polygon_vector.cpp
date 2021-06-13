#include "polygon_vector.h"

#include "polygon_3D.h"
#include "define_float_type.h"
#include <thread>

namespace smoke_simulation {
polygon_vector::polygon_vector(): num_threads(std::thread::hardware_concurrency()){
    data.resize(std::thread::hardware_concurrency());
}
}// namespace smoke_simulation
