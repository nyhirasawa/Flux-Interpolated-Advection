#pragma once

#include <vector>
#include "cell_face.h"
#include "grid.h"

namespace smoke_simulation {
    std::vector<cell_face> split_face(
        const Grid& all_grid,
        const cell_face& face,
        const std::string split_method
    );
}
