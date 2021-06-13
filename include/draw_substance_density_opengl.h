#ifndef DRAW_SUBSTANCE_DENSITY_OPENGL_H
#define DRAW_SUBSTANCE_DENSITY_OPENGL_H

#include "grid.h"
#include "physical_const.h"
#include <GLFW/glfw3.h>
#include <GL/glut.h>

#include "define_float_type.h"

namespace smoke_simulation{
    void draw_substance_density_opengl(
        GLFWwindow* window,
        std::string img_file_path,
        Grid all_grid,
        const int i_frame,
        int& i_movie_frame,
        const bool write_frame_to_movie_file,
        const MY_FLOAT_TYPE time_step_length,
        const bool color_negative_density
    );

    void draw_test(
        GLFWwindow* window,
        std::string img_file_path,
        int i_frame
    );

}//namespace smoke_simulation

#endif //DRAW_SUBSTANCE_DENSITY_OPENGL_H
