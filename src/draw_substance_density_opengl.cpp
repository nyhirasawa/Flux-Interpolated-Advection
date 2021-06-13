#include "draw_substance_density_opengl.h"

#include "grid.h"
#include "physical_const.h"
#include "utils.h"
#include <iostream>
#include <GLFW/glfw3.h>
#include <GL/glut.h>
#include <save_image.h>
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
    ){
        //密度場の描画
        // OpenGL animation calls
        glClear(GL_COLOR_BUFFER_BIT);
        for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
            for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                const MY_FLOAT_TYPE cell_face_length_on_window_x = 2.0 / all_grid.Grid_num_x;
                const MY_FLOAT_TYPE cell_face_length_on_window_y = 2.0 / all_grid.Grid_num_y;

                //セルの 4頂点の位置
                std::vector<VEC3_TYPE> vertex_pos_list(4);
                vertex_pos_list[0] = VEC3_TYPE(-1.0 + cell_face_length_on_window_x * ix, -1.0 + cell_face_length_on_window_y * iy, 0.0);
                vertex_pos_list[1] = VEC3_TYPE(-1.0 + cell_face_length_on_window_x * (ix + 1), -1.0 + cell_face_length_on_window_y * iy, 0.0);
                vertex_pos_list[2] = VEC3_TYPE(-1.0 + cell_face_length_on_window_x * (ix + 1), -1.0 + cell_face_length_on_window_y * (iy + 1), 0.0);
                vertex_pos_list[3] = VEC3_TYPE(-1.0 + cell_face_length_on_window_x * ix, -1.0 + cell_face_length_on_window_y * (iy + 1), 0.0);

                MY_FLOAT_TYPE density_value = all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
//                density_value *= all_grid._cell_volume / all_grid.cell_volume_cell_center[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
                //密度場を強調
//                density_value *= 100.0;
                // [0, 1]に正規化されたセルの色
                VEC3_TYPE cell_color;
                //セルの色を計算
                if (density_value > 0.0) {
                    if (density_value > 1.0) {
                        cell_color[0]=1.0;
                        cell_color[1]=1.0;
                        cell_color[2]=1.0;
                    }
                    else {
                        cell_color[0]=density_value;
                        cell_color[1]=density_value;
                        cell_color[2]=density_value;
                    }
                }
                // density が負の場合
                else {
                    if (color_negative_density) {
                        //src.at<unsigned char>(scale * (all_grid.Grid_num_y - 1 - iy) + j, scale * (ix)+k) = (unsigned char)255;
//                        if (density_value < -0.0) {
                            if(density_value < -1.0){
                                cell_color[0]=1.0;
                                cell_color[1]=0.0;
                                cell_color[2]=0.0;
                            }
                            else{
                                cell_color[0]=-density_value;
                                cell_color[1]=0.0;
                                cell_color[2]=0.0;
                            }
                        }
                        else {
                            cell_color[0]=0.0;
                            cell_color[1]=0.0;
                            cell_color[2]=0.0;
                        }
//                    }
                }
                //セルの色
                glColor4f(cell_color[0], cell_color[1], cell_color[2], 1.0);
                //セルを表す四角形を描画
                glBegin(GL_QUADS);
                glVertex2d(vertex_pos_list[0][0], vertex_pos_list[0][1]);
                glVertex2d(vertex_pos_list[1][0], vertex_pos_list[1][1]);
                glVertex2d(vertex_pos_list[2][0], vertex_pos_list[2][1]);
                glVertex2d(vertex_pos_list[3][0], vertex_pos_list[3][1]);
                glEnd();
            }
        }

        ////// バックトレース後のグリッドを描画 //////
        if (physical_const::kDraw_backtrace_grid) {
            //draw_backtrace_grid(src, all_grid, scale, time_step_length);
        }
        ////// バックトレース後のグリッドの描画終了 //////

        //cv::imshow(" ", src);
        if (i_frame == 0 || write_frame_to_movie_file) {
            //std::cout << "write mov file:  " << i_movie_frame << "/" << physical_const::kMax_num_movie_frames << std::endl;;
            std::cout << "write mov file:  " << i_movie_frame <<"  frame"<< std::endl;;
            //画像を保存
            save_image(window, img_file_path+"/screenshot_"+std::to_string(i_movie_frame)+".png");

            //writer << src;
            //画像を書き出す
            /*
            if (physical_const::kWrite_image_files) {
                std::stringstream ss;
                ss << "img/result" << i_movie_frame << ".jpg";
                cv::imwrite(ss.str(), src);
            }
            */
            i_movie_frame += 1;
        }
        ///*
        if (physical_const::kWrite_image_files) {
            std::stringstream ss;
            ss << "../../flux-advection_data/img/result" << i_frame << ".jpg";
            //cv::imwrite(ss.str(), src);
        }
        //*/
        //else if (write_frame_to_movie_file) {
        //	std::cout << "write mov file:  " << i_movie_frame << "/" << physical_const::kMax_num_movie_frames << std::endl;;
        //	writer << src;
        //	i_movie_frame += 1;
        //	//画像を書き出す
        //	if (physical_const::kWrite_image_files) {
        //		std::stringstream ss;
        //		ss << "img/result" << i_movie_frame << ".jpg";
        //		cv::imwrite(ss.str(), src);
        //	}
        //}

        //cv::waitKey(1);

        //バッファの情報を描画
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}//namespace smoke_simulation
