CXX = g++
CXXFLAGS = -O3 -MMD -MP -Wno-unused-result
BUILDDIR = ./build
OBJDIR = $(BUILDDIR)/obj
TARGET = $(BUILDDIR)/a.out
#OBJS := main.o sginterp2.o sginterp3.o sginterp2_ispc_float.o sginterp2_ispc_double.o sginterp3_ispc_float.o sginterp3_ispc_double.o
OBJS := $(OBJDIR)/calc_psi.o $(OBJDIR)/main.o $(OBJDIR)/sginterp2.o $(OBJDIR)/sginterp3.o $(OBJDIR)/sginterp2_ispc_float.o $(OBJDIR)/sginterp2_ispc_double.o $(OBJDIR)/sginterp3_ispc_float.o $(OBJDIR)/sginterp3_ispc_double.o $(OBJDIR)/calc_time_step_length_from_CFL_number.o $(OBJDIR)/cell_face.o $(OBJDIR)/cell_face_3D.o $(OBJDIR)/cell_vertex.o $(OBJDIR)/cell_vertex_3D.o $(OBJDIR)/correct_volume_concentration_error.o $(OBJDIR)/draw_substance_density_opengl.o $(OBJDIR)/gauss_quadrature_points.o $(OBJDIR)/grid.o $(OBJDIR)/grid_3d.o $(OBJDIR)/initialize_grid.o $(OBJDIR)/linear_solver.o $(OBJDIR)/lodepng.o $(OBJDIR)/move_substances.o $(OBJDIR)/move_substances_3d.o $(OBJDIR)/physical_const.o $(OBJDIR)/polygon_3D.o $(OBJDIR)/polygon_vector.o $(OBJDIR)/save_image.o $(OBJDIR)/sparse_matrix.o $(OBJDIR)/split_face_3D.o $(OBJDIR)/triangle_3D.o $(OBJDIR)/update_fluid_velocity.o $(OBJDIR)/update_fluid_velocity_3d.o $(OBJDIR)/utils.o $(OBJDIR)/vdbexport.o $(OBJDIR)/write_substance_density_data.o  $(OBJDIR)/calc_backtraced_face_3D.o $(OBJDIR)/backtrace_and_calc_all_quadrature_point_list.o $(OBJDIR)/calc_density_in_cell_from_quadrature_point_list.o $(OBJDIR)/clamping_in_MacCormack_scheme.o $(OBJDIR)/linear_interpolation_3d.o $(OBJDIR)/WENO4_interpolation_3d.o $(OBJDIR)/WENO6_interpolation_3d.o $(OBJDIR)/linear_interpolation_1d.o $(OBJDIR)/WENO6_interpolation_1d.o $(OBJDIR)/backtrace_and_calc_all_quadrature_point_list_for_integrate_density.o $(OBJDIR)/calc_density_in_cell_from_quadrature_point_list_for_integrate_density.o $(OBJDIR)/calc_cell_volumes.o $(OBJDIR)/calc_backtraced_face.o $(OBJDIR)/split_face.o $(OBJDIR)/linear_interpolation_2d.o $(OBJDIR)/WENO6_interpolation_2d.o $(OBJDIR)/calc_pressure.o $(OBJDIR)/calc_psi_3D.o $(OBJDIR)/calc_velocity_from_quadrature_point_list.o $(OBJDIR)/calc_velocity_from_quadrature_point_list_use_integral.o $(OBJDIR)/advect_velocity_semi_lagrangian_3d.o $(OBJDIR)/calc_pressure_3d.o $(OBJDIR)/advect_velocity_flux_advection_3d.o $(OBJDIR)/parallelize_functions.o $(OBJDIR)/grid_1d.o $(OBJDIR)/move_substances_1d.o $(OBJDIR)/calc_psi_1D.o
LDFLAGS = -lglut -lGLU -lGL -lglfw /nfsdata/opt/lib/libopenvdb.a /nfsdata/opt/lib/libblosc.a -ltbb -lHalf -no-pie -pthread -fopenmp -lgomp -lboost_iostreams -lz -lboost_system -lboost_thread
LIBS = -L/usr/local/lib -L/lib/x86_64-linux-gnu
INCLUDE = -I./include -I/nfsdata/opt/include -I./include/eigen-3.3.7 -I./include/amgcl

main: $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/main.o: ./src/main.cpp
	$(CXX) $(CXXFLAGS) -c ./src/main.cpp -o $(OBJDIR)/main.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/main.cpp -o $(OBJDIR)/main.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/sginterp2.o: ./src/sginterp2.cpp ./include/sginterp2.h
	$(CXX) $(CXXFLAGS) -c ./src/sginterp2.cpp -o $(OBJDIR)/sginterp2.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/sginterp2.cpp -o $(OBJDIR)/sginterp2.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/sginterp3.o: ./src/sginterp3.cpp ./include/sginterp3.h
	$(CXX) $(CXXFLAGS) -c ./src/sginterp3.cpp -o $(OBJDIR)/sginterp3.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/sginterp3.cpp -o $(OBJDIR)/sginterp3.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/sginterp2_ispc_float.o: ./src/sginterp2.ispc ./src/weno.ispc
	ispc -O3 ./src/sginterp2.ispc -o $(OBJDIR)/sginterp2_ispc_float.o -DT=float -DUSE_FLOAT --target=avx2
	#
$(OBJDIR)/sginterp2_ispc_double.o: ./src/sginterp2.ispc ./src/weno.ispc
	ispc -O3 ./src/sginterp2.ispc -o $(OBJDIR)/sginterp2_ispc_double.o -DT=double -DUSE_DOUBLE --target=avx2
	#
$(OBJDIR)/sginterp3_ispc_float.o: ./src/sginterp3.ispc ./src/weno.ispc
	ispc -O3 ./src/sginterp3.ispc -o $(OBJDIR)/sginterp3_ispc_float.o -DT=float -DUSE_FLOAT --target=avx2
	#
$(OBJDIR)/sginterp3_ispc_double.o: ./src/sginterp3.ispc ./src/weno.ispc
	ispc -O3 ./src/sginterp3.ispc -o $(OBJDIR)/sginterp3_ispc_double.o -DT=double -DUSE_DOUBLE --target=avx2

$(OBJDIR)/calc_psi.o: ./src/calc_psi.cpp ./include/calc_psi.h
	$(CXX) $(CXXFLAGS) -c ./src/calc_psi.cpp -o $(OBJDIR)/calc_psi.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/calc_psi.cpp -o $(OBJDIR)/calc_psi.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/calc_time_step_length_from_CFL_number.o: ./src/calc_time_step_length_from_CFL_number.cpp ./include/calc_time_step_length_from_CFL_number.h
	$(CXX) $(CXXFLAGS) -c ./src/calc_time_step_length_from_CFL_number.cpp -o $(OBJDIR)/calc_time_step_length_from_CFL_number.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/calc_time_step_length_from_CFL_number.cpp -o $(OBJDIR)/calc_time_step_length_from_CFL_number.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/cell_face.o: ./src/cell_face.cpp ./include/cell_face.h
	$(CXX) $(CXXFLAGS) -c ./src/cell_face.cpp -o $(OBJDIR)/cell_face.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/cell_face.cpp -o $(OBJDIR)/cell_face.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/cell_face_3D.o: ./src/cell_face_3D.cpp ./include/cell_face_3D.h
	$(CXX) $(CXXFLAGS) -c ./src/cell_face_3D.cpp -o $(OBJDIR)/cell_face_3D.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/cell_face_3D.cpp -o $(OBJDIR)/cell_face_3D.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/cell_vertex.o: ./src/cell_vertex.cpp ./include/cell_vertex.h
	$(CXX) $(CXXFLAGS) -c ./src/cell_vertex.cpp -o $(OBJDIR)/cell_vertex.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/cell_vertex.cpp -o $(OBJDIR)/cell_vertex.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/cell_vertex_3D.o: ./src/cell_vertex_3D.cpp ./include/cell_vertex_3D.h
	$(CXX) $(CXXFLAGS) -c ./src/cell_vertex_3D.cpp -o $(OBJDIR)/cell_vertex_3D.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/cell_vertex_3D.cpp -o $(OBJDIR)/cell_vertex_3D.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/correct_volume_concentration_error.o: ./src/correct_volume_concentration_error.cpp ./include/correct_volume_concentration_error.h
	$(CXX) $(CXXFLAGS) -c ./src/correct_volume_concentration_error.cpp -o $(OBJDIR)/correct_volume_concentration_error.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/correct_volume_concentration_error.cpp -o $(OBJDIR)/correct_volume_concentration_error.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/draw_substance_density_opengl.o: ./src/draw_substance_density_opengl.cpp ./include/draw_substance_density_opengl.h
	$(CXX) $(CXXFLAGS) -c ./src/draw_substance_density_opengl.cpp -o $(OBJDIR)/draw_substance_density_opengl.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/draw_substance_density_opengl.cpp -o $(OBJDIR)/draw_substance_density_opengl.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/gauss_quadrature_points.o: ./src/gauss_quadrature_points.cpp ./include/gauss_quadrature_points.h
	$(CXX) $(CXXFLAGS) -c ./src/gauss_quadrature_points.cpp -o $(OBJDIR)/gauss_quadrature_points.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/gauss_quadrature_points.cpp -o $(OBJDIR)/gauss_quadrature_points.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/grid.o: ./src/grid.cpp ./include/grid.h
	$(CXX) $(CXXFLAGS) -c ./src/grid.cpp -o $(OBJDIR)/grid.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/grid.cpp -o $(OBJDIR)/grid.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/grid_3d.o: ./src/grid_3d.cpp ./include/grid_3d.h
	$(CXX) $(CXXFLAGS) -c ./src/grid_3d.cpp -o $(OBJDIR)/grid_3d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/grid_3d.cpp -o $(OBJDIR)/grid_3d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/initialize_grid.o: ./src/initialize_grid.cpp ./include/initialize_grid.h
	$(CXX) $(CXXFLAGS) -c ./src/initialize_grid.cpp -o $(OBJDIR)/initialize_grid.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/initialize_grid.cpp -o $(OBJDIR)/initialize_grid.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/linear_solver.o: ./src/linear_solver.cpp ./include/linear_solver.h
	$(CXX) $(CXXFLAGS) -c ./src/linear_solver.cpp -o $(OBJDIR)/linear_solver.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/linear_solver.cpp -o $(OBJDIR)/linear_solver.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/lodepng.o: ./src/lodepng.cpp ./include/lodepng.h
	$(CXX) $(CXXFLAGS) -c ./src/lodepng.cpp -o $(OBJDIR)/lodepng.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/lodepng.cpp -o $(OBJDIR)/lodepng.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/move_substances.o: ./src/move_substances.cpp ./include/move_substances.h
	$(CXX) $(CXXFLAGS) -c ./src/move_substances.cpp -o $(OBJDIR)/move_substances.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/move_substances.cpp -o $(OBJDIR)/move_substances.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/move_substances_3d.o: ./src/move_substances_3d.cpp ./include/move_substances_3d.h
	$(CXX) $(CXXFLAGS) -c ./src/move_substances_3d.cpp -o $(OBJDIR)/move_substances_3d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/move_substances_3d.cpp -o $(OBJDIR)/move_substances_3d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/physical_const.o: ./src/physical_const.cpp ./include/physical_const.h
	$(CXX) $(CXXFLAGS) -c ./src/physical_const.cpp -o $(OBJDIR)/physical_const.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/physical_const.cpp -o $(OBJDIR)/physical_const.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/polygon_3D.o: ./src/polygon_3D.cpp ./include/polygon_3D.h
	$(CXX) $(CXXFLAGS) -c ./src/polygon_3D.cpp -o $(OBJDIR)/polygon_3D.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/polygon_3D.cpp -o $(OBJDIR)/polygon_3D.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/polygon_vector.o: ./src/polygon_vector.cpp ./include/polygon_vector.h
	$(CXX) $(CXXFLAGS) -c ./src/polygon_vector.cpp -o $(OBJDIR)/polygon_vector.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/polygon_vector.cpp -o $(OBJDIR)/polygon_vector.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/save_image.o: ./src/save_image.cpp ./include/save_image.h
	$(CXX) $(CXXFLAGS) -c ./src/save_image.cpp -o $(OBJDIR)/save_image.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/save_image.cpp -o $(OBJDIR)/save_image.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/sparse_matrix.o: ./src/sparse_matrix.cpp ./include/sparse_matrix.h
	$(CXX) $(CXXFLAGS) -c ./src/sparse_matrix.cpp -o $(OBJDIR)/sparse_matrix.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/sparse_matrix.cpp -o $(OBJDIR)/sparse_matrix.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/split_face_3D.o: ./src/split_face_3D.cpp ./include/split_face_3D.h ./src/facecutter.h
	$(CXX) $(CXXFLAGS) -c ./src/split_face_3D.cpp -o $(OBJDIR)/split_face_3D.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/split_face_3D.cpp -o $(OBJDIR)/split_face_3D.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/triangle_3D.o: ./src/triangle_3D.cpp ./include/triangle_3D.h
	$(CXX) $(CXXFLAGS) -c ./src/triangle_3D.cpp -o $(OBJDIR)/triangle_3D.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/triangle_3D.cpp -o $(OBJDIR)/triangle_3D.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/update_fluid_velocity.o: ./src/update_fluid_velocity.cpp ./include/update_fluid_velocity.h
	$(CXX) $(CXXFLAGS) -c ./src/update_fluid_velocity.cpp -o $(OBJDIR)/update_fluid_velocity.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/update_fluid_velocity.cpp -o $(OBJDIR)/update_fluid_velocity.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/update_fluid_velocity_3d.o: ./src/update_fluid_velocity_3d.cpp ./include/update_fluid_velocity_3d.h
	$(CXX) $(CXXFLAGS) -c ./src/update_fluid_velocity_3d.cpp -o $(OBJDIR)/update_fluid_velocity_3d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/update_fluid_velocity_3d.cpp -o $(OBJDIR)/update_fluid_velocity_3d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/utils.o: ./src/utils.cpp ./include/utils.h
	$(CXX) $(CXXFLAGS) -c ./src/utils.cpp -o $(OBJDIR)/utils.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/utils.cpp -o $(OBJDIR)/utils.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/vdbexport.o: ./src/vdbexport.cpp ./include/vdbexport.h
	$(CXX) $(CXXFLAGS) -c ./src/vdbexport.cpp -o $(OBJDIR)/vdbexport.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/vdbexport.cpp -o $(OBJDIR)/vdbexport.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/write_substance_density_data.o: ./src/write_substance_density_data.cpp ./include/write_substance_density_data.h
	$(CXX) $(CXXFLAGS) -c ./src/write_substance_density_data.cpp -o $(OBJDIR)/write_substance_density_data.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/write_substance_density_data.cpp -o $(OBJDIR)/write_substance_density_data.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/backtrace_and_calc_all_quadrature_point_list.o: ./src/backtrace_and_calc_all_quadrature_point_list.cpp ./include/backtrace_and_calc_all_quadrature_point_list.h
	$(CXX) $(CXXFLAGS) -c ./src/backtrace_and_calc_all_quadrature_point_list.cpp -o $(OBJDIR)/backtrace_and_calc_all_quadrature_point_list.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/backtrace_and_calc_all_quadrature_point_list.cpp -o $(OBJDIR)/backtrace_and_calc_all_quadrature_point_list.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/calc_backtraced_face_3D.o: ./src/calc_backtraced_face_3D.cpp ./include/calc_backtraced_face_3D.h
	$(CXX) $(CXXFLAGS) -c ./src/calc_backtraced_face_3D.cpp -o $(OBJDIR)/calc_backtraced_face_3D.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/calc_backtraced_face_3D.cpp -o $(OBJDIR)/calc_backtraced_face_3D.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/calc_density_in_cell_from_quadrature_point_list.o: ./src/calc_density_in_cell_from_quadrature_point_list.cpp ./include/calc_density_in_cell_from_quadrature_point_list.h
	$(CXX) $(CXXFLAGS) -c ./src/calc_density_in_cell_from_quadrature_point_list.cpp -o $(OBJDIR)/calc_density_in_cell_from_quadrature_point_list.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/calc_density_in_cell_from_quadrature_point_list.cpp -o $(OBJDIR)/calc_density_in_cell_from_quadrature_point_list.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/clamping_in_MacCormack_scheme.o: ./src/clamping_in_MacCormack_scheme.cpp ./include/clamping_in_MacCormack_scheme.h
	$(CXX) $(CXXFLAGS) -c ./src/clamping_in_MacCormack_scheme.cpp -o $(OBJDIR)/clamping_in_MacCormack_scheme.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/clamping_in_MacCormack_scheme.cpp -o $(OBJDIR)/clamping_in_MacCormack_scheme.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/linear_interpolation_3d.o: ./src/linear_interpolation_3d.cpp ./include/linear_interpolation_3d.h
	$(CXX) $(CXXFLAGS) -c ./src/linear_interpolation_3d.cpp -o $(OBJDIR)/linear_interpolation_3d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/linear_interpolation_3d.cpp -o $(OBJDIR)/linear_interpolation_3d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/WENO4_interpolation_3d.o: ./src/WENO4_interpolation_3d.cpp ./include/WENO4_interpolation_3d.h
	$(CXX) $(CXXFLAGS) -c ./src/WENO4_interpolation_3d.cpp -o $(OBJDIR)/WENO4_interpolation_3d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/WENO4_interpolation_3d.cpp -o $(OBJDIR)/WENO4_interpolation_3d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/WENO6_interpolation_3d.o: ./src/WENO6_interpolation_3d.cpp ./include/WENO6_interpolation_3d.h
	$(CXX) $(CXXFLAGS) -c ./src/WENO6_interpolation_3d.cpp -o $(OBJDIR)/WENO6_interpolation_3d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/WENO6_interpolation_3d.cpp -o $(OBJDIR)/WENO6_interpolation_3d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/linear_interpolation_1d.o: ./src/linear_interpolation_1d.cpp ./include/linear_interpolation_1d.h
	$(CXX) $(CXXFLAGS) -c ./src/linear_interpolation_1d.cpp -o $(OBJDIR)/linear_interpolation_1d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/linear_interpolation_1d.cpp -o $(OBJDIR)/linear_interpolation_1d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/WENO6_interpolation_1d.o: ./src/WENO6_interpolation_1d.cpp ./include/WENO6_interpolation_1d.h
	$(CXX) $(CXXFLAGS) -c ./src/WENO6_interpolation_1d.cpp -o $(OBJDIR)/WENO6_interpolation_1d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/WENO6_interpolation_1d.cpp -o $(OBJDIR)/WENO6_interpolation_1d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/backtrace_and_calc_all_quadrature_point_list_for_integrate_density.o: ./src/backtrace_and_calc_all_quadrature_point_list_for_integrate_density.cpp ./include/backtrace_and_calc_all_quadrature_point_list_for_integrate_density.h
	$(CXX) $(CXXFLAGS) -c ./src/backtrace_and_calc_all_quadrature_point_list_for_integrate_density.cpp -o $(OBJDIR)/backtrace_and_calc_all_quadrature_point_list_for_integrate_density.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/backtrace_and_calc_all_quadrature_point_list_for_integrate_density.cpp -o $(OBJDIR)/backtrace_and_calc_all_quadrature_point_list_for_integrate_density.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/calc_density_in_cell_from_quadrature_point_list_for_integrate_density.o: ./src/calc_density_in_cell_from_quadrature_point_list_for_integrate_density.cpp ./include/calc_density_in_cell_from_quadrature_point_list_for_integrate_density.h
	$(CXX) $(CXXFLAGS) -c ./src/calc_density_in_cell_from_quadrature_point_list_for_integrate_density.cpp -o $(OBJDIR)/calc_density_in_cell_from_quadrature_point_list_for_integrate_density.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/calc_density_in_cell_from_quadrature_point_list_for_integrate_density.cpp -o $(OBJDIR)/calc_density_in_cell_from_quadrature_point_list_for_integrate_density.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/calc_cell_volumes.o: ./src/calc_cell_volumes.cpp ./include/calc_cell_volumes.h
	$(CXX) $(CXXFLAGS) -c ./src/calc_cell_volumes.cpp -o $(OBJDIR)/calc_cell_volumes.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/calc_cell_volumes.cpp -o $(OBJDIR)/calc_cell_volumes.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/calc_backtraced_face.o: ./src/calc_backtraced_face.cpp ./include/calc_backtraced_face.h
	$(CXX) $(CXXFLAGS) -c ./src/calc_backtraced_face.cpp -o $(OBJDIR)/calc_backtraced_face.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/calc_backtraced_face.cpp -o $(OBJDIR)/calc_backtraced_face.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/split_face.o: ./src/split_face.cpp ./include/split_face.h
	$(CXX) $(CXXFLAGS) -c ./src/split_face.cpp -o $(OBJDIR)/split_face.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/split_face.cpp -o $(OBJDIR)/split_face.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/linear_interpolation_2d.o: ./src/linear_interpolation_2d.cpp ./include/linear_interpolation_2d.h
	$(CXX) $(CXXFLAGS) -c ./src/linear_interpolation_2d.cpp -o $(OBJDIR)/linear_interpolation_2d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/linear_interpolation_2d.cpp -o $(OBJDIR)/linear_interpolation_2d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/WENO6_interpolation_2d.o: ./src/WENO6_interpolation_2d.cpp ./include/WENO6_interpolation_2d.h
	$(CXX) $(CXXFLAGS) -c ./src/WENO6_interpolation_2d.cpp -o $(OBJDIR)/WENO6_interpolation_2d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/WENO6_interpolation_2d.cpp -o $(OBJDIR)/WENO6_interpolation_2d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/calc_pressure.o: ./src/calc_pressure.cpp ./include/calc_pressure.h
	$(CXX) $(CXXFLAGS) -c ./src/calc_pressure.cpp -o $(OBJDIR)/calc_pressure.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/calc_pressure.cpp -o $(OBJDIR)/calc_pressure.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/calc_psi_3D.o: ./src/calc_psi_3D.cpp ./include/calc_psi_3D.h
	$(CXX) $(CXXFLAGS) -c ./src/calc_psi_3D.cpp -o $(OBJDIR)/calc_psi_3D.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/calc_psi_3D.cpp -o $(OBJDIR)/calc_psi_3D.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/calc_velocity_from_quadrature_point_list.o: ./src/calc_velocity_from_quadrature_point_list.cpp ./include/calc_velocity_from_quadrature_point_list.h
	$(CXX) $(CXXFLAGS) -c ./src/calc_velocity_from_quadrature_point_list.cpp -o $(OBJDIR)/calc_velocity_from_quadrature_point_list.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/calc_velocity_from_quadrature_point_list.cpp -o $(OBJDIR)/calc_velocity_from_quadrature_point_list.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/calc_velocity_from_quadrature_point_list_use_integral.o: ./src/calc_velocity_from_quadrature_point_list_use_integral.cpp ./include/calc_velocity_from_quadrature_point_list_use_integral.h
	$(CXX) $(CXXFLAGS) -c ./src/calc_velocity_from_quadrature_point_list_use_integral.cpp -o $(OBJDIR)/calc_velocity_from_quadrature_point_list_use_integral.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/calc_velocity_from_quadrature_point_list_use_integral.cpp -o $(OBJDIR)/calc_velocity_from_quadrature_point_list_use_integral.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/advect_velocity_semi_lagrangian_3d.o: ./src/advect_velocity_semi_lagrangian_3d.cpp ./include/advect_velocity_semi_lagrangian_3d.h
	$(CXX) $(CXXFLAGS) -c ./src/advect_velocity_semi_lagrangian_3d.cpp -o $(OBJDIR)/advect_velocity_semi_lagrangian_3d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/advect_velocity_semi_lagrangian_3d.cpp -o $(OBJDIR)/advect_velocity_semi_lagrangian_3d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/calc_pressure_3d.o: ./src/calc_pressure_3d.cpp ./include/calc_pressure_3d.h
	$(CXX) $(CXXFLAGS) -c ./src/calc_pressure_3d.cpp -o $(OBJDIR)/calc_pressure_3d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/calc_pressure_3d.cpp -o $(OBJDIR)/calc_pressure_3d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/advect_velocity_flux_advection_3d.o: ./src/advect_velocity_flux_advection_3d.cpp ./include/advect_velocity_flux_advection_3d.h
	$(CXX) $(CXXFLAGS) -c ./src/advect_velocity_flux_advection_3d.cpp -o $(OBJDIR)/advect_velocity_flux_advection_3d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/advect_velocity_flux_advection_3d.cpp -o $(OBJDIR)/advect_velocity_flux_advection_3d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/parallelize_functions.o: ./src/parallelize_functions.cpp ./include/parallelize_functions.h
	$(CXX) $(CXXFLAGS) -c ./src/parallelize_functions.cpp -o $(OBJDIR)/parallelize_functions.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/parallelize_functions.cpp -o $(OBJDIR)/parallelize_functions.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/grid_1d.o: ./src/grid_1d.cpp ./include/grid_1d.h
	$(CXX) $(CXXFLAGS) -c ./src/grid_1d.cpp -o $(OBJDIR)/grid_1d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/grid_1d.cpp -o $(OBJDIR)/grid_1d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/move_substances_1d.o: ./src/move_substances_1d.cpp ./include/move_substances_1d.h
	$(CXX) $(CXXFLAGS) -c ./src/move_substances_1d.cpp -o $(OBJDIR)/move_substances_1d.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/move_substances_1d.cpp -o $(OBJDIR)/move_substances_1d.o $(LDFLAGS) $(LIBS) $(INCLUDE)

$(OBJDIR)/calc_psi_1D.o: ./src/calc_psi_1D.cpp ./include/calc_psi_1D.h
	$(CXX) $(CXXFLAGS) -c ./src/calc_psi_1D.cpp -o $(OBJDIR)/calc_psi_1D.o $(LIBS) $(INCLUDE)
#	$(CXX) -c ./src/calc_psi_1D.cpp -o $(OBJDIR)/calc_psi_1D.o $(LDFLAGS) $(LIBS) $(INCLUDE)

clean:
	rm -f $(OBJDIR)/*.o $(TARGET)
