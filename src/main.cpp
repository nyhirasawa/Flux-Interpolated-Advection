// xvfb-run --auto-servernum --server-args="-screen 0 1920x1080x24" ./build/a.out --dimension=2 --grid_num_x=128 --grid_num_y=128 --num_movie_frames=60 --cfl=3.0 --interpolation_method=linear --num_gauss_quad_boundary=1 --num_gauss_quad_bulk=1 --split_method=y-CellFaceAligned --integral_method=gauss

//#include <opencv2/opencv.hpp>
#include <GLFW/glfw3.h>
#include <GL/glut.h>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <save_image.h>
#include <cxxopts.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

#include "calc_time_step_length_from_CFL_number.h"
#include "grid.h"
#include "grid_3d.h"
#include "grid_1d.h"
#include "initialize_grid.h"
#include "physical_const.h"
#include "draw_substance_density_opengl.h"
#include "update_fluid_velocity.h"
#include "update_fluid_velocity_3d.h"
#include "move_substances.h"
#include "move_substances_3d.h"
#include "move_substances_1d.h"
#include "utils.h"
#include "gauss_quadrature_points.h"
#include "vdbexport.h"
#include "write_substance_density_data.h"
#include "split_face_3D.h"
#include "define_float_type.h"
#include "sginterp3.h"
#include "calc_cell_volumes.h"

int main(int argc, const char* argv[]) {
	//
	// Define arguments
	cxxopts::Options options("SampleProgram", "An example of OpenGL program");
	options.add_options()
		("d,dimension", "dimension", cxxopts::value<int>()->default_value("3"))
		("use_flux_advection", "use flux advection(density)", cxxopts::value<bool>()->default_value("true"))
		("use_flux_advection_velocity", "use flux advection(velocity)", cxxopts::value<bool>()->default_value("true"))
		("use_lentines_advection", "use Lentines advection(density)", cxxopts::value<bool>()->default_value("false"))
		("c,cfl", "cfl number for determine timestep length", cxxopts::value<MY_FLOAT_TYPE>()->default_value("1.0"))
		("grid_num_x", "number of grid cell (x direction)", cxxopts::value<int>()->default_value("64"))
		("grid_num_y", "number of grid cell (y direction)", cxxopts::value<int>()->default_value("64"))
		("grid_num_z", "number of grid cell (z direction)", cxxopts::value<int>()->default_value("64"))
		("f,fps", "fps of movie", cxxopts::value<MY_FLOAT_TYPE>()->default_value("60.0"))
		("num_movie_frames", "number of movie frames", cxxopts::value<int>()->default_value("600"))
		("i,interpolation_method", "interpolation method(linear or WENO6)", cxxopts::value<std::string>())
		("color_negative_density", "color negative density(true of false)", cxxopts::value<bool>()->default_value("true"))
		("set_zero_normal_velocity_at_boundary", "set zero normal velocity at boundary(false: allow normal velocity at boundary, true: zero normal velocity at bondary)", cxxopts::value<bool>()->default_value("true"))
		("fix_velocity", "fix all velocity field(false: evolve velocity field, true: ignore e.o.m(use initial velocity field at all frame))", cxxopts::value<bool>()->default_value("false"))
		("use_integral", "use integral instead interpolation for calc psi", cxxopts::value<bool>()->default_value("true"))
		("use_MacCormack_scheme", "use MacCormack method(false: use single semi-lagrangian, true: use MacCormack method)", cxxopts::value<bool>()->default_value("false"))
		("use_clamping_in_MacCormack_scheme", "use clamping in MacCormack method", cxxopts::value<bool>()->default_value("true"))
		("num_gauss_quad_boundary", "number of gauss quadrature points on boundary", cxxopts::value<int>()->default_value("1"))
		("num_gauss_quad_bulk", "number of gauss quadrature points on bulk", cxxopts::value<int>()->default_value("1"))
		("num_gauss_quadrature_point_for_integrate_density", "number of gauss quadrature points for integrate density", cxxopts::value<int>()->default_value("1"))
		("enable_eliminate_zero_velocity_false_diffusion", "enable eliminate zero velocity false diffusion", cxxopts::value<bool>()->default_value("true"))
		("enable_cell_volume_correction", "enable cell volume correction", cxxopts::value<bool>()->default_value("false"))
		("o,output", "Output directory path", cxxopts::value<std::string>()->default_value("img"))
		("total_mass_file_path", "total_mass_file_path", cxxopts::value<std::string>()->default_value("./dat/"))
//		("total_mass_filename", "total_mass_filename", cxxopts::value<std::string>()->default_value("total_mass.dat"))
//		("total_momentum_x_filename", "total_momentum_x_filename", cxxopts::value<std::string>()->default_value("total_momentum_x.dat"))
//		("total_momentum_y_filename", "total_momentum_y_filename", cxxopts::value<std::string>()->default_value("total_momentum_y.dat"))
//		("total_momentum_z_filename", "total_momentum_z_filename", cxxopts::value<std::string>()->default_value("total_momentum_z.dat"))
		("density_data_path", "density_data_path", cxxopts::value<std::string>()->default_value("./density_data"))
		("split_method", "face split method", cxxopts::value<std::string>()->default_value("y-CellFaceAligned"))
		("integral_method", "integral method(analytical or gauss)", cxxopts::value<std::string>()->default_value("gauss"))
		("prefix_of_data_files", "prefix_of_data_files", cxxopts::value<std::string>()->default_value("data"))
		("minimal_area", "minimal_area of triangle split", cxxopts::value<MY_FLOAT_TYPE>()->default_value("1.0"))
		("h,help", "Print usage");
	//

	// Parse arguments
	auto params = options.parse(argc,argv);
	if (params.count("help")) {
		std::cout << options.help() << std::endl;
		return 0;
	}

	////パースしたコマンドライン引数を変数に代入していく
	const int dimension = params["dimension"].as<int>();
	const bool use_flux_advection = params["use_flux_advection"].as<bool>();
	const bool use_flux_advection_velocity = params["use_flux_advection_velocity"].as<bool>();
	const bool use_lentines_advection = params["use_lentines_advection"].as<bool>();
	const MY_FLOAT_TYPE movie_fps = params["fps"].as<MY_FLOAT_TYPE>();
	const std::string img_file_path = params["output"].as<std::string>();
    const MY_FLOAT_TYPE CFL_number = params["cfl"].as<MY_FLOAT_TYPE>();
	const int grid_num_x = params["grid_num_x"].as<int>();
	const int grid_num_y = params["grid_num_y"].as<int>();
	const int grid_num_z = params["grid_num_z"].as<int>();
	//セルの一辺の長さ
	const MY_FLOAT_TYPE cell_length = 1.0 / grid_num_x;
	//セルの体積
	const MY_FLOAT_TYPE cell_volume = cell_length * cell_length;
	//補間の方法
	const std::string interpolation_method = params["interpolation_method"].as<std::string>();
	// 動画のフレーム数
	const int num_movie_frames = params["num_movie_frames"].as<int>();
	// 負の密度を色付けするかどうか
	const bool color_negative_density = params["color_negative_density"].as<bool>();
	// 境界上で速度の法線成分を0と置くかどうか
	const bool set_zero_normal_velocity_at_boundary = params["set_zero_normal_velocity_at_boundary"].as<bool>();
	// 固定した速度場を使うかどうか(trueの場合速度場は時間発展しない。また、速度場は void initialize_grid(Grid& all_grid) 関数内で定義する)
	const bool fix_velocity = params["fix_velocity"].as<bool>();
	// use integral instead interpolation for calc psi
	const bool use_integral = params["use_integral"].as<bool>();
	// advect に MacCormack を使うかどうか
	const bool use_MacCormack_scheme = params["use_MacCormack_scheme"].as<bool>();
	//MacCormack 法でclamping を行うかどうか
	const bool use_clamping_in_MacCormack_scheme = params["use_clamping_in_MacCormack_scheme"].as<bool>();
	// 境界でのガウス求積の評価点数
	const int num_gauss_quad_boundary = params["num_gauss_quad_boundary"].as<int>();
	// バルクでのガウス求積の評価点数
	const int num_gauss_quad_bulk = params["num_gauss_quad_bulk"].as<int>();
	// number of gauss quadrature points for integrate density
	const int num_gauss_quadrature_point_for_integrate_density = params["num_gauss_quadrature_point_for_integrate_density"].as<int>();
	// zero velocity false diffusion を除去するかどうか
	const bool enable_eliminate_zero_velocity_false_diffusion = params["enable_eliminate_zero_velocity_false_diffusion"].as<bool>();
	// volume correction を行うかどうか
	const bool enable_cell_volume_correction = params["enable_cell_volume_correction"].as<bool>();
	//total mass を記録するファイルを格納するパス
	const std::string total_mass_file_path = params["total_mass_file_path"].as<std::string>();
	//total mass を記録するファイル名
//	const std::string total_mass_filename = params["total_mass_filename"].as<std::string>();
	//total mass を記録するファイル名
//	const std::string total_momentum_x_filename = params["total_momentum_x_filename"].as<std::string>();
	//total mass を記録するファイル名
//	const std::string total_momentum_y_filename = params["total_momentum_y_filename"].as<std::string>();
	//total mass を記録するファイル名
//	const std::string total_momentum_z_filename = params["total_momentum_z_filename"].as<std::string>();
	//density を記録するファイル
	const std::string density_data_path = params["density_data_path"].as<std::string>();
	//面の切断の方法
	const std::string split_method = params["split_method"].as<std::string>();
	//積分の方法(analytical or gauss)
	const std::string integral_method = params["integral_method"].as<std::string>();
	//データを記録するファイルの接頭辞
	const std::string prefix_of_data_files = params["prefix_of_data_files"].as<std::string>();
	const MY_FLOAT_TYPE minimal_area = params["minimal_area"].as<MY_FLOAT_TYPE>();

/////////////////////

	//total mass を記録するファイル
	//	std::string total_mass_filename = "./dat/total_mass_tmp.dat";
//	system(("mkdir -p "+total_mass_file_path+"; rm "+total_mass_file_path+total_mass_filename).c_str());
	system(("mkdir -p "+total_mass_file_path+"; rm "+total_mass_file_path+prefix_of_data_files+"total_mass.dat").c_str());
	system(("touch "+total_mass_file_path+prefix_of_data_files+"total_mass.dat"+";").c_str());
	std::ofstream total_mass_writing_file;
	total_mass_writing_file.open(total_mass_file_path+prefix_of_data_files+"total_mass.dat", std::ios::out);

	//total momentum x を記録するファイル
	system(("mkdir -p "+total_mass_file_path+"; rm "+total_mass_file_path+prefix_of_data_files+"total_momentum_x.dat").c_str());
	system(("touch "+total_mass_file_path+prefix_of_data_files+"total_momentum_x.dat"+";").c_str());
	std::ofstream total_momentum_x_writing_file;
	total_momentum_x_writing_file.open(total_mass_file_path+prefix_of_data_files+"total_momentum_x.dat", std::ios::out);
	//total momentum y を記録するファイル
	system(("mkdir -p "+total_mass_file_path+"; rm "+total_mass_file_path+prefix_of_data_files+"total_momentum_y.dat").c_str());
	system(("touch "+total_mass_file_path+prefix_of_data_files+"total_momentum_y.dat"+";").c_str());
	std::ofstream total_momentum_y_writing_file;
	total_momentum_y_writing_file.open(total_mass_file_path+prefix_of_data_files+"total_momentum_y.dat", std::ios::out);
	//total momentum z を記録するファイル
	system(("mkdir -p "+total_mass_file_path+"; rm "+total_mass_file_path+prefix_of_data_files+"total_momentum_z.dat").c_str());
	system(("touch "+total_mass_file_path+prefix_of_data_files+"total_momentum_z.dat"+";").c_str());
	std::ofstream total_momentum_z_writing_file;
	total_momentum_z_writing_file.open(total_mass_file_path+prefix_of_data_files+"total_momentum_z.dat", std::ios::out);

	//total_loss_of_momentum_x を記録するファイル
	system(("mkdir -p "+total_mass_file_path+"; rm "+total_mass_file_path+prefix_of_data_files+"total_loss_of_momentum_x.dat").c_str());
	system(("touch "+total_mass_file_path+prefix_of_data_files+"total_loss_of_momentum_x.dat"+";").c_str());
	std::ofstream total_loss_of_momentum_x_writing_file;
	total_loss_of_momentum_x_writing_file.open(total_mass_file_path+prefix_of_data_files+"total_loss_of_momentum_x.dat", std::ios::out);
	//total_loss_of_momentum_y を記録するファイル
	system(("mkdir -p "+total_mass_file_path+"; rm "+total_mass_file_path+prefix_of_data_files+"total_loss_of_momentum_y.dat").c_str());
	system(("touch "+total_mass_file_path+prefix_of_data_files+"total_loss_of_momentum_y.dat"+";").c_str());
	std::ofstream total_loss_of_momentum_y_writing_file;
	total_loss_of_momentum_y_writing_file.open(total_mass_file_path+prefix_of_data_files+"total_loss_of_momentum_y.dat", std::ios::out);
	//total_loss_of_momentum_z を記録するファイル
	system(("mkdir -p "+total_mass_file_path+"; rm "+total_mass_file_path+prefix_of_data_files+"total_loss_of_momentum_z.dat").c_str());
	system(("touch "+total_mass_file_path+prefix_of_data_files+"total_loss_of_momentum_z.dat"+";").c_str());
	std::ofstream total_loss_of_momentum_z_writing_file;
	total_loss_of_momentum_z_writing_file.open(total_mass_file_path+prefix_of_data_files+"total_loss_of_momentum_z.dat", std::ios::out);

	//momentum_x + total_loss_of_momentum_x を記録するファイル
	system(("mkdir -p "+total_mass_file_path+"; rm "+total_mass_file_path+prefix_of_data_files+"momentum_x_plus_total_loss_of_momentum_x.dat").c_str());
	system(("touch "+total_mass_file_path+prefix_of_data_files+"momentum_x_plus_total_loss_of_momentum_x.dat"+";").c_str());
	std::ofstream momentum_x_plus_total_loss_of_momentum_x_writing_file;
	momentum_x_plus_total_loss_of_momentum_x_writing_file.open(total_mass_file_path+prefix_of_data_files+"momentum_x_plus_total_loss_of_momentum_x.dat", std::ios::out);
	//momentum_y + total_loss_of_momentum_y を記録するファイル
	system(("mkdir -p "+total_mass_file_path+"; rm "+total_mass_file_path+prefix_of_data_files+"momentum_y_plus_total_loss_of_momentum_y.dat").c_str());
	system(("touch "+total_mass_file_path+prefix_of_data_files+"momentum_y_plus_total_loss_of_momentum_y.dat"+";").c_str());
	std::ofstream momentum_y_plus_total_loss_of_momentum_y_writing_file;
	momentum_y_plus_total_loss_of_momentum_y_writing_file.open(total_mass_file_path+prefix_of_data_files+"momentum_y_plus_total_loss_of_momentum_y.dat", std::ios::out);
	//momentum_z + total_loss_of_momentum_z を記録するファイル
	system(("mkdir -p "+total_mass_file_path+"; rm "+total_mass_file_path+prefix_of_data_files+"momentum_z_plus_total_loss_of_momentum_z.dat").c_str());
	system(("touch "+total_mass_file_path+prefix_of_data_files+"momentum_z_plus_total_loss_of_momentum_z.dat"+";").c_str());
	std::ofstream momentum_z_plus_total_loss_of_momentum_z_writing_file;
	momentum_z_plus_total_loss_of_momentum_z_writing_file.open(total_mass_file_path+prefix_of_data_files+"momentum_z_plus_total_loss_of_momentum_z.dat", std::ios::out);

	if(dimension==1){
		//系の情報が乗ったグリッド
		smoke_simulation::Grid_1D all_grid(grid_num_y, cell_length);
		//グリッドを初期化
		smoke_simulation::initialize_grid(all_grid);
		////メインの計算部分
		//動画ファイルの何フレーム目までを描画したか数える変数
		int i_movie_frame = 0;
		//経過時刻(現実時刻ではなく、シミュレーション中の時刻)を記録する変数
		MY_FLOAT_TYPE now_time = 0.0;

		int i_frame = 0;

		while (true) {
			std::cout << "simulaition frame: " << i_frame << std::endl;
			//time step の長さを計算
			MY_FLOAT_TYPE time_step_length;
			bool write_frame_to_movie_file = false;
//			std::cout<<"calc_time_step_length_from_CFL_number_1D() begin"<<std::endl;
//			getchar();
			time_step_length
				= smoke_simulation::calc_time_step_length_from_CFL_number_1D(
					all_grid,
					now_time,
					i_movie_frame,
					write_frame_to_movie_file,
					CFL_number,
					movie_fps
				);
//			std::cout<<"calc_time_step_length_from_CFL_number_1D() end"<<std::endl;
			std::cout<<"time_step_length: "<<time_step_length<<std::endl;

			//動画に必要なフレーム数シミュレーションが走ったら終了
			if (i_movie_frame > num_movie_frames) {
				break;
			}

//			std::cout<<"write density data begin"<<std::endl;
			if(write_frame_to_movie_file || i_frame == 0){
				std::stringstream ss;
				ss << "density_" << i_movie_frame << ".dat";
				std::ofstream density_writing_file;
				density_writing_file.open(ss.str(), std::ios::out);
				for(int iy = 0; iy < all_grid.Grid_num_y; ++iy){
					density_writing_file << iy * all_grid._cell_length << " " << all_grid.substance_density[iy] << std::endl;
				}
				++i_movie_frame;
			}
//			std::cout<<"write density data end"<<std::endl;

			//substance densityの更新
//			std::cout<<"move_substances_1D() begin"<<std::endl;
			smoke_simulation::move_substances_1D(
				all_grid,
				i_frame,
				time_step_length,
				use_flux_advection,
				use_MacCormack_scheme,
				use_clamping_in_MacCormack_scheme,
				num_gauss_quad_boundary,
				num_gauss_quad_bulk,
				enable_eliminate_zero_velocity_false_diffusion,
				split_method,
				integral_method,
				interpolation_method,
				use_integral,
				num_gauss_quadrature_point_for_integrate_density,
				enable_cell_volume_correction,
				use_lentines_advection
			);
//			std::cout<<"move_substances end"<<std::endl;

			//シミュレーション中の時刻を発展させる
			now_time += time_step_length;

			//フレームのカウントを増やす
			i_frame++;
		}
	}
	else if(dimension==2){
		//系の情報が乗ったグリッド
		smoke_simulation::Grid all_grid(grid_num_x, grid_num_y, cell_length, interpolation_method);
		//グリッドを初期化
		smoke_simulation::initialize_grid(all_grid, fix_velocity, true, true);

		//
		// Initialize GLFW
		assert(glfwInit());
		//
		// Create window
		const int scale = 8;
		int width = all_grid.Grid_num_x * scale;
		int height = all_grid.Grid_num_y * scale;
		GLFWwindow* window = glfwCreateWindow(width, height, "Simple example", nullptr, nullptr);
		assert(window);
		glfwMakeContextCurrent(window);
		//
		// Set clear color
		glClearColor(0.0,0.0,0.0,1.0);
		//
		// Make screenshot directory and remove existing images/video if exist
		system(("mkdir -p "+img_file_path+"; rm "+img_file_path+"/*.png "+img_file_path+"/out.mp4").c_str());

		////メインの計算部分
		//動画ファイルの何フレーム目までを描画したか数える変数
		int i_movie_frame = 0;
		//経過時刻(現実時刻ではなく、シミュレーション中の時刻)を記録する変数
		MY_FLOAT_TYPE now_time = 0.0;
		//0フレーム目でのpsiの最大値 (heat map描画用の変数)
		MY_FLOAT_TYPE max_psi_at_initial_frame = 0.0;
		int i_frame = 0;
		while (!glfwWindowShouldClose(window)) {
			//質量の総和を計算
			MY_FLOAT_TYPE total_sum = 0.0;
			MY_FLOAT_TYPE max_density = 0.0;
			int max_density_ix, max_density_iy;
			for (int ix = 0; ix < all_grid.Grid_num_x; ++ix) {
				for (int iy = 0; iy < all_grid.Grid_num_y; ++iy) {
					total_sum
						+= all_grid.substance_density[smoke_simulation::get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
//						* all_grid.cell_volume_cell_center[smoke_simulation::get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
						* all_grid._cell_volume;
/*
					if (max_density < all_grid.substance_density[smoke_simulation::get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]) {
						max_density = all_grid.substance_density[smoke_simulation::get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
						max_density_ix = ix;
						max_density_iy = iy;
					}
*/
				}
			}
			//ファイル書き込み
			std::cout << "writing " << prefix_of_data_files+"total_mass.dat" << "..." << std::endl;
			total_mass_writing_file << std::setprecision(16) << now_time << " " << total_sum << std::endl;
			//total_mass_writing_file << now_time << " " << max_density << " " << max_density_ix << " " << max_density_iy << std::endl;

			std::cout << "simulaition frame: " << i_frame << std::endl;
			//time step の長さを計算
			MY_FLOAT_TYPE time_step_length;
			bool write_frame_to_movie_file = false;

			if (smoke_simulation::physical_const::kUse_variable_time_step_length) {
				time_step_length
				= smoke_simulation::calc_time_step_length_from_CFL_number(
					all_grid,
					now_time,
					i_movie_frame,
					write_frame_to_movie_file,
					CFL_number,
					movie_fps
				);
			}
			else {
				time_step_length
				= smoke_simulation::calc_time_step_length_from_fixed_time_step_length(
					now_time,
					i_movie_frame,
					write_frame_to_movie_file,
					movie_fps
				);
			}
			std::cout<<"time_step_length: "<<time_step_length<<std::endl;

			//substance densityの描画
			smoke_simulation::draw_substance_density_opengl(
				window,
				img_file_path,
				all_grid,
				i_frame,
				i_movie_frame,
				write_frame_to_movie_file,
				time_step_length,
				color_negative_density
			);

			//smoke_simulation::draw_substance_density(all_grid, scale, writer, i_frame, i_movie_frame, write_frame_to_movie_file, time_step_length);
			//smoke_simulation::draw_test(window, img_file_path, i_frame);

			////psi density のヒートマップの描画
	/*
			//0フレーム目でのpsiの最大値を計算
			if (i_frame == 0) {
				//cell face での psiを計算する
				all_grid.calc_psi_substance_density_cell_face();
				for (int ix = 0; ix < all_grid.Grid_num_x; ++ix) {
					for (int iy = 0; iy < all_grid.Grid_num_y + 1; ++iy) {
						max_psi_at_initial_frame
							= std::max(max_psi_at_initial_frame, all_grid.psi_substance_density_cell_face_y[smoke_simulation::get_voxel_face_index_y(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]);
					}
				}
			}
			//psi density のヒートマップの描画
			smoke_simulation::draw_psi_density_heatmap(all_grid, scale, writer, i_frame, i_movie_frame, write_frame_to_movie_file, max_psi_at_initial_frame);
	*/

			//動画に必要なフレーム数シミュレーションが走ったら終了
			if (i_movie_frame > num_movie_frames) {
				break;
			}

			//流体の速度場の更新
//			if (!fix_velocity) {
				smoke_simulation::update_fluid_velocity(
					all_grid,
					i_frame,
					time_step_length,
					use_flux_advection_velocity,
					use_MacCormack_scheme,
					use_clamping_in_MacCormack_scheme,
					num_gauss_quad_boundary,
					num_gauss_quad_bulk,
					enable_eliminate_zero_velocity_false_diffusion,
					split_method,
					integral_method,
					interpolation_method,
					enable_cell_volume_correction,
					fix_velocity,
					use_integral,
					num_gauss_quadrature_point_for_integrate_density
				);
				std::cout<<"update_fluid_velocity end"<<std::endl;
				std::cout<<"getchar()"<<std::endl;
//				getchar();

//			}

//			std::cout<<"move_substances() begin"<<std::endl;
			//substance densityの更新
			smoke_simulation::move_substances(
				all_grid,
				i_frame,
				time_step_length,
				use_flux_advection,
//				set_zero_normal_velocity_at_boundary,
				use_MacCormack_scheme,
				use_clamping_in_MacCormack_scheme,
				num_gauss_quad_boundary,
				num_gauss_quad_bulk,
				enable_eliminate_zero_velocity_false_diffusion,
				split_method,
				integral_method,
				interpolation_method,
				use_integral,
				num_gauss_quadrature_point_for_integrate_density,
				enable_cell_volume_correction,
				use_lentines_advection
			);
			std::cout<<"move_substances end"<<std::endl;

			std::cout<<"getchar()"<<std::endl;
//			getchar();


//			std::cout<<"move_substances() end"<<std::endl;
/*
			smoke_simulation::correct_volume_concentration_error(
				all_grid,
				time_step_length
			);
*/
			//シミュレーション中の時刻を発展させる
			now_time += time_step_length;

			//フレームのカウントを増やす
			i_frame++;
			//getchar();
		}
		system(("ffmpeg -y -r "+std::to_string(movie_fps)+" -i "+img_file_path+"/screenshot_%d.png -pix_fmt yuv420p -b:v 12000k "+img_file_path+"/out.mp4").c_str());
	}
	else if(dimension==3){
		////時間計測開始( http://vivi.dyndns.org/tech/cpp/timeMeasurement.html 参照 )
		auto start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
		//系の情報が乗ったグリッド
		smoke_simulation::Grid_3D all_grid(grid_num_x, grid_num_y, grid_num_z, cell_length);
//		smoke_simulation::Grid_3D all_grid(grid_num_x, grid_num_y, grid_num_z, cell_length, interpolation_method);
		////時間計測終了
		auto end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
	    auto duration_time = end_time - start_time;        // 要した時間を計算
	    auto duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
	    // 要した時間をミリ秒（1/1000秒）に変換して表示
	    std::cout <<"Grid_3D constructor: "<< duration_time_msec << " milli sec \n";

		////時間計測開始
		start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
		//グリッドを初期化
		smoke_simulation::initialize_Grid_3D(all_grid, fix_velocity);

		////時間計測終了
		end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
	    duration_time = end_time - start_time;        // 要した時間を計算
	    duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
	    // 要した時間をミリ秒（1/1000秒）に変換して表示
	    std::cout <<"initialize_Grid_3D: "<< duration_time_msec << " milli sec \n";

		// Make density directory and remove existing density data if exist
	    system(("mkdir -p "+density_data_path+"; rm "+density_data_path+"/*.dat "+"; rm "+density_data_path+"/*.vdb ").c_str());
		////メインの計算部分
		//動画ファイルの何フレーム目までを描画したか数える変数
		int i_movie_frame = 0;
		//経過時刻(現実時刻ではなく、シミュレーション中の時刻)を記録する変数
		MY_FLOAT_TYPE now_time = 0.0;
		//シミュレーションのフレーム数
		int i_frame = 0;
		//シミュレーションのメインループを開始した時刻
		auto start_main_loop_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
		while (true) {
			std::cout << "simulaition frame: " << i_frame << std::endl;
			//質量の総和を計算
			MY_FLOAT_TYPE total_sum = 0.0;
			MY_FLOAT_TYPE max_density = 0.0;
			int max_density_ix, max_density_iy;
			for (int ix = 0; ix < all_grid.Grid_num_x; ++ix) {
				for (int iy = 0; iy < all_grid.Grid_num_y; ++iy) {
					for (int iz = 0; iz < all_grid.Grid_num_z; ++iz) {
						total_sum
							+= all_grid.substance_density[smoke_simulation::get_voxel_center_index_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
							* all_grid._cell_volume;
					}
				}
			}
			//ファイル書き込み
			std::cout << "writing " << prefix_of_data_files+"total_mass.dat" << "..." << std::endl;
			total_mass_writing_file << std::setprecision(16) << now_time << " " << total_sum << std::endl;

			//// momentum_x の総和を計算
			MY_FLOAT_TYPE total_momentum_x = 0.0;
			for (int ix = 0; ix < all_grid.Grid_num_x + 1; ++ix) {
				for (int iy = 0; iy < all_grid.Grid_num_y; ++iy) {
					for (int iz = 0; iz < all_grid.Grid_num_z; ++iz) {
						total_momentum_x
							+= all_grid.velocity_in_voxel_face_x[smoke_simulation::get_voxel_face_index_x_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
							* all_grid._cell_volume;
					}
				}
			}
			// momentum_x のファイル書き込み
			std::cout << "writing " << prefix_of_data_files+"total_momentum_x.dat" << "..." << std::endl;
			total_momentum_x_writing_file << std::setprecision(16) << now_time << " " << total_momentum_x << std::endl;
			// total_loss_of_momentum_x のファイル書き込み
			std::cout << "writing " << prefix_of_data_files+"total_loss_of_momentum_x.dat" << "..." << std::endl;
			total_loss_of_momentum_x_writing_file << std::setprecision(16) << now_time << " " << all_grid.total_loss_of_momentum_x << std::endl;
			//momentum_x + total_loss_of_momentum_x のファイル書き込み
			std::cout << "writing " << prefix_of_data_files+"momentum_x_plus_total_loss_of_momentum_x.dat" << "..." << std::endl;
			momentum_x_plus_total_loss_of_momentum_x_writing_file << std::setprecision(16) << now_time << " " << total_momentum_x + all_grid.total_loss_of_momentum_x << std::endl;

			//// momentum_y の総和を計算
			MY_FLOAT_TYPE total_momentum_y = 0.0;
			for (int ix = 0; ix < all_grid.Grid_num_x; ++ix) {
				for (int iy = 0; iy < all_grid.Grid_num_y + 1; ++iy) {
					for (int iz = 0; iz < all_grid.Grid_num_z; ++iz) {
						total_momentum_y
							+= all_grid.velocity_in_voxel_face_y[smoke_simulation::get_voxel_face_index_y_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
							* all_grid._cell_volume;
					}
				}
			}
			// momentum_y のファイル書き込み
			std::cout << "writing " << prefix_of_data_files+"total_momentum_y.dat" << "..." << std::endl;
			total_momentum_y_writing_file << std::setprecision(16) << now_time << " " << total_momentum_y << std::endl;
			// total_loss_of_momentum_y のファイル書き込み
			std::cout << "writing " << prefix_of_data_files+"total_loss_of_momentum_y.dat" << "..." << std::endl;
			total_loss_of_momentum_y_writing_file << std::setprecision(16) << now_time << " " << all_grid.total_loss_of_momentum_y << std::endl;
			//momentum_y + total_loss_of_momentum_y のファイル書き込み
			std::cout << "writing " << prefix_of_data_files+"momentum_y_plus_total_loss_of_momentum_y.dat" << "..." << std::endl;
			momentum_y_plus_total_loss_of_momentum_y_writing_file << std::setprecision(16) << now_time << " " << total_momentum_y + all_grid.total_loss_of_momentum_y << std::endl;

			//// momentum_z の総和を計算
			MY_FLOAT_TYPE total_momentum_z = 0.0;
			for (int ix = 0; ix < all_grid.Grid_num_x; ++ix) {
				for (int iy = 0; iy < all_grid.Grid_num_y; ++iy) {
					for (int iz = 0; iz < all_grid.Grid_num_z + 1; ++iz) {
						total_momentum_z
							+= all_grid.velocity_in_voxel_face_z[smoke_simulation::get_voxel_face_index_z_3D(ix, iy, iz, all_grid.Grid_num_x, all_grid.Grid_num_y, all_grid.Grid_num_z)]
							* all_grid._cell_volume;
					}
				}
			}
			// momentum_z のファイル書き込み
			std::cout << "writing " << prefix_of_data_files+"total_momentum_z.dat" << "..." << std::endl;
			total_momentum_z_writing_file << std::setprecision(16) << now_time << " " << total_momentum_z << std::endl;
			// total_loss_of_momentum_z のファイル書き込み
			std::cout << "writing " << prefix_of_data_files+"total_loss_of_momentum_z.dat" << "..." << std::endl;
			total_loss_of_momentum_z_writing_file << std::setprecision(16) << now_time << " " << all_grid.total_loss_of_momentum_z << std::endl;
			//momentum_z + total_loss_of_momentum_z のファイル書き込み
			std::cout << "writing " << prefix_of_data_files+"momentum_z_plus_total_loss_of_momentum_z.dat" << "..." << std::endl;
			momentum_z_plus_total_loss_of_momentum_z_writing_file << std::setprecision(16) << now_time << " " << total_momentum_z + all_grid.total_loss_of_momentum_z << std::endl;

			//time step の長さを計算
			MY_FLOAT_TYPE time_step_length;
			bool write_frame_to_movie_file = false;
			////時間計測開始
			start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
			if (smoke_simulation::physical_const::kUse_variable_time_step_length) {
				time_step_length
				= smoke_simulation::calc_time_step_length_from_CFL_number_3D(
					all_grid,
					now_time,
					i_movie_frame,
					write_frame_to_movie_file,
					CFL_number,
					movie_fps
				);
			}
			else {
				time_step_length
				= smoke_simulation::calc_time_step_length_from_fixed_time_step_length(
					now_time,
					i_movie_frame,
					write_frame_to_movie_file,
					movie_fps
				);
			}
			std::cout<<"time_step_length: "<<time_step_length<<std::endl;

			////時間計測終了
			end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
			duration_time = end_time - start_time;        // 要した時間を計算
			duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
			std::cout <<"calc_time_step_length_from_CFL_number_3D: "<< duration_time_msec << " milli sec \n";

			//substance densityをファイルに書き出す
			if(i_frame == 0 || write_frame_to_movie_file){
				std::cout << "write density data:  " << i_movie_frame <<"  frame"<< std::endl;;
				//自分のレイマーチングの実装
//				smoke_simulation::write_substance_density_data(all_grid, i_movie_frame, density_data_path);

				// vdb ファイル名
				std::ostringstream vdb_filename;
				vdb_filename<<density_data_path<<"/substance_density_"<<std::setfill('0') << std::right << std::setw(4)<<i_movie_frame<<".vdb"<<std::flush;
				// vdb ファイルの出力
				export2vdb(
					vdb_filename.str().c_str(),
					all_grid.Grid_num_x,
					all_grid.Grid_num_y,
					all_grid.Grid_num_z,
					-0.5 * all_grid._cell_length * all_grid.Grid_num_x,
					-0.5 * all_grid._cell_length * all_grid.Grid_num_y,
					-0.5 * all_grid._cell_length * all_grid.Grid_num_z,
					all_grid._cell_length,
					all_grid.substance_density.data()
				);
				++i_movie_frame;
			}
			//動画に必要なフレーム数シミュレーションが走ったら終了
			if (i_movie_frame > num_movie_frames) {
				break;
			}
			////時間計測開始
			start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
			//流体の速度場の更新
			if (!fix_velocity) {
				smoke_simulation::update_fluid_velocity_3D(
					all_grid,
					i_frame,
					time_step_length,
					use_flux_advection_velocity,
				    use_MacCormack_scheme,
				    use_clamping_in_MacCormack_scheme,
				    num_gauss_quad_boundary,
				    num_gauss_quad_bulk,
				    enable_eliminate_zero_velocity_false_diffusion,
				    split_method,
				    integral_method,
				    interpolation_method,
				    enable_cell_volume_correction,
				    fix_velocity,
					num_gauss_quadrature_point_for_integrate_density,
				    use_integral,
			        minimal_area
				);
			}
			////時間計測終了
			end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
			duration_time = end_time - start_time;        // 要した時間を計算
			duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
			std::cout <<"update_fluid_velocity_3D: "<< duration_time_msec << " milli sec \n";

			////時間計測開始
			start_time = std::chrono::system_clock::now();      // 計測スタート時刻を保存
			//substance densityの更新
			smoke_simulation::move_substances_3D(
				all_grid,
				i_frame,
				time_step_length,
				use_flux_advection,
//				set_zero_normal_velocity_at_boundary,
				use_MacCormack_scheme,
				use_clamping_in_MacCormack_scheme,
				num_gauss_quad_boundary,
				num_gauss_quad_bulk,
				enable_eliminate_zero_velocity_false_diffusion,
			    split_method,
				integral_method,
				num_gauss_quad_bulk,
				interpolation_method,
				use_integral,
			    num_gauss_quadrature_point_for_integrate_density,
				enable_cell_volume_correction,
				minimal_area
			);
			////時間計測終了
			end_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
			duration_time = end_time - start_time;        // 要した時間を計算
			duration_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_time).count();
			std::cout <<"move_substances_3D: "<< duration_time_msec << " milli sec \n";

			//シミュレーション中の時刻を発展させる
			now_time += time_step_length;
			//フレームのカウントを増やす
			i_frame++;
		}
		////時間計測終了
		auto end_main_loop_time = std::chrono::system_clock::now();       // 計測終了時刻を保存
		auto duration_main_loop_time = end_main_loop_time - start_main_loop_time;        // 要した時間を計算
		auto duration_main_loop_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(duration_main_loop_time).count();
		// 要した時間をミリ秒（1/1000秒）に変換して表示
		std::cout <<"duration_main_loop_time("<<i_frame<<" time steps): "<< duration_main_loop_time_msec << " milli sec \n";
	}
	return 0;
}
