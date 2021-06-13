#ifndef PHYSICAL_CONST_H
#define PHYSICAL_CONST_H

#include <string>

#include "define_float_type.h"

namespace smoke_simulation {
	struct physical_const {
		static constexpr MY_FLOAT_TYPE kPI = 3.14159265358979323846264338327950288;
		//x, y方向のグリッドセルの個数(壁は含まない)
//		static const int kGrid_num_x = 128;
//		static const int kGrid_num_y = 128;
		//static const int kGrid_num=(kGrid_num_x+2)*(kGrid_num_y+2);
		//セルの一辺の長さ
//		static constexpr MY_FLOAT_TYPE kCell_length = 1.0 / kGrid_num_x;
		//セルの体積
//		static constexpr MY_FLOAT_TYPE kCell_volume = kCell_length * kCell_length;
		//1時間ステップの長さ
		static constexpr MY_FLOAT_TYPE kDt_fixed = 1.0 / 100.0;
		//static constexpr MY_FLOAT_TYPE kDt = 1.0 / 50;
		//最大動画フレーム数
//		static const int kMax_num_movie_frames = 240;
		//最大シミュレーションフレーム数
//		static const int kMax_num_simulation_frames = 200000000;
		// frame per second of movie
//		static constexpr MY_FLOAT_TYPE kFps_of_movie = 24.0;

		//cubic interpolationを使うか
		//static const bool kUse_cubic_interpolation = true;
		//どの interpolation method を使うか("nearest" or "linear" or "cubic" or "monotone cubic")
//		static const std::string kInterpolation_method;

		//flux advection を使うか
//		static const bool kUse_flux_advection = false;

		// psi のどの定義を使うか("x" or "y" or "xy")(cpp ファイルの方に定義あり)
		static const std::string kPsi_definition;
		//負の質量密度を色付けするか
//		static const bool kColor_negative_density = true;
		//可変な時間ステップの大きさを使うか(true の場合、毎フレーム CFL 数から時間ステップの大きさを計算する)
		static const bool kUse_variable_time_step_length = true;
		//可変な時間ステップの大きさを使うときに、許容する最大のCFL数
//		static constexpr MY_FLOAT_TYPE kMax_CFL_number = 1.0;
		// 画像ファイルを書き出すかどうか
		static const bool kWrite_image_files = true;
		//境界を固定するかどうか(trueの場合、境界に接する面はバックトレースされない)
//		static const bool kFix_boundary_face = true;

		// 初期状態での密度の形状( "square" or "circle" )
		static const std::string kInitial_density_shape;
		//固定した速度場を使うかどうか(trueの場合速度場は時間発展しない。また、速度場は void initialize_grid(Grid& all_grid) 関数内で定義する)
//		static const bool kUse_fixed_velocity_field = false;

		//faceを等間隔に切断する際、いくつにfaceを分割するか
		static const int kSplit_num_uniformly = 64;
		//ガウス求積法の評価点数
//		static const int kNum_gauss_quadrature_points = 11;
		//バックトレースした後のグリッドを描画するかどうか(trueの場合描画するが、バグを取り切れている自信なし)
		static const bool kDraw_backtrace_grid = false;

		// advect に MacCormack を使うかどうか
//		static const bool kUse_MacCormack_scheme = false;
		//MacCormack 法でclamping を行うかどうか
//		static const bool kUse_clamping_in_MacCormack_scheme = true;
		//zero velocity false diffusion を除去するかどうか
//		static const bool kEliminate_zero_velocity_false_diffusion = true;
	};
} // namespace smoke_simulation

#endif
