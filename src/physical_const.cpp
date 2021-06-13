#include "physical_const.h"

#include "define_float_type.h"

namespace smoke_simulation {
	//どの interpolation method を使うか("nearest" or "linear" or "cubic" or "monotone cubic" or "WENO6")
//	const std::string physical_const::kInterpolation_method = "WENO6";
	// psi のどの定義を使うか(cpp ファイルの方に定義あり)
	const std::string physical_const::kPsi_definition = "y";
	// 初期状態での密度の形状( "square" or "circle" )
	const std::string physical_const::kInitial_density_shape = "square";
//	const std::string physical_const::kInitial_density_shape = "Tylor-Green";
//	const std::string physical_const::kInitial_density_shape = "vortex_sheet_2d";

}//namespace smoke_simulation
