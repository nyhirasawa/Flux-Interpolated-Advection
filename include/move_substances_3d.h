#ifndef MOVING_SUBSTANCES_3D_H
#define MOVING_SUBSTANCES_3D_H
#include <vector>
#include "grid_3d.h"
#include "cell_face_3D.h"
#include "cell_vertex_3D.h"
#include "physical_const.h"
#include "polygon_vector.h"
#include "define_float_type.h"
#include "gauss_quadrature_points.h"

namespace smoke_simulation{
    //流体の速度場によってsubstanceが運ばれる項(流体のadvect項に相当)
    //速度場を時間 -dt だけバックトレースして計算する
    void advect_density_semi_lagrangian_3D(
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const bool use_MacCormack_scheme,
        const bool use_clamping_in_MacCormack_scheme,
        const int num_gauss_quad_boundary,
        const int num_gauss_quad_bulk,
        const std::string interpolation_method
    );
    //
    //密度場の移流項を flux advection で計算
    void advect_density_flux_advection_3D_std_thread(
        Grid_3D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
//        const bool set_zero_normal_velocity_at_boundary,
        const bool use_MacCormack_scheme,
        const bool use_clamping_in_MacCormack_scheme,
        const int num_gauss_quad_boundary,
        const int num_gauss_quad_bulk,
        const bool eliminate_zero_velocity_false_diffusion,
        const std::string split_method,
        const std::string integral_method,
        const int num_gauss_quadrature_points,
        std::vector<MY_FLOAT_TYPE> &advected_values,
        const std::string interpolation_method,
        const std::string interpolation_method_in_calc_psi,
        const bool use_zero_velocity_for_backtrace,
        const bool use_integral,
	    const int num_gauss_quadrature_point_for_integrate_density,
		const bool enable_cell_volume_correction,
		const MY_FLOAT_TYPE minimal_area
    );
    //上の4ステップをまとめただけの関数(substance densityの1時間ステップ分の更新に相当)
    void move_substances_3D(
        Grid_3D& all_grid,
    	int i_frame,
    	const MY_FLOAT_TYPE time_step_length,
    	const bool use_flux_advection,
//    	const bool set_zero_normal_velocity_at_boundary,
    	const bool use_MacCormack_scheme,
    	const bool use_clamping_in_MacCormack_scheme,
    	const int num_gauss_quad_boundary,
    	const int num_gauss_quad_bulk,
    	const bool enable_eliminate_zero_velocity_false_diffusion,
        const std::string split_method,
        const std::string integral_method,
        const int num_gauss_quadrature_points,
        const std::string interpolation_method,
        bool use_integral,
	    const int num_gauss_quadrature_point_for_integrate_density,
		const bool enable_cell_volume_correction,
		const MY_FLOAT_TYPE minimal_area
    );
    // セルの体積によって質量密度を補正
    void correct_mass_and_volume_by_pressure_solve_3D(
        std::vector<MY_FLOAT_TYPE> &corrected_values,
        const std::vector<MY_FLOAT_TYPE> &cell_volume,
        const Grid_3D& all_grid
    );
//------------------------------------------------------------------------------
// 並列化関連
//------------------------------------------------------------------------------
//backtrace_and_split_cell_faces()をstd::thread で並列に実行する関数
void parallelize_backtrace_and_calc_all_quadrature_point_list(
    auto func,
    const int dim,
    size_t size,
    quadrature_point_vector &quadrature_point_list,
    bool use_MacCormack_scheme,
    const Grid_3D &all_grid,
    MY_FLOAT_TYPE time_step_length,
//    bool set_zero_normal_velocity_at_boundary,
    std::string split_method,
//    bool is_zero_velocity_false_diffusion_correction_mode,
    const std::string integral_method,
    const int num_gauss_quadrature_points,
    const std::string interpolation_method,
    const bool use_zero_velocity_for_backtrace,
    const MY_FLOAT_TYPE minimal_area
);
//calc_density_in_cell_from_quadrature_point_list()をstd::thread で並列に実行する関数
void parallelize_calc_density_in_cell_from_quadrature_point_list (
    auto func,
    const quadrature_point_vector &quadrature_point_list,
    const Grid_3D &all_grid,
    std::vector<MY_FLOAT_TYPE> &substance_density_after_advect,
    const std::string interpolation_method
);
//calc_cell_volumes_from_quadrature_point_list()をstd::thread で並列に実行する関数
void parallelize_calc_cell_volumes_from_quadrature_point_list (
    auto func,
    const quadrature_point_vector &quadrature_point_list,
    const Grid_3D &all_grid,
    std::vector<MY_FLOAT_TYPE> &cell_volume_after_advect
);
//backtrace_and_split_cell_faces()をstd::thread で並列に実行する関数
void parallelize_backtrace_and_calc_all_quadrature_point_list_for_integrate_density(
    auto func,
    const int dim,
    size_t size,
    quadrature_point_vector &quadrature_point_list_for_interpolate_density,
    quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
    quadrature_point_vector &quadrature_point_list_for_calc_cell_volume,
    bool use_MacCormack_scheme,
    const Grid_3D &all_grid,
    MY_FLOAT_TYPE time_step_length,
//    bool set_zero_normal_velocity_at_boundary,
    std::string split_method,
//    bool is_zero_velocity_false_diffusion_correction_mode,
    const std::string integral_method,
    const int num_gauss_quadrature_points,
    const std::string interpolation_method,
    const bool use_zero_velocity_for_backtrace,
    const int num_gauss_quadrature_point_for_integrate_density,
    const bool enable_cell_volume_correction,
    const std::string backtrace_usage,
    const MY_FLOAT_TYPE minimal_area
);
//calc_density_in_cell_from_quadrature_point_list()をstd::thread で並列に実行する関数
void parallelize_calc_density_in_cell_from_quadrature_point_list_for_integrate_density (
    auto func,
    const quadrature_point_vector &quadrature_point_list_for_interpolate_density,
    const quadrature_point_vector &quadrature_point_list_for_interpolate_psi,
    const Grid_3D &all_grid,
    std::vector<MY_FLOAT_TYPE> &substance_density_after_advect,
    const std::string interpolation_method,
    const std::vector<MY_FLOAT_TYPE> &advected_values
);

}//namespace smoke_simulation

#endif //MOVING_SUBSTANCES_3D_H
