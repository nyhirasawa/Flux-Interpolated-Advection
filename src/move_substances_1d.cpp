#include "move_substances_1d.h"

#include "calc_psi_1D.h"
#include <iostream>


namespace smoke_simulation {
    void move_substances_1D(
        Grid_1D& all_grid,
        int i_frame,
        const MY_FLOAT_TYPE time_step_length,
        const bool use_flux_advection,
    //		const bool set_zero_normal_velocity_at_boundary,
        const bool use_MacCormack_scheme,
        const bool use_clamping_in_MacCormack_scheme,
        const int num_gauss_quad_boundary,
        const int num_gauss_quad_bulk,
        const bool enable_eliminate_zero_velocity_false_diffusion,
        const std::string split_method,
        const std::string integral_method,
        const std::string interpolation_method,
        const bool use_integral,
        const int num_gauss_quadrature_point_for_integrate_density,
        const bool enable_cell_volume_correction,
        const bool use_lentines_advection
    ) {
        std::vector<MY_FLOAT_TYPE> before_advect_density = all_grid.substance_density;
        //移流項の計算
//        if (use_flux_advection) {
            std::string interpolation_method_in_calc_psi;
            if(interpolation_method =="linear"){
                interpolation_method_in_calc_psi = "linear";
            }
            else if(interpolation_method=="1Dy_linear"){
                interpolation_method_in_calc_psi = "const";
            }
            else{
                interpolation_method_in_calc_psi = "const";
            }
            advect_density_flux_advection_1D(
                all_grid,
				time_step_length,
				num_gauss_quad_boundary,
				num_gauss_quad_bulk,
				enable_eliminate_zero_velocity_false_diffusion,
				split_method,
			    integral_method,
				all_grid.substance_density,
				interpolation_method,
				interpolation_method_in_calc_psi,
				use_integral,
				num_gauss_quadrature_point_for_integrate_density
            );
/*
            // セルの体積によって質量密度を補正
            if(enable_cell_volume_correction){
                correct_mass_and_volume_by_pressure_solve(
                    all_grid.substance_density,
                    all_grid.cell_volume_cell_center,
                    all_grid
                );
            }
*/
/*
            if (enable_eliminate_zero_velocity_false_diffusion) {
                // 各セル中心での密度の補正量
                std::vector<MY_FLOAT_TYPE> zero_velocity_correction_values(all_grid.Grid_num_x * all_grid.Grid_num_y);
                ////// "バックトレース前の面で定義されたpsiの離散値"を足す
                for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
                    for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                        zero_velocity_correction_values[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
                            = before_advect_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
                    }
                }
                for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
                    for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                        //(ix, iy)番目のセルのface を走るループ
                        int dx[4], dy[4];
                        dx[0] = 0; dx[1] = 1; dx[2] = 1; dx[3] = 0;
                        dy[0] = 0; dy[1] = 0; dy[2] = 1; dy[3] = 1;
                        for (int i_face = 0; i_face < 4; ++i_face) {
                            // 考えるcell face を構成するvertex のindex
                            std::vector<int> vertex_index_1{ ix + dx[i_face], iy + dy[i_face] }, vertex_index_2{ ix + dx[(i_face + 1) % 4], iy + dy[(i_face + 1) % 4] };
                            // バックトレース前の cell face を構成するvertexを定義する
                            cell_vertex before_backtrace_vertex_1(vertex_index_1, all_grid._cell_length);
                            cell_vertex before_backtrace_vertex_2(vertex_index_2, all_grid._cell_length);
                            // バックトレース前の頂点から面を定義
                            cell_face before_backtrace_face;
                            before_backtrace_face._vertex_list[0] = before_backtrace_vertex_1;
                            before_backtrace_face._vertex_list[1] = before_backtrace_vertex_2;
                            ////// エラーの補正
                            ////// "バックトレース前の面上で積分した値" と "バックトレース前の面で定義されたpsiの離散値"の差を引いていく
                            zero_velocity_correction_values[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
                                -= integrate_normal_component_of_psi_on_face(
                                    all_grid,
                                    before_backtrace_face,
                                    num_gauss_quad_boundary,
                                    num_gauss_quad_bulk,
                                    split_method,
                                    interpolation_method,
                                    use_integral,
                                    num_gauss_quadrature_point_for_integrate_density
                                    ) / (all_grid._cell_length * all_grid._cell_length);
                        }
                    }
                }
                std::string interpolation_method_of_eccor_advection;
                if(interpolation_method == "linear" || interpolation_method == "1Dy_linear"){
                    interpolation_method_of_eccor_advection = "1Dy_linear";
                }
                else if(interpolation_method == "WENO6"||interpolation_method == "WENO6-optimized"){
                    interpolation_method_of_eccor_advection = "1Dy_WENO6";
                }
                std::string split_method_of_eccor_advection = "y-AxisAligned";
                std::string interpolation_method_in_calc_psi = "const";
                // エラーを移流させる
                advect_density_flux_advection(
                    all_grid,
                    time_step_length,
//					set_zero_normal_velocity_at_boundary,
                    use_MacCormack_scheme,
                    use_clamping_in_MacCormack_scheme,
                    num_gauss_quad_boundary,
                    num_gauss_quad_bulk,
                    enable_eliminate_zero_velocity_false_diffusion,
                    split_method_of_eccor_advection,
                    integral_method,
                    zero_velocity_correction_values,
                    interpolation_method_of_eccor_advection,
                    interpolation_method_in_calc_psi,
                    use_integral,
                    num_gauss_quadrature_point_for_integrate_density,
                    false
                );

                //密度を補正
                for (int ix = 0; ix < all_grid.Grid_num_x; ix++) {
                    for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
                        all_grid.substance_density[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)]
                            +=zero_velocity_correction_values[get_voxel_center_index(ix, iy, all_grid.Grid_num_x, all_grid.Grid_num_y)];
                    }
                }
            }
*/
//        }
/*
        else {
            if(use_lentines_advection){
                advect_density_Lentine(
                    all_grid,
                    time_step_length,
                    enable_cell_volume_correction
                );
            }
            else{
                advect_density_semi_lagrangian(
                    all_grid,
                    time_step_length,
                    use_MacCormack_scheme,
                    use_clamping_in_MacCormack_scheme,
                    num_gauss_quad_boundary,
                    num_gauss_quad_bulk,
                    interpolation_method
                );
            }
        }
    }
*/
}

    //密度場の移流項を flux advection で計算
    void advect_density_flux_advection_1D(
        Grid_1D& all_grid,
        const MY_FLOAT_TYPE time_step_length,
        const int num_gauss_quad_boundary,
        const int num_gauss_quad_bulk,
        const bool enable_eliminate_zero_velocity_false_diffusion,
        const std::string split_method,
        const std::string integral_method,
        std::vector<MY_FLOAT_TYPE> &advected_values,
        const std::string interpolation_method,
        const std::string interpolation_method_in_calc_psi,
        const bool use_integral,
        const int num_gauss_quadrature_point_for_integrate_density
    ) {
    // psi の離散値を計算
    //OK
	calc_psi_on_cell_face_from_density_on_cell_center_1D(
		all_grid,
		all_grid.psi_substance_density_cell_face_y,
		advected_values,
		all_grid.Grid_num_y,
		all_grid._cell_length,
		interpolation_method_in_calc_psi
	);

    std::cout << "calc_interpolated_psi: " << all_grid.psi_substance_density_cell_face_y[26] <<std::endl;
    std::cout << "calc_interpolated_psi: " << all_grid.calc_interpolated_psi(26.5 * all_grid._cell_length) <<std::endl;
    std::cout << "calc_interpolated_psi: " << all_grid.psi_substance_density_cell_face_y[27] <<std::endl;
    getchar();
	//計算結果を格納する変数
	std::vector<MY_FLOAT_TYPE> mass_after_advect(all_grid.Grid_num_y);

//    std::cout << "getchar()" <<std::endl;
//    getchar();
    //全セルを走査するループ
	for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
		MY_FLOAT_TYPE mass_in_the_cell = 0.0;

        ////左側の点からの寄与を計算
        MY_FLOAT_TYPE before_backtrace_pos = iy * all_grid.Grid_num_y;
//        std::cout << "calc_interpolated_velocity() begin" <<std::endl;
        MY_FLOAT_TYPE velocity_at_vertex = all_grid.calc_interpolated_velocity(before_backtrace_pos);
//        std::cout << "calc_interpolated_velocity() end" <<std::endl;
//        std::cout << "calc_interpolated_psi() begin" <<std::endl;
        MY_FLOAT_TYPE psi_at_vertex = all_grid.calc_interpolated_psi(before_backtrace_pos - velocity_at_vertex * time_step_length);
//        std::cout << "calc_interpolated_psi() begin" <<std::endl;
        mass_in_the_cell -= psi_at_vertex;

        ////右側の点からの寄与を計算
        before_backtrace_pos = (iy + 1) * all_grid.Grid_num_y;
//        std::cout << "calc_interpolated_velocity() begin" <<std::endl;
        velocity_at_vertex = all_grid.calc_interpolated_velocity(before_backtrace_pos);
//        std::cout << "calc_interpolated_velocity() end" <<std::endl;
//        std::cout << "calc_interpolated_psi() begin" <<std::endl;
        psi_at_vertex = all_grid.calc_interpolated_psi(before_backtrace_pos - velocity_at_vertex * time_step_length);
//        std::cout << "calc_interpolated_psi() begin" <<std::endl;
        mass_in_the_cell += psi_at_vertex;

        //質量の計算結果を格納
    	mass_after_advect[iy]
			= mass_in_the_cell;
    }
    //計算結果をコピー
	for (int iy = 0; iy < all_grid.Grid_num_y; iy++) {
		//密度の計算結果
		advected_values[iy]
			= mass_after_advect[iy]
			/ all_grid._cell_length;
	}
}

} //namespace smoke_simulation
