#ifndef __CALC_RESPONSE__
#define __CALC_RESPONSE__

bool calc_response(
		const std::vector<double>& m,
		const std::vector<std::vector<double>>& k,
		const std::vector<std::vector<double>>& d_yield,
		double c,
		double dt,
		const std::vector<double>& input_a_x,
		const std::vector<double>& input_a_y,
		double input_dt,
		double scale_factor,
		std::vector<double>& output_a_x,
		std::vector<double>& output_a_y,
		double output_dt,
		double output_extent,
		double beta,
		bool damping_tangent,
		std::vector<double>& d_r_max);

#endif // __CALC_RESPONSE__
