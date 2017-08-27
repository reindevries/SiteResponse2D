#include "pch.hpp"
#include "calc_response.hpp"
#include "elasto_plastic.hpp"
#include "interpolate.hpp"

//
// (surface)                         (bottom)
//
//                                       F = c v(t)
//                                      -->
//  ------     k0    ------     k1    ------    c    |/
//  | m0 |---/\/\/---| m1 |---/\/\/---| m2 |----[----|/
//  ------           ------           ------         |/
//    =>               =>               =>
//    d0               d1               d2
//

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
		std::vector<double>& d_r_max)
{
	using namespace std;

	// check paramaters

	if (k.size() < 1 || k.size() != m.size()-1 || input_a_x.size() < 1)
	{
		output_a_x.resize(0);
		output_a_y.resize(0);
		return false;
	}

	// determine the number of output samples and adjust dt

	double input_duration = (input_a_x.size() - 1) * input_dt;
	double output_duration = input_duration + output_extent;

	size_t output_samples = (size_t)round(output_duration / output_dt) + 1;
	size_t samples_per_output_dt = (size_t)ceil(output_dt / dt);
	dt = output_dt / samples_per_output_dt;

	output_a_x.resize(output_samples, 0.0);
	output_a_y.resize(output_samples, 0.0);

	// initialise vectors needed in the calculation

	size_t nodes = m.size();
	size_t elements = k.size();
	size_t layers = k[0].size();

	vector<double> d_x_n(nodes, 0.0);        // node displacements at t_n
	vector<double> d_y_n(nodes, 0.0);	     
										     
	vector<double> d_x_np1(nodes, 0.0);      // node displacements at t_{n+1}
	vector<double> d_y_np1(nodes, 0.0);	     
										     
	vector<double> d_x_nm1(nodes, 0.0);      // node displacements at t_{n-1}
	vector<double> d_y_nm1(nodes, 0.0);	     
										     
	vector<double> v_x_n(nodes);             // node velocities at t_n
	vector<double> v_y_n(nodes);		     
										     
	vector<double> F_s_x_n(elements);        // element spring force
	vector<double> F_s_y_n(elements);	     
										     
	vector<double> F_d_x_n(elements);        // element damper force (Rayleigh, stiffness proporional)
	vector<double> F_d_y_n(elements);		 
											 
	vector<vector<double>> d_o_x(elements);  // origin of elasto-plastic layers
	vector<vector<double>> d_o_y(elements);

	for (size_t i = 0; i < elements; ++i)
	{
		d_o_x[i].resize(layers, 0.0);
		d_o_y[i].resize(layers, 0.0);
	}

	d_r_max.resize(elements, 0.0);

	double input_a_x_nm1 = -interpolate(input_a_x, input_dt, 0.0);  // cancel first contribution if a(0) != 0
	double input_a_y_nm1 = -interpolate(input_a_y, input_dt, 0.0);

	double input_v_x_nm1 = 0.0;
	double input_v_y_nm1 = 0.0;

	// perform the calculation

	for (size_t i_output = 0; i_output < output_samples; ++i_output)
	{
		for (size_t i_per_output = 0; i_per_output < samples_per_output_dt; ++i_per_output)
		{
			double t_n = i_output * output_dt + i_per_output * dt;

			// integrate ground motion acceleration to velocity

			double input_a_x_n = interpolate(input_a_x, input_dt, t_n);
			double input_a_y_n = interpolate(input_a_y, input_dt, t_n);

			double input_v_x_n = input_v_x_nm1 + (input_a_x_nm1 + input_a_x_n) / 2.0 * dt;
			double input_v_y_n = input_v_y_nm1 + (input_a_y_nm1 + input_a_y_n) / 2.0 * dt;

			// store previous displacements

			for (size_t i = 0; i < nodes; ++i)
			{
				d_x_nm1[i] = d_x_n[i];
				d_y_nm1[i] = d_y_n[i];

				d_x_n[i] = d_x_np1[i];
				d_y_n[i] = d_y_np1[i];
			}

			// calculate velocities at time t_n

			for (size_t i = 0; i < nodes; ++i)
			{
				v_x_n[i] = (d_x_n[i] - d_x_nm1[i]) / dt;
				v_y_n[i] = (d_y_n[i] - d_y_nm1[i]) / dt;
			}

			// calculate internal forces using displacements at time t_n

			for (size_t i = 0; i < elements; ++i)
			{
				double d_x_rel = d_x_n[i+1] - d_x_n[i];
				double d_y_rel = d_y_n[i+1] - d_y_n[i];

				double v_x_rel = v_x_n[i+1] - v_x_n[i];
				double v_y_rel = v_y_n[i+1] - v_y_n[i];

				F_s_x_n[i] = 0.0;
				F_s_y_n[i] = 0.0;

				F_d_x_n[i] = 0.0;
				F_d_y_n[i] = 0.0;

				for (size_t j = 0; j < layers; ++j)
				{
					// calculate forces and tangential stiffness for the 2D elasto-plasic layer

					double F_x, F_y, k_tan;
					elasto_plastic(k[i][j], d_yield[i][j], d_o_x[i][j], d_o_y[i][j], 
							d_x_rel, d_y_rel, F_x, F_y, k_tan);

					// add them to the total element spring force

					F_s_x_n[i] += F_x;
					F_s_y_n[i] += F_y;

					// and add the corresponding damper forces (stiffness proportional)

					if (damping_tangent)
					{
						F_d_x_n[i] += beta * k_tan * v_x_rel;
						F_d_y_n[i] += beta * k_tan * v_y_rel;
					}
					else
					{
						F_d_x_n[i] += beta * k[i][j] * v_x_rel;
						F_d_y_n[i] += beta * k[i][j] * v_y_rel;
					}
				}

				// store max relative displacement (radius)

				double d_r = sqrt(sq(d_x_rel) + sq(d_y_rel));
				d_r_max[i] = max(d_r, d_r_max[i]);
			}

			// left node (surface)

			{
				size_t i = 0;

				double F_e_x_n = 0.0;                    // external force acting on node
				double F_e_y_n = 0.0;

				double F_s_x_n_node = 0.0 - F_s_x_n[i];  // spring force acting on node
				double F_s_y_n_node = 0.0 - F_s_y_n[i];

				double F_d_x_n_node = 0.0 - F_d_x_n[i];  // damper force acting on node
				double F_d_y_n_node = 0.0 - F_d_y_n[i];

				d_x_np1[i] = sq(dt) / m[i] * (F_e_x_n - F_s_x_n_node - F_d_x_n_node - 
						1.0 / sq(dt) * m[i] * d_x_nm1[i] + 2.0 / sq(dt) * m[i] * d_x_n[i]);

				d_y_np1[i] = sq(dt) / m[i] * (F_e_y_n - F_s_y_n_node - F_d_y_n_node - 
						1.0 / sq(dt) * m[i] * d_y_nm1[i] + 2.0 / sq(dt) * m[i] * d_y_n[i]);
			}

			// middle nodes

			for (size_t i = 1; i < nodes-1; ++i)
			{
				double F_e_x_n = 0.0;                             // external force acting on node
				double F_e_y_n = 0.0;

				double F_s_x_n_node = F_s_x_n[i-1] - F_s_x_n[i];  // spring force acting on node
				double F_s_y_n_node = F_s_y_n[i-1] - F_s_y_n[i];

				double F_d_x_n_node = F_d_x_n[i-1] - F_d_x_n[i];  // damper force acting on node
				double F_d_y_n_node = F_d_y_n[i-1] - F_d_y_n[i];

				d_x_np1[i] = sq(dt) / m[i] * (F_e_x_n - F_s_x_n_node - F_d_x_n_node - 
						1.0 / sq(dt) * m[i] * d_x_nm1[i] + 2.0 / sq(dt) * m[i] * d_x_n[i]);

				d_y_np1[i] = sq(dt) / m[i] * (F_e_y_n - F_s_y_n_node - F_d_y_n_node - 
						1.0 / sq(dt) * m[i] * d_y_nm1[i] + 2.0 / sq(dt) * m[i] * d_y_n[i]);
			}

			// right node (bottom)

			{
				size_t i = nodes-1;

				double F_e_x_n = c * scale_factor * input_v_x_n;    // external force acting on node
				double F_e_y_n = c * scale_factor * input_v_y_n;

				double F_s_x_n_node = F_s_x_n[i-1];                 // spring force acting on node
				double F_s_y_n_node = F_s_y_n[i-1];

				double F_d_x_n_node = F_d_x_n[i-1] + c * v_x_n[i];  // damper force acting on node
				double F_d_y_n_node = F_d_y_n[i-1] + c * v_y_n[i];

				d_x_np1[i] = sq(dt) / m[i] * (F_e_x_n - F_s_x_n_node - F_d_x_n_node - 
						1.0 / sq(dt) * m[i] * d_x_nm1[i] + 2.0 / sq(dt) * m[i] * d_x_n[i]);

				d_y_np1[i] = sq(dt) / m[i] * (F_e_y_n - F_s_y_n_node - F_d_y_n_node - 
						1.0 / sq(dt) * m[i] * d_y_nm1[i] + 2.0 / sq(dt) * m[i] * d_y_n[i]);
			}

			// output accelerations at the surface node

			if (i_per_output == 0)
			{
				output_a_x[i_output] = (d_x_np1[0] - 2.0*d_x_n[0] + d_x_nm1[0]) / sq(dt);
				output_a_y[i_output] = (d_y_np1[0] - 2.0*d_y_n[0] + d_y_nm1[0]) / sq(dt);
			}

			// store ground motion acceleration and velocity for next timestep

			input_a_x_nm1 = input_a_x_n;
			input_a_y_nm1 = input_a_y_n;

			input_v_x_nm1 = input_v_x_n;
			input_v_y_nm1 = input_v_y_n;
		}
	}

	return true;
}
