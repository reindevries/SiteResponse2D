#include "pch.hpp"
#include "read_csv.hpp"
#include "read_gm.hpp"
#include "calc_response.hpp"
#include "interpolate.hpp"

#if defined(_WIN32) && defined(_DEBUG)
	#define RETURN {cout << endl; system("pause"); return 0;}
	#define RETURN_ERROR(error_msg) {cout << error_msg << endl << endl; system("pause"); return -1;}
#else
	#define RETURN {cout << endl; return 0;}
	#define RETURN_ERROR(error_msg) {cout << error_msg << endl << endl; return -1;}
#endif

int main(int argc, const char* argv[])
{
	using namespace std;

	const double pi = acos(-1.0);

	//
	// Open and read the input file ----------------------------------------------------------------
	//

	cout << "Input parameters" << endl;
	cout << "----------------" << endl;

	if (argc == 1)
		RETURN_ERROR("Please specify input file");

	string input_file_name(argv[1]);
	ifstream input_file(input_file_name);

	if (!input_file.is_open())
		RETURN_ERROR("Unable to open input file");

	string gm_x_file_name;
	string gm_y_file_name;
	string gm_z_file_name;
	string output_gm_x_file_name;
	string output_gm_y_file_name;
	string output_gm_z_file_name;
	string soil_elements_file_name;
	string max_strains_file_name;
	double scale_factor = 1.0;
	double damping_ratio = 0.01;
	double damping_period = 0.1;
	bool damping_tangent = true;
	double p_wave_vel_z = 1500.0;
	double timestep = 5e-4;
	double output_timestep = 5e-3;
	double output_extent = 0.0;
	int verbosity = 1;

	string line;

	while (getline(input_file, line))
	{
		while (replace(line, " =", "="));
		while (replace(line, "= ", "="));
		replace(line, "\r", "");

		if (line.compare(0, 1, "#") == 0)
			continue;
		else if (line.compare(0, 16, "ground_motion_x=") == 0)
			gm_x_file_name = line.substr(16);
		else if (line.compare(0, 16, "ground_motion_y=") == 0)
			gm_y_file_name = line.substr(16);
		else if (line.compare(0, 16, "ground_motion_z=") == 0)
			gm_z_file_name = line.substr(16);
		else if (line.compare(0, 23, "output_ground_motion_x=") == 0)
			output_gm_x_file_name = line.substr(23);
		else if (line.compare(0, 23, "output_ground_motion_y=") == 0)
			output_gm_y_file_name = line.substr(23);
		else if (line.compare(0, 23, "output_ground_motion_z=") == 0)
			output_gm_z_file_name = line.substr(23);
		else if (line.compare(0, 14, "soil_elements=") == 0)
			soil_elements_file_name = line.substr(14);
		else if (line.compare(0, 12, "max_strains=") == 0)
			max_strains_file_name = line.substr(12);
		else if (line.compare(0, 13, "scale_factor=") == 0)
			scale_factor = stod(line.substr(13));
		else if (line.compare(0, 14, "damping_ratio=") == 0)
			damping_ratio = stod(line.substr(14));
		else if (line.compare(0, 15, "damping_period=") == 0)
			damping_period = stod(line.substr(15));
		else if (line.compare(0, 16, "damping_tangent=") == 0)
			damping_tangent = (line.compare(16, 4, "true") == 0);
		else if (line.compare(0, 13, "p_wave_vel_z=") == 0)
			p_wave_vel_z = stod(line.substr(13));
		else if (line.compare(0, 9, "timestep=") == 0)
			timestep = stod(line.substr(9));
		else if (line.compare(0, 16, "output_timestep=") == 0)
			output_timestep = stod(line.substr(16));
		else if (line.compare(0, 14, "output_extent=") == 0)
			output_extent = stod(line.substr(14));
		else if (line.compare(0, 10, "verbosity=") == 0)
			verbosity = stoi(line.substr(10));
		else
			cout << "Unknown parameter: " << line << endl;
	}

	cout << "ground_motion_x = " << gm_x_file_name << endl;
	cout << "ground_motion_y = " << gm_y_file_name << endl;
	cout << "ground_motion_z = " << gm_z_file_name << endl;
	cout << "output_ground_motion_x = " << output_gm_x_file_name << endl;
	cout << "output_ground_motion_y = " << output_gm_y_file_name << endl;
	cout << "output_ground_motion_z = " << output_gm_z_file_name << endl;
	cout << "soil_elements = " << soil_elements_file_name << endl;
	cout << "max_strains = " << max_strains_file_name << endl;
	cout << "scale_factor = " << scale_factor << endl;
	cout << "damping_ratio = " << damping_ratio << endl;
	cout << "damping_period = " << damping_period << endl;
	cout << "damping_tangent = " << (damping_tangent ? "true" : "false") << endl;
	cout << "p_wave_vel_z = " << p_wave_vel_z << endl;
	cout << "timestep = " << timestep << endl;
	cout << "output_timestep = " << output_timestep << endl;
	cout << "output_extent = " << output_extent << endl;
	cout << "verbosity = " << verbosity << endl;
	cout << endl;

	bool analyse_2d = !gm_y_file_name.empty();
	bool process_z = !gm_z_file_name.empty();
	bool output_max_strains = !max_strains_file_name.empty();

	//
	// Extract head and put in front of given file names -------------------------------------------
	//

	size_t head_pos = input_file_name.find_last_of("\\/");

	string head;
	if (head_pos != string::npos)
		head = input_file_name.substr(0, head_pos + 1);

	gm_x_file_name          = head + gm_x_file_name;
	gm_y_file_name          = head + gm_y_file_name;
	gm_z_file_name          = head + gm_z_file_name;
	output_gm_x_file_name   = head + output_gm_x_file_name;
	output_gm_y_file_name   = head + output_gm_y_file_name;
	output_gm_z_file_name   = head + output_gm_z_file_name;
	soil_elements_file_name = head + soil_elements_file_name;
	max_strains_file_name   = head + max_strains_file_name;

	//
	// Read ground motions -------------------------------------------------------------------------
	//

	cout << "Ground motions" << endl;
	cout << "--------------" << endl;

	double gm_dt_x, gm_dt_y, gm_dt_z;
	vector<double> gm_a_x, gm_a_y, gm_a_z;

	if (!read_gm(gm_x_file_name, gm_a_x, gm_dt_x))
		RETURN_ERROR("Failure reading x-direction data");

	if (analyse_2d)
	{
		if (!read_gm(gm_y_file_name, gm_a_y, gm_dt_y))
			RETURN_ERROR("Failure reading y-direction data");

		if (gm_a_x.size() != gm_a_y.size() || abs(gm_dt_x / gm_dt_y - 1.0) > 1e-6)
			RETURN_ERROR("Ground motion pair timestep/duration does not match");
	}
	else
	{
		gm_a_y.resize(gm_a_x.size(), 0.0);
		gm_dt_y = gm_dt_x;
	}

	if (process_z)
	{
		if (!read_gm(gm_z_file_name, gm_a_z, gm_dt_z))
			RETURN_ERROR("Failure reading z-direction data");

		if (gm_a_x.size() != gm_a_z.size() || abs(gm_dt_x / gm_dt_z - 1.0) > 1e-6)
			RETURN_ERROR("Ground motion pair timestep/duration does not match");
	}

	double duration = (gm_a_x.size() - 1) * gm_dt_x;

	double PGA_x = 0.0;
	for (auto&& a : gm_a_x)
		PGA_x = max(abs(a), PGA_x);

	double PGA_y = 0.0;
	for (auto&& a : gm_a_y)
		PGA_y = max(abs(a), PGA_y);

	double PGA_z = 0.0;
	for (auto&& a : gm_a_z)
		PGA_z = max(abs(a), PGA_z);

	cout << "Read " << gm_a_x.size() << " points" << endl;
	cout << "Timestep is " << gm_dt_x << " s" << endl;
	cout << "Duration is " << duration << " s" << endl;
	cout << "Unscaled PGA is " << PGA_x << " m/s^2 in x-direction" << endl;
	cout << "Unscaled PGA is " << PGA_y << " m/s^2 in y-direction" << endl;
	cout << "Unscaled PGA is " << PGA_z << " m/s^2 in z-direction" << endl;
	cout << endl;

	//
	// Read soil elements --------------------------------------------------------------------------
	//

	cout << "Soil elements" << endl;
	cout << "-------------" << endl;

	vector<vector<string>> values;
	if (!read_csv(soil_elements_file_name, values))
		RETURN_ERROR("Unable to open soil elements file");

	if (values.size() <= 2)
		RETURN_ERROR("No elements defined");

	vector<double> gamma;
	for (string& value : values[0])
	{
		stringstream ss(value);
		double d;
		ss >> d;
		if (!ss.fail())
			gamma.push_back(d);
	}

	struct soil_element
	{
		string name;
		double z, rho, Vs;
		vector<double> G_sec_over_G_0;
	};

	vector<soil_element> soil_elements;
	for (size_t i = 1; i < values.size(); ++i)
	{
		soil_element se;
		se.z = se.rho = se.Vs = nan("");
		for (size_t j = 0; j < values[i].size(); ++j)
		{
			if (values[0][j].compare("name") == 0)
				se.name = values[i][j];
			else if (values[0][j].compare("z") == 0)
				se.z = stod(values[i][j]);
			else if (values[0][j].compare("rho") == 0)
				se.rho = stod(values[i][j]);
			else if (values[0][j].compare("Vs") == 0)
				se.Vs = stod(values[i][j]);
			else if (!values[i][j].empty())
				se.G_sec_over_G_0.push_back(stod(values[i][j]));
		}

		if (std::isnan(se.z) || std::isnan(se.rho) || std::isnan(se.Vs))
			RETURN_ERROR("Invalid element properties");

		soil_elements.push_back(se);
	}

	size_t n = soil_elements.size();  // number of nodes (includes boundary)

	cout << "Read " << n-1 << " elements" << endl;
	cout << "G_sec/G_0 curves discretised into " << gamma.size() << " points" << endl;

	vector<double> G_0(n-1);
	for (size_t i = 0; i < n-1; ++i)
		G_0[i] = soil_elements[i].rho * sq(soil_elements[i].Vs);

	vector<double> h(n-1);
	for (size_t i = 0; i < n-1; ++i)
		h[i] = abs(soil_elements[i+1].z) - abs(soil_elements[i].z);

	if (verbosity > 1)
		for (size_t i = 0; i < n-1; ++i)
			cout << "Element " << i+1 << ": h = " << h[i] <<
					" m, G_0 = " << G_0[i] << " N/m^2" << endl;

	cout << "Elastic half-space boundary defined at " << soil_elements[n-1].z << " m" << endl;
	cout << "with mass density " << soil_elements[n-1].rho << " kg/m^3" << endl;
	cout << "and shear wave velocity " << soil_elements[n-1].Vs << " m/s" << endl;

	cout << "Difference in elastic impedance of boundary is " <<
			(int)round(100.0 * ((soil_elements[n-1].rho * soil_elements[n-1].Vs) /
			(soil_elements[n-2].rho * soil_elements[n-2].Vs) - 1.0)) << "% " << endl;
	cout << endl;

	//
	// Perform site response analysis --------------------------------------------------------------
	//

	cout << "Site response analysis" << endl;
	cout << "----------------------" << endl;

	double A = 1.0;
	if (verbosity > 1)
		cout << "Using unit area A = 1 m^2" << endl;

	vector<double> m(n, 0.0);
	for (size_t i = 0; i < n-1; ++i)
	{
		double M = soil_elements[i].rho * A * h[i];
		m[i] += M / 2.0;
		m[i + 1] += M / 2.0;
	}

	if (verbosity > 1)
		for (size_t i = 0; i < n; ++i)
			cout << "Node " << i+1 << ": m = " << m[i] << " kg" << endl;

	vector<vector<double>> G(n-1);
	for (size_t i = 0; i < n-1; ++i)
	{
		G[i].resize(gamma.size());

		for (size_t j = gamma.size()-1; j > 0; --j)  // work backwards
		{
			// calculate tangential stiffness

			G[i][j] = (soil_elements[i].G_sec_over_G_0[j] * G_0[i] * gamma[j] -
					soil_elements[i].G_sec_over_G_0[j-1] * G_0[i] * gamma[j-1]) /
					(gamma[j] - gamma[j-1]);

			// subtract all that follow to obtain individual stiffness

			for (size_t k = j+1; k < gamma.size(); ++k)
				G[i][j] -= G[i][k];
		}

		G[i][0] = G_0[i];

		for (size_t j = 1; j < gamma.size(); ++j)
			G[i][0] -= G[i][j];
	}

	vector<double> k_0(n-1);
	for (size_t i = 0; i < n-1; ++i)
		k_0[i] = A * G_0[i] / h[i];

	vector<vector<double>> k(n-1), d_yield(n-1);
	for (size_t i = 0; i < n-1; ++i)
	{
		k[i].resize(gamma.size());
		d_yield[i].resize(gamma.size());

		for (size_t j = 0; j < gamma.size(); ++j)
		{
			k[i][j] = A * G[i][j] / h[i];
			d_yield[i][j] = gamma[j] * h[i];
		}
	}

	// calculate critical timestep

	vector<double> dt_crit(n-1, 0.0);
	for (size_t i = 0; i < n-1; ++i)
		dt_crit[i] = 2.0 * sqrt((m[i] * m[i+1]) / (m[i] + m[i+1]) * 1.0 / k_0[i]);

	double min_dt_crit = 1e10;
	for (size_t i = 0; i < n-1; ++i)
		min_dt_crit = min(dt_crit[i], min_dt_crit);

	if (verbosity > 1)
		for (size_t i = 0; i < n-1; ++i)
			cout << "Element " << i+1 << ": k_0 = " << k_0[i]
					<< " N/m, dt_crit = " << dt_crit[i] << " s" << endl;
	cout << "Critical timestep is 0.9 x " << min_dt_crit << " = " <<
			0.9 * min_dt_crit << " s" << endl;

	timestep = min(timestep, min_dt_crit);

	cout << "Using timestep " << timestep << " s" << endl;

	// calculate the damping value (Rayleigh, only stiffness proportional)

	double avg_Vs = 0.0;
	for (size_t i = 0; i < n-1; ++i)
		avg_Vs += h[i] / soil_elements[i].Vs;
	double total_height = abs(soil_elements.back().z);
	avg_Vs = total_height / avg_Vs;
	cout << "avg_Vs = " << avg_Vs << " m/s" << endl;

	double T_1 = 4.0 * total_height / avg_Vs;
	cout << "T_1 = " << T_1 << " s (estimated)" << endl;

	if (abs(damping_period) < 1e-6)
	{
		cout << "Using T_1 for damping_period" << endl;
		damping_period = T_1;
	}
	
	double beta = damping_ratio * damping_period / pi;
	cout << "beta = " << beta << endl;

	double c = soil_elements.back().rho * soil_elements.back().Vs;
	cout << "c = " << c << endl;

	// calculate the response

	vector<double> output_a_x, output_a_y, d_r_max;

	if (!calc_response(m, k, d_yield, c, timestep, gm_a_x, gm_a_y, gm_dt_x, scale_factor,
			output_a_x, output_a_y, output_timestep, output_extent, beta, damping_tangent, d_r_max))
		RETURN_ERROR("Error in calculating the response");

	cout << "Calculation done" << endl;

	vector<double> max_strains(d_r_max.size());
	for (size_t i = 0; i < max_strains.size(); ++i)
		max_strains[i] = d_r_max[i] / h[i];

	if (verbosity > 1)
	{
		cout << "Maximum strains:" << endl;
		for (size_t i = 0; i < max_strains.size(); ++i)
			cout << soil_elements[i].z << "\t" << max_strains[i] << endl;
	}

	//
	// Apply scale factor and delay the z-component ------------------------------------------------
	//

	vector<double> output_a_z;

	if (process_z)
	{
		double total_height = abs(soil_elements[n-1].z - soil_elements[0].z);
		double delay = total_height / p_wave_vel_z;
		cout << "Delay for z-component is " << delay << " s" << endl;

		gm_a_z.insert(gm_a_z.begin(), 0.0);
		gm_a_z.push_back(0.0);
		
		output_a_z.resize(output_a_x.size());
		for (size_t i = 0; i < output_a_z.size(); ++i)
			output_a_z[i] = scale_factor * interpolate(gm_a_z, 
					gm_dt_x, i * output_timestep - delay, -gm_dt_x);
	}

	//
	// Write surface acceleration to file ----------------------------------------------------------
	//

	ofstream output_gm_x_file(output_gm_x_file_name);
	if (!output_gm_x_file.is_open())
		RETURN_ERROR("Error openening x-direction output file");

	output_gm_x_file << "t [s],a_x [m/s^2]" << endl;

	for (size_t i = 0; i < output_a_x.size(); ++i)
		output_gm_x_file << i * output_timestep << "," << 
				(abs(output_a_x[i]) > 1e-100 ? output_a_x[i] : 0.0) << endl;

	if (analyse_2d)
	{
		ofstream output_gm_y_file(output_gm_y_file_name);
		if (!output_gm_y_file.is_open())
			RETURN_ERROR("Error openening y-direction output file");

		output_gm_y_file << "t [s],a_y [m/s^2]" << endl;

		for (size_t i = 0; i < output_a_y.size(); ++i)
			output_gm_y_file << i * output_timestep << "," << 
					(abs(output_a_y[i]) > 1e-100 ? output_a_y[i] : 0.0) << endl;
	}

	if (process_z)
	{
		ofstream output_gm_z_file(output_gm_z_file_name);
		if (!output_gm_z_file.is_open())
			RETURN_ERROR("Error openening z-direction output file");

		output_gm_z_file << "t [s],a_z [m/s^2]" << endl;

		for (size_t i = 0; i < output_a_z.size(); ++i)
			output_gm_z_file << i * output_timestep << "," << 
					(abs(output_a_z[i]) > 1e-100 ? output_a_z[i] : 0.0) << endl;
	}

	cout << "Written surface acceleration to file" << endl;

	if (output_max_strains)
	{
		ofstream max_strains_file(max_strains_file_name);
		if (!max_strains_file.is_open())
			RETURN_ERROR("Error opening max_strains file");

		max_strains_file << "z,gamma_max" << endl;
		for (size_t i = 0; i < max_strains.size(); ++i)
			max_strains_file << soil_elements[i].z << "," << max_strains[i] << endl;

		cout << "Written maximum strains to file" << endl;		
	}

	RETURN;
}
