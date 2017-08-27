#include "pch.hpp"
#include "read_gm.hpp"

bool read_gm(
		const std::string& file_name, 
		std::vector<double>& acc_data, 
		double& timestep)
{
	using namespace std;

	acc_data.resize(0);

	ifstream file(file_name);
	
	if (!file.is_open())
		return false;

	string line;
	double prev_val1 = 0.0;

	for (int i = 0; getline(file, line); ++i)
	{
#ifndef _WIN32
		replace(line, "\r", "");
#endif

		size_t pos = line.find_first_of(",; \t");

		if (pos == string::npos)
		{
			if (i == 0)
				continue;
			else
				return false;
		}

		stringstream ss_val1(line.substr(0, pos));
		stringstream ss_val2(line.substr(pos + 1));

		double val1, val2;
		ss_val1 >> val1;
		ss_val2 >> val2;

		if (ss_val1.fail() || ss_val2.fail())
		{
			if (i == 0)
				continue;
			else
				return false;
		}

		if (i == 2)
			timestep = val1 - prev_val1;

		acc_data.push_back(val2);

		prev_val1 = val1;
	}

	return (acc_data.size() > 0);
}
