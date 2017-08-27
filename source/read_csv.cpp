#include "pch.hpp"
#include "read_csv.hpp"

bool read_csv(
		std::string file_name, 
		std::vector<std::vector<std::string>>& values,
		const char delimiter)
{
	// [0][0], [0][1], [0][2]
	// [1][0], [1][1], [1][2]
	// etc.

	using namespace std;

	ifstream file(file_name);
	if (!file.is_open())
		return false;

	string line, token;

	values.resize(0);

	while (getline(file, line))
	{
#ifndef _WIN32
		replace(line, "\r", "");
#endif
		
		values.resize(values.size() + 1);
		istringstream ss(line);
		while (getline(ss, token, delimiter))
			values.back().push_back(token);
	}

	return true;
}
