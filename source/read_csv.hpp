#ifndef __READ_CSV__
#define __READ_CSV__

bool read_csv(
		std::string file_name, 
		std::vector<std::vector<std::string>>& values,
		const char delimiter = ',');

#endif // __READ_CSV__
