#pragma once
#include "Matrix.h"
#include <string>
#include <vector>

#define MAX_CITIES 128

namespace data
{
struct InputData
{
	std::string city_name;
	int city_id = 0;
	int lat_deg = 0;
	int lat_min = 0;
	int long_deg = 0;
	int long_min = 0;
};
typedef std::vector<InputData> InputDataBuffer;

struct OutputData
{
	std::vector<int> path;
};

struct System
{
	int num_cities;
	InputDataBuffer input_data_buffer;
	OutputData output_data;
	Matrix distance_matrix;
};
}; // namespace data

typedef std::vector<int> Path;
