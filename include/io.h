#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdint>
#include <algorithm>
bool save_vti_scalar_double(
	const std::string& filename,
	const double* phi,
	int nx, int ny, int nz,
	int ghost,
	double dx, double dy, double dz);
bool write_pvd(
	const std::string& pvd_name,
	const std::vector<std::pair<double, std::string>>& files);
