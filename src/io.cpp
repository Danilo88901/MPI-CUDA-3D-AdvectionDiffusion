#include"io.h"
bool save_vti_scalar_double(
	const std::string& filename,
	const double* phi,
	int nx, int ny, int nz,
	int ghost,
	double dx, double dy, double dz)
{
	const int px = nx - 2 * ghost;
	const int py = ny - 2 * ghost;
	const int pz = nz - 2 * ghost;

	if (px <= 0 || py <= 0 || pz <= 0) {
		std::cerr << "invalid sizes\n";
		return false;
	}

	const size_t npts = static_cast<size_t>(px) * py * pz;
	std::vector<double> buf(npts);

	// VTK ordering: x fastest, then y, then z
	for (int z = 0; z < pz; ++z) {
		for (int y = 0; y < py; ++y) {
			for (int x = 0; x < px; ++x) {

				const size_t src =
					static_cast<size_t>(z + ghost) * ny * nx +
					static_cast<size_t>(y + ghost) * nx +
					static_cast<size_t>(x + ghost);

				const size_t dst =
					static_cast<size_t>(z) * py * px +
					static_cast<size_t>(y) * px +
					static_cast<size_t>(x);

				buf[dst] = phi[src];
			}
		}
	}

	std::ofstream f(filename, std::ios::binary | std::ios::trunc);
	if (!f) {
		std::cerr << "cannot open file\n";
		return false;
	}

	// ---------------- XML HEADER ----------------
	f << "<?xml version=\"1.0\"?>\n";
	f << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";

	f << "  <ImageData WholeExtent=\"0 " << (px - 1)
		<< " 0 " << (py - 1)
		<< " 0 " << (pz - 1)
		<< "\" Origin=\"0 0 0\" Spacing=\""
		<< dx << " " << dy << " " << dz << "\">\n";

	f << "    <Piece Extent=\"0 " << (px - 1)
		<< " 0 " << (py - 1)
		<< " 0 " << (pz - 1) << "\">\n";

	f << "      <PointData Scalars=\"phi\">\n";

	f << "        <DataArray type=\"Float64\" "
		"Name=\"phi\" "
		"NumberOfComponents=\"1\" "
		"format=\"appended\" offset=\"0\"/>\n";

	f << "      </PointData>\n";
	f << "      <CellData/>\n";

	f << "    </Piece>\n";
	f << "  </ImageData>\n";

	// ---------------- APPENDED DATA ----------------
	f << "  <AppendedData encoding=\"raw\">_";

	const uint32_t byte_count = static_cast<uint32_t>(npts * sizeof(double));

	f.write(reinterpret_cast<const char*>(&byte_count), sizeof(uint32_t));
	f.write(reinterpret_cast<const char*>(buf.data()), byte_count);

	f << "\n</AppendedData>\n";
	f << "</VTKFile>\n";

	f.flush();
	f.close();

	return true;
}
bool write_pvd(
	const std::string& pvd_name,
	const std::vector<std::pair<double, std::string>>& files)
{
	std::ofstream pvd(pvd_name);
	if (!pvd) {
		std::cerr << "write_pvd: cannot open file " << pvd_name << "\n";
		return false;
	}

	pvd << "<?xml version=\"1.0\"?>\n";
	pvd << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	pvd << "  <Collection>\n";

	for (size_t i = 0; i < files.size(); ++i) {
		const double t = files[i].first;
		const std::string& file = files[i].second;

		pvd << "    <DataSet timestep=\"" << std::setprecision(17) << t
			<< "\" group=\"\" part=\"0\" file=\"" << file << "\"/>\n";
	}

	pvd << "  </Collection>\n";
	pvd << "</VTKFile>\n";
	return static_cast<bool>(pvd);
}
