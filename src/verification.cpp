#include "verification.h"
#include<math.h>
#include<iostream>
constexpr double M_PI = 3.14159265358979323846;
void fill_velocity(double* u, double* v, double* w,
    int nx, int ny, int nz,
    double dx, double dy, double dz,
    double /*xc*/, double /*yc*/,
    double u0, double v0, double w0)
{
    (void)dx; (void)dy; (void)dz;

    for (int i = 0; i < nz; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nx; ++k) {
                int idx = i * ny * nx + j * nx + k;
                u[idx] = u0;
                v[idx] = v0;
                w[idx] = w0;
            }
        }
    }
}
void init(double* phi, int nx, int ny, int nz,
    double dx, double dy, double dz,
    double xc, double yc, double zc, double /*sigma0*/)
{
    const double Lx = (nx - 4) * dx;
    const double Ly = (ny - 4) * dy;
    const double Lz = (nz - 4) * dz;

    // ????????? ???, ????? ???? ???? “???????”
    const double a0 = 0.15;

    const int m1x = 4, m1y = 3, m1z = 2;
    const int m2x = 5, m2y = 2, m2z = 4;
    const int m3x = 7, m3y = 5, m3z = 1;
    const int m4x = 3, m4y = 4, m4z = 6;

    const double k1x = 2.0 * M_PI * m1x / Lx;
    const double k1y = 2.0 * M_PI * m1y / Ly;
    const double k1z = 2.0 * M_PI * m1z / Lz;

    const double k2x = 2.0 * M_PI * m2x / Lx;
    const double k2y = 2.0 * M_PI * m2y / Ly;
    const double k2z = 2.0 * M_PI * m2z / Lz;

    const double k3x = 2.0 * M_PI * m3x / Lx;
    const double k3y = 2.0 * M_PI * m3y / Ly;
    const double k3z = 2.0 * M_PI * m3z / Lz;

    const double k4x = 2.0 * M_PI * m4x / Lx;
    const double k4y = 2.0 * M_PI * m4y / Ly;
    const double k4z = 2.0 * M_PI * m4z / Lz;

    for (int i = 0; i < nz; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nx; ++k) {
                const int idx = i * ny * nx + j * nx + k;

                const double x = (k - 2) * dx;
                const double y = (j - 2) * dy;
                const double z = (i - 2) * dz;

                const double x0 = x - xc;
                const double y0 = y - yc;
                const double z0 = z - zc;

                phi[idx] =
                    a0
                    + 0.18 * std::cos(k1x * x0) * std::cos(k1y * y0) * std::cos(k1z * z0)
                    + 0.14 * std::sin(k2x * x0 + 0.30) * std::cos(k2y * y0 - 0.20) * std::sin(k2z * z0 + 0.40)
                    + 0.12 * std::cos(k3x * x0 - 0.50) * std::sin(k3y * y0 + 0.10) * std::cos(k3z * z0)
                    + 0.10 * std::sin(k4x * x0) * std::sin(k4y * y0 + 0.70) * std::cos(k4z * z0 - 0.20);
            }
        }
    }
}
void exact_solution(double* phi_exact, int nx, int ny, int nz,
    double dx, double dy, double dz,
    double t,
    double xc, double yc, double zc,
    double u0, double v0, double w0,
    double nu)
{
    const double Lx = (nx - 4) * dx;
    const double Ly = (ny - 4) * dy;
    const double Lz = (nz - 4) * dz;

    const double a0 = 0.15;

    const int m1x = 4, m1y = 3, m1z = 2;
    const int m2x = 5, m2y = 2, m2z = 4;
    const int m3x = 7, m3y = 5, m3z = 1;
    const int m4x = 3, m4y = 4, m4z = 6;

    const double k1x = 2.0 * M_PI * m1x / Lx;
    const double k1y = 2.0 * M_PI * m1y / Ly;
    const double k1z = 2.0 * M_PI * m1z / Lz;
    const double k1_2 = k1x * k1x + k1y * k1y + k1z * k1z;

    const double k2x = 2.0 * M_PI * m2x / Lx;
    const double k2y = 2.0 * M_PI * m2y / Ly;
    const double k2z = 2.0 * M_PI * m2z / Lz;
    const double k2_2 = k2x * k2x + k2y * k2y + k2z * k2z;

    const double k3x = 2.0 * M_PI * m3x / Lx;
    const double k3y = 2.0 * M_PI * m3y / Ly;
    const double k3z = 2.0 * M_PI * m3z / Lz;
    const double k3_2 = k3x * k3x + k3y * k3y + k3z * k3z;

    const double k4x = 2.0 * M_PI * m4x / Lx;
    const double k4y = 2.0 * M_PI * m4y / Ly;
    const double k4z = 2.0 * M_PI * m4z / Lz;
    const double k4_2 = k4x * k4x + k4y * k4y + k4z * k4z;

    const double d1 = std::exp(-nu * k1_2 * t);
    const double d2 = std::exp(-nu * k2_2 * t);
    const double d3 = std::exp(-nu * k3_2 * t);
    const double d4 = std::exp(-nu * k4_2 * t);

    for (int i = 0; i < nz; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nx; ++k) {
                const int idx = i * ny * nx + j * nx + k;

                const double x = (k - 2) * dx;
                const double y = (j - 2) * dy;
                const double z = (i - 2) * dz;

                const double x0 = x - xc - u0 * t;
                const double y0 = y - yc - v0 * t;
                const double z0 = z - zc - w0 * t;

                phi_exact[idx] =
                    a0
                    + 0.18 * d1 * std::cos(k1x * x0) * std::cos(k1y * y0) * std::cos(k1z * z0)
                    + 0.14 * d2 * std::sin(k2x * x0 + 0.30) * std::cos(k2y * y0 - 0.20) * std::sin(k2z * z0 + 0.40)
                    + 0.12 * d3 * std::cos(k3x * x0 - 0.50) * std::sin(k3y * y0 + 0.10) * std::cos(k3z * z0)
                    + 0.10 * d4 * std::sin(k4x * x0) * std::sin(k4y * y0 + 0.70) * std::cos(k4z * z0 - 0.20);
            }
        }
    }
}
void compute_error(double* phi, double* phi_exact,
	int nx, int ny, int nz,
	double* l2_rel, double* linf)
{
	double num = 0.0;
	double den = 0.0;
	double maxerr = 0.0;

	for (int i = 2; i < nz - 2; i++) {
		for (int j = 2; j < ny - 2; j++) {
			for (int k = 2; k < nx - 2; k++) {
				int idx = i * ny * nx + j * nx + k;
				double d = phi[idx] - phi_exact[idx];

				num += d * d;
				den += phi_exact[idx] * phi_exact[idx];

				double ad = fabs(d);
				if (ad > maxerr) maxerr = ad;
			}
		}
	}

	*l2_rel = sqrt(num / den);
	*linf = maxerr;
}


