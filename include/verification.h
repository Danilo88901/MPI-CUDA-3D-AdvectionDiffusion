#pragma once
void fill_velocity(double* u, double* v, double* w,
	int nx, int ny, int nz,
	double dx, double dy, double dz,
	double /*xc*/, double /*yc*/,
	double u0, double v0, double w0);
void init(double* phi, int nx, int ny, int nz,
		double dx, double dy, double dz,
		double xc, double yc, double zc, double /*sigma0*/);
void exact_solution(double* phi_exact, int nx, int ny, int nz,
	double dx, double dy, double dz,
	double t,
	double xc, double yc, double zc,
	double u0, double v0, double w0,
	double nu);
void compute_error(double* phi, double* phi_exact,
	int nx, int ny, int nz,
	double* l2_rel, double* linf);
