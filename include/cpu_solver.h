#pragma once
double vanleer(double a, double b);
void advection_diffusion(int height, int width, int lentgh, double* phi, const double* u,
	const double* v, const double* w,
	double* rhs_phi, double nu, double dx_reverse, double dy_reverse, double dz_reverse,
	int lz, int ly, int  lx);
