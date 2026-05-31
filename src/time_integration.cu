#include"time_integration.cuh"
#include<cmath>
#include "device_launch_parameters.h"
__global__ void ssprk2_stage1(int height, int width, int length, double* phi_tmp, double* rhs_phi, double* phi, double dt) {
	int i = blockDim.z * blockIdx.z + threadIdx.z;
	int j = blockDim.y * blockIdx.y + threadIdx.y;
	int k = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < 2 || i >= height - 2 || j < 2 || j >= width - 2 || k < 2 || k >= length - 2)return;
	int idx = i * width * length + j * length + k;
	phi_tmp[idx] = phi[idx] + dt * rhs_phi[idx];
}
__global__ void ssprk2_stage2(int height, int width, int length, double* phi_tmp, double* rhs_tmp, double* phi, double dt) {
	int i = blockDim.z * blockIdx.z + threadIdx.z;
	int j = blockDim.y * blockIdx.y + threadIdx.y;
	int k = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < 2 || i >= height - 2 || j < 2 || j >= width - 2 || k < 2 || k >= length - 2)return;
	int idx = i * width * length + j * length + k;
	phi[idx] =
		0.5 * phi[idx]
		+ 0.5 * (phi_tmp[idx] + dt * rhs_tmp[idx]);
}
double compute_dt(int height, int width, int length,
	double* u, double* v, double* w,
	double dx, double dy, double dz,
	double nu,
	double cfl)
{
	double max_adv = 0.0;

	// ==========================================
	// Advection constraint
	// ==========================================
	for (int i = 2; i < height - 2; i++) {
		for (int j = 2; j < width - 2; j++) {
			for (int k = 2; k < length - 2; k++) {

				int idx = i * width * length + j * length + k;

				double adv =
					fabs(u[idx]) / dx +
					fabs(v[idx]) / dy +
					fabs(w[idx]) / dz;

				if (adv > max_adv)
					max_adv = adv;
			}
		}
	}

	double dt_adv = 1e30;

	if (max_adv > 1e-14)
		dt_adv = cfl / max_adv;

	// ==========================================
	// Diffusion constraint
	// ==========================================
	double inv_dx2 = 1.0 / (dx * dx);
	double inv_dy2 = 1.0 / (dy * dy);
	double inv_dz2 = 1.0 / (dz * dz);

	double dt_diff = 1e30;

	if (nu > 0.0) {
		dt_diff =
			0.5 /
			(nu * (inv_dx2 + inv_dy2 + inv_dz2));
	}

	// ==========================================
	// Final timestep
	// ==========================================
	double dt = std::min(dt_adv, dt_diff);

	return dt;
}
