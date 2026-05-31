#pragma once
#include <cuda_runtime.h>
#include <cuda.h>
__global__ void ssprk2_stage1(int height, int width, int length, double* phi_tmp, double* rhs_phi, double* phi, double dt);
__global__ void ssprk2_stage2(int height, int width, int length, double* phi_tmp, double* rhs_tmp, double* phi, double dt);
void ssprk2_stage1_cpu(int height, int width, int lentgh, int lz, int ly, int  lx, double* phi_tmp, 
	double* rhs_phi, double* phi, double dt);
void ssprk2_stage2_cpu(int height, int width, int length, int lz, int ly, int lx, double* phi_tmp,
	double* rhs_tmp, double* phi, double dt);

double compute_dt(int height, int width, int length,
	double* u, double* v, double* w,
	double dx, double dy, double dz,
	double nu,
	double cfl);
