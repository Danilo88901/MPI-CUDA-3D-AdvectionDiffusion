#pragma once
#include <cuda_runtime.h>
#include <cuda.h>
__device__ double vanleer_gpu(double a, double b);
__global__ void optimized_gpu_kernel(int height, int width, int lentgh, const double* __restrict__ phi, const double* __restrict__ u,
	const double* __restrict__ v, const double* __restrict__ w,
	double* __restrict__ rhs_phi, double nu, double dx_reverse, double dy_reverse, double dz_reverse);
