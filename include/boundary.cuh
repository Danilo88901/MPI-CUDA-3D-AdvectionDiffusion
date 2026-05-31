#pragma once
#include <cuda_runtime.h>
#include <cuda.h>
__global__ void boundary_conditions_x(int height, int width, int length, double* phi);
__global__ void boundary_conditions_y(int height, int width, int length, double* phi);
__global__ void boundary_conditions_z(int height, int width, int length, double* phi);
void bc_periodic(int height, int width, int length, double* phi);
