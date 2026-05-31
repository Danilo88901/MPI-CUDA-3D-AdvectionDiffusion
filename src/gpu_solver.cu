#include"gpu_solver.cuh"
#include "device_launch_parameters.h"
__device__ double vanleer_gpu(double a, double b) {
	if (a * b <= 0)return 0.0;
	return 2 * a * b / (a + b);
}
__global__ void optimized_gpu_kernel(int height, int width, int lentgh, const double* __restrict__ phi, const double* __restrict__ u,
	const double* __restrict__ v, const double* __restrict__ w,
	double* __restrict__ rhs_phi, double nu, double dx_reverse, double dy_reverse, double dz_reverse) {
	extern __shared__ double data[];
	int i = blockDim.z * blockIdx.z + threadIdx.z;
	int j = blockDim.y * blockIdx.y + threadIdx.y;
	int k = blockDim.x * blockIdx.x + threadIdx.x;
	int tid = threadIdx.z * (blockDim.x * blockDim.y) + threadIdx.y * blockDim.x + threadIdx.x;
	int bsize = blockDim.x * blockDim.y * blockDim.z;
	int s_height = (blockDim.z + 4);
	int s_width = (blockDim.y + 4);
	int s_length = (blockDim.x + 4);
	//int slice_area = s_width * s_length;
	//int volume=s_height * s_width * s_length


	for (int idx = tid; idx < s_height * s_width * s_length; idx += bsize) {
		int s_i = idx / (s_width * s_length);   // Сколько целых слоев по Z прошли
		int rem = idx % (s_width * s_length);   // Остаток на плоскость Y-X

		// ВАЖНО: делим остаток на длину строки (X), чтобы узнать строку (Y)
		int s_j = rem / s_length;     // Было rem / s_width -> ИСПРАВЛЕНО
		int s_k = rem % s_length;     // Было rem % s_width -> ИСПРАВЛЕНО

		// Вычисляем глобальные координаты для этой ячейки shared memory
		int raw_i = blockDim.z * blockIdx.z + s_i - 2;
		int raw_j = blockDim.y * blockIdx.y + s_j - 2;
		int raw_k = blockDim.x * blockIdx.x + s_k - 2;

		// Честная периодика (добавляем размер, чтобы не уйти в отрицательный остаток)
		int g_i = (raw_i + height) % height;
		int g_j = (raw_j + width) % width;
		int g_k = (raw_k + lentgh) % lentgh;
		data[idx] = phi[g_i * width * lentgh + g_j * lentgh + g_k];


	}
	__syncthreads();
	if (i < 2 || i >= height - 2 || j < 2 || j >= width - 2 || k < 2 || k >= lentgh - 2)return;
	int idx_main = i * width * lentgh + j * lentgh + k;


	int idx = (threadIdx.z + 2) * s_width * s_length + (threadIdx.y + 2) * s_length + threadIdx.x + 2;
	//int e = idx + 1;
	//int west = idx - 1;
	int n = idx - s_length;
	int s = idx + s_length;
	int top = idx - (s_width * s_length);
	int bot = idx + (s_width * s_length);

	double rhs_main = 0.0;
	rhs_main += nu * ((data[idx + 1] - 2 * data[idx] + data[idx - 1]) * (dx_reverse * dx_reverse) +
		(data[n] - 2 * data[idx] + data[s]) * (dy_reverse * dy_reverse) +
		(data[bot] - 2 * data[idx] + data[top]) * (dz_reverse * dz_reverse));

	double ue = 0.5 * (u[idx_main] + u[idx_main + 1]);
	double uw = 0.5 * (u[idx_main] + u[idx_main - 1]);

	double vn = 0.5 * (v[idx_main] + v[idx_main - lentgh]);
	double vs = 0.5 * (v[idx_main] + v[idx_main + lentgh]);

	double wtop = 0.5 * (w[idx_main] + w[idx_main - (width * lentgh)]);
	double wbot = 0.5 * (w[idx_main] + w[idx_main + (width * lentgh)]);

	double vanleer_x = vanleer_gpu(data[idx] - data[idx - 1], data[idx + 1] - data[idx]);
	//EAST
	double phiL_east = data[idx] + 0.5 * vanleer_x;
	double phiR_east = data[idx + 1] - 0.5 * vanleer_gpu(data[idx + 1] - data[idx], data[idx + 2] - data[idx + 1]);
	//double ae = fabs(ue);
	double fe = 0.5 * (ue * phiL_east + ue * phiR_east) - 0.5 * fabs(ue) * (phiR_east - phiL_east);
	//WEST
	double phiL_west = data[idx - 1] + 0.5 * vanleer_gpu(data[idx - 1] - data[idx - 2], data[idx] - data[idx - 1]);
	double phiR_west = data[idx] - 0.5 * vanleer_x;

	//double aw = fabs(uw);
	double fw = 0.5 * (uw * phiL_west + uw * phiR_west) - 0.5 * fabs(uw) * (phiR_west - phiL_west);

	//y part
//north
	double vanleer_y = vanleer_gpu(data[idx] - data[n], data[s] - data[idx]);

	double phiL_north = data[n] + 0.5 * vanleer_gpu(data[n] - data[n - s_length], data[idx] - data[n]);
	double phiR_north = data[idx] - 0.5 * vanleer_y;

	//double an = fabs(vn);
	double fn = 0.5 * (vn * phiL_north + vn * phiR_north) - 0.5 * fabs(vn) * (phiR_north - phiL_north);

	//south
	//double as = fabs(vs);
	double phiL_south = data[idx] + 0.5 * vanleer_y;
	double phiR_south = data[s] - 0.5 * vanleer_gpu(data[s] - data[idx], data[s + s_length] - data[s]);

	double fs = 0.5 * (vs * phiL_south + vs * phiR_south) - 0.5 * fabs(vs) * (phiR_south - phiL_south);

	//z part
	double vanleer_z = vanleer_gpu(data[idx] - data[top], data[bot] - data[idx]);
	double phiL_top = data[top] + 0.5 * vanleer_gpu(data[top] - data[top - s_width * s_length], data[idx] - data[top]);
	double phiR_top = data[idx] - 0.5 * vanleer_z;

	//double at = fabs(wtop);
	double ft = 0.5 * (wtop * phiL_top + wtop * phiR_top) - 0.5 * fabs(wtop) * (phiR_top - phiL_top);

	//double ab = fabs(wbot);
	double phiL_bot = data[idx] + 0.5 * vanleer_z;
	double phiR_bot = data[bot] - 0.5 * vanleer_gpu(data[bot] - data[idx], data[bot + s_width * s_length] - data[bot]);

	double fb = 0.5 * (wbot * phiL_bot + wbot * phiR_bot) - 0.5 * fabs(wbot) * (phiR_bot - phiL_bot);

	rhs_main += -((fe - fw) * dx_reverse + (fs - fn) * dy_reverse + (fb - ft) * dz_reverse);
	rhs_phi[idx_main] = rhs_main;
}
