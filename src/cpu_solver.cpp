#include"cpu_solver.h"
#include<cmath>
double vanleer(double a, double b) {
	if (a * b <= 0)return 0.0;
	return 2 * a * b / (a + b);
}
void advection_diffusion(int height, int width, int lentgh, double* phi, const double* u,
	const double* v, const double* w,
	double* rhs_phi, double nu, double dx_reverse, double dy_reverse, double dz_reverse,
	int lz, int ly, int  lx) {
	for (int i = 2; i < height + 2; i++) {
		for (int j = 2; j < width + 2; j++) {
			for (int k = 2; k < lentgh + 2; k++) {
				int idx = i * ly * lx + j * lx + k;
				int e = idx + 1;
				int west = idx - 1;
				int n = i * ly * lx + (j - 1) * lx + k;
				int s = i * ly * lx + (j + 1) * lx + k;;
				int top = (i - 1) * ly * lx + j * lx + k;
				int bot = (i + 1) * ly * lx + j * lx + k;

				double rhs_main = 0.0;
				rhs_main += nu * ((phi[e] - 2 * phi[idx] + phi[west]) * (dx_reverse * dx_reverse) +
					(phi[n] - 2 * phi[idx] + phi[s]) * (dy_reverse * dy_reverse) +
					(phi[bot] - 2 * phi[idx] + phi[top]) * (dz_reverse * dz_reverse));

				double ue = 0.5 * (u[idx] + u[e]);
				double uw = 0.5 * (u[idx] + u[west]);

				double vn = 0.5 * (v[idx] + v[n]);
				double vs = 0.5 * (v[idx] + v[s]);

				double wtop = 0.5 * (w[idx] + w[top]);
				double wbot = 0.5 * (w[idx] + w[bot]);

				//x part
				double vanleer_x = vanleer(phi[idx] - phi[west], phi[e] - phi[idx]);
				//EAST
				double phiL_east = phi[idx] + 0.5 * vanleer_x;
				double phiR_east = phi[e] - 0.5 * vanleer(phi[e] - phi[idx], phi[e + 1] - phi[e]);
				double ae = fabs(ue);
				double fe = 0.5 * (ue * phiL_east + ue * phiR_east) - 0.5 * ae * (phiR_east - phiL_east);
				//WEST
				double phiL_west = phi[west] + 0.5 * vanleer(phi[west] - phi[west - 1], phi[idx] - phi[west]);
				double phiR_west = phi[idx] - 0.5 * vanleer_x;

				double aw = fabs(uw);
				double fw = 0.5 * (uw * phiL_west + uw * phiR_west) - 0.5 * aw * (phiR_west - phiL_west);

				//y part
				//north
				double vanleer_y = vanleer(phi[idx] - phi[n], phi[s] - phi[idx]);

				double phiL_north = phi[n] + 0.5 * vanleer(phi[n] - phi[n - lx], phi[idx] - phi[n]);
				double phiR_north = phi[idx] - 0.5 * vanleer_y;

				double an = fabs(vn);
				double fn = 0.5 * (vn * phiL_north + vn * phiR_north) - 0.5 * an * (phiR_north - phiL_north);

				//south
				double as = fabs(vs);
				double phiL_south = phi[idx] + 0.5 * vanleer_y;
				double phiR_south = phi[s] - 0.5 * vanleer(phi[s] - phi[idx], phi[s + lx] - phi[s]);

				double fs = 0.5 * (vs * phiL_south + vs * phiR_south) - 0.5 * as * (phiR_south - phiL_south);

				//z part
				double vanleer_z = vanleer(phi[idx] - phi[top], phi[bot] - phi[idx]);
				double phiL_top = phi[top] + 0.5 * vanleer(phi[top] - phi[top - ly * lx], phi[idx] - phi[top]);
				double phiR_top = phi[idx] - 0.5 * vanleer_z;

				double at = fabs(wtop);
				double ft = 0.5 * (wtop * phiL_top + wtop * phiR_top) - 0.5 * at * (phiR_top - phiL_top);

				double ab = fabs(wbot);
				double phiL_bot = phi[idx] + 0.5 * vanleer_z;
				double phiR_bot = phi[bot] - 0.5 * vanleer(phi[bot] - phi[idx], phi[bot + lx * ly] - phi[bot]);

				double fb = 0.5 * (wbot * phiL_bot + wbot * phiR_bot) - 0.5 * ab * (phiR_bot - phiL_bot);

				rhs_main += -((fe - fw) * dx_reverse + (fs - fn) * dy_reverse + (fb - ft) * dz_reverse);
				rhs_phi[idx] = rhs_main;
			}
		}
	}
}
