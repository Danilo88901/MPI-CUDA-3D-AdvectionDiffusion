# 3D Advection–Diffusion Solver (CPU / MPI / CUDA)

A 3D solver for the advection–diffusion equation with periodic boundary conditions, explicit second-order SSPRK2 time integration, and numerical verification against an analytical reference solution.

The repository contains:

- a serial CPU implementation,
- an MPI-based CPU implementation with 3D Cartesian domain decomposition,
- a CUDA implementation,
- an optimized CUDA kernel using shared memory and `__restrict__`,
- verification utilities,
- ParaView output generation (`.vti` / `.pvd`),
- runtime and error measurements.

---

## Problem Statement

This project solves the 3D advection–diffusion equation on a structured Cartesian grid with periodic boundary conditions.

The numerical solution is compared against a manufactured analytical solution built from a periodic trigonometric field with diffusion-advection evolution.

---

## Numerical Method

The solver uses:

- a second-order explicit SSPRK2 time integrator,
- a Rusanov / Van Leer-type flux reconstruction for advection,
- central finite differences for diffusion,
- periodic boundary conditions in all three spatial directions.

The solution is verified using:

- relative L2 error,
- Linf error,
- mass conservation checks.

---

## Implementations

### Serial CPU baseline
A straightforward reference implementation without MPI.

### MPI CPU version
The MPI version uses:

- `MPI_Dims_create` to build a 3D process grid,
- `MPI_Cart_create` for Cartesian communication topology,
- halo exchange between neighboring subdomains,
- `MPI_Type_create_subarray` for sending boundary slabs,
- `MPI_Gather` for global field reconstruction during verification/output.

This implementation uses a **static block decomposition**.

**Important limitation:**  
The current MPI decomposition assumes that each physical dimension is divisible by the corresponding process-grid dimension. If this condition is not satisfied, remainder cells are not redistributed, so the current MPI implementation is only correct for compatible grid/process configurations.

### CUDA version
The CUDA version includes:

- a shared-memory optimized stencil kernel,
- precomputed inverse grid spacing,
- `__restrict__` pointers for improved compiler optimization,
- periodic boundary handling,
- separate kernels for SSPRK2 stages and boundary updates.

---

## Verification

The code includes a manufactured analytical solution, which makes it possible to verify correctness by comparing:

- numerical solution,
- exact solution,
- L2 relative error,
- Linf error.

The CUDA and MPI implementations are checked against the same reference solution. For compatible decompositions, the numerical error metrics agree.

---

## Performance Notes

Benchmarks are reported for cubic grids:

- 32³
- 64³
- 128³
- 256³

The performance comparison includes:

- serial CPU,
- MPI CPU,
- CUDA,
- optimized CUDA with shared memory.

Timing should be interpreted carefully:

- some runs include output generation,
- some runs include verification,
- some runs include profiling-oriented configurations.

For pure compute comparison, benchmark runs should be separated from I/O and postprocessing.

---

## Output

The project can export results for visualization in ParaView:

- `.vti` files for scalar fields,
- `.pvd` collection files for time series.

This is used to generate snapshots and animations of the evolving field.

---

## Repository Structure

├── include/
│   ├── cpu_solver.hpp
│   ├── gpu_solver.cuh
│   ├── verification.hpp
│   ├── boundary.hpp
│   └── io.hpp
│
├── src/
│   ├── cpu_solver.cpp
│   ├── gpu_solver.cu
│   ├── verification.cpp
│   ├── boundary.cpp
│   └── io.cpp
│
├── results/
│   ├── figures/
│   ├── plots/
│   ├── animations/
│   └── profiling/
│
├── CMakeLists.txt
└── README.md
