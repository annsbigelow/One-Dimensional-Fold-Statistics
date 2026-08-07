#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>

#include "mesh.hh"

int main() {
	// Create mesh and initialize acceleration
	mesh_param par(0.05, 0.03, 0.001, false, false);
	mesh_rk4 mp(par, "sh48_3x3.bin");

	// Centralize and scale the mesh
	double wx, wy, wz;
	mp.centralize(wx, wy, wz);
	printf("Centralized\n");

	// Setup the output directory and allocate memory for integrator.
	mp.setup_output_dir("Argyris_Test.odr");
	double start_time = omp_get_wtime();
	mp.setup_springs();
	mp.build_matrices();
	//printf("Mass matrix and stiffness matrix setup time: %g seconds\n",omp_get_wtime()-start_time);

	// Apply perturbation in z-direction
	mp.linear_gradient();

	mp.debug();

	// Solve!
	//mp.solve_fixed(1e-2, 2, true);
	//mp.solve_adaptive(1, 1e-4, 1e-4, false, 1);
	//printf("Elapsed solution time: %g seconds\n", omp_get_wtime() - start_time);
}