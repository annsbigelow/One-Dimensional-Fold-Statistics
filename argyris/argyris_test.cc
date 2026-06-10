#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>

#include "mesh.hh"

int main() {

	// Create mesh and initialize acceleration
	mesh_param par(0.05, 0.03, 0.001, false, false);
	mesh_rk4 mp(par, "../crumple/sh48_43x43.bin");

	// Centralize and scale the mesh
	//TODO: the dofs in centralize() are off.
	double wx, wy, wz;
	mp.centralize(wx, wy, wz);

	mp.PCG = false; 
	mp.CG = false;
	double start_time = omp_get_wtime();
	mp.setup_springs();
	printf("Triangle table and mass matrix setup time: %g seconds\n",omp_get_wtime()-start_time);

	// Apply perturbation in z-direction
	...
	// Setup the output directory and allocate memory for integrator.
	mp.setup_output_dir("Argyris_Test.odr");

	// Solve!
	//mp.solve_adaptive(75, 1e-4, 1e-4, false, 100);
	//printf("Elapsed solution time: %g seconds\n", omp_get_wtime() - start_time);
}
