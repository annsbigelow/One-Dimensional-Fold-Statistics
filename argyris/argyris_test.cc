#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>

#include "mesh.hh"

int main() {

	// Create mesh and initialize acceleration
	mesh_param par(0.05, 0.03, 0.001, false, false);
	mesh_rk4 mp(par, "sh48_6x6.bin");

	// Centralize and scale the mesh
	double wx, wy, wz;
	mp.centralize(wx, wy, wz);

	mp.PCG = false; 
	mp.CG = false;
	double start_time = omp_get_wtime();
	mp.setup_springs();
	printf("Triangle table and mass matrix setup time: %g seconds\n",omp_get_wtime()-start_time);

	// Apply perturbation in z-direction
	// Add Gaussian displacement around the central location
    for(double *p=mp.pts,*xy=mp.xyz; p<mp.pts+6*mp.n; p+=6,xy+=3) {
        *p+=-0.1+0.02*exp(-0.1*(*xy*(*xy)+xy[1]*xy[1]));
	}

	// Setup the output directory and allocate memory for integrator.
	mp.setup_output_dir("Argyris_Test.odr");

	// Solve!
	mp.solve_adaptive(2, 1e-4, 1e-4, false, 2);
	printf("Elapsed solution time: %g seconds\n", omp_get_wtime() - start_time);
}