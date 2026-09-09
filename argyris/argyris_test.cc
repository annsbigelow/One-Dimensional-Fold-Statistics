#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>

#include "mesh.hh"

int main() {
	// Create mesh and initialize acceleration
	const bool r_all_dofs=false;
	const bool wr_all_dofs=true;
	const double kb=.0001;
	mesh_param par(0.05, 0.02, kb, false,true,r_all_dofs,wr_all_dofs);
	mesh_rk4 mp(par, "sh48_7x7.bin");
	if(mp.fix_boundary==true) printf("sup\n");

	float s=1; int nx=7;
	printf("The side length is set to %f and the number of nodes in one direction is %d.\n",s,nx);

	// Centralize and scale the mesh
	double wx, wy, wz;
	mp.centralize(wx, wy, wz);

	// Setup the output directory and allocate memory for integrator.
	mp.setup_output_dir("Argyris_Test.odr");
	double start_time = omp_get_wtime();
	mp.setup_springs(); mp.setup_fem();
	mp.build_matrices();
	printf("Mass matrix and stiffness matrix setup time: %g seconds\n",omp_get_wtime()-start_time);


	// Apply perturbation in z-direction
	mp.Gauss_displacement(s,nx);
	//mp.mesh_print_last_step();

	// Solve!
	mp.solve_fixed(30,300000, true); printf("Elapsed solution time: %g seconds\n",
										omp_get_wtime() - start_time);
	//mp.solve_adaptive(1, 1e-4, 1e-4, false, 1);
}