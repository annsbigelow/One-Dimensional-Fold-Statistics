#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>

#include "mesh.hh"

int main() {
	// Set parameters
	bool adaptive=true, fix_boundary=true;
	float s=.8; int nx=7;
	char buf[50], buf1[50];
	const bool r_all_dofs=false, wr_all_dofs;
	if(adaptive) wr_all_dofs=false; else wr_all_dofs=true;
	const double K=.05,drag=.02,kb=.00001;
	printf("The side length of the coarse mesh is set to %f"
			"and the number of nodes in one direction is %d.\n",s,nx);

	// Create mesh
	mesh_param par(K,drag,kb,false,fix_boundary,r_all_dofs,wr_all_dofs,s);
	sprintf(buf,"./sheet_gen rec48 %f %d %d",s,nx,nx);
	std::system(buf);
	sprintf(buf1,"sh48_%dx%d.bin",nx,nx);
	mesh_rk4 mp(par,buf1);

	// Centralize and scale the mesh
	double wx, wy, wz;
	mp.centralize(wx, wy, wz);

	// Setup the output directory and allocate memory for integrator.
	// Build FEM matrices.
	if(adaptive) {
		mp.np=6; mp.output_refined=true;
	}
	else mp.output_refined=false;
	mp.setup_output_dir("dense_run.odr");
	double start_time = omp_get_wtime();
	mp.setup_springs(); mp.setup_fem(); mp.build_matrices();
	printf("Mass matrix and stiffness matrix setup time: %g seconds\n",omp_get_wtime()-start_time);

	// Apply perturbation in z-direction
	mp.Gauss_displacement(s,nx);
	//mp.mesh_print_last_step(); // Switch on/off to observe initial displacement

	// Solve!
	if(adaptive) {
		printf("Number of partitions for bendy triangle visualization is set to %d.\n",mp.np);
		mp.solve_adaptive(100, 1e-4, 1e-4, false, 50);
	}
	else {
		double t=1; int steps=1000;
		printf("The fixed timestep is equal to %g.\n",t/steps);
		mp.solve_fixed(t,steps,true); 
	}
	printf("Elapsed solution time: %g seconds\n",omp_get_wtime() - start_time);
}