#include "mesh.hh"
#include <omp.h>

int main() {
	
	// Create meshes and initialize acceleration
	mesh_param par(0.05, 0.03, 0.001, false, false);
	mesh_rk4 mp(par, "sh48_211x211.bin");
	mp.lump=false;
	mp.CG=false;
	mp.PCG=true;

	// Centralize and scale the mesh
	double wx, wy, wz;
	mp.centralize(wx, wy, wz);

	double start_time = omp_get_wtime();

	/* Set up springs. Assign lengths manually, since initializing out of
	equilibrium configuration. */
	mp.setup_springs();

	printf("Triangle table and mass matrix setup time: %g seconds\n",omp_get_wtime()-start_time);
	start_time = omp_get_wtime();
	
	// Stretch perturbation in-plane
	double eps = .4;
	for(double *p=mp.pts,*pe=p+3*mp.n;p<pe;p+=3) {
		p[0]=(1+eps)*p[0]; p[1]=(1+eps)*p[1];
	}
	
	// Pure shear perturbation
	/*for (double* p = mp.pts, *pe = p + 3 * mp.n; p < pe; p += 3) {
		p[0] = (1 - eps) * p[0]; p[1] = (1 + eps) * p[1];
	}*/

	// Setup the output directory and allocate memory for integrator.
	mp.setup_output_dir("PCG_Test.odr");

	mp.solve_adaptive(75, 1e-4, 1e-4, false, 100); 

	printf("Elapsed solution time: %g seconds\n",omp_get_wtime()-start_time);
}