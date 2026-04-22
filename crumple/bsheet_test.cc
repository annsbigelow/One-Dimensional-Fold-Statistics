#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "mesh.hh"

int main() {

    // Create mesh and initialize acceleration
    mesh_param par(0.05,0.03,0.001,false,false);
    mesh_rk4 mp(par,"sh48_101x101.bin");

    // Centralize and scale the mesh
    double wx,wy,wz;
    mp.centralize(wx,wy,wz);

    /* Set up springs. Assign lengths manually, since initializing out of
    equilibrium configuration. */
	mp.lump=true;
    mp.setup_springs();

    // XXX - Diagnostic routine to see triangles
    //mp.print_triangle_table();
    //return 1;

	// XXX - ref, reg are no longer setup
    //for(double *p=mp.reg;p<mp.reg+mp.ns;p++) *p=1.;
    //for(double *p=mp.ref;p<mp.ref+mp.nh;p++) *p=2./sqrt(3);
	
    // Stretch perturbation in-plane
	double eps=.4;
    for(double *p=mp.pts,*pe=p+3*mp.n;p<pe;p+=3) {
		p[0]=(1+eps)*p[0]; p[1]=(1+eps)*p[1];
	} 
	
	// Pure shear perturbation
	/*for (double *p=mp.pts, *pe=p+3*mp.n;p<pe;p+=3) {
		p[0]=(1-eps)*p[0]; p[1]=(1+eps)*p[1];
	}*/
	
	
    // Add centering potential
	// I don't think this does anything...
    ep_centering epc(0.001);
    mp.add(&epc);
	
    // Setup the output directory and allocate memory for integrator.
    mp.setup_output_dir("brun_d.out");
    //mp.allocate(6*mp.n); // XXX - Mesh constructor already allocates mem

    // Carry out the simulation, reporting the time to compute the timesteps
    // between each output frame
    mp.solve_adaptive(50,1e-4,1e-4,false,75); // this is equivalent to running section below:

    /*
    mp.output_positions(0);
    puts("Frame 0");
    double t0=wtime(),t1;
    for(int i=1;i<=50;i++) {
        for(int j=0;j<2000;j++) mp.step(mp.pts,5);
        mp.output_positions(i);
        printf("Frame %d [%8.5f s]\n",i,(t1=wtime())-t0);
        t0=t1;
    }
    */

	// XXX
    // Now do some minimization steps
    /*for(int i=51;i<=70;i++) {
        mp.minimize_energy(mp.t);
        mp.output_positions(i);
    }
	*/
}
