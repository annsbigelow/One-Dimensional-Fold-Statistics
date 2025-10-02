#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "mesh.hh"

int main() {

    // Read in the mesh
    int len=100;
    int inc=20;
    char buf[128];
    sprintf(buf,"sheet_%dx%d.bin",len,len);
    //sprintf(buf,"rsheet_2500_2.bin");
    mesh_param par(0.5,0.01,0,0.2,false,true,0.001,1.3);
    mesh_rk4 mp(par,buf);

    // Centralize and scale the mesh
    double wx,wy,wz,w;
    mp.centralize(wx,wy,wz);
    w=1./sqrt(wx+wy);
    for(double *p=mp.pts,*pe=p+3*mp.n;p<pe;p+=3) {
        *p*=w;p[1]*=w;
        *p*=50;p[1]*=50;
    }

    // Set up springs.
    mp.setup_springs();

    // Add Gaussian displacement around the central location
    for(double *p=mp.pts,*pe=p+3*mp.n;p<pe;p+=3) {
        p[2]+=-0.1+0.02*exp(-0.1*(*p*(*p)+p[1]*p[1]));
	}

	// Copy initial positions and randomize shrink rates.
	if (mp.shrink) mp.init_shrink(.0005,.002);

    // Add external potential.
    //ep_spherical eps(80,10,5000,0.0002);
    //mp.add(&eps);

    /* Increment spring rest lengths by random amount in a chosen interval
    (Oppenheimer paper). */
    //srand(12);
    //mp.perturb_springs(1,1.+0.01*inc);

    // Setup the output directory.
    mp.setup_output_dir("srun_h.odr");

    // Evolve in time with equally spaced output
	mp.solve_adaptive(1000, 1e-3, 1e-3, false, 250);
}
