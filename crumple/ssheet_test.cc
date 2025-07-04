#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "mesh.hh"

int main() {

    // Read in the mesh
    int len=201;
    int inc=20;
    char buf[128];
    sprintf(buf,"sheet_%dx%d.bin",len,len);
    //sprintf(buf,"rsheet_2500_2.bin");
    mesh_param par(0.05,0.02,0.005,false,false);
    mesh_rk4 mp(par,buf); // use an RK4 integrator
    //mesh_rkfsal mp(par,buf); // use an RK-FSAL adaptive integrator

    // Centralize and scale the mesh
    double wx,wy,wz,w;
    mp.centralize(wx,wy,wz);
    w=1./sqrt(wx+wy);
    for(double *p=mp.pts,*pe=p+3*mp.n;p<pe;p+=3) {
        //*p*=w;p[1]*=w;
        //*p*=50;p[1]*=50;
    }

    // Set up springs.
    mp.setup_springs();
    mp.allocate(6*mp.n);

    // Add Gaussian displacement around the central location
    for(double *p=mp.pts,*pe=p+3*mp.n;p<pe;p+=3)
    p[2]+=-0.1+0.0001*exp(-10*(*p*(*p)+p[1]*p[1]));

    // Add external potential.
    ep_spherical eps(40,0.000001);
    mp.add(&eps);

    /* Increment spring rest lengths by random amount in a chosen interval
    (Oppenheimer paper). */
    //srand(12);
    //mp.perturb_springs(1,1.+0.01*inc);

    // Setup the output directory.
    mp.setup_output_dir("srun_b.out");

    // Evolve in time with equally spaced output
    mp.solve_adaptive(1,1e-4,1e-4,false,101);
}
