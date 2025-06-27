#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "mesh.hh"

int main() {

    // Create mesh and initialize acceleration
    mesh_param par(0.05,0.02,0.001,false,false);
    mesh_rk4 mp(par,"disloc_30.bin");

    // Centralize and scale the mesh
    double wx,wy,wz;
    mp.centralize(wx,wy,wz);

    /* Set up springs. Assign lengths manually, since initializing out of
    equilibrium configuration. */
    mp.setup_springs();
    for(double *p=mp.reg;p<mp.reg+mp.ns;p++) *p=1.;
    for(double *p=mp.ref;p<mp.ref+mp.nh;p++) *p=2./sqrt(3);

    // Perturb the z positions
    for(double *p=mp.pts,*pe=p+3*mp.n;p<pe;p+=3)
       p[2]+=1e-3*sin(*p*(*p)-p[1]*p[1])*exp(-(*p*(*p)+p[1]*p[1]));

    // Add centering potential
    ep_centering epc(0.001);
    mp.add(&epc);

    // Setup the output directory and allocate memory for integrator.
    mp.setup_output_dir("brun_d.out");
    mp.allocate(6*mp.n);

    // Carry out the simulation, reporting the time to compute the timesteps
    // between each output frame
    mp.solve(mp.pts,0,500000,5,51); // this is equivalent to running section below:
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

    // Now do some minimization steps
    for(int i=51;i<=70;i++) {
        mp.minimize_energy(mp.t);
        mp.output_positions(i);
    }
}
