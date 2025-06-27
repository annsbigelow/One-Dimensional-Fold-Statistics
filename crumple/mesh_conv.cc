#include <cstdio>
#include <cmath>
#include <cstring>

#include "mesh.hh"

int main() {

    // Maximum time to integrate to
    double tmax=10,d,h,nor;

    // Initial number of integration steps
    int n=10,nref=2000,i;
    mesh_param par(0.05,0.02,false);
    mesh_rk4 mp(par,"sheet_101x101.bin");
    mp.setup_springs();
    mp.allocate(6*mp.n);
    size_t sz=sizeof(double)*mp.dof;

    // Displace a section of the mesh
    double mid=50,dx,dy;
    for(double *p=mp.pts,*pe=p+3*mp.n;p<pe;p+=3) {

        // Compute distance to central location
        dx=*p-mid;
        dy=p[1]-mid;

        // Add Gaussian displacement around the central location
        p[2]+=0.5*exp(-0.1*(dx*dx+dy*dy));
    }

    // Store this initial state
    double *ipts=new double[mp.dof];
    memcpy(ipts,mp.pts,sz);

    // Compute a reference solution and store it
    h=tmax/nref;
    for(i=0;i<nref;i++) mp.fixed_step(h);
    double *refsol=new double[mp.dof];
    memcpy(refsol,mp.pts,sz);

    // Try simulating with progressively higher resolution
    while(n<1000) {

        // Compute timestep, initialize ODE system, and integrate
        h=tmax/n;
        memcpy(mp.pts,ipts,sz);
        mp.t=0;
        for(i=0;i<n;i++) mp.fixed_step(h);

        // Print difference between numerical and exact solution
        nor=0.;
        for(int i=0;i<mp.dof;i++) {
            d=refsol[i]-mp.pts[i];
            nor+=d*d;
        }
        printf("%.14g %.14g\n",h,sqrt(nor));

        // Increase resolution
        n+=n/2;
    }

    // Free the dynamically allocated memory
    delete [] refsol;
    delete [] ipts;
}
