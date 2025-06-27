#include "brusselator.hh"

#include <cstdio>
#include <cmath>

#include "omp.h"

// The duration for the ODE convergence test
const double dur=20.;

int main() {

    // Compute reference solution
    brus_rk4 br;
    br.init();
    br.solve_adaptive(dur,3e-15,3e-15);
    double ref0=br.q[0],ref1=br.q[1];

    // Allocate space for storing accuracy
    const int N=101;
    int fcount[2*N];
    double err[2*N],t0=omp_get_wtime();

    // Loop over a variety of grid sizes (given by i index) and integration
    // schemes (given by j index)
#pragma omp parallel for schedule(dynamic) collapse(2)
    for(int i=0;i<N;i++) for(int j=0;j<2;j++) {

        // Compute the number of steps according to a power law. Minimum
        // steps of 100, and maximum steps of 100000.
        brus_rk4 *sb=new brus_rk4();
        sb->init();

        // Perform integration and compute the difference to the reference
        // solution
        if(j==0) {
            double tol=0.003*pow(0.1,0.1*i);
            sb->solve_adaptive(dur,tol,tol);
        } else sb->solve_fixed(dur,int(10*pow(10,0.03*i)));
        double dy0=ref0-sb->q[0],
               dy1=ref1-sb->q[1];
        err[j*N+i]=sqrt(dy0*dy0+dy1*dy1);
        fcount[j*N+i]=sb->fcount;
        delete sb;
    }

    // Output the results in a Gnuplot-readable format
    printf("# Time taken : %g ms\n",1e3*(omp_get_wtime()-t0));
    for(int j=0;j<2;j++) {
        FILE *fp=fopen(j==0?"conv_adapt.dat":"conv_fixed.dat","w");
        if(fp==NULL) {
            fputs("Error opening file",stderr);
            return 1;
        }
        for(int i=0;i<N;i++) fprintf(fp,"%d %g\n",fcount[j*N+i],err[j*N+i]);
        fclose(fp);
    }
}
