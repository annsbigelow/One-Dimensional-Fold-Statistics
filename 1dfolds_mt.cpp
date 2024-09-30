#include <cstdio>
#include <iostream>
using namespace std;
#include "declarations.h"
#include <gsl/gsl_rng.h>

#include <omp.h>

int main() {
    double t0 = omp_get_wtime();
    
    //fs.log_fixedn();
    int n = 160;
#pragma omp parallel for
    for(int j=0;j<16;j++) {
        int max_t=omp_get_max_threads(),th=omp_get_thread_num();
        printf("j=%d is processed by %d of %d threads\n",th,max_t);

        //folds_stats fs; 
        //fs.segdens(n);
    }


    printf("Time elapsed is %g s\n", omp_get_wtime() - t0);
    return 0;
}


