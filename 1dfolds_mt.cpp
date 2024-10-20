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
#pragma omp parallel
    {
        int th=omp_get_thread_num();
        folds_stats fs(th);
        fs.segdens(n);
    }

    printf("Time elapsed is %g s\n", omp_get_wtime() - t0);
    return 0;
}
