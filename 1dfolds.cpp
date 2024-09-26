#include <iostream>
using namespace std;
#include "declarations.h"
#include <omp.h>
#include <gsl/gsl_rng.h>


int main() {
    double t0 = omp_get_wtime();
    folds_stats fs; 
    
    fs.logavg();
    printf("Time elapsed is %g s\n", omp_get_wtime() - t0);

    return 0;
}


