#include <iostream>
using namespace std;
#include "declarations.h"
#include <omp.h>
#include "gsl/gsl_rng.h"


int main() {
    double t0 = omp_get_wtime();
    
    folds_stats fs; 

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus2); // Create instance of the Tausworthe generator
    gsl_rng_set(rng, 1729); // Seed the above generator
    
    fs.logavg();
    printf("Time elapsed is %g s\n", omp_get_wtime() - t0);

    gsl_rng_free(rng);

    return 0;
}


