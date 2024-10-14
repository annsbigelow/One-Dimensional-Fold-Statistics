#include <iostream>
using namespace std;
#include "declarations.h"
#include <omp.h>
#include <gsl/gsl_rng.h>


int main() {
    double t0 = omp_get_wtime();
    
  //  omp_set_num_threads(4);
//#pragma omp parallel
    //{
        int s_divide = 200;
    //    //int th = omp_get_thread_num();
        folds_stats fs(15); // if you want new sequence of random nums, change fs() input.
    //    fs.segdens(s_divide);

   // }
        fs.altFold_segdens(s_divide);

    printf("Time elapsed is %g s\n", omp_get_wtime() - t0);

    return 0;
}


