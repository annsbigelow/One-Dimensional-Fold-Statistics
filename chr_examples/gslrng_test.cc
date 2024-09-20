#include <cstdio>

#include "gsl/gsl_rng.h"

int main() {

    gsl_rng* rng=gsl_rng_alloc(gsl_rng_taus2);

    gsl_rng_set(rng,1729);

    for(int i=0;i<10;i++) {
        printf("%g\n",gsl_rng_uniform(rng));
    }

    gsl_rng_free(rng);
}
