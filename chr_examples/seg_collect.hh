#ifndef SEG_COLLECT_HH
#define SEG_COLLECT_HH

#include "gsl/gsl_rng.h"

/** The initial memory allocation for the seg_collect class. */
const int init_mem=4;

/** The maximum allowed memory allocation, to prevent a runaway number of requests. */
const int max_mem=16777216;

class seg_collect {
    public:
        /** The number of memory slots used. */
        int n;
        /** The total number of allocated memory slots. */
        int mem;
        /** A pointer to the allocated memory. */
        double* s;
        /** The class constructor dynamically allocates memory for the initial
         * number of slots. */
        seg_collect() : n(0), mem(init_mem),
            s(new double[2*mem]), rng(gsl_rng_alloc(gsl_rng_taus2)) {}
        /** The class destructor frees the dynamically allocated memory. */
        ~seg_collect() {
            gsl_rng_free(rng);
            delete [] s;
        }
        inline void seed(int k) {
            gsl_rng_set(rng,k);
        }
        void print_segs();
        void add(double x,double y);
        void add_memory();
    private:
        gsl_rng *rng;
};

#endif
