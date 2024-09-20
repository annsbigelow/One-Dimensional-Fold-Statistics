#include "seg_collect.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>

/** Prints a list of the entries in the allocated slots. */
void seg_collect::print_segs() {
    for(int i=0;i<n;i++)
        printf("%g %g\n",s[2*i],s[2*i+1]);
}

/** Adds an entry into the first available slot, allocating more slots if none
 * are currently available.
 * \param[in] (x,y) the entry values. */
void seg_collect::add(double x,double y) {
    if(n==mem) add_memory();
    s[2*n]=x;
    s[2*n++]=y;
}

/** Doubles the number of available slots. */
void seg_collect::add_memory() {

    // Calculate the new number of slots, and ensure that it doesn't exceed the
    // safe limit
    int nmem=mem*2;
    if(nmem>max_mem) {
        fputs("Maximum memory allocation exceeded\n",stderr);
        exit(1);
    }

    // Allocate the new array and copy the contents of the old one into it
    double *ss=new double[2*nmem];
    memcpy(ss,s,2*mem*sizeof(double));
    
    // Delete the old one and update the class constants
    delete [] s;
    s=ss;
    mem=nmem;
    printf("# Memory scaled up to %d\n",nmem);
}

