#ifndef COMMON_HH
#define COMMON_HH

#include <cstdio>
#include <cstdlib>
#include <limits>

#ifdef _OPENMP
#include "omp.h"
inline double wtime() {return omp_get_wtime();}
#else
#include <ctime>
inline double wtime() {return static_cast<double>(clock())*(1./CLOCKS_PER_SEC);}
#endif

void fatal_error(const char *p,int status);
FILE* safe_fopen(const char *filename,const char* mode);
void safe_fread(void *ptr,size_t size,size_t count,FILE *fp,const char* p);

#endif
