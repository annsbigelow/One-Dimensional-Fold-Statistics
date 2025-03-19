# Compiler and compilation flags
cc=gcc
cxx=g++
cflags=-Wall -ansi -pedantic -march=native -O3 -fopenmp

# MPI compiler
mpicxx=mpicxx -Wno-long-long

# Flags for linking to PNG library
png_iflags=-DHAS_PNG
png_lflags=-lpng

# Flags for FFTW library
fftw_iflags=
fftw_lflags=-lfftw3

# Flags for GSL
gsl_iflags=`gsl-config --cflags`
gsl_lflags=`gsl-config --libs`

# Flags for Voro++
voropp_iflags=-I/usr/local/include/voro++
voropp_lflags=-L/usr/local/lib -lvoro++

# LAPACK flags for dense linear algebra
lp_lflags=-llapack -lblas
