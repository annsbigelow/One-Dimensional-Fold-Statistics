#ifndef BRUSSELATOR_HH
#define BRUSSELATOR_HH

#include "rk4.hh"

#include <cstdio>

/** This class has functions to specify the test Brusselator problem. */
class brus_rk4 : public rk4 {
    public:
        brus_rk4() : rk4(2) {}
        /** Evaluates the function f(x,y) on the RHS of the ODE.
         * \param[in] t the dependent x variable in Hairer et al.'s notation.
         * \param[in] in the array containing y variable.
         * \param[in] out the function f(x,y). */
        virtual void ff(double t_,double *in,double *out) {
            double &y1=*in,&y2=in[1];
            *out=1+y1*(y1*y2-4);
            out[1]=y1*(3-y1*y2);
        }
        /** Sets up the initial conditions for the ODE. */
        virtual void init() {
            *q=1.5;
            q[1]=3.;
        }
        virtual void print_step() {
            printf("%g %g %g\n",t,*q,q[1]);
        }
        virtual void print_dense(int fr,double t_,double *in) {
            printf("%g %g %g\n",t_,*in,in[1]);
        }
};

#endif
