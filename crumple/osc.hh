#ifndef OSC_HH
#define OSC_HH

#include "rk4.hh"

#include <cstdio>
#include <cmath>

/** This class has functions to specify the test oscillator problem. */
class osc_rk4 : public rk4 {
    public:
        osc_rk4() : rk4(2) {}
        /** Evaluates the function f(x,y) on the RHS of the ODE.
         * \param[in] t the dependent x variable in Hairer et al.'s notation.
         * \param[in] in the array containing y variable.
         * \param[in] out the function f(x,y). */
        virtual void ff(double t_,double *in,double *out) {
            *out=-in[1]*t_;
            out[1]=*in*t_;
        }
        /** Sets up the initial conditions for the ODE. */
        virtual void init() {
            *q=1;
            q[1]=0;
        }
        virtual void print_step() {
            printf("%g %g %g %g %g\n",t,*q,q[1],*q-sol0(t),q[1]-sol1(t));
        }
        virtual void print_dense(double t_,double *in) {
            printf("%g %g %g %g %g\n",t_,*in,in[1],*in-sol0(t_),in[1]-sol1(t_));
        }
        inline double sol0(double t) {return cos(0.5*t*t);}
        inline double sol1(double t) {return sin(0.5*t*t);}
};

#endif
