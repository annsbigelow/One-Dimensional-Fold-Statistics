#include "osc.hh"
#include <cstdio>

int main() {

    // Solve the oscillator problem using the sixth-order
    // Richardson-extrapolated Cash-Karp method
    osc_rk4 bh;
    //brus_ckr bh;
    double t=8,lambda=1e-3;
    int d_steps=600;

    // Output the integration time points
    bh.init();
    bh.solve_adaptive(t,lambda,lambda,true);
    puts("\n");

    // Do dense output
    bh.t=0;
    bh.init();
    bh.solve_adaptive(t,lambda,lambda,false,d_steps);
}
