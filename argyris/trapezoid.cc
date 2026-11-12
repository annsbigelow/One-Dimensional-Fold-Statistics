#include "trapezoid.hh"

#include <cstdio>
#include <cmath>

/** The class constructor initializes the class constants and dynamically
 * allocates memory for the solution.
 * \param[in] dof_ the number of degrees of freedom. */
trapezoid::trapezoid(int dof_,bool linear_) : dof(dof_), linear(linear_), fcount(0),
    i_alloc(true), t(0.), dq(new double[dof]) {}

trapezoid::trapezoid() : dof(0), linear(true), fcount(0), i_alloc(false), t(0.) {}

/** Allocates memory for the solution. */
void trapezoid::allocate(double *q_) {
    q=q_;
    dq=new double[dof];
}

/** The class constructor frees the dynamically allocated memory. */
trapezoid::~trapezoid() {
    delete [] dq;
    if(i_alloc) delete [] q;
}

/** Solves the ODE problem with a fixed integration step.
 * \param[in] duration the time to solve for.
 * \param[in] steps the number of integration steps
 * \param[in] output whether to print each integration step. */
void trapezoid::solve_fixed_linear(double duration,int steps,bool output) {

    // Set up initial condition and compute timestep
    double dt=duration/steps;

    // Assemble the matrix involved in the solve
    mat_assemble(dt);

    // Perform integration steps
    if(output) print_step();
    for(int i=0;i<steps;i++) {
        if(linear) fixed_step_linear(dt);
        else fixed_Newton_step(dt);
        if(output) print_step();
        printf("Outputted solution at time %g.\n",t);
    }
}

/** Performs an integration step via the trapezoid method for a linear problem.
 * \param[in] dt the integration step. */
void trapezoid::fixed_step_linear(double dt) {
    // Solve for the solution step
    fcount+=1;
    t+=dt;
    FF(dt,q,dq);
    for(int i=0;i<dof;i++) q[i]+=dq[i];
}

/** Prints the current state of the solution. */
void trapezoid::print_step() {
    for(int i=0;i<dof;i++) printf(" %g",q[i]);
    puts("");
}

void trapezoid::fixed_Newton_step(double dt) {
    // TODO : implement the Newton step.
    printf("Hello!\n");
}