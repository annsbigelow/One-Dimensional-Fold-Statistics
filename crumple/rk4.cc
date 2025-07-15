#include "rk4.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

/** The class constructor initializes the class constants and dynamically
 * allocates memory for the solution and intermediate steps.
 * \param[in] dof_ the number of degrees of freedom. */
rk4::rk4(int dof_) : dof(dof_), fcount(0), num_acc(0), num_tot(0),
    t(0.), q(new double[dof]), dq(new double[dof]), k1(new double[dof]),
    k2(new double[dof]), k3(new double[dof]), k4(new double[dof]) {}

/** This alternative class constructor is used for cases where the number of
 * degrees of freedom is not known ahead of time. In this case, no memory is
 * allocated. */
rk4::rk4() : dof(0), fcount(0), num_acc(0), num_tot(0), t(0.) {}

/** Allocates memory for the solution and intermediate steps. */
void rk4::allocate(double *q_) {
    q=q_;
    dq=new double[dof];
    k1=new double[dof];
    k2=new double[dof];
    k3=new double[dof];
    k4=new double[dof];
}

/** The class constructor frees the dynamically allocated memory. */
rk4::~rk4() {
    delete [] k4;
    delete [] k3;
    delete [] k2;
    delete [] k1;
    delete [] dq;
}

/** Solves the ODE problem with a fixed integration step.
 * \param[in] duration the time to solve for.
 * \param[in] steps the number of integration steps
 * \param[in] output whether to print each integration step. */
void rk4::solve_fixed(double duration,int steps,bool output) {

    // Set up initial condition and compute timestep
    double dt=duration/steps;

    // Perform integration steps
    if(output) print_step();
    for(int i=0;i<steps;i++) {
        fixed_step(dt);
        if(output) print_step();
    }
}

/** Performs an integration step with the fourth-order Runge-Kutta solver.
 * \param[in] dt the integration step. */
void rk4::fixed_step(double dt) {

    // First RK step
    ff(t,q,k1);

    // Second RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+0.5*dt*k1[i];
    ff(t+0.5*dt,dq,k2);

    // Third RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+0.5*dt*k2[i];
    ff(t+0.5*dt,dq,k3);

    // Fourth RK step
    t+=dt;fcount+=4;
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*k3[i];
    ff(t,dq,k4);

    // Complete solution
    for(int i=0;i<dof;i++) q[i]+=dt*(1/6.)*(k1[i]+2*(k2[i]+k3[i])+k4[i]);
}

/** Prints the current state of the solution. */
void rk4::print_step() {
    for(int i=0;i<dof;i++) printf(" %g",q[i]);
    puts("");
}

/** Prints the current state of the solution. */
void rk4::print_dense(int fr,double t_,double *in) {
    for(int i=0;i<dof;i++) printf(" %g",in[i]);
    puts("");
}

/** Performs a time integration of the ODE problem, using adaptive step size
 * selection.
 * \param[in] duration the end point of the integration.
 * \param[in] atol the absolute tolerance.
 * \param[in] rtol the relative tolerance.
 * \param[in] output whether to output the integration steps.
 * \param[in] d_steps the number of dense output intervals, set to zero if
 *                    dense output is not required. */
void rk4::solve_adaptive(double duration,double atol,double rtol,bool output,int d_steps) {
    double t_start=t,t_final=t+duration;

    // Parameters for dense output
    do_count=0;int new_do;
    double sf=d_steps>0?duration/d_steps:0,isf=d_steps>0?1.0/sf:0,t_den;

    // Set up initial condition and compute timestep
    double dt=initial_step_size(atol,rtol),err;
    bool last=false;
    if(t+dt>t_final) {last=true;dt=t_final-t;}

    // Do any required output of the initial step
    if(output) print_step();
    if(d_steps>0) {
        print_dense(0,t,q);
    }

    num_acc=0,num_tot=0;
    while(true) {

        num_tot++;
        // Perform integration step and check if the scaled error term is
        // acceptable
        if((err=step_and_error(dt,atol,rtol))<1.0) {
            t+=dt;
            num_acc++;

            // Do any dense output interpolation
            if(d_steps>0) {
                new_do=int((t-t_start)*isf);
                if(new_do>=d_steps) new_do=d_steps-1;
                while(do_count<new_do) {
                    do_count++;
                    t_den=t_start+do_count*sf;
                    dense_output(1.+(t_den-t)/dt,dt);
                    print_dense(do_count,t_start+do_count*sf,k3);

                    num_acc=num_tot=0;
                }
            }

            // Copy the solution and first RK step into the correct arrays
            memcpy(k1,k4,dof*sizeof(double));
            memcpy(q,dq,dof*sizeof(double));

            // Print solution and check for the termination condition
            if(output) print_step();
            if(last) {
                if(d_steps>0) {
                    do_count++;print_dense(do_count,t_final,q);
                    num_acc=num_tot=0;
                }
                return;
            }
        }

        // Compute the new timestep. If it exceeds the end of the integration
        // interval, then truncate the timestep and mark this as the last step.
        dt*=min(facmax,max(facmin,safe_fac*pow(err,-1./4.)));
        if(t+dt>t_final) {last=true;dt=t_final-t;} else last=false;
    }
}

/** Estimates an initial step size choice.
 * \param[in] atol the absolute tolerance.
 * \param[in] rtol the relative tolerance. */
double rk4::initial_step_size(double atol,double rtol) {

    // Compute d_0 and d_1
    ff(t,q,k1);
    double d0=scaled_norm(q,atol,rtol),
           d1=scaled_norm(k1,atol,rtol),d2=0,h0,h1,md12,o;

    // Make initial step size guess
    h0=fabs(d0)<1e-5||fabs(d1)<1e-5?1e-6:0.01*(d0/d1);

    // Perform one explicit step with initial step size
    for(int i=0;i<dof;i++) k2[i]=q[i]+h0*k1[i];
    ff(t+h0,k2,k3);fcount+=2;

    // Estimate second derivative of the solution
    for(int i=0;i<dof;i++) {
        o=(k1[i]-k3[i])/(atol+rtol*max(fabs(q[i]),fabs(k2[i])));
        d2+=o*o;
    }
    d2=sqrt(d2/dof)/h0;

    // Compute second step size estimate
    md12=max(d1,d2);
    h1=md12<=1e-15?max(1e-6,h0*1e-3)
                  :pow(0.01/md12,1./5.);

    // Choose the minimum of the two possibilities
    return min(100*h0,h1);
}

/** Computes the scaled norm used in initial step size selection.
 * \param[in] atol the absolute tolerance.
 * \param[in] rtol the relative tolerance. */
double rk4::scaled_norm(double *in,double atol,double rtol) {
    double no=0.,o;
    for(int i=0;i<dof;i++) {
        o=in[i]/(atol+rtol*fabs(q[i]));
        no+=o*o;
    }
    return sqrt(no/dof);
}

/** Computes a Hermite interpolation of the solution, for dense output. The
 * result is stored into the k3 array.
 * \param[in] th the fraction of the timestep at which to evaluate the Hermite
 *               interpolant.
 * \param[in] dt the length of the current timestep. */
void rk4::dense_output(double th,double dt) {
    double mth=1-th;

    // The function assumes that the current solution is in q, the new solution
    // is in dq, the current derivative is in k1, and the new derivative is in
    // k4
    for(int i=0;i<dof;i++)
        k3[i]=mth*q[i]+th*dq[i]
             -th*mth*((1-2*th)*(dq[i]-q[i])+dt*(th*k4[i]-mth*k1[i]));
}

/** Performs a Runge-Kutta step. It assumes k1 is already available.
 * \param[in] dt the step size.
 * \param[in] atol the absolute tolerance.
 * \param[in] rtol the relative tolerance.
 * \return A scaled estimate of error over the current integration step. */
double rk4::step_and_error(double dt,double atol,double rtol) {

    // Second RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+(1/3.)*dt*k1[i];
    ff(t+(1/3.)*dt,dq,k2);

    // Third RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*(-(1/3.)*k1[i]+k2[i]);
    ff(t+(2/3.)*dt,dq,k3);

    // Fourth RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*(k1[i]-k2[i]+k3[i]);
    ff(t+dt,dq,k4);

    // Complete solution that is fourth-order accurate
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*0.125*(k1[i]+3*(k2[i]+k3[i])+k4[i]);

    // Perform FSAL step needed for error estimation, reusing k4 since it is no
    // longer needed
    ff(t+dt,dq,k4);fcount+=4;

    // Compute normalized error estimate
    double err=0.,qhat,o;
    for(int i=0;i<dof;i++) {

        // Compute third-order solution for error
        qhat=q[i]+dt*((1/12.)*k1[i]+0.5*k2[i]+0.25*k3[i]+(1/6.)*k4[i]);

        // Compute scaled error contribution
        o=(dq[i]-qhat)/(atol+rtol*max(fabs(dq[i]),fabs(qhat)));
        err+=o*o;
    }
    return sqrt(err/dof);
}
