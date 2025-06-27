#ifndef RK4_HH
#define RK4_HH

/** The minimum amount that the timestep can be reduced by in one step. */
const double facmin=1/3.;

/** The maximum amount that the timestep can be enlarged by in one step. */
const double facmax=3.;

/** A safety factor to scale the optimal timestep choice by, so that the next
 * step will be accepted with high probability. */
const double safe_fac=0.9;

/** Class for solving an ODE IVP using the fourth-order Runge-Kutta method. */
class rk4 {
    public:
        /** The total number of degrees of freedom in the ODE system. */
        int dof;
        /** A counter for the number of function evaluations. */
        int fcount;
        /** A counter for the dense outputs. */
        int do_count;
        /** The counter of accepted timesteps between dense outputs */
        int num_acc;
        /** The counter of total timesteps between dense outputs. */
        int num_tot;
        /** The current time. */
        double t;
        /** The solution vector. */
        double *q;
        double *dq;
        double *k1;
        double *k2;
        double *k3;
        double *k4;
        rk4(int dof_);
        rk4();
        virtual ~rk4();
        void allocate();
        /** Allocates memory for the solution and intermediate steps. */
        inline void allocate(int dof_) {
            dof=dof_;
            allocate();
        }
        void solve_fixed(double duration,int iters,bool output=false);
        void solve_adaptive(double duration,double atol,double rtol,bool output=false,int d_steps=0);
        double initial_step_size(double atol,double rtol);
        void dense_output(double theta,double dt);
        void fixed_step(double dt);
        double step_and_error(double dt,double atol,double rtol);
        virtual void print_step();
        virtual void print_dense(double t_,double *in);
        virtual void ff(double t_,double *in,double *out) = 0;
    private:
        double scaled_norm(double *in,double atol,double rtol);
        inline double min(double a,double b) {return a<b?a:b;}
        inline double max(double a,double b) {return a>b?a:b;}
};

#endif
