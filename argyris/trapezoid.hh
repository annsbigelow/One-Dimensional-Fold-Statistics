#ifndef TRAPEZOID_HH
#define TRAPEZOID_HH

class trapezoid {
    public:
        /** The total number of degrees of freedom in the ODE system. */
        int dof;
        /** The current time. */
        double t;
        /** The solution vector. */
        double *q;
        double *dq;
        /** Whether the underlying problem is linear or not. */
        bool linear;
        /** A counter for the number of function evaluations. */
        int fcount;
        
        trapezoid(int dof_,bool linear_);
        trapezoid();
        virtual ~trapezoid();
        void allocate(double *q_);
        /** Allocates memory for the solution. */
        inline void allocate(int dof_,double *q_) {
            dof=dof_;
            allocate(q_);
        }
        void solve_fixed_linear(double duration,int iters,bool output=false);
        virtual void mat_assemble(double dt) = 0;
        void fixed_step_linear(double dt);
        void fixed_Newton_step(double dt);
        virtual void print_step();
        virtual void FF(double dt,double t_,double *in,double *out) = 0;
};

#endif