#ifndef EXT_POTENTIAL_HH
#define EXT_POTENTIAL_HH

class ext_potential {
    public:
        ext_potential() {}
        virtual ~ext_potential() {}
        virtual void accel(double t,int n,double *in,double *acc) = 0;
        virtual double energy(double t,int n,double *in) = 0;
};

class ep_centering : public ext_potential {
    public:
        /** Strength of the centering force. */
        const double cforce;
        ep_centering(double cforce_) : cforce(cforce_) {}
        virtual ~ep_centering() {}
        virtual void accel(double t,int n,double *in,double *acc);
        virtual double energy(double t,int n,double *in);
};

class ep_quadratic : public ext_potential {
    public:
        /** Rate at which the quadratic potential increases. */
        const double lambda;
        ep_quadratic(double lambda_) : lambda(lambda_) {}
        virtual ~ep_quadratic() {}
        virtual void accel(double t,int n,double *in,double *acc);
        virtual double energy(double t,int n,double *in);
};

class ep_spherical : public ext_potential {
    public:
        /** The initial cutoff radius for the radial potential. */
        double r_start;
        /** The final cutoff radius for the radial potential. */
        double r_end;
        /** The ramping time to reach the final cutoff radius. */
        double ramp;
        /** The strength of the potential. */
        double C;
        /** The class constructor for the case when the parameters are not
         * initialized ahead of time. */
        ep_spherical() {}
        /** The class constructor when parameters can be initialized right
         * away. */
        ep_spherical(double r_start_,double r_end_,double ramp_,double C_)
            : r_start(r_start_), r_end(r_end_), ramp(ramp_), C(C_),
            fac((r_end-r_start)/ramp) {}
        virtual ~ep_spherical() {}
        virtual void accel(double t,int n,double *in,double *acc);
        virtual double energy(double t,int n,double *in);
        void set_params(double r_start_,double r_end_,double ramp_,double C_);
    private:
        /** A multiplicative factor used in computing the current cutoff radius. */
        double fac;
};

class ep_piston : public ext_potential {
    public:
        /** The radius of the piston. */
        double r_piston;
        /** The radius of the piston squared. */
        double r_piston2;
        /** The z coordinate of the piston top. */
        double z_top;
        /** The z coordinate of the piston bottom (symmetric about xy plane). */
        double z_bot;
        /** The small contact region for compression. */
        double delta;
        /** The frequency of compression. */
        double freq;
        /** The class constructor for the case when the piston parameters
         * are not known ahead of time. */
        ep_piston() {}
        /** The class constructor when the piston parameters are known. */
        ep_piston(double r_piston_,double z_piston,double delta_,double freq_) :
            r_piston(r_piston_), r_piston2(r_piston*r_piston), z_top(z_piston),
            z_bot(-z_piston), delta(delta_), freq(freq_) {}
        virtual ~ep_piston() {}
        virtual void accel(double t,int n,double *in,double *acc);
        virtual double energy(double t,int n,double *in);
        void set_params(double r_piston_,double z_piston_,double delta_,double freq_);
};

class ep_shell : public ext_potential {
    public:
        /** Constant velocity at which shell radius decreases. */
        const double gamma;
        /** Strength of potential. */
        const double C;
        /** Initial radius of the shell. */
        double r0;
        /** Initial time. */
        double t0;
        ep_shell(double gamma_,double C_) : gamma(gamma_), C(C_) {}
        virtual ~ep_shell() {}
        virtual void accel(double t,int n,double *in,double *acc);
        virtual double energy(double t,int n,double *in);
        void set_params(double r_,double t_);
};

#endif
