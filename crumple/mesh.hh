#ifndef MESH_HH
#define MESH_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "common.hh"
#include "rk4.hh"
#include "mesh_param.hh"
#include "ext_potential.hh"
#include "p_grid.hh"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

double mesh_f(const gsl_vector *v,void *params);
void mesh_df(const gsl_vector *v,void *params,gsl_vector *df);
void mesh_fdf(const gsl_vector *x,void *params,double *f,gsl_vector *df);

// Bit flag to signify a mesh point without a full loop of triangles
const unsigned int bflag=1<<31;

// Maximum number of allowed external potentials
const int max_ep=4;

class mesh : public mesh_param {
    public:
        /** The total number of nodes. */
        int n;
        /** The total number of connections. */
        int nc;
        /** The total number of springs. */
        int ns;
        /** The total number of "hinge" edges. */
        int nh;
        /** The total number of external potentials. */
        int n_ep;
        /** Pointers to external potentials. */
        ext_potential* ex_pot[max_ep];
        /** The node positions. */
        double *pts;
		/** Whether to use random shrink spring constants or not. */
		bool rand_sh;
		/** Whether to use random bend spring constants or not. */
		bool rand_b;
		/** Whether to use random stretch spring constants or not. */
		bool rand_st;
		/** The contracting node positions. */
		double *sh_pts;
		/** The depth of the sheet's edges to ignore. */
		int R;
		/** The scale of random features on the shrinking sheet. */
		double set_scale;
		/** The shrink strengths for each node. */
		double *shs;
		/** The bend spring constants for each node. */
		double *kappas;
		/** The stretch spring constants for each node. */
		double *kss;
		/** GSL RNG */
		gsl_rng* rng;
        /** The node velocities. */
        double *vel;
        /** The number of connections by node. */
        unsigned int *ncn;
        /** Pointers to the general connection information. */
        int **ed;
        /** Memory for the general connection information. */
        int *edm;
        /** Pointers to the one-sided connection information. */
        int **eo;
        /** Memory for the one-sided connection information. */
        int *eom;
        /** Relaxed edge lengths. */
        double *reg;
        /** Relaxed edge factors for bending rigidity. */
        double *ref;
        /** Pointers to the one-sided triangle information. */
        int **to;
        /** Memory for one-sided triangle information. */
        int *tom;
        mesh(mesh_param &mp,const char* filename);
        mesh(mesh_param &mp,const char* f_topo,const char* f_pts);
        virtual ~mesh();
		void init_shrink(bool shflag, bool bendflag, bool stflag,int nx,int ny);
        void mesh_ff(double t_,double *in,double *out);
        void mesh_init() {};
        void mesh_print_dense(int fr,double t_,double *in);
        void centralize(double &wx,double &wy,double &wz);
        void read_topology(FILE *fp);
        void read_positions(FILE *fp);
        void setup_springs();
        void perturb_springs(double min_fac,double max_fac);
        void reset_relaxed();
        void contact_forces(double *in,double *out);
        void acceleration(double t_,double *in,double *acc);
        void accel_springs(double *in,double *acc);
        void accel_bsheet(double *in,double *acc);
        void accel_rbsheet(double t_,double *in,double *acc);
        double energy(double t_,double *in);
        double energy_springs(double *in);
        double energy_bsheet(double *in);
        int bandwidth();
        void add(ext_potential *ep);
		double sdev(int nx,int ny);
		double tot_area(int nx,int ny);
		int find_pos(int &i, int &j,int nx);
		bool inside(int i, int j, int nt, int ny);
		void gen_spring_params(double* out,double mu,double sig,int nx,int ny);
        //void accel_repulsive(double *in,double *acc);
        void check_deriv(double t_);
        inline void draw_nodes(const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_nodes(fp);
            fclose(fp);
        }
        void draw_nodes(FILE *fp=stdout);
        inline void draw_nodes_pov(const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_nodes_pov(fp);
            fclose(fp);
        }
        void draw_nodes_pov(FILE *fp=stdout);
        inline void draw_mesh_gnuplot(const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_mesh_gnuplot(fp);
            fclose(fp);
        }
        void draw_mesh_gnuplot(FILE *fp=stdout);
        inline void draw_mesh_pov(const char *filename,bool normals=true) {
            FILE *fp=safe_fopen(filename,"w");
            draw_mesh_pov(fp);
            fclose(fp);
        }
        void draw_mesh_pov(FILE *fp=stdout,bool normals=true);
        inline void draw_cylinders_pov(const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_cylinders_pov(fp);
            fclose(fp);
        }
        void draw_cylinders_pov(FILE *fp=stdout);
        inline void output_topology(const char *filename) {
            FILE *fp=safe_fopen(filename,"wb");
            output_topology(fp);
            fclose(fp);
        }
        void output_topology(FILE *fp);
        inline void edge_diagnostic(const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            output_positions(fp);
            fclose(fp);
        }
        void edge_diagnostic(FILE *fp=stdout);
        inline void output_positions(const char *filename) {
            FILE *fp=safe_fopen(filename,"wb");
            output_positions(fp);
            fclose(fp);
        }
        void output_positions(FILE *fp);
        void output_positions(int l);
        inline void output_velocities(const char *filename) {
            FILE *fp=safe_fopen(filename,"wb");
            output_velocities(fp);
            fclose(fp);
        }
        void output_velocities(FILE *fp);
        void output_velocities(int l);
        void edge_diagnostic();
        void setup_output_dir(const char* odir);
        void minimize_energy(double t_);
        inline double fun(double *in) {
            return energy(fzt,in);
        }
        inline void dfun(double *in,double *out) {
            for(int i=0;i<3*n;i++) out[i]=0.;
            acceleration(fzt,in,out);
            for(int i=0;i<3*n;i++) out[i]=-out[i];
        }
    private:
        /** The output directory filename, if present. */
        char *odir;
        /** The buffer for assembling output filenames, if present. */
        char *obuf;
        /** A reference to the proximity grid data structure. */
        proximity_grid pg;
        /** "Frozen" time point, needed when halting simulation to perform
        derivative checks or energy minimization. */
        double fzt;
        void normal_vector(int i,double &x,double &y,double &z);
        void n_vec_contrib(double &x,double &y,double &z,int i,int j,int k);
        int edge_mark(int i,int *edp);
        void cyl_print(FILE *fp,double *p,double *q);
        void edge_force(double *in,double *acc,int i,int k);
        void stretch_force(double *in,double *acc,int i,int k,double sf);
		void shrink_force(double fac,double* in, double* acc, int i);
        void damp_force(double *in,double *acc,int i,int k);
        double edge_energy(double *in,int i,int k);
        void triangle_force(double *in,double *acc,int i,int j,int k,int l);
        double triangle_energy(double *in,int i,int j,int k,int l);
        double edge_factor(double *in,int i,int j,int k,int l);
        void bend_force(double *in,double *acc,int i,int j,int k,int l,double ef);
        void repulsive_force(double *in,double *acc,int i,int k);
        void min_message(int &piter,int iter,gsl_multimin_fdfminimizer *s,double &t0);
        inline bool find_unmarked(int j,int *&p) {
            for(p=ed[j];p<ed[j+1];p++) if((*p&bflag)==0) return true;
            return false;
        }
        inline int ind(int i,int j,int ldab,int kd) {
            return((i+kd-j)+j*ldab);
        }
        void spatial_sort(int *s,int N,int offset);
        void width(int *s,int N,double &wx,double &wy,double &wz);
};

/** Class to set up a mesh with the RK4 integrator. */
class mesh_rk4 : public mesh, public rk4 {
    public:
        mesh_rk4(mesh_param &mp,const char* filename) : mesh(mp,filename), rk4() {
            allocate(6*n,pts);
        }
        mesh_rk4(mesh_param &mp,const char* f_topo,const char* f_pts) : mesh(mp,f_topo,f_pts), rk4() {
            allocate(6*n,pts);
        }
        virtual void ff(double t_,double *in,double *out) {mesh_ff(t_,in,out);}
        virtual void init(double *q) {mesh_init();}
        virtual void print_dense(int fr,double t_,double *in) {mesh_print_dense(fr,t_,in);}
};

#endif
