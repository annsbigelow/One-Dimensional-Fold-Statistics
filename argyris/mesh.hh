#ifndef MESH_HH
#define MESH_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../crumple/common.hh"
#include "../crumple/rk4.hh"
#include "../crumple/mesh_param.hh"
#include "../crumple/ext_potential.hh"
#include "../crumple/p_grid.hh"
#include "pre_comps.hh"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>

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
        /** The total number of triangles. */
        int ntri;
        /** The total number of "hinge" edges. */
        int nh;
        /** The total number of external potentials. */
        int n_ep;
        /** Pointers to external potentials. */
        ext_potential* ex_pot[max_ep];
        /** The reference mesh. */
        double *xyz;
        /** The node positions. */
        double *pts;
		/** Whether to use random shrink spring constants or not. */
		bool rand_sh;
		/** Whether to use random bend spring constants or not. */
		bool rand_b;
		/** Whether to use random stretch spring constants or not. */
		bool rand_st;
		/** The depth of the sheet's edges to ignore. */
		int R;
		/** The scale of random features on the shrinking sheet. */
		double set_scale;
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
        /** Pointers to the triangular elements of the mesh. */
        int **to;
        /** Memory for the triangular elements of the mesh. */
        int *tom;
		/** Argyris degrees of freedom */
		int Adof;
		/** Half the Argyris degrees of freedom */
		int Adof2;
		/** FEM Sparse mass matrix */
		Eigen::SparseMatrix<double, Eigen::RowMajor> M_sp;
		/** FEM Mass matrix Sparse Cholesky solver */
		Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > Msolver;
		/** Argyris stiffness matrix */
		Eigen::SparseMatrix<double, Eigen::RowMajor> Kd;
		/** FEM Biharmonic term */
		Eigen::VectorXd Kq;
		/** Argyris change of bases matrix */
		double *C_glob;
		/** Global normal vectors to edges */
		double *normals;
		/** Whether to use quadrature or exact integration for the FEM. */
		bool quadrature;
		/** Gauss-Legendre quadrature points */
		double *xi;
		/** Gauss-Legendre quadrature weights */
		double *w; 

        mesh(mesh_param &mp,const char* filename);
        mesh(mesh_param &mp,const char* f_topo,const char* f_pts);
        virtual ~mesh();
		void init_shrink(bool shflag,bool bendflag,bool stflag,double shm,double shv,double bm,double bv,
			double ksm,double ksv,int nx,int ny);
        void mesh_ff(double t_,double *in,double *out);
        void mesh_init() {};
        void mesh_print_dense(int fr,double t_,double *in);
		void mesh_print_last_step();
        void centralize(double &wx,double &wy,double &wz);
        void read_topology(FILE *fp);
        void read_positions(FILE *fp);
        void setup_springs();
        void reset_relaxed();
        //void contact_forces(double *in,double *out);
        void print_triangle_table();
        double energy(double t_,double *in);
        int bandwidth();
        void add(ext_potential *ep);
		void gen_spring_params_rec(double* out,double mu,double sig,int nx,int ny);
		/** Roughness (deformation) functions */
		void Sq_Sa(double& Sq,double& Sa,int nx,int ny);
		void Sdq_Sdr(double &Sdq,double &Sdr,int nx,int ny,double sed);
		double tot_area_rec(int nx,int ny);
		int find_pos_rec(int &i, int &j,int nx);
		bool inside(int i,int j,int nt,int ny,int sub);
		/** Small FEM helper functions */
		void print_pts(double* pt_array);
		void arr_zeros(double* A, int size);
		void Kq_multiply(double* in);
		void local_Kq_multiply(int tri, int Ti);
		void get_argv(int* argv, int v[3], int ed[3]);
		void tri_geo(int v[3], double* vb, double* l, double* na);
		/** Initial displacement functions */
		void Gauss_displacement();
		void linear_gradient();
		void const_displacement();
		/** FEM Matrix Assembly functions */
		void setup_fem_matrices();
		void buildC();
		void assemble_K();
		void assemble_M();
		void global_normals();
		
		/** Quadrature helper functions */
		void setup_quad_matrices();
		/** Mass matrix assembly (with quadrature) functions */
		void assembleM_quad();
			void monomials(double x, double y, double z[21]);
			void arg_z(double x, double y, double phi[21]);
			double mass_integral(double GL[21][6][6], int tri, int I, int J);
		/** Stiffness matrix assembly (with quadrature) functions */
		void assembleK_quad();
			void ders(double x, double y, double mxx[21],double mxy[21],double myy[21]);
			void arg_ders(double x, double y, double xx[21], double xy[21], double yy[21]);
			double stiff_integral(double xxGL[21][6][6], double xyGL[21][6][6], double yyGL[21][6][6],
							double The_inv[3][3], int tri,
							int I, int J);
			double laplace(double xxGL[21][6][6],double xyGL[21][6][6],double yyGL[21][6][6],
							int j,int i, 
							double The_inv[3][3],int tri,
							int I,int a);

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
            for(int i=0;i<Adof2;i++) out[i]=0.;
            Kq_multiply(in);
            for(int i=0;i<Adof2;i++) out[i]=-out[i];
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
		void shrink_force(double fac,double* in, double* acc, int i);
        //void damp_force(double *in,double *acc,int i,int k);
        //void repulsive_force(double *in,double *acc,int i,int k);
        int edge_lookup(int i,int j);
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
            allocate(Adof,pts);
        }
        mesh_rk4(mesh_param &mp,const char* f_topo,const char* f_pts) : mesh(mp,f_topo,f_pts), rk4() {
            allocate(Adof,pts);
        }
        virtual void ff(double t_,double *in,double *out) {mesh_ff(t_,in,out);}
        virtual void init(double *q) {mesh_init();}
        virtual void print_dense(int fr,double t_,double *in) {mesh_print_dense(fr,t_,in);}
		virtual void print_step() {mesh_print_last_step();}
};

#endif
