#include <cstring>

#include "mesh.hh"
#include "vec3.hh"

#include <Eigen/SparseCholesky>
#include <gsl/gsl_randist.h>

/** The constructor reads in a mesh from a file, and sets up the vertices and
 * edge tables.
 * \param[in] mp a mesh_param structure containing simulation constants.
 * \param[in] filename the file to read from. */
mesh::mesh(mesh_param &mp,const char* filename) : mesh_param(mp),
    n_ep(0), odir(NULL) {
    FILE *fp=safe_fopen(filename,"rb");
    read_topology(fp);
    read_positions(fp);
    fclose(fp);

	// Seed the RNG
	if(shrink) { rng=gsl_rng_alloc(gsl_rng_taus2); gsl_rng_set(rng,22); }
}

/** The constructor reads in mesh topology and vertex positions from separate
 * files.
 * \param[in] mp a mesh_param structure containing simulation constants.
 * \param[in] f_topo the file to read the mesh topology from.
 * \param[in] f_pts the file to read the vertex positions from. */
mesh::mesh(mesh_param &mp,const char* f_topo,const char* f_pts) :
    mesh_param(mp), n_ep(0), odir(NULL) {

    // Read in the mesh topology
    FILE *fp=safe_fopen(f_topo,"rb");
    read_topology(fp);
    fclose(fp);

    // Read in the vertex positions
    fp=safe_fopen(f_pts,"rb");
    read_positions(fp);
    fclose(fp);

	// Seed the RNG
	if(shrink) {
		rng=gsl_rng_alloc(gsl_rng_taus2); gsl_rng_set(rng,22);
	}
}

/** The class destructor frees the dynamically allocated memory. */
mesh::~mesh() {
    if(odir!=NULL) delete [] odir;
    if(bsheet_model) { delete [] tom; delete [] to;}
     delete [] eom; delete [] eo;
    delete [] edm;delete [] ed;
    delete [] ncn;
    delete [] pts;

	if(shrink){ 
		delete[] shs; delete[] kappas; delete[] kss;
		gsl_rng_free(rng);
	}

	// FEM terms
	delete[] M_lump; delete[] P;
}

/** Sets up the spring network table and initializes the spring rest lengths to
 * be fully relaxed. Sets up the boundaries of the subsheet, if applicable.
 */
void mesh::setup_springs() {
	// Copy initial positions to define the reference configuration
	sh_pts = new double[3*n];
	std::memcpy(sh_pts,pts,3*n*sizeof(double));

    // For this case, the number of springs will just be half the number of
    // connections
    ns=nc>>1;

    // Scan the triangle table, and store an edge (i,j) if i<j to ensure that
    // only one copy of each is kept
    eo=new int*[n+1];
    eom=new int[ns];
    int i,*eop=eom,*edp=edm;
    for(i=0;i<n;i++) {
        eo[i]=eop;
        while(edp<ed[i+1]) {
            if(*edp>i) *(eop++)=*edp;
            edp++;
        }
    }
    eo[n]=eop;

    // Compute the number of triangles
    ntri=0;
    for(int i=0;i<n;i++)
        ntri+=(ncn[i]&bflag)==0?ncn[i]:ncn[i]-1;
    if(ntri%3!=0) {
        fputs("Triangle count should be divisible by 3",stderr);
        exit(1);
    }
    ntri/=3;

    // Compute the triangle table
    to=new int*[n+1];
    tom=new int[2*ntri];
    int *top=tom;edp=edm;
    for(int i=0;i<n;i++) {
        to[i]=top;

        // Enumerate triangles between successive edge pairs
        while(edp+1<ed[i+1]) {
            if(*edp>i&&edp[1]>i) {
                *(top++)=*edp;
                *(top++)=edp[1];
            }
            edp++;
        }

        // If this isn't a boundary case, then enumerate an
        // additional triangle
        if((ncn[i]&bflag)==0) {
            if(*edp>i && *ed[i]>i) {
                *(top++)=*edp;
                *(top++)=*ed[i];
            }
        }
        edp++;
    }
	
	// Initialize force vector in FEM computations
	P=new double[3*n];
	top=tom;

	// WITH MASS LUMPING
	if (lump) {
		double a=rho/6;
		M_lump=new double[3*n]; arr_zeros(M_lump,3*n);
		// Loop through generating indices
		for(int Ti=0;Ti<n;Ti++) {
			while(top<to[Ti+1]) {
				// Get coordinates from ref. domain
				int v[3]={Ti,*top,top[1]};
				double *v1=sh_pts+3*v[0], x1=*v1, y1=v1[1],
						*v2=sh_pts+3*v[1], x2=*v2, y2=v2[1],
						*v3=sh_pts+3*v[2], x3=*v3, y3=v3[1];
				// Reference mapping
				double detF=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
				double mass=detF*a;
				if (detF<1e-16) fprintf(stderr, "Reference mapping matrix is singular.\n"); // DIAGNOSTIC
				// Loop through nodes and vector components of nodes to assemble M
				for (int i=0;i<3;i++)
				for (int k=0;k<3;k++) M_lump[3*v[i]+k]+=mass;
				top+=2;
			}
		}
		// DIAGNOSTIC: Check whether M_lump is singular, and print elements. 
		/*printf("Lumped mass matrix:\n");
		for (int i=0;i<3*n;i++) {
			if (M_lump[i]<1e-16) fprintf(stderr, "Lumped mass matrix is singular.\n");
			printf("%g\n",M[i]);
		}
		printf("End mass matrix diagnostic.\n");*/
	}

	// W/O MASS LUMPING
	if (!lump) {
		// Build a list of triplets to store mass matrix in compressed form
		M_sp.resize(3*n,3*n);
		typedef Eigen::Triplet<double> T;
		std::vector<T> triplets;
		// Get (estimated) memory up front for performance
		triplets.reserve(18*n);

		int S[9]={2,1,1, 1,2,1, 1,1,2};
		double a=rho/24;
		for(int Ti=0;Ti<n;Ti++) {
			while(top<to[Ti+1]) {
				// Get coordinates from ref. domain
				int v[3]={Ti,*top,top[1]};
				double *v1=sh_pts+3*v[0], x1=*v1, y1=v1[1],
						*v2=sh_pts+3*v[1], x2=*v2, y2=v2[1],
						*v3=sh_pts+3*v[2], x3=*v3, y3=v3[1];
				// Reference mapping
				double detF=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
				double mass=detF*a;
				if (detF<1e-16) fprintf(stderr, "Reference mapping matrix is singular.\n"); // DIAGNOSTIC
				// Loop through nodes and vector components of nodes to assemble M
				for (int i=0;i<3;i++)
				for (int j=0;j<3;j++)
				for (int k=0;k<3;k++) { 
					int row=3*v[i]+k, col=3*v[j]+k;
					triplets.push_back(Eigen::Triplet<double>(row,col,mass*S[3*i+j]));
				}
				top+=2;
			}
		}
		
		// Convert triplets list to SparseMatrix.
		// Contributions in same row,col are summed automatically.
		M_sp.setFromTriplets(triplets.begin(),triplets.end());
		printf("Sparse mass matrix assembled\n");

		// DIAGNOSTIC: look at nonzero entries of mass matrix
		/*for (int i=0;i<M_sp.outerSize();++i)
		for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(M_sp, i);it;++it)
				printf("(%ld,%ld)=%g\n",it.row(),it.col(),it.value());*/

		// Factorize sparse matrix using LLT Cholesky factorization
		solver.analyzePattern(M_sp);
		solver.factorize(M_sp);
		if (solver.info()!=Eigen::Success) printf("Matrix factorization failed\n");
		printf("Finished sparse matrix factorization.\n");
	}
	
    //if(shrink) {
		set_scale=1.;
		double h=static_cast<double>(sed);
		R = static_cast<int>(std::ceil(set_scale / h));
	//}
}

void mesh::arr_zeros(double *A,int size) {
	for (int i=0;i<size;i++) A[i]=0.;
}

void mesh::print_triangle_table() {
    int *top=tom;
    for(int i=0;i<n;i++) {
        while(top<to[i+1]) {
            printf("(%d,%d,%d)\n",i,*top,top[1]);
            top+=2;
        }
    }
}

void mesh::print_pts(double *pt_array) {
	for (int i = 0; i < n; i++) {
		double* pt = pt_array + 3 * i;
		printf("(%g,%g,%g)", pt[0], pt[1], pt[2]);
		printf(" ");
	}
	printf("\n");
}

/** Copies initial node positions in the presence of a shrinking substrate and applies
*	a random perturbation to the rate of each contracting node.
*	\param[in] shflag, bendflag, stflag: Flags setting the choices to use random spring constants.
*	param[in] shm,shv,...,ksv: Mean and variance for each set of springs.
*/
void mesh::init_shrink(bool shflag,bool bendflag,bool stflag,double shm,double shv,double bm,double bv,
	double ksm,double ksv,int nx,int ny) {
	shs=new double[n]; kappas=new double[n]; kss=new double[n];
	rand_sh=shflag;rand_b=bendflag;rand_st=stflag;  

	if(rand_sh) gen_spring_params_rec(shs,shm,shv,nx,ny);
	if(rand_b) gen_spring_params_rec(kappas,bm,bv,nx,ny);
	if(rand_st) gen_spring_params_rec(kss,ksm,ksv,nx,ny);
}

void mesh::mesh_ff(double t_,double *in,double *out) {
    double *acc=out+3*n;
    int i;
	// Add forces coming from finite-element (FEM) computations
	arr_zeros(P,3*n);
	fem_forces(t_,in);

	if (lump) {
		// Air-resistance-type drag
		// (acceleration) = - (const.)*(velocity)
		for(double *ap=acc,*vp=in+3*n,*Mp=M_lump;ap<acc+3*n;) 
			*(ap++)=-*(vp++)*drag / *(Mp++);
	
		// XXX - we can ignore this for now
		// contact_forces(in,out);
		for(double *ap=acc,*Fp=P,*Mp=M_lump;ap<acc+3*n;) *(ap++) += *(Fp++) / *(Mp++);
	}
	else {
		Eigen::VectorXd f_sum(3*n), av(3*n);
		double *vp=in+3*n;
		for (int i=0;i<3*n;i++) f_sum[i] = P[i]-drag*vp[i];
		av=solver.solve(f_sum);
		if (solver.info()!=Eigen::Success) printf("Matrix solving failed\n");
		std::memcpy(acc,av.data(),3*n*sizeof(double));
	}
	
    // Assemble the velocities in the first part of the out array. In addition,
    // zero out the forces for nodes on the boundary, if required.
    if(fix_boundary) {
        for(i=0;i<n;i++) {
            if(ncn[i]&bflag) {
                out[3*i]=out[3*i+1]=out[3*i+2]=0;
                acc[3*i]=acc[3*i+1]=acc[3*i+2]=0;
            } else {
                out[3*i]=in[3*n+3*i];
                out[3*i+1]=in[3*n+3*i+1];
                out[3*i+2]=in[3*n+3*i+2];
            }
        }
    } else for(i=0;i<n;i++) {
        out[3*i]=in[3*n+3*i];
        out[3*i+1]=in[3*n+3*i+1];
        out[3*i+2]=in[3*n+3*i+2];
    }
}

/** Adds in the contact forces between nodes in the mesh.
 * \param[in] in the mesh point positions.
 * \param[in] out the mesh point accelerations (cumulative). */
void mesh::contact_forces(double *in,double *out) {
    double *acc=out+3*n,K=100;
    const double diamsq=diam*diam,
                 screen=6,screensq=screen*screen;

    // Build the proximity grid data structure
    pg.setup(in,n);
    pg.populate(in,n);

    // Loop over all of the balls in the simulation
    for(int i=0;i<n;i++) {

        // Determine the grid subregion that could possibly interact with this
        // ball
        int li,ui,lj,uj,lk,uk;
        pg.subregion(in+3*i,diam,li,ui,lj,uj,lk,uk);

        for(int ck=lk;ck<=uk;ck++) for(int cj=lj;cj<=uj;cj++) {
            for(int ci=li;ci<=ui;ci++) {
                int ijk=ci+pg.m*(cj+pg.n*ck);

                // Loop over all mesh points in this block
                point_info *pip=pg.p[ijk],*pie=pip+pg.co[ijk];
                for(;pip<pie;pip++) if(pip->id>i) {
                    double dx=pip->x-in[3*i],dy=pip->y-in[3*i+1],dz=pip->z-in[3*i+2],
                           rsq=dx*dx+dy*dy+dz*dz;

                    // If the mesh point is in contact then compute the acceleration
                    // contribution
                    if(rsq<diamsq) {

                        // Rule out nearby points on the mesh using initial positions
                        int i2=pip->id;
                        double *q=sh_pts+3*i,
                               *r=sh_pts+3*i2,
                               ex=q[0]-r[0],
                               ey=q[1]-r[1],
                               ez=q[2]-r[2];

                        // If the points are far away from each other in the
                        // mesh coordinates, then apply a contact force
                        if(ex*ex+ey*ey+ez*ez>screensq) {
                            double fac=K*(1-diam/sqrt(rsq));
                            dx*=fac;dy*=fac;dz*=fac;
                            acc[3*i]+=dx;
                            acc[3*i+1]+=dy;
                            acc[3*i+2]+=dz;
                            acc[3*i2]-=dx;
                            acc[3*i2+1]-=dy;
                            acc[3*i2+2]-=dz;
                        }
                    }
                }
            }
        }
    }
}

/** Computes the finite-elements forces from sheet mechanics.
 * \param[in] t_ the time at which to evaluate the acceleration.
 * \param[in] in the mesh point positions. */
void mesh::fem_forces(double t_,double *in) {
	// Gradients of basis functions listed as dPsi/dX, dPsi/dY
	int dPdX[6]={-1,1,0, -1,0,1};

    int *top=tom;
    for(int Ti=0;Ti<n;Ti++) {
        while(top<to[Ti+1]) {
			int v[3]={Ti,*top,top[1]};
			double *v1=sh_pts+3*v[0], x1=*v1, y1=v1[1],
					*v2=sh_pts+3*v[1], x2=*v2, y2=v2[1],
					*v3=sh_pts+3*v[2], x3=*v3, y3=v3[1];
			double F[4]={x2-x1,x3-x1,y2-y1,y3-y1};
			double detF=F[0]*F[3]-F[1]*F[2];
			if (detF < 1e-12) { //DIAGNOSTIC 
				fprintf(stderr,"Error: detF is small or negative. A triangle's vertices may be stored clockwise.\n");
				break;
			}

			// Build Pk (independent of psi_i)
			double Pk[3][2];
			double *qT[3]={in+3*v[0], in+3*v[1], in+3*v[2]};
			for (int k=0;k<3;k++) get_hatP(Pk[k],k,qT,dPdX,F,detF);

			// Assemble fem force vector
			for (int i=0;i<3;i++) // Loop through triangle vertices
			for (int k=0;k<3;k++) { // Loop through vector components
				double Ak=Pk[k][0]*FdPI(F,detF,dPdX,i,0) + Pk[k][1]*FdPI(F,detF,dPdX,i,1);
				P[3*v[i]+k]-=detF*Ak/2;
			}
            top+=2;
        }
    }
    // Add accelerations from the external potentials
    //for(int i=0;i<n_ep;i++) ex_pot[i]->accel(t_,n,in,acc);
}

/* Compute a row of the block matrix P_hat used FEM force calculations. 
*	\param[in] k the row of P to return
*/
void mesh::get_hatP(double(&hatP_k)[2],int k,double* qT[3],int dPdX[6],double F[4],double detF) {
	hatP_k[0]=0.;
	hatP_k[1]=0.;

	// grad q 
	double G[3][2] = {0.};
	for (int p=0;p<3;p++)
	for (int m=0;m<2;m++)
	for (int i=0;i<3;i++) {
		double *qTi=qT[i]; // 3 components of i-th node
		G[p][m] += qTi[p]*FdPI(F,detF,dPdX,i,m);
	}
	// trace(G^T G)
	double C=0.;
	for (int p=0;p<3;p++)
	for (int m=0;m<2;m++) C+=G[p][m]*G[p][m];
	// G^T G
	double D[2][2] = {0};
	for (int m=0;m<2;m++)
	for (int l=0;l<2;l++)
	for (int p=0;p<3;p++) D[m][l]+=G[p][m]*G[p][l];

	for (int l=0;l<2;l++)
	for (int m=0;m<2;m++) hatP_k[l] += G[k][m]*( lambda*.5*(C - 2.)*delta(m,l) + mu_fem*(D[m][l]-delta(m,l)) );
}

int mesh::delta(int m, int l) {
	return (l==m?1:0);
}

/** A component of the inverse-transpose of the reference mapping matrix
	multiplied by the gradient of a basis function. */
double mesh::FdPI(double F[4],double detF,int dPdX[6],int i,int m) {
	if (m==0) return (dPdX[i]*F[3] - dPdX[3+i]*F[2])/detF;
	else return (dPdX[3+i]*F[0] - dPdX[i]*F[1])/detF;
}

/** Computes the energy.
 * \param[in] in the mesh point positions. */
double mesh::energy(double t_,double *in) {

    // Compute energies due to edge springs
    /*double en=bsheet_model?energy_bsheet(in):energy_springs(in);

    // Add energies due to external potentials
    for(int i=0;i<n_ep;i++) en+=ex_pot[i]->energy(t_,n,in);
    return en;*/
    return 0;
}

/** Centralizes the mesh, and calculates its square width in each of the three
 * coordinate directions.
 * \param[out] (wx,wy,wz) the square widths in the three coordinate directions. */
void mesh::centralize(double &wx,double &wy,double &wz) {
    double sx=0.,sy=0.,sz=0.,fac=1./static_cast<double>(n);

    // Compute the centroid
    for(double *p=pts;p<pts+3*n;p+=3) {
        sx+=*p;sy+=p[1];sz+=p[2];
    }
    sx*=fac;sy*=fac;sz*=fac;

    // Displace the mesh so its centroid is at the origin, and compute the
    // variance of the vertices in each coordinate
    wx=wy=wz=0;
    for(double *p=pts;p<pts+3*n;p+=3) {
        *p-=sx;p[1]-=sy;p[2]-=sz;
        wx+=*p*(*p);wy+=p[1]*p[1];wz+=p[2]*p[2];
    }
    wx*=fac;wy*=fac;wz*=fac;
}

/** Calculates the acceleration due to the contraction of nodes
 * towards the centroid. */
void mesh::shrink_force(double fac,double *in,double *acc,int i) {
	// Define "springs" between shrinking points and current nodes
	double *is=sh_pts+3*i, *ip=in+3*i,
			dx=*is*fac-*ip, dy=is[1]*fac-ip[1], dz=is[2]*fac-ip[2],
			*ai=acc+3*i;

	// Add the force contributions to the vertex
	if(rand_sh) {
		dx*=shs[i]; dy*=shs[i]; dz*=shs[i];
	}
	else {
		dx*=ks;dy*=ks;dz*=ks;
	}
	*ai+=dx;ai[1]+=dy;ai[2]+=dz;
}

void mesh::damp_force(double *in,double *acc,int i,int k) {
    double *ii=in+3*(i+n),*ik=in+3*(k+n),*ai=acc+3*i,*ak=acc+3*k,
           dx=*ii-*ik,dy=ii[1]-ik[1],dz=ii[2]-ik[2];

    // Add the force contributions to the two vertices
    *ai-=dx*B;ai[1]-=dy*B;ai[2]-=dz*B;
    *ak+=dx*B;ak[1]+=dy*B;ak[2]+=dz*B;
}

void mesh::repulsive_force(double *in,double *acc,int i,int k) {
    double *ii=in+3*i,*ik=in+3*k,*ai=acc+3*i,*ak=acc+3*k,
           dx=*ii-*ik,dy=ii[1]-ik[1],dz=ii[2]-ik[2],
           r=sqrt(dx*dx+dy*dy+dz*dz);

    if(r<sigma) {
        //printf("%d %d %.4f %.4f\n",i,k,r,sigma);
        //double rs=sigma/r-1; // linear spring
        double xinv=sigma/(r-sigma),rs=exp(xinv)*xinv*xinv*(sigma/r); // differentiable test function

        // Add the force contributions to the two vertices
        dx*=rs*K;dy*=rs*K;dz*=rs*K;
        *ai+=dx;
        ai[1]+=dy;
        ai[2]+=dz;
        *ak-=dx;
        ak[1]-=dy;
        ak[2]-=dz;
    }
}

/** Adds an external potential to the class.
 * \param[in] ep a pointer to the potential to add. */
void mesh::add(ext_potential *ep) {
    if(n_ep==max_ep) {
        fputs("Too many external potentials\n",stderr);
        exit(1);
    }
    ex_pot[n_ep++]=ep;
}

/** Checks whether a node lies within a given boundary.
*	\param[in] row,col the location of the node
*	\param[in] nt,ny the x- and y-dimensions of the current row
*	\param[in] sub the layers of edge nodes to ignore
*/
bool mesh::inside(int row,int col,int nt,int ny,int sub) {
	if (col>sub-1&&col<nt-sub&&row>sub-1&&row<ny-sub) return true;
	else return false;
}

/** Fills an array with n log-normally distributed values after using a Gaussian filter.
*	Valid for a mesh with regular hexagonal topology.
*	\param[in] m the mean of the log-normal distribution.
*	\param[in] s the standard deviation of the log-normal distribution.
*	\param[in] nx,ny the dimensions of the mesh.
*	\return the array of random, filtered values.
*/
void mesh::gen_spring_params_rec(double* out,double m,double s,int nx,int ny) {
	double mm=m*m,ss=s*s;
	double mu=0.,sig=0.,R2=static_cast<double>(R*R),val;
	int i,g=0;

	// Calculate filter weights. Each weight is applied to an entire support layer.
	double* weights=new double[R+1];
	double norm=0.,sum=0.,sum2=0.,A,B;
	// Loop through neighbor layers, starting with closest hexagonal layer
	for (i=0;i<=R;i++) {
		weights[i]=exp(-(static_cast<double>(i)*i)/(2*R2));
		norm+=weights[i];
	}
	// Normalize the weights and sum the squares
	for (i=0;i<=R;i++) {
		weights[i]/=norm;
		sum+=i*weights[i];
		sum2+=i*weights[i]*weights[i];
	}
	A=weights[0]+(6*sum);
	B=(weights[0]*weights[0])+(6*sum2);
	mu=log(mm/sqrt(ss+mm))/A;
	sig=sqrt(log(ss/mm+1)/B);

	// Generate grid of normal values
	double *in=new double[n];
	for (i=0;i<n;i++) in[i]=mu+gsl_ran_gaussian_ziggurat(rng,sig);

	// Apply the Gaussian filter to interior nodes
	std::memcpy(out,in,n*sizeof(double));
	int* seen=new int[n];
	for (i=0;i<n;i++) seen[i]=-1;
	int visit=0;
	int nn=3*R*(R+1)+1,cct,nct,nt; 
	int* cpts=new int[nn]; int* npts=new int[nn]; // Current and next layers
	for (int j=0;j<ny;j++) {
		nt=nx+(j&1);
		for (i=0;i<nt;i++,g++) {
			if(inside(i,j,nt,ny,R)){
				// Loop through layers of neighbors
				cct=0; nct=0;
				// Layer 0: focus node
				visit++;
				seen[g]=visit; cpts[cct++]=g;
				val=weights[0]*in[g];	
				for (int layer=1;layer<=R;layer++) {
					nct=0;
					// Current layer
					for (int ci=0;ci<cct;ci++) {
						int v=cpts[ci];
						// Neighbors of points in current layer
						for (int* p=ed[v];p<ed[v+1];p++) {
							int u=*p;
							// Only count contributions of unvisited neighbors
							if (seen[u]!=visit) {
								seen[u]=visit; npts[nct++]=u;
								val+=weights[layer]*in[u];
							}
						}
					}
					cct=nct;
					int* tmp=cpts; cpts=npts; npts=tmp;
					// If there are no new neighbors
					if (cct==0) break;
				}
				out[g]=val;
			}
			// The edges are unfiltered but will still be Lognormal
			else out[g]=mu*A+gsl_ran_gaussian_ziggurat(rng,s*sqrt(B));
		}
	}
	for (i=0;i<n;i++) {
		out[i]=exp(out[i]);
		//printf("%g\n",out[i]);
	}
	//printf("\n");
	delete [] cpts; delete[] npts;
	delete [] seen;
	delete [] weights;
	delete [] in;
}
