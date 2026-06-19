#include <cstring>
#include <fstream> //DEBUG
#include <iostream>//DEBUG

#include "mesh.hh"
#include "../crumple/vec3.hh"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>
#include <gsl/gsl_randist.h>

// TODO - update output routines to print the smooth Argyris triangles

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
    delete [] xyz;
    delete [] pts;

	if(shrink){
		delete[] kappas; delete[] kss;
		gsl_rng_free(rng);
	}

	// FEM terms
	delete[] M_lump; delete[] P;
	delete[] C_inv;
}

/** Sets up the spring network table and initializes the spring rest lengths to
 * be fully relaxed. Sets up the boundaries of the subsheet, if applicable.
 */
void mesh::setup_springs() {
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
    tom=new int[5*ntri];
    int *top=tom;edp=edm;
    for(int i=0;i<n;i++) {
        to[i]=top;

        // Enumerate triangles between successive edge pairs
        while(edp+1<ed[i+1]) {
            if(*edp>i&&edp[1]>i) {
                *(top++)=*edp;
                *(top++)=edp[1];
                *(top++)=edge_lookup(i,*edp);
				*(top++)=edge_lookup(edp[1],i);
                *(top++)=edge_lookup(*edp,edp[1]);
            }
            edp++;
        }

        // If this isn't a boundary case, then enumerate an
        // additional triangle
        if((ncn[i]&bflag)==0) {
            if(*edp>i && *ed[i]>i) {
                *(top++)=*edp;
                *(top++)=*ed[i];
                *(top++)=edge_lookup(i,*edp);
				*(top++)=edge_lookup(*ed[i],i);
                *(top++)=edge_lookup(*edp,*ed[i]);
            }
        }
        edp++;
    }

	// Build change of bases matrices
	buildC();
	printf("Change of bases matrices have been built.\n");

	// Initialize force vector in FEM computations
	P=new double[Adof2];

	// Build the mass matrix
	// Use a list of triplets for fast computation
	M_sp.resize(Adof2,Adof2);
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	// Get (estimated) memory up front for performance
	triplets.reserve(6*Adof2);

	top=tom; 
	int tri=0;
	for(int Ti=0;Ti<n;Ti++) { // Loop through generating indices
		while(top<to[Ti+1]) { // Loop through triangles
			// Get coordinates from ref. domain
			int v[3]={Ti,*top,top[1]}; // Vertices 1,2,3
			int ed[3]={top[2],top[3],top[4]}; // Edges 1,2,3
			double *v1=xyz+3*v[0], x1=*v1, y1=v1[1],
					*v2=xyz+3*v[1], x2=*v2, y2=v2[1],
					*v3=xyz+3*v[2], x3=*v3, y3=v3[1];
			
			// Reference map Jacobian
			double detF=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
			double mass = detF*rho;

			int argv[21]; // Global triangle dofs indices
			int j=3,j1=9,k;
			for (k=0;k<3;k++) argv[k]=6*v[k]; // Function values
			for (int i=0;i<3;i++) {
				argv[18+i]=6*n+ed[i]; // Normal derivatives
				for (k=1;k<3;k++) {
					argv[j]=6*v[i]+k; // First derivatives
					j++;
				}
				for (k=3;k<6;k++) { // Second derivatives
					argv[j1]=6*v[i]+k;
					j1++;
				}
			}
			/*for (int i=0;i<3;i++) { // TODO - Is this right? M, F, q ordered as they were?
				// Edges dofs
				argv[18+i]=6*n+ed[i];
				// Nodal dofs
				for (int k=0;k<6;k++) argv[6*i+k]=6*v[i]+k;
			}*/

			for (int I=0;I<21;I++)
			for (int J=0;J<21;J++) {
				double sum=0.;
				for (int a=0;a<21;a++)
				for (int b=0;b<21;b++) {
					sum += C_inv[441*tri+21*a+I]*C_inv[441*tri+21*b+J]*S[21*a+b];
				}
				triplets.push_back(Eigen::Triplet<double>(argv[I],argv[J],mass*sum));
			}

			top+=5; tri+=1;
		}
	}

	// Convert triplets list to SparseMatrix.
	// Contributions in same row,col are summed automatically.
	M_sp.setFromTriplets(triplets.begin(),triplets.end());
	printf("Sparse mass matrix assembled\n");
	printf("Is the sparse mass matrix compressed? %d\n",M_sp.isCompressed());

	//TODO: Delete: BEGIN DEBUG
	printf("Nonzeros: %ld\n",M_sp.nonZeros());
	printf("Size of M: %d\n",Adof2*Adof2);
	/*Eigen::MatrixXd M_d(M_sp);
	std::ofstream outputFile("Sparse matrix test 6x6.csv");
	for (int i=0;i<M_d.rows();i++) {
		for (int j=0;j<M_d.cols();j++) {
			outputFile << M_d(i,j) << " ";
		}
		outputFile << std::endl; 
	}
	outputFile.close();*/
	// END DEBUG

	// Factorize sparse matrix using LLT Cholesky factorization
	solver.analyzePattern(M_sp);
	solver.factorize(M_sp);
	if (solver.info()!=Eigen::Success) {
		printf("Matrix factorization failed\n");
		exit(1);
	}
	printf("Finished sparse matrix factorization.\n");
}

int mesh::edge_lookup(int i,int j) {
    if(i>j) {int k=j;j=i;i=k;}
    for(int *eop=eo[i];eop<eo[i+1];eop++) {
        if(*eop==j) return int(eop-eom);
    }
    fputs("Can't find edge index\n",stderr);
    exit(1);
}

// Build change of bases matrices for each triangle
void mesh::buildC() {
	C_inv=new double[ntri*441];

	int *top=tom;
	double D[504], E[504];
	int tri=0,i,j;
	for (int Ti=0;Ti<n;Ti++) // Loop through triangles
	while (top<to[Ti+1]) {
		// Vertices 1,2,3
		int v[3]={Ti,*top,top[1]};
		double *v1=xyz+3*v[0], x1=*v1, y1=v1[1],
				*v2=xyz+3*v[1], x2=*v2, y2=v2[1],
				*v3=xyz+3*v[2], x3=*v3, y3=v3[1];
		// Sides of the triangle
		double vb[6]={	x2-x1, y2-y1,
						x3-x1, y3-y1,
						x3-x2, y3-y2	}; 
		// Affine transformation Jacobian
		double B[4]={ vb[0],vb[2], vb[1],vb[3] };
		// Hessian
		double H[9]={	B[0]*B[0], 2*B[0]*B[2], B[2]*B[2],
						B[1]*B[0], B[1]*B[2]+B[0]*B[3], B[2]*B[3],
						B[1]*B[1], 2*B[3]*B[1], B[3]*B[3]	};
		// Lengths of each side, squared 
		double l2[3]={	vb[0]*vb[0]+vb[1]*vb[1],
						vb[2]*vb[2]+vb[3]*vb[3],
						vb[4]*vb[4]+vb[5]*vb[5]	};
		// Build the change of bases for the normal derivatives
		double a[6]={ vb[2]/l2[0], vb[3]/l2[0],
					-vb[0]/l2[1], -vb[1]/l2[1],
			-(vb[0]+vb[2])/(sqrt(2)*l2[2]), -(vb[1]+vb[3])/(sqrt(2)*l2[2]) };
		// Normal vectors
		double Rv[6]={	-vb[1], vb[0],
						-vb[3], vb[2],
						-vb[5], vb[4]	};
		double f[3], g[3];
		for (i=0;i<3;i++) {
			f[i] = a[2*i]*Rv[2*i] + a[2*i+1]*Rv[2*i+1];
			g[i] = a[2*i]*vb[2*i] + a[2*i+1]*vb[2*i+1];
		}

		// Build D, the left rectangular block matrix
		arr_zeros(D,504);
		for (i=0;i<3;i++) {
			D[24*i+i]=1;
			for (int I=0;I<2;I++) for (int J=0;J<2;J++)
				D[24*(3+2*i+I)+3+2*i+J]=B[2*J+I];
			for (int I=0;I<3;I++) for (int J=0;J<3;J++)
				D[24*(9+3*i+I)+9+3*i+J]=H[3*I+J];
			D[24*(18+i) + 18+i]=f[i];
			D[24*(18+i) + 18+i + 3]=g[i];
		}

		// Build E, the right rectangular block matrix
		arr_zeros(E,504);
		for (i=0;i<18;i++) 
			E[21*i+i]=1;
		for (i=0;i<3;i++)
			E[21*(18+i)+18+i]=sqrt(l2[i]);
		// Final "T" block of E
		double T1[9]={-1,1,0, -1,0,1, 0,-1,1};
		double T2[18]={	vb[0],vb[1], vb[0],vb[1], 0,0,
						vb[2],vb[3], 0,0, vb[2],vb[3],
						0,0, vb[4],vb[5], vb[4],vb[5] };
		double wa[9]={	vb[0]*vb[0], 2*vb[0]*vb[1], vb[1]*vb[1],
						vb[2]*vb[2], 2*vb[2]*vb[3], vb[3]*vb[3],
						vb[4]*vb[4], 2*vb[4]*vb[5], vb[5]*vb[5] };
		double T3[27]={	-wa[0],-wa[1],-wa[2], wa[0],wa[1],wa[2], 0,0,0,
						-wa[3],-wa[4],-wa[5], 0,0,0, wa[3],wa[4],wa[5],
						0,0,0, -wa[6],-wa[7],-wa[8], wa[6],wa[7],wa[8] };
		for (i=0;i<3;i++) {
			for (j=0;j<3;j++) E[21*(21+i)+j] = 15*T1[3*i+j]/18;
			for (j=0;j<6;j++) E[21*(21+i)+j+3] = -7*T2[6*i+j]/16;
			for (j=0;j<9;j++) E[21*(21+i)+j+9] = T3[9*i+j]/32;
		}

		// Multiply D and E to make C
		double C[441];
		arr_zeros(C,441);
		for (i=0;i<21;i++)
		for (j=0;j<21;j++)
		for (int k=0;k<24;k++)
			C[21*i+j] += D[24*i+k]*E[21*k+j];

		// Invert C using Eigen
		Eigen::Map<Eigen::Matrix<double,21,21,Eigen::RowMajor> > C_eig(C);
		Eigen::FullPivLU<Eigen::Matrix<double,21,21,Eigen::RowMajor> > lu(C_eig);
		if (!lu.isInvertible()) {
			printf("Error: Change of bases matrix is not invertible.\n");
			exit(1);
		}
		Eigen::Matrix<double,21,21,Eigen::RowMajor> Ce_inv = C_eig.inverse();
		double C_loc_inv[441];
		std::memcpy(C_loc_inv, Ce_inv.data(), 441*sizeof(double));

		// Copy local change of basis matrix into global C_inv table
		for (int i=0;i<441;i++) C_inv[441*tri+i]=C_loc_inv[i];

		top+=5; tri+=1;
	}
}

void mesh::Gauss_displacement() {
	const double eps=0.1;
	// Initialize function values and gradients at all nodes
	for(int i=0;i<n;i++) {
		double x=xyz[3*i], y=xyz[3*i+1];
		double power = exp(-eps*(x*x+y*y));
		// Displace the z-component
		pts[6*i]+=-eps+0.02*power;
		// First derivatives
		pts[6*i+1]+=-eps*.04*x*power;
		pts[6*i+2]+=-eps*.04*y*power;
		// Second derivativs
		pts[6*i+3]+=eps*(-.04 + eps*.08*x*x)*power;
		pts[6*i+4]+=eps*eps*.08*x*y*power;
		pts[6*i+5]+=eps*(-.04 + eps*.08*y*y)*power;
	}

	// Initialize normal derivatives
	int *top=tom;
	for (int Ti=0;Ti<n;Ti++)
	while (top<to[Ti+1]) {
		int v[3]={Ti,*top,top[1]};
		double *v1=xyz+3*v[0], x1=*v1, y1=v1[1],
				*v2=xyz+3*v[1], x2=*v2, y2=v2[1],
				*v3=xyz+3*v[2], x3=*v3, y3=v3[1];
		// Sides of the triangle
		double vb[6]={	x2-x1, y2-y1,
						x3-x1, y3-y1,
						x3-x2, y3-y2	};
		// Lengths of each side
		double l[3]={	sqrt(vb[0]*vb[0]+vb[1]*vb[1]),
						sqrt(vb[2]*vb[2]+vb[3]*vb[3]),
						sqrt(vb[4]*vb[4]+vb[5]*vb[5])	};
		// Midpoints of sides
		double m[6]={	(x2+x1)/2, (y2+y1)/2,
						(x3+x1)/2, (y3+y1)/2,
						(x3+x2)/2, (y3+y2)/2 };
		// Normal vectors
		double na[6]={	-vb[1]/l[0], vb[0]/l[0],
						-vb[3]/l[1], vb[2]/l[1],
						-vb[5]/l[2], vb[4]/l[2] };

		pts[6*n+top[2]] += -eps*.04*m[0]*exp(-eps*(m[0]*m[0]+m[1]*m[1]))*na[0]
							-eps*.04*m[1]*exp(-eps*(m[0]*m[0]+m[1]*m[1]))*na[1] ;
		pts[6*n+top[3]] += -eps*.04*m[2]*exp(-eps*(m[2]*m[2]+m[3]*m[3]))*na[2]
							-eps*.04*m[3]*exp(-eps*(m[2]*m[2]+m[3]*m[3]))*na[3] ;
		pts[6*n+top[4]] += -eps*.04*m[4]*exp(-eps*(m[4]*m[4]+m[5]*m[5]))*na[4]
							-eps*.04*m[5]*exp(-eps*(m[4]*m[4]+m[5]*m[5]))*na[5] ;

		// Normal derivatives should cancel?
		for (int i=2;i<5;i++) if (pts[6*n+top[i]]<1e-16) pts[6*n+top[i]]=0.;

		top+=5;
	}
	// TODO - delete
	/*printf("Nodes\n");
	for (int i=0;i<n;i++) {
		for (int k=0;k<6;k++) printf("%g ",pts[6*i+k]);
		printf("\n");
	}
	printf("Edges\n");
	for (int i=0;i<ns;i++) {
		printf("%g ",pts[6*n+i]);
	}
	printf("\n");*/

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
	for (int k=0; k<6;k++) {
		double* pt = pt_array + 6 * i;
		printf("%g ",pt[k]);
	}
	printf("\n");
	}
}

void mesh::mesh_ff(double t_,double *in,double *out) {
    double *acc=out+Adof2;
    int i;
	// Add forces coming from finite-element (FEM) computations
	arr_zeros(P,Adof2);
	fem_forces(t_,in);

	Eigen::VectorXd f_sum(Adof2), av(Adof2);
	double *vp=in+Adof2;
	for (i=0;i<Adof2;i++) f_sum[i] = P[i]-drag*vp[i];
	// Cholesky direct solver, a = f/m
	av=solver.solve(f_sum);
	if (solver.info()!=Eigen::Success) printf("Matrix solving failed\n");
	std::memcpy(acc,av.data(),Adof2*sizeof(double));

    // Assemble the velocities in the first part of the out array. In addition,
    // zero out the forces for nodes on the boundary, if required.
    if(fix_boundary) {
        for(i=0;i<n;i++)
		for(int k=0;k<6;k++){
            if(ncn[i]&bflag) {
                out[6*i+k]=0;
                acc[6*i+k]=0;
            } else out[6*i+k]=in[Adof2+6*i+k];
        }
		// Edge dofs
		for(i=0;i<ns;i++) {
			if(ncn[i]&bflag) { // TODO - this is wrong. ncn has n elements. Loop through triangles.
				out[6*n+i]=0; acc[6*n+i]=0;
			} else out[6*n+i]=in[Adof2+6*n+i];
		}
    } else {
		for(i=0;i<n;i++)
		for(int k=0;k<6;k++) {
        out[6*i+k]=in[Adof2+6*i+k];
		}
		// Edge dofs
		for (i=0;i<ns;i++) out[6*n+i]=in[Adof2+6*n+i];
	}
}

/** Adds in the contact forces between nodes in the mesh.
 * \param[in] in the mesh point positions.
 * \param[in] out the mesh point accelerations (cumulative). */
/*void mesh::contact_forces(double *in,double *out) {
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
}*/

/** Computes the finite-elements forces from sheet mechanics.
 * \param[in] t_ the time at which to evaluate the acceleration.
 * \param[in] in the mesh point positions. */
void mesh::fem_forces(double t_,double *in) {
    int *top=tom, tri=0;
    for(int Ti=0;Ti<n;Ti++) {
        while(top<to[Ti+1]) {
			int v[3]={Ti,*top,top[1]};
			int ed[3]={top[2],top[3],top[4]};
			double *v1=xyz+3*v[0], x1=*v1, y1=v1[1],
					*v2=xyz+3*v[1], x2=*v2, y2=v2[1],
					*v3=xyz+3*v[2], x3=*v3, y3=v3[1];
			double F_T[4]={x2-x1,x3-x1,y2-y1,y3-y1};
			double detF=F_T[0]*F_T[3]-F_T[1]*F_T[2];
			double F_inv[4] = {F_T[3], -F_T[1], 
								-F_T[2], F_T[0]};
			double fac=1/(detF*detF*detF*detF);
			double prefac=kappa*detF; // TODO: Use the same bending modulus as before?

			int argv[21]; // Global triangle dofs indices
			int j=3,j1=9,k;
			for (k=0;k<3;k++) argv[k]=6*v[k]; // Function values
			for (int i=0;i<3;i++) {
				argv[18+i]=6*n+ed[i]; // Normal derivatives
				for (k=1;k<3;k++) {
					argv[j]=6*v[i]+k; // First derivatives
					j++;
				}
				for (k=3;k<6;k++) { // Second derivatives
					argv[j1]=6*v[i]+k;
					j1++;
				}
			}
			/*for (int i=0;i<3;i++) {
				// Edges dofs
				argv[18+i]=6*n+ed[i];
				// Nodal dofs
				for (int k=0;k<6;k++) argv[6*i+k]=6*v[i]+k;
			}*/
			
			for (int J=0;J<21;J++) { // TODO - this is a lot of loops and is super slow. Consolidate? Precompute C_inv*C_inv?
			double w=0.;
			for (int I=0;I<21;I++) {
				double GH = 0.;
				for (int a=0;a<21;a++)
				for (int b=0;b<21;b++) {
					double C_prod = C_inv[441*tri+21*a+I]*C_inv[441*tri+21*b+J];
					for (int alpha=0;alpha<2;alpha++)
					for (int beta=0;beta<2;beta++) {
						int r = dmap(alpha,beta);
						for (int del=0;del<2;del++)
						for (int gam=0;gam<2;gam++) {
							int s = dmap(del,gam);
							for (int k1=0;k1<2;k1++)
							for (int k2=0;k2<2;k2++)
								GH += fac*F_inv[2*alpha+k1]*F_inv[2*beta+k1]*
										F_inv[2*del+k2]*F_inv[2*gam+k2]*
										F[9*(21*I+J)+3*r+s]*C_prod;
						}
					}
				}
				double w_i = *(in+argv[I]);
				w += w_i*GH;
			}
			P[argv[J]] += prefac*w;
			}
            top+=5; tri+=1;
        }
    }
}

int mesh::dmap(int alpha, int beta) {
	int r;
	if (alpha==0 && beta==0) r=0;
	else if (alpha==1 && beta==1) r=2;
	else r=1;
	return r;
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

// OLD CENTRALIZE()
/** Centralizes the mesh, and calculates its square width in each of the three
 * coordinate directions.
 * \param[out] (wx,wy,wz) the square widths in the three coordinate directions. */
 /*void mesh::centralize(double &wx,double &wy,double &wz) {
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
}*/


// NEW centralize()
void mesh::centralize(double &wx,double &wy,double &wz) {
	// TODO: unsure: I think all z-positions and gradients are zero initially
	arr_zeros(pts, Adof2);

    double sx=0.,sy=0.,sz=0.,fac=1./static_cast<double>(n);

    // Compute the centroid
    for(double *p=xyz;p<xyz+3*n;p+=3) {
        sx+=*p;sy+=p[1];sz+=p[2];
    }
    sx*=fac;sy*=fac;sz*=fac;

    // Displace the mesh so its centroid is at the origin, and compute the
    // variance of the vertices in each coordinate
    wx=wy=wz=0;
	// Nodes
    for(double *p=pts,*xy=xyz;p<pts+6*n;p+=6,xy+=3) {
        *xy-=sx;xy[1]-=sy;xy[2]-=sz;*p-=sz;
        wx+=*xy*(*xy);wy+=xy[1]*xy[1];wz+=*p*(*p);
    }
    wx*=fac;wy*=fac;wz*=fac;
	// Edges
	for(double *p=pts+6*n;p<pts+Adof2;p++) *p-=sz;
}

/*void mesh::damp_force(double *in,double *acc,int i,int k) {
    double *ii=in+3*(i+n),*ik=in+3*(k+n),*ai=acc+3*i,*ak=acc+3*k,
           dx=*ii-*ik,dy=ii[1]-ik[1],dz=ii[2]-ik[2];

    // Add the force contributions to the two vertices
    *ai-=dx*B;ai[1]-=dy*B;ai[2]-=dz*B;
    *ak+=dx*B;ak[1]+=dy*B;ak[2]+=dz*B;
}*/

/*void mesh::repulsive_force(double *in,double *acc,int i,int k) {
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
}*/

/** Adds an external potential to the class.
 * \param[in] ep a pointer to the potential to add. */
void mesh::add(ext_potential *ep) {
    if(n_ep==max_ep) {
        fputs("Too many external potentials\n",stderr);
        exit(1);
    }
    ex_pot[n_ep++]=ep;
}