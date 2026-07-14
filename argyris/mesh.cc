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
	delete[] normals; delete[] C_glob;
}

/** Sets up the spring network table and builds FEM matrices. */
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
                *(top++)=edge_lookup(*edp,edp[1]);
				*(top++)=edge_lookup(edp[1],i);
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
				*(top++)=edge_lookup(*edp,*ed[i]);
				*(top++)=edge_lookup(*ed[i],i);
            }
        }
        edp++;
    }

	if (quadrature)
		setup_quad_matrices();
	else
		setup_fem_matrices();
}

int mesh::edge_lookup(int i,int j) {
    if(i>j) {int k=j;j=i;i=k;}
    for(int *eop=eo[i];eop<eo[i+1];eop++) {
        if(*eop==j) return int(eop-eom);
    }
    fputs("Can't find edge index\n",stderr);
    exit(1);
}

/** Use quadrature to build the FEM change of bases, mass, and 
	stiffness matrices. */
void mesh::setup_quad_matrices() {
	global_normals();
	buildC();
	printf("Change of bases matrices have been built.\n");

	// Evaluate second derivatives of Argyris basis functions on reference triangle
	// Define the Gauss-Legendre points 
	double tmp_x[6] = { -0.932469514203152, -0.661209386466265, -0.238619186083197,
					0.238619186083197,  0.661209386466265, 0.932469514203152 };
	double tmp_w[6] = { 0.171324492379170, 0.360761573048139, 0.467913934572691,
					0.467913934572691, 0.360761573048139, 0.171324492379170 };
	//std::memcpy(xi,tmp_x,6*sizeof(double));
	//std::memcpy(w,tmp_w,6*sizeof(double));
	xi = tmp_x; w = tmp_w;

	// Build the mass matrix and global stiffness matrix
	assembleM_quad();
	assembleK_quad();
}

/** Assembles the mass matrix using 6-point Gauss-Legendre quadrature. */
void mesh::assembleM_quad() {
	// Evaluate each basis function at each Gauss-Legendre point
	double GL[21][6][6];
	double phi[21];
	for (int j=0;j<6;j++) {
		for (int i=0;i<6;i++) {
			double x = (xi[j] + 1)/2;
			double y = (1-x)*(xi[i] + 1)/2;
			arg_z(x,y,phi);
			for (int k=0;k<21;k++) {
				GL[k][j][i] = phi[k];
			}
		}
	}

	M_sp.resize(Adof2,Adof2);
	typedef Eigen::Triplet<double> T; // Use a list of triplets for fast computation
	std::vector<T> Mtriplets;
	Mtriplets.reserve(6*Adof2); // Get (estimated) memory up front for performance

	int *top=tom, tri=0, argv[21];
	double vb[6], l[3], na[6];
	for(int Ti=0;Ti<n;Ti++) { // Loop through generating indices
		while(top<to[Ti+1]) { // Loop through triangles
			// Get coordinates from reference domain
			int v[3]={Ti,*top,top[1]};
			int ed[3]={top[2],top[4],top[3]}; // Edges 1,2,3
			tri_geo(v,vb,l,na);

			// Reference map Jacobian
			double detF=vb[0]*vb[3]-vb[2]*vb[1];
			double mass = detF*rho;

			get_argv(argv,v,ed);
	
			for (int I=0;I<21;I++)
			for (int J=0;J<21;J++) { // TODO - incorporate signs[]?
				double z = mass_integral(GL,tri, I,J);
				Mtriplets.push_back(Eigen::Triplet<double>(argv[I],argv[J],mass*z));
			}

			top+=5; tri+=1;
		}
	}

	// Convert triplets list to SparseMatrix.
	// Contributions in same row,col are summed automatically.
	M_sp.setFromTriplets(Mtriplets.begin(),Mtriplets.end());
	//printf("Is the sparse mass matrix compressed? %d\n",M_sp.isCompressed());

	// Debug
	Eigen::MatrixXd M_d(M_sp);
	std::ofstream outputFile("Mass matrix.csv");
	for (int i=0;i<M_d.rows();i++) {
		for (int j=0;j<M_d.cols();j++) {
			outputFile << M_d(i,j) << " ";
		}
		outputFile << std::endl; 
	}
	outputFile.close();

	// Factorize sparse matrix using LLT Cholesky factorization
	Msolver.analyzePattern(M_sp);
	Msolver.factorize(M_sp);
	if (Msolver.info()!=Eigen::Success) {
		printf("Mass matrix factorization failed\n");
		exit(1);
	}
	printf("Mass matrix factorization finished.\n");
}

/** Performs an integration of two basis functions multiplied together
*	using 6-point Gauss-Legendre quadrature. 
*	\param[in] GL each of the 21 Argyris basis functions, evaluated 
				at the 36 quadrature points.
*	\param[in] tri a particular triangle in the "physical" space. 
	\param[in] I,J the indices of the basis functions in the integrand. */
double mesh::mass_integral(double GL[21][6][6], int tri, int I, int J) {
	double out = 0.,prefac,x;
	for (int a=0;a<21;a++)
	for (int b=0;b<21;b++){
		// Integrate
		double integral_val = 0.;
		for (int j=0;j<6;j++) {
			x = (xi[j] + 1)/2;
			prefac = w[j]*(1-x)/2;
			double sumi=0.,phi_a,phi_b;
			for (int i=0;i<6;i++) {
				// First function
				phi_a = C_glob[441*tri+21*a+I]*GL[a][j][i];
				// Second function
				phi_b = C_glob[441*tri+21*b+J]*GL[b][j][i];
				sumi += w[i]*phi_a*phi_b;
			}
			integral_val += prefac*sumi;
		}	
		integral_val /= 2; 
		out += integral_val; // TODO - Check that this is the proper implementation of the double-sum
	}
	return out;
}

/** Evaluates the monomial basis at a point (x,y).
*	\return z[], a 21-array storing the evaluated monomials. */
void mesh::monomials(double x,double y,double z[21]) {
	double tmp[21] = { 1, y, y*y, y*y*y, y*y*y*y, y*y*y*y*y,
		 x, x*y, x*y*y, x*y*y*y, x*y*y*y*y, 
		 x*x, x*x*y, x*x*y*y, x*x*y*y*y,
		 x*x*x, x*x*x*y, x*x*x*y*y,
		 x*x*x*x, x*x*x*x*y, x*x*x*x*x };
	std::memcpy(z,tmp,21*sizeof(double));
}

/** Evaluates each Argyris basis function at a point (x,y) 
	on the reference triangle. 
*	\return phi[] */
void mesh::arg_z(double x, double y, double phi[21]) {
	double z[21];
	monomials(x,y,z);

	for (int i=0;i<21;i++) {
		phi[i]=0.;
		for (int j=0;j<21;j++) {
			phi[i] += M[21*i+j]*z[j];
		}
	}
}

/** Assembles the K matrix using 6-point Gauss-Legendre quadrature. */
void mesh::assembleK_quad() {
	// Evaluate the second derivatives of each basis function at each Gauss-Legendre point
	double xxGL[21][6][6], xyGL[21][6][6], yyGL[21][6][6]; // Indexing: basis func., then x, then y
	double xx[21], xy[21], yy[21];
	for (int j=0;j<6;j++) {
		for (int i=0;i<6;i++) {
			double x = (xi[j] + 1)/2;
			double y = (1-x)*(xi[i] + 1)/2;
			arg_ders(x,y,xx,xy,yy);
			for (int k=0;k<21;k++) {
				xxGL[k][j][i] = xx[k]; xyGL[k][j][i] = xy[k]; yyGL[k][j][i] = yy[k];
			}
		}
	}

	Kd.resize(Adof2, Adof2);
	typedef Eigen::Triplet<double> T;
	std::vector<T> Ktriplets;
	// Get (estimated) memory up front for performance
	Ktriplets.reserve(6*Adof2);

	int *top=tom, tri=0, argv[21];
	double vb[6], l[3], na[6];
	for (int Ti=0;Ti<n;Ti++)
		while (top<to[Ti+1]) {
			int v[3]={ Ti,*top,top[1] };
			int ed[3]={ top[2],top[4],top[3] };
			tri_geo(v, vb, l, na);
			double B[4]={ vb[0],vb[2],vb[1],vb[3] };
			double detF=vb[0]*vb[3] - vb[2]*vb[1];
			double fac = 1/(detF*detF);
			double prefac = kappa*detF; // TODO: Use the same bending modulus as before?

			get_argv(argv, v, ed);

			double The_inv[3][3] = {{B[3]*B[3]*fac,-2*B[2]*B[3]*fac,B[2]*B[2]*fac},
									{-B[1]*B[3]*fac,B[0]*B[3]*fac+B[1]*B[2]*fac,-B[0]*B[2]*fac},
									{B[1]*B[1]*fac,-2*B[0]*B[1]*fac,B[0]*B[0]*fac}};
			for (int I=0;I<21;I++)
			for (int J=0;J<21;J++) {
				// Compute the integral using Gauss-Legendre quadrature (loop through j,i) for these particular (I,J) functions
				double z=stiff_integral(xxGL, xyGL, yyGL, The_inv, tri,I,J);
				Ktriplets.push_back(Eigen::Triplet<double>(argv[I], argv[J], prefac*z)); // TODO - incorporate signs?
			}

			top+=5; tri+=1;
		}
	// Convert triplets list to SparseMatrix.
	Kd.setFromTriplets(Ktriplets.begin(), Ktriplets.end());

	// Debug
	std::ofstream outputFile("Stiffness matrix.csv");
	Eigen::MatrixXd K_dense(Kd);
	for (int i = 0; i < K_dense.rows(); i++) {
		for (int j = 0; j < K_dense.cols(); j++) {
			outputFile << K_dense(i, j) << " ";
			double diff = abs(K_dense(i, j) - K_dense(j, i));
			if (diff > 1e-13)
				printf("Stiffness symmetry break. diff=%g\n", diff);
		}
		outputFile << std::endl;
	}
	outputFile.close();

	Eigen::SimplicialLLT<Eigen::SparseMatrix< double, Eigen::RowMajor> > llt(Kd);
	if (llt.info() == Eigen::NumericalIssue) {
		printf("Error: Stiffness matrix is not symmetric positive definite.\n");
		//exit(1);
	}
}

/** Performs an integration of the Laplacians of two basis functions 
	multiplied together using 6-point Gauss-Legendre quadrature.
*	\param[in] xxGL the xx-derivatives of each of the 21 Argyris basis 
				functions, evaluated at the 36 quadrature points.
*	\param[in] The_inv[][] the 3x3 Theta^-1 matrix for the physical-to-
*				reference transformation.
*	\param[in] tri a particular triangle in the "physical" space.
	\param[in] I,J the indices of the basis functions in the integrand. */
double mesh::stiff_integral(double xxGL[21][6][6],double xyGL[21][6][6],double yyGL[21][6][6],
							double The_inv[3][3],int tri,
							int I,int J) {
	double out = 0.,prefac,x;
	for (int a=0;a<21;a++) {
		for (int b=0;b<21;b++){
			// Integrate
			double integral_val = 0.;
			for (int j=0;j<6;j++) {
				x = (xi[j] + 1)/2;
				prefac = w[j]*(1-x)/2;
				double sumi=0.,Da,Db;
				for (int i=0;i<6;i++) {
					// First Laplacian
					Da = laplace(xxGL,xyGL,yyGL, j,i, The_inv, tri,I,a);
					// Second Laplacian
					Db = laplace(xxGL,xyGL,yyGL, j,i, The_inv, tri,J,b);
					sumi += w[i]*Da*Db;
				}
				integral_val += prefac*sumi;	
			}	
			integral_val /= 2; 
			out += integral_val; // TODO - I think is the proper implementation of the double-sum
		}
	}
	return out;
}

/** Evaluates the second derivatives of the monomial basis at a point (x,y).
*	\return z, a 21-array storing the evaluated monomial derivatives. */
void mesh::ders(double x,double y,double mxx[21],double mxy[21],double myy[21]) {
	double tmp[21] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
		2, 2*y, 2*y*y, 2*y*y*y, 6*x, 6*x*y, 6*x*y*y,
		12*x*x, 12*x*x*y, 20*x*x*x};
	double tmp1[21] = {0, 0, 0, 0, 0, 0, 0, 1, 2*y, 3*y*y, 
			4*y*y*y, 0, 2*x, 4*x*y, 6*x*y*y, 0, 3*x*x,
			6*x*x*y, 0, 4*x*x*x, 0};
	double tmp2[21] = {0, 0, 2, 6*y, 12*y*y, 20*y*y*y,
            0, 0, 2*x, 6*x*y, 12*x*y*y,
            0, 0, 2*x*x, 6*x*x*y,
			0, 0, 2*x*x*x, 0, 0, 0};
	std::memcpy(mxx,tmp,21*sizeof(double));
	std::memcpy(mxy,tmp1,21*sizeof(double));
	std::memcpy(myy,tmp2,21*sizeof(double));
}

/** Evaluates the second derivatives of each Argyris basis function at a point (x,y)
	on the reference triangle.
*	\return xx,xy,yy */
void mesh::arg_ders(double x, double y, double xx[21],double xy[21],double yy[21]) {
	// Evaluate second derivatives in monomial basis
	double mxx[21],mxy[21],myy[21];
	ders(x,y,mxx,mxy,myy);

	for (int i=0;i<21;i++) {
		xx[i]=0.; xy[i]=0.; yy[i]=0.;
		for (int j=0;j<21;j++) {
			xx[i] += M[21*i+j]*mxx[j];
			xy[i] += M[21*i+j]*mxy[j];
			yy[i] += M[21*i+j]*myy[j];
		}
	}
}

/** Calculates the Laplacian of a basis function in reference coordinates, 
	evaluated at a quadrature point. */
double mesh::laplace(double xxGL[21][6][6],double xyGL[21][6][6],double yyGL[21][6][6],
						int j,int i, 
						double The_inv[3][3],int tri,
						int I,int a) {
	double dxx, dyy;
	dxx = xxGL[a][j][i]*The_inv[0][0] + xyGL[a][j][i]*The_inv[1][0] + yyGL[a][j][i]*The_inv[2][0];
	dyy = xxGL[a][j][i]*The_inv[0][2] + xyGL[a][j][i]*The_inv[1][2] + yyGL[a][j][i]*The_inv[2][2];
	double out = C_glob[441*tri+21*a+I]*(dxx + dyy);

	return out;
}

/** Build the FEM change of bases, mass, and stiffness matrices. */
void mesh::setup_fem_matrices() {
	// Define the global normal directions
	global_normals();

	// Build change of bases matrices
	buildC();
	printf("Change of bases matrices have been built.\n");

	// Build the mass matrix and global stiffness matrix
	assemble_M();
	assemble_K();
}

/** Assembles the mass matrix using pre-computed integral values. */
void mesh::assemble_M() {
	M_sp.resize(Adof2,Adof2);
	typedef Eigen::Triplet<double> T; // Use a list of triplets for fast computation
	std::vector<T> triplets;
	triplets.reserve(6*Adof2); // Get (estimated) memory up front for performance

	int *top=tom, tri=0, argv[21];
	double vb[6], l[3], na[6];
	for(int Ti=0;Ti<n;Ti++) { // Loop through generating indices
		while(top<to[Ti+1]) { // Loop through triangles
			// Get coordinates from reference domain
			int v[3]={Ti,*top,top[1]};
			int ed[3]={top[2],top[4],top[3]}; // Edges 1,2,3
			tri_geo(v,vb,l,na);

			// Reference map Jacobian
			double detF=vb[0]*vb[3]-vb[2]*vb[1];
			double mass = detF*rho;
			
			// Check the local normal direction against global
			float signs[21];
			for (int i=0;i<18;i++) signs[i]=1;
			for (int i=0;i<3;i++) {
				signs[18+i] = na[2*i]*normals[2*ed[i]]+na[2*i+1]*normals[2*ed[i]+1];
			}

			get_argv(argv,v,ed);
	
			for (int I=0;I<21;I++)
			for (int J=0;J<21;J++) {
				double sum=0.;
				for (int a=0;a<21;a++)
				for (int b=0;b<21;b++) {
					//sum += C_glob[441*tri+21*a+I]*C_glob[441*tri+21*b+J]*S[21*a+b];
					sum += signs[a]*signs[b]*C_glob[441*tri+21*a+I]*C_glob[441*tri+21*b+J]*S[21*a+b];
				}
				triplets.push_back(Eigen::Triplet<double>(argv[I],argv[J],mass*sum));
			}

			top+=5; tri+=1;
		}
	}

	// Convert triplets list to SparseMatrix.
	// Contributions in same row,col are summed automatically.
	M_sp.setFromTriplets(triplets.begin(),triplets.end());
	//printf("Is the sparse mass matrix compressed? %d\n",M_sp.isCompressed());

	// Debug
	Eigen::MatrixXd M_d(M_sp);
	std::ofstream outputFile("Mass matrix.csv");
	for (int i=0;i<M_d.rows();i++) {
		for (int j=0;j<M_d.cols();j++) {
			outputFile << M_d(i,j) << " ";
		}
		outputFile << std::endl; 
	}
	outputFile.close();

	// Factorize sparse matrix using LLT Cholesky factorization
	Msolver.analyzePattern(M_sp);
	Msolver.factorize(M_sp);
	if (Msolver.info()!=Eigen::Success) {
		printf("Mass matrix factorization failed\n");
		exit(1);
	}
	printf("Mass matrix factorization finished.\n");
}

/** Assembles the stiffness matrix from the biharmonic term in the FEM computations.
*	Also checks that the stiffness matrix is positive definite. */
void mesh::assemble_K() {
	Kd.resize(Adof2, Adof2);
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	// Get (estimated) memory up front for performance
	triplets.reserve(6*Adof2);

	int *top=tom, tri=0, argv[21];
	double vb[6],l[3],na[6];
	for (int Ti=0;Ti<n;Ti++)
		while (top<to[Ti+1]) {
			int v[3]={ Ti,*top,top[1] };
			int ed[3]={ top[2],top[4],top[3] };
			tri_geo(v,vb,l,na);
			double B[4]={ vb[0],vb[2],vb[1],vb[3] };
			double detF = vb[0]*vb[3] - vb[2]*vb[1];
			double fac = 1/(detF*detF);
			double prefac = kappa*detF; // TODO: Use the same bending modulus as before?

			get_argv(argv,v,ed);

			float signs[21];
			for (int i=0;i<18;i++) signs[i]=1;
			for (int i=0;i<3;i++)
				signs[18+i] = na[2*i]*normals[2*ed[i]]+na[2*i+1]*normals[2*ed[i]+1];

			double The[3] = { (B[3]*B[3]+B[1]*B[1])*fac,
				-2*(B[2]*B[3]+B[0]*B[1])*fac, (B[2]*B[2]+B[0]*B[0])*fac };

			for (int I=0;I<21;I++)
				for (int J=0;J<21;J++) {
					double HaHb=0.;
					for (int a=0;a<21;a++)
						for (int b=0;b<21;b++) {
							//double C_prod = C_glob[441*tri+21*a+I]*C_glob[441*tri+21*b+J];
							double C_prod=signs[a]*signs[b]*C_glob[441*tri+21*a+I]*C_glob[441*tri+21*b+J];
							for (int r=0;r<3;r++)
								for (int s=0;s<3;s++)
									HaHb += C_prod*The[r]*The[s]*F[9*(21*a+b)+3*r+s];
						}
					triplets.push_back(Eigen::Triplet<double>(argv[I],argv[J],prefac*HaHb));
				}
			top+=5; tri+=1;
		}
	// Convert triplets list to SparseMatrix.
	Kd.setFromTriplets(triplets.begin(), triplets.end());

	// Debug
	std::ofstream outputFile("Stiffness matrix.csv");
	Eigen::MatrixXd K_dense(Kd);
	for (int i=0;i<K_dense.rows();i++) {
		for (int j=0;j<K_dense.cols();j++) {
			outputFile << K_dense(i,j) << " ";
			double diff = abs(K_dense(i,j) - K_dense(j,i));
			if (diff > 1e-13)
				printf("Stiffness symmetry break. diff=%g\n", diff);
		}
		outputFile << std::endl;
	}
	outputFile.close();

	Eigen::SimplicialLLT<Eigen::SparseMatrix< double, Eigen::RowMajor> > llt(Kd);
	if (llt.info()==Eigen::NumericalIssue) {
		printf("Error: Stiffness matrix is not symmetric positive definite.\n");
		//exit(1);
	}
}


/** Define the global normal orientations */
void mesh::global_normals() {
	normals = new double[2*ns];

	// Keep track of edges which have been seen already
	int *seen = new int[ns];
	for (int i=0;i<ns;i++) seen[i]=0;

	int *top=tom;
	double vb[6], l[3], na[6];
	// Loop through triangles
	for (int Ti=0;Ti<n;Ti++) 
	while (top<to[Ti+1]) {
		int ed[3]={top[2],top[4],top[3]};
		int v[3]={Ti,*top,top[1]};
		tri_geo(v, vb, l, na);

		for (int i=0;i<3;i++) 
		if (!seen[ed[i]]) {
			normals[2*ed[i]]=na[2*i];
			normals[2*ed[i]+1]=na[2*i+1];
			// Mark this edge as "seen"
			seen[ed[i]]=1;
		}
		top+=5;
	}
	delete[] seen;
}

/** Collects geometry for a given triangle.
*	\param[in] v the indices for the vertices of the triangle.
*	\return vb the vectors defining the sides of the triangle.
*	\return l the lengths of each side.
*	\return na the normal vectors at each side, found by rotating each
*		element of vb by pi/2.
*/
void mesh::tri_geo(int v[3], double* vb, double* l, double* na) {
	double  *v1=xyz+3*v[0], x1=*v1, y1=v1[1],
			*v2=xyz+3*v[1], x2=*v2, y2=v2[1],
			*v3=xyz+3*v[2], x3=*v3, y3=v3[1];
	// Sides of the triangle
	vb[0]=x2-x1; vb[1]=y2-y1;
	vb[2]=x3-x1; vb[3]=y3-y1;
	vb[4]=x3-x2; vb[5]=y3-y2;

	// Lengths of each side
	l[0]=sqrt(vb[0]*vb[0]+vb[1]*vb[1]);
	l[1]=sqrt(vb[2]*vb[2]+vb[3]*vb[3]);
	l[2]=sqrt(vb[4]*vb[4]+vb[5]*vb[5]);

	// Normal vectors
	na[0]=-vb[1]/l[0]; na[1]=vb[0]/l[0];
	na[2]=-vb[3]/l[1]; na[3]=vb[2]/l[1];
	na[4]=-vb[5]/l[2]; na[5]=vb[4]/l[2];
}

/** Get the global indices for a given triangle.
*	\param[in] v the vertices of the triangle.
	\param[in] ed the edges of the triangle. 
	\return argv the global indexing. */
void mesh::get_argv(int* argv, int v[3], int ed[3]) {
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
}

// Build change of bases matrices for each triangle
void mesh::buildC() {
	C_glob = new double[ntri*441];

	int *top=tom;
	double D[504], E[504];
	int tri=0,i,j;
	double vb[6], l[3], na[6];
	for (int Ti=0;Ti<n;Ti++) // Loop through triangles
	while (top<to[Ti+1]) {
		// Vertices 1,2,3
		int v[3]={Ti,*top,top[1]};
		tri_geo(v,vb,l,na);
		// Lengths of each side, squared 
		double l2[3];
		for (i=0;i<3;i++) l2[i]=l[i]*l[i];
		// Affine transformation Jacobian
		double B[4]={ vb[0],vb[2], vb[1],vb[3] };
		// Hessian
		double H[9]={	B[0]*B[0], 2*B[0]*B[2], B[2]*B[2],
						B[1]*B[0], B[1]*B[2]+B[0]*B[3], B[2]*B[3],
						B[1]*B[1], 2*B[3]*B[1], B[3]*B[3]	};
		// Normal vectors (un-normalized)
		double Rv[6];
		for (i=0;i<3;i++) {
			Rv[2*i]=na[2*i]*l[i];
			Rv[2*i+1]=na[2*i+1]*l[i];
		}
		// Change of bases for the normal derivatives
		double a[6]={ vb[2]/l2[0], vb[3]/l2[0],
					-vb[0]/l2[1], -vb[1]/l2[1],
			-(vb[0]+vb[2])/(sqrt(2)*l2[2]), -(vb[1]+vb[3])/(sqrt(2)*l2[2]) };

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

		for (int i=0;i<441;i++) C_glob[441*tri+i]=C[i];
		top+=5; tri+=1;
	}
}

/** Displace the points in the "z" direction with a Gaussian */
void mesh::Gauss_displacement() {
	const double eps=0.001;
	// Initialize function values and gradients at all nodes
	for(int i=0;i<n;i++) {
		double x=xyz[3*i], y=xyz[3*i+1];
		double power = exp(-eps*(x*x+y*y));
		// Displace the z-component
		pts[6*i]+=-eps+0.02*power;
		// First derivatives
		pts[6*i+1]+=-eps*.04*x*power;
		pts[6*i+2]+=-eps*.04*y*power;
		// Second derivatives
		pts[6*i+3]+=eps*(-.04 + eps*.08*x*x)*power;
		pts[6*i+4]+=eps*eps*.08*x*y*power;
		pts[6*i+5]+=eps*(-.04 + eps*.08*y*y)*power;
	}

	// Keep track of edges which have been seen already
	int *seen = new int[ns];
	for (int i=0;i<ns;i++) seen[i]=0;

	// Initialize normal derivatives
	int *top=tom;
	for (int Ti=0;Ti<n;Ti++)
	while (top<to[Ti+1]) {
		int v[3]={Ti,*top,top[1]};
		int ed[3]={top[2],top[4],top[3]};
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

		int j=0;
		for (int i=0;i<3;i++) {
			if (!seen[ed[i]]) {
				pts[6*n+ed[i]] += -eps*.04*exp(-eps*(m[j]*m[j]+m[j+1]*m[j+1]))*(m[j]*na[j]+m[j+1]*na[j+1]);
				seen[ed[i]]=1;
			}
			j+=2;
		}
		top+=5;
	}
}

void mesh::linear_gradient() {
	// Slope of the gradient
	const double eps=.001;
	// Initialize function values and gradients at all nodes
	for(int i=0;i<n;i++) {
		double x=xyz[3*i], y=xyz[3*i+1];
		// Displace the z-component
		pts[6*i]+=eps*(x+y);
		// First derivatives
		pts[6*i+1]+=eps;
		pts[6*i+2]+=eps;
	}

	// Keep track of edges which have been seen already
	int *seen = new int[ns];
	for (int i=0;i<ns;i++) seen[i]=0;

	// Initialize normal derivatives
	int *top=tom;
	double vb[6], l[3], na[6];
	for (int Ti=0;Ti<n;Ti++)
	while (top<to[Ti+1]) {
		int v[3]={Ti,*top,top[1]};
		int ed[3]={top[2],top[4],top[3]};
		tri_geo(v,vb,l,na);

		int j=0;
		for (int i=0;i<3;i++) {
			if (!seen[ed[i]]) {
				pts[6*n+ed[i]] += eps*(na[j]+na[j+1]);
				seen[ed[i]]=1;
			}
			j+=2;
		}
		top+=5;
	}
}

void mesh::const_displacement() {
	const float eps = .001;
	for (int i=0;i<n;i++) pts[6*i]+=eps;
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
	printf("Edges\n");
	for (int i=0;i<ns;i++) printf("%g ",pt_array[6*n+i]);
	printf("\n\n");
}

void mesh::mesh_ff(double t_,double *in,double *out) {
    double *acc=out+Adof2;
    int i;
	// Add biharmonic term from FEM computations
	Kq_multiply(in);

	Eigen::VectorXd f_sum(Adof2), av(Adof2);
	double *vp=in+Adof2;
	for (i=0;i<Adof2;i++)
		f_sum[i] = Kq[i]-drag*vp[i];
	// Cholesky direct solver, a = M^-1*f
	av=Msolver.solve(f_sum);
	if (Msolver.info()!=Eigen::Success) printf("Matrix solving failed\n");
	std::memcpy(acc,av.data(),Adof2*sizeof(double));

	// Debug: Calculate the residual
	/*Eigen::MatrixXd M_d(M_sp);
	for (int i=0;i<Adof2;i++) {
		double Mx = 0.;
		for (int j=0;j<Adof2;j++) {
			Mx += M_d(i,j)*acc[j];
		}
		double residual = Mx-f_sum[i];
		printf("Residual: %g\n",residual);
	}
	printf("\n");*/

    // Assemble the velocities in the first part of the out array. 
	for(i=0;i<n;i++)
		for(int k=0;k<6;k++) {
			out[6*i+k]=in[Adof2+6*i+k];
		}
	for (i=0;i<ns;i++) out[6*n+i]=in[Adof2+6*n+i];
}

/** Computes the biharmonic term from sheet mechanics.
 * \param[in] in the mesh degrees of freedom. */
void mesh::Kq_multiply(double *in) {
	// Copy "in" into Eigen vector
	// .eval() creates a copy so that "in" is not modified
	Eigen::VectorXd in_e = Eigen::Map<Eigen::VectorXd>(in, Adof2).eval();
	Kq.setZero();

	// Do a matrix-vector product with K
	//printf("Before multiply:\n"); // Debug
	//print_pts(in);
	Kq = Kd*in_e;
	/*double *tmp = new double[Adof2]; // Debug
	std::memcpy(tmp,Kq.data(),Adof2*sizeof(double));
	printf("Kq:\n");
	print_pts(tmp);
	delete[] tmp;*/
}

/** Test a multiplication by a local stiffness matrix on one triangle. 
*	\param[in] tri the triangle to test 
	\param[in] Ti the generating node of the triangle. */
void mesh::local_Kq_multiply(int tri, int Ti) {
	Eigen::MatrixXd K_dense(Kd);

	int *top = tom + 5*tri;
	int v[3] = {Ti, top[0],top[1]};
	int ed[3] = {top[2],top[4],top[3]};
	int argv[21];
	get_argv(argv,v,ed);

	double loc_Kq[21];
	printf("local multiply:\n");
	for (int I=0;I<21;I++) {
		loc_Kq[I]=0.;
		for (int J=0;J<21;J++) {
			if (argv[I]==11) {
				//printf("K val: %g, pts val: %g, test: %g\n",K_dense(argv[I],argv[J]),pts[argv[J]],K_dense(argv[I],argv[J])*pts[argv[J]]);
			}
			loc_Kq[I] += K_dense(argv[I],argv[J])*pts[argv[J]];
		}
		printf("index: %d, val: %g\n",argv[I],loc_Kq[I]);
	}
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

void mesh::centralize(double &wx,double &wy,double &wz) {
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