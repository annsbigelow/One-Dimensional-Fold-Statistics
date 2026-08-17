#include <cstring>
#include <fstream>

#include "mesh.hh"

#include <Eigen/Dense>

/** Check that the physical Argyris basis functionals are "nodal" -- in the sense that they are
	1 and 0 when they should be. */
void mesh::check_dofs() {
	double bigL[ntri][21][21];
	std::ofstream outputFile("L.csv");

	// Loop through triangles
	int *top=tom,tri=0;
	for(int Ti=0;Ti<n;Ti++) {
		while(top<to[Ti+1]) {
			// L[i][j] = L_i(N_j) for the current triangle.
			double L[21][21];
			for(int a=0;a<21;a++) for(int b=0;b<21;b++) L[a][b]=0;

			// Invert C using Eigen
			double C[441];
			for(int i=0;i<21;i++){
				for(int j=0;j<21;j++) {
					C[21*i+j]=C_glob[441*tri+21*i+j];
				}
			}
			Eigen::Map<Eigen::Matrix<double,21,21,Eigen::RowMajor> > C_eig(C);
			Eigen::FullPivLU<Eigen::Matrix<double,21,21,Eigen::RowMajor> > lu(C_eig);
			if (!lu.isInvertible()) {
				printf("Error: Triangle %d change of bases matrix is not invertible.\n",tri);
				exit(1);
			}
			Eigen::Matrix<double,21,21,Eigen::RowMajor> Ce_inv = C_eig.inverse();
			double C_inv[441];
			std::memcpy(C_inv, Ce_inv.data(), 441*sizeof(double));
			
			// Fill L_i(N_j)
			int m,l,mode,d,k;
			for(int i=0;i<21;i++) {
				for(int j=0;j<21;j++) {
					k=0;
					// Vertex-value functionals
					mode=0; d=0;
					for(m=0;m<3;m++) { // Loop through vertices
						double sum=0.;
						for(l=0;l<21;l++) { // Loop through reference basis functions
							sum+=C[21*l+j]*hatL(m,mode,d,l);
						}
						L[i][j]+=C_inv[21*i+k]*sum;
						k++;
					}
					// Vertex-gradient functionals
					mode=1;
					for(m=0;m<3;m++) { 
						double sum=0.;
						d=0; // x-partial derivative
						for(l=0;l<21;l++) sum+=C[21*l+j]*hatL(m,mode,d,l);
						L[i][j]+=C_inv[21*i+k]*sum;
						k++;
						sum=0.;
						d=1; // y-partial derivative
						for(l=0;l<21;l++) sum+=C[21*l+j]*hatL(m,mode,d,l);
						L[i][j]+=C_inv[21*i+k]*sum;
						k++;
					}
					// Vertex-Hessian functionals
					mode=2;
					for(m=0;m<3;m++) { 
						double sum=0.;
						d=0; // xx-derivative
						for(l=0;l<21;l++) sum+=C[21*l+j]*hatL(m,mode,d,l);
						L[i][j]+=C_inv[21*i+k]*sum;
						k++;
						sum=0.;
						d=1; // xy-derivative
						for(l=0;l<21;l++) sum+=C[21*l+j]*hatL(m,mode,d,l);
						L[i][j]+=C_inv[21*i+k]*sum;
						k++;
						sum=0.;
						d=2; // yy-derivative
						for(l=0;l<21;l++) sum+=C[21*l+j]*hatL(m,mode,d,l);
						L[i][j]+=C_inv[21*i+k]*sum;
						k++;
					}
					// Edge-normal functionals 
					mode=3; d=0;
					for(m=0;m<3;m++) { // Loop through edges
						double sum=0.;
						for(l=0;l<21;l++) sum+=C[21*l+j]*hatL(m,mode,d,l);
						L[i][j]+=C_inv[21*i+k]*sum;
						k++;
					}
				}
			}

			for(int a=0;a<21;a++) for(int b=0;b<21;b++) bigL[tri][a][b]=L[a][b];
			top+=5; tri++;
		}
	}
	double val;
	for(int T=0;T<ntri;T++) {
		for(int a=0;a<21;a++){
			for(int b=0;b<21;b++) {
				val=bigL[T][a][b];
				if(std::abs(val)<1e-13) val=0.;
				outputFile << val << " ";
			}
			outputFile << std::endl;
		}
		outputFile << std::endl;
	}
	outputFile.close();
	printf("DOF matrix outputted to .csv file.\n");

	// Output the same matrix for the reference triangle
	double L_ref[21][21];
	double x[3]={0,1,0}, y[3]={0,0,1};
	double ptx,pty,phi_ref[21],phi_refx[21],phi_refy[21];
	double phi_refxx[21], phi_refxy[21], phi_refyy[21];
	// Vertex-valued DOFs
	for(int k=0;k<3;k++) {
		ptx=x[k]; pty=y[k];
		arg_z(ptx,pty,phi_ref);
		arg_grads(ptx,pty,phi_refx,phi_refy);
		arg_ders(ptx,pty,phi_refxx,phi_refxy,phi_refyy);
		for(int j=0;j<21;j++) {
			L_ref[k][j]=phi_ref[j];
			L_ref[3+2*k][j]=phi_refx[j]; L_ref[3+2*k+1][j]=phi_refy[j];
			L_ref[9+3*k][j]=phi_refxx[j]; 
			L_ref[9+3*k+1][j]=phi_refxy[j]; L_ref[9+3*k+2][j]=phi_refyy[j];
		}
	}
	// Edge-normal DOFs
	double mx[3]={.5,0,.5}, my[3]={0,.5,.5};
	double n_ref[6]={ 0,1,-1,0,-1/sqrt(2),-1/sqrt(2) },dn;
	for(int k=0;k<3;k++) {
		ptx=mx[k]; pty=my[k];
		arg_grads(ptx,pty,phi_refx,phi_refy);
		for(int j=0;j<21;j++) {
			dn=n_ref[2*k]*phi_refx[j]+n_ref[2*k+1]*phi_refy[j];
			L_ref[18+k][j]=dn;
		}
	}
	std::ofstream outputFile1("L ref.csv");

	for(int a=0;a<21;a++){
		for(int b=0;b<21;b++) {
			val=L_ref[a][b];
			if(std::abs(val)<1e-13) val=0.;
			outputFile1 << val << " ";
		}
		outputFile1 << std::endl;
	}
	outputFile1.close();
	printf("Reference DOF matrix outputted to .csv file.\n");
}

/** Computes \hat L(\hat N_l) for a specified reference functional \hat L.
*	\param[in] k=0,1,2 the index of the vertex or the edge.
*	\param[in] mode vertex value/gradient/Hessian or edge-normals
*	\param[in] d: 0=(dx or dxx), 1=(dy or dxy), 2=dyy
*	\param[in] l=0,...,20 the input reference basis function. */
double mesh::hatL(int k,int mode,int d,int l) {
	double refx[3]={0,1,0},refy[3]={0,0,1};
	double ptx = refx[k], pty = refy[k];
	// Vertex-value functionals
	if(mode==0) {
		double phi_ref[21];
		arg_z(ptx,pty,phi_ref);
		return phi_ref[l];
	}
	// Derivatives
	if(mode==1) {
		double phi_refx[21],phi_refy[21];
		arg_grads(ptx,pty,phi_refx,phi_refy);
		if(d==0) { // x-partial derivative
			return phi_refx[l];
		}
		if(d==1) {
			return phi_refy[l];
		}
	}
	// Second derivatives
	if(mode==2) {
		double phi_refxx[21],phi_refxy[21],phi_refyy[21];
		arg_ders(ptx,pty,phi_refxx,phi_refxy,phi_refyy);
		if(d==0) return phi_refxx[l];
		if(d==1) return phi_refxy[l];
		if(d==2) return phi_refyy[l];
	}
	// Edge-normal functionals
	if(mode==3) {
		double mrefx[3]={.5,0,.5}, mrefy[3]={0,.5,.5};
		double nrefx[3]={0,-1,-1/sqrt(2)}, nrefy[3]={1,0,-1/sqrt(2)};
		ptx=mrefx[k]; pty=mrefy[k];
		double phi_refx[21],phi_refy[21];
		arg_grads(ptx,pty,phi_refx,phi_refy);
		return nrefx[k]*phi_refx[l]+nrefy[k]*phi_refy[l];
	}
	printf("Mode unknown.\n"); return 0.;
}

/** Map a point (xhat,yhat) on the reference element to an arbitrary triangle 
	under the affine reference mapping. */
void mesh::khat2k(double B[4],double x1,double y1,double xhat,double yhat,double &x,double &y) {
	x = xhat*B[0] + yhat*B[1] + x1;
	y = xhat*B[2] + yhat*B[3] + y1;
}

/** A deluxe version of draw_mesh_gnuplot(), allowing for 
*	curvy triangle visualization by finding the solution along the edges of the triangles. */
void mesh::draw_mesh_gnuplot_deluxe(FILE *fp) {
	// Set up tables to be able to access eo[], used in edge_lookup()
	setup_springs();
	// Calculate the solution at each new edge point
	interpolate();

	int i,j,*edp=edm,*edp2,l,idx,p;
	double *pt;

    for(i=0;i<n;i++) {
        while(edp<ed[i+1]) {
            // If this edge hasn't been marked, then trace a path starting from
            // this edge
			if((*edp&bflag)==0) {
                // Mark and print this edge
                j=edge_mark(i,edp); 
				if(i>j) {printf("Error: i>j.\n"); exit(1);}
				
				// Write the generating node
				fprintf(fp,"%g %g %g\n",xyz[3*i],xyz[3*i+1],pts[6*i]);
				
				// Select the points on this edge
				l=edge_lookup(i,j);
				for(int k=0;k<np-1;k++) {
					idx=(np-1)*l+k; // Global indexing of new points
					pt=edp_vals[idx];
					// Write the new points
					fprintf(fp,"%g %g %g\n",pt[0],pt[1],pt[2]);
				}
				// Write the connecting node
                fprintf(fp,"%g %g %g\n",xyz[3*j],xyz[3*j+1],pts[6*j]);

                // Follow and print as many unmarked edges as possible
                while(find_unmarked(j,edp2)) {
					p=j;
                    j=edge_mark(j,edp2);
					if(p>j) { // Keep contiguous ordering of new points
						// Get the index of the current edge
						l=edge_lookup(p,j);
						for(int k=np-2;k>=0;k--) {
							idx=(np-1)*l+k;
							pt=edp_vals[idx];
							fprintf(fp,"%g %g %g\n",pt[0],pt[1],pt[2]);
						}
					}
					else {
						l=edge_lookup(p,j);
						for(int k=0;k<np-1;k++){
							idx=(np-1)*l+k;
							pt=edp_vals[idx];
							fprintf(fp,"%g %g %g\n",pt[0],pt[1],pt[2]);
						}
					}
                    fprintf(fp,"%g %g %g\n",xyz[3*j],xyz[3*j+1],pts[6*j]);
                }
                fputs("\n\n",fp);
			}
            edp++;
        }
    }

    // Clear the markers
    for(edp=edm;edp<edm+nc;) *(edp++)&=~bflag;
}
/** Solves for the solution values at the new edge points. */
void mesh::interpolate() {
	int *new_pts=new int[3*(np-1)]; // Contiguous storage of new points on one triangle

	edp_vals=new double*[(np-1)*ns]; // Pointers to new (x,y,z) points
	ed_vals=new double[3*(np-1)*ns]; // Memory for new points

	int *top=tom,tri=0;
	for (int Ti=0;Ti<n;Ti++)
		while (top<to[Ti+1]) {
			// Get physical triangle global indexing
			int v[3]={ Ti,*top,top[1] },
				ed[3]={ top[2],top[4],top[3] };
			
			int l1=edge_lookup(v[0],v[1]),
				l2=edge_lookup(v[1],v[2]),
				l3=edge_lookup(v[2],v[0]);

			for (int k=0;k<np-1;k++) { 
				new_pts[k]=(np-1)*l1+k; // First edge
				if(v[1]>v[2]) new_pts[2*(np-1)-1-k]=(np-1)*l2+k; // Second edge
				else new_pts[np-1+k]=(np-1)*l2+k; 
				if(v[2]>v[0]) new_pts[3*(np-1)-1-k]=(np-1)*l3+k; // Third edge
				else new_pts[2*(np-1)+k]=(np-1)*l3+k;
			}
			triangle_interpolate(new_pts,tri,v,ed);
			
			top+=5; tri++;
		}
	delete[] new_pts;
}

/** Calculates the curvy solution at each new point along
	the edges of a given triangle.
*	\param[in] new_pts an array of global indices for each new point,
				based on the indices of the edges.
*	\param[in] tri the current triangle.
*	\param[in] v the vertices' indices.
*	\param[in] ed the edges' indices.
*/
void mesh::triangle_interpolate(int* new_pts,int tri,int v[3],int ed[3]) {
	// Traverse the edges of the triangle
	int i,j;
	// Bottom edge
	j=0;
	for (i=1;i<np;i++) {
		// Point to the correct place
		edp_vals[*new_pts]=ed_vals + 3*(*new_pts);
		fill_ed_vals(edp_vals[*new_pts],tri,v,ed,i,j);
		new_pts++;
	}
	// Hypotenuse (of the reference triangle)
	for (i=np-1,j=1;j<np;j++,i--) {
		edp_vals[*new_pts]=ed_vals + 3*(*new_pts);
		fill_ed_vals(edp_vals[*new_pts],tri,v,ed,i,j);
		new_pts++;
	}
	// Left edge
	i=0;
	for (j=np-1;j>0;j--)  {
		edp_vals[*new_pts]=ed_vals + 3*(*new_pts);
		fill_ed_vals(edp_vals[*new_pts],tri,v,ed,i,j);
		new_pts++;
	}
}

/** Calculates (x,y,z) values of the solution using Argyris basis functions 
*	and inserts them into the ed_vals[] array at a specified location.
*	\param[in] tri the current triangle.
*	\param[in] v the vertices of the current triangle.
*	\param[in] ed the edges of the current triangle.
*	\param[in] i,j the local indexing of the specified point on the reference triangle.
*	\return pt the (x,y,z) values.
*/
void mesh::fill_ed_vals(double *pt,int tri,int v[3],int ed[3],int i,int j) {
	double physx, physy, soln;
	// Calculate the solution at the new point
	calculate_q(v,ed,tri,i,j,physx,physy,soln);
	// (x,y,z) values for added point
	pt[0]=physx; pt[1]=physy; pt[2]=soln;
}

/** Uses Argyris basis functions to calculate the solution 
*	at some point on a triangle. 
*	\param[in] v the vertices of the triangle.
*	\param[in] ed the edges of the triangle.
*	\param[in] tri the current triangle.
*	\param[in] i,j the local indexing of the specified point on the triangle.
*	\return (physx,physy,soln) the (x,y,z) values of the solution. */
void mesh::calculate_q(int v[3],int ed[3],int tri,int i, int j,
						double &physx,double &physy,double &soln) {
	int argv[21];
	get_argv(argv, v, ed);
	double  *v1=xyz+3*v[0], x1=*v1, y1=v1[1],
			*v2=xyz+3*v[1], x2=*v2, y2=v2[1],
			*v3=xyz+3*v[2], x3=*v3, y3=v3[1];
	double B[4]={ x2-x1,x3-x1,y2-y1,y3-y1 };
	double detF = B[0]*B[3] - B[2]*B[1];
	// The (x,y) coordinates on the physical triangle
	physx = ((np-i-j)*x1 + i*x2 + j*x3)/np;
	physy = ((np-i-j)*y1 + i*y2 + j*y3)/np;
	// Transform the physical coordinates to reference
	double refx,refy;
	k2khat(B,detF,x1,y1,refx,refy,physx,physy);

	// Evaluate the reference Argyris basis functions at the calculated point
	double phi_ref[21];
	arg_z(refx,refy,phi_ref);

	// Use the correct global normal directions
	double na[6];
	double nref[6] = {0,1,-1,0,-1,-1};
	for(int i=0;i<3;i++) {
		nhat2n(B,detF,nref[2*i],nref[2*i+1],na[2*i],na[2*i+1]);
	}
	
	float signs[21];
	for (int i=0;i<18;i++) signs[i]=1;
	for (int i=0;i<3;i++) {
		signs[18+i] = na[2*i]*normals[2*ed[i]]+na[2*i+1]*normals[2*ed[i]+1];
	}

	// Compute the solution value at the specified point 
	//		using the physical basis functions and the solution info at 
	//		the Argyris DOFs
	soln = 0.;
	double phys_phi;
	for (int k=0;k<21;k++) {
		// Map the basis function back to the physical triangle
		phys_phi = signs[k]*phys_phi_eval(tri,k,phi_ref);
		soln += pts[argv[k]]*phys_phi;
	}
}

/** Evaluates one of the physical basis functions at a given point. 
*	\param[in] tri the physical triangle we're dealing with.
*	\param[in] phi[21] the 21 Argyris basis functions evaluated at a point on 
*				the reference triangle.
*	\param[in] j the index of the basis function, from 0-20 */
double mesh::phys_phi_eval(int tri,int j,double phi[21]) {
	double out = 0.;
	for (int a=0;a<21;a++) {
		out += C_glob[441*tri+21*a+j]*phi[a];
	}
	return out;
}

/** Maps a point on an arbitrary triangle to the reference element under
	the affine reference mapping. */
void mesh::k2khat(double B[4],double detF,double x1,double y1,
			double &xhat,double &yhat,double x,double y) {
	double	diffx = x - x1, diffy = y - y1;
	xhat = ( diffx*B[3] - diffy*B[1] )/detF;
	yhat = ( -diffx*B[2] + diffy*B[0] )/detF;
}

/** Maps a normal vector of an edge of the reference element to an edge-normal on 
	an arbitrary triangle. */
void mesh::nhat2n(double B[4],double detF,double nhatx,double nhaty,double &nx,double &ny) {
	nx=(B[3]*nhatx-B[2]*nhaty)/detF;
	ny=(-B[1]*nhatx+B[0]*nhaty)/detF;
	double norm=sqrt(nx*nx+ny*ny);
	nx/=norm; ny/=norm;
}