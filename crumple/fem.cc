
void matrix_inversion() {

}


void assemble_MP(q,lam,mu,rho) {
	// M will be sparse
	// Initialize M,P full of zeros, size 3nx3n and 3nx1 respectively, elsewhere
	// Maybe make lam, mu global params and dPdX global

	// Gradients of basis functions listed as dPsi/dX, dPsi/dY
	int dPdX[6]={-1,1,0, 
			-1,0,1};

	// Assemble rows of P one at a time
	// P = zeros(3n,1)

	int S[9]={2,1,1, 1,2,1, 1,1,2};

	// Assemble the mass matrix M and "load" P
	// Loop through elements 
	for (int T=0; T<num_elements; T++) {
		// Select the vertices' indices
		int* v=elements_array+3*T;
		int d=*v, e=v[1], f=v[2];
		// Extract the coordinates from the reference domain
		double *v1=sh_pts+3*d, x1=*v1, y1=v1[1],
				*v2=sh_pts+3*e, x2=*v2, y2=v2[1],
				*v3=sh_pts+3*f, x3=*v3, y3=v3[1];

		// Reference mapping matrix 
		double F[4]={x2-x1,x3-x1,
					y2-y1,y3-y1};
		double detF=F[0]*F[3]-F[1]*F[2];

		// Compute the entries of M corresponding to this element since all other element contributions will be zero
		// Loop through nodes and vector components of nodes
		int vd[3]={d,e,f};
		for (int i=0;i<3;i++) 
			for (int j=0;j<3;j++) {

				int I=vd[i], J=vd[j];
				double* qi=q+3*I;
				// Not sure about casting dPdX elements here
				double FdPJ[2]={ ((double)dPdX(j)*F[3]-(double)dPdX(3+j)*F[2])/detF, ((double)dPdX(3+j)*F[0]-(double)dPdX(j)*F[1])/detF };

				for (int k=0;k<3;k++) {
					int row=3*I+k;

					// Assemble P blocks
					double hatP_k[2];
					get_hatP_k(hatP_k,k,lam,mu);
					double Ak=detF*(hatP_k[0]*FdPJ[0]+hatP_k[1]*FdPJ[1]);
					int tmp=(detF>0?1:-1);
					P[row]+=tmp*Ak/2;

					for (int ii=0; ii<3;ii++) {
						// Assemble M blocks
						int col=3*J+ii;
						M[3*n*row+col]+=std::abs(detF)*rho*S[3*i+j]/24;
					}
				}
			}
	
	}
}


/* k-th row of \hat P */
void get_hatP_k(double (&hatP_k)[2],int k,double lambda,double mu) {
	double C=0.;
	double D[2]={0.,0.};
	double E[2]={0.,0.};
	for (int p=0;p<3;p++) {
		double A=qSum(p,dPdX,0), B=qSum(p,dPdX,1);
		C+=A*A+B*B;
		D[0]+=A*A; D[1]+=A*B;
		E[0]+=B*A; E[1]+=B*B;
	}

	for (int l=0;l<2;l++) hatP_k[l]=qSum(k,dPdX,0)*(lambda*(C-2)*(l==0?1:0)/2 + mu*(D[l]-(l==0?1:0))) +
							qSum(k,dPdX,1)*(lambda*(C-2)*(l==1?1:0)/2 + mu*(E[l]-(l==1?1:0)));

}

/** Sums q_i[k]\nabla\psi over all elements. Used in computing the stress \hat P.
*	\param[in] k the component of q
*	\param[in] m the component of X to differentiate with respect to
*/
double qSum(int k,int dPdX,int m) {
	double result=0.;
	for (int T=0; T<num_elements; T++) {
		int* v=elements_array+3*T;
		// Vertices
		int d=*v, e=v[1], f=v[2];
		// Vertex coordinates
		double *v1=sh_pts+3*d, x1=*v1, y1=v1[1],
				*v2=sh_pts+3*e, x2=*v2, y2=v2[1],
				*v3=sh_pts+3*f, x3=*v3, y3=v3[1];

		// Reference mapping matrix 
		double F[4]={x2-x1,x3-x1, y2-y1,y3-y1};
		double detF=F[0]*F[3]-F[1]*F[2];

		// Loop through vertices
		int vd[3]={d,e,f};
		for (int i=0;i<3;i++) {
			int I=vd[i];
			double FdPI[2]={ ((double)dPdX(i)*F[3]-(double)dPdX(3+i)*F[2])/detF, ((double)dPdX(3+i)*F[0]-(double)dPdX(i)*F[1])/detF };
			// Choose the current I-th node and sum all contributions
			double* qi=q+3*I;
			result+=qi[k]*FdPI[m];
		}
	}
	return result;
}