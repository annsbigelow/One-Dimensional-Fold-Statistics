#include <cstring>
#include "mesh.hh"
#include "vec3.hh"

/** Finds the row and column position of a node, given its global index.
	Valid for a regular hexagonal topology.
*	\param[in] nx The horizontal dimension of the sheet.
*	Input row as 0 and col as the global index.
*	\param[out] row,col the row/col indices
*/
int mesh::find_pos_rec(int& row,int& col,int nx) {
	int nt;
	while (true) {
		nt=nx+(row&1);
		if(col<nt) break;
		// Move up a row
		col-=nt;
		row++;
	}
	return nt;
}

/** Calculates the Developed Interfacial Area Ratio as a roughness measure.
*	(UNFINISHED)
*/
/*double mesh::Sdr(int nx,int ny) {
	// Loop through all points 
	for(int i=0;i<n;i++) {
		int row=0, col=i;
		// Check if the current point is within the subsheet, minus another layer (due to the gradient calculation)
		int nt=find_pos_rec(row,col,nx);
		if(inside(row,col,nt,ny,R+1)) {
			pt=pts+3*i; 
			right=pts+3*(i+1); left=pts+3*(i-1); 
			//up=pts+3*(i+nt)
		}
	}
	A=tot_area_rec(nx,ny);
}*/

/** Calculates the statistical standard deviation of the height and the arithmetical mean height of the sheet.
*	\param[in] nx,ny The dimensions of the sheet.
*/
void mesh::Sq_Sa(double& Sq,double& Sa,int nx, int ny) {
	double* pt, mean=0., sd=0.; Sa=0.;
	double* z=new double[n];
	bool* intflag=new bool[n];
	int intn=0;

	// Loop through each node
	for (int i=0;i<n;i++) {
		int row=0, col=i;
		int nt=find_pos_rec(row,col,nx);
		// Collect statistics only for interior nodes.
		intflag[i]=inside(row,col,nt,ny,R);
		if(intflag[i]) {
			pt=pts+3*i; z[i]=pt[2];
			mean+=z[i]; intn+=1;
			Sa+=std::abs(z[i]);
		}
	}
	mean/=intn; Sa/=intn;
	for (int i=0;i<n;i++) {
		if(intflag[i]) sd+=(z[i]-mean)*(z[i]-mean);
	}
	sd/=(intn-1);

	delete[] z; delete[] intflag;
	Sq=sqrt(sd);
}

/** Calculates the surface area of the sheet.
*	\param[in] nx,ny The dimensions of the sheet.
*/
double mesh::tot_area_rec(int nx, int ny) {
	int* tp = to[0], i, k, l, ci, cj, di, dj, ei, ej, nt, nt1, nt2;
	double A=0., magn;
	// Make sure all triangles are counted 
	//FILE* fp2 = safe_fopen("tri_int.gnu", "wb");

	// Loop through all of the triangles
	for (i=0;i<n;i++) while (tp<to[i+1]) {
		k=tp[1], l=tp[2];
		// Check if i,k,l are ALL within the interior. If not, don't add the area to the total.
		ci=0, cj=i, di=0, dj=k, ei=0, ej=l;
		nt=find_pos_rec(ci,cj,nx); nt1=find_pos_rec(di,dj,nx); nt2=find_pos_rec(ei,ej,nx);
		if (inside(ci,cj,nt,ny,R) && inside(di,dj,nt1,ny,R) && inside(ei,ej,nt2,ny,R)) {
			double* ii=pts+3*i, *ik=pts+3*k, *il=pts+3*l;
			vec3 c(*ik-*ii, ik[1]-ii[1], ik[2]-ii[2]),
				d(*il-*ii, il[1]-ii[1], il[2]-ii[2]), e;

			//fprintf(fp2, "%g %g %g\n%g %g %g\n%g %g %g\n", *il, il[1], il[2], *ii, ii[1], ii[2], *ik, ik[1], ik[2]);

			// Compute cross product of two edges for one triangle out of the pair
			e=d*c; magn=e.magnitude();
			A+=magn/2;
		}
		tp+=3;
	}
	//fclose(fp2);
	return A;
}