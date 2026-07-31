#include <cstring>

#include "mesh.hh"

/** A deluxe version of draw_mesh_gnuplot(), allowing for 
*	curvy triangle visualization by finding the solution along the edges
*	of the triangles. */
void mesh::draw_mesh_gnuplot_deluxe(FILE *fp) {
	// Set the number of partitions along the edges.
	np=6;
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
	// Hypotenuse
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
	// Transform the physical coordinates to reference triangle
	double	diffx = physx - x1,
			diffy = physy - y1;
	double	refx = ( diffx*B[3] - diffy*B[1] )/detF,
			refy = ( -diffx*B[2] + diffy*B[0] )/detF;
	// Evaluate the reference Argyris basis functions at the calculated point
	double phi_ref[21];
	arg_z(refx,refy,phi_ref);

	// Compute the solution value at the specified point 
	//		using the physical basis functions and the solution info at 
	//		the Argyris DOFs
	soln = 0.;
	double phys_phi;
	for (int k=0;k<21;k++) {
		// Map the basis function back to the physical triangle
		phys_phi = phys_phi_eval(tri,k,phi_ref);
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