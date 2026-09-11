#include "mesh.hh"

#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>

/** Reads in the integer part of an input file, containing the topological
 * information about the mesh.
 * \param[in] fp a file handle to read from. */
void mesh::read_topology(FILE *fp) {

    // Find the number of nodes and the number of edges
    safe_fread(&n,sizeof(int),2,fp,"mesh size information");
    if(n<=0) fatal_error("mesh size must be positive",1);
    if(nc<=0) fatal_error("the number of connections must be positive",1);
    //printf("# nodes: %d, total connections: %d\n",n,nc);

    // Read in connection numbers
    ncn=new unsigned int[n];
    safe_fread(ncn,sizeof(int),n,fp,"connection numbers");

    // Allocate memory and set up pointers
    ed=new int*[n+1];
    edm=new int[nc];
    *ed=edm;
    for(int i=1;i<=n;i++) ed[i]=ed[i-1]+(ncn[i-1]&~bflag);

    // Read in connections
    safe_fread(edm,sizeof(int),nc,fp,"connections");
}

/** Reads the vertex positions from an input file.
 * \param[in] fp a file handle to read from. */
void mesh::read_positions(FILE *fp) {

    // Determine the total number of edges
    ns=nc>>1;

	// Set the number of degrees of freedom
	Adof=2*(6*n+ns);
	Adof2=6*n+ns;

    // Set up the memory for the integration routines
    pts=new double[Adof];vel=pts+Adof2;

	// Read in the points and set velocities to zero
	if(r_all_dofs) {
		double *tmp=new double[8*n+ns];
		safe_fread(tmp,sizeof(double),8*n+ns,fp,"Argyris positions");

		for(int i=0;i<Adof2;i++) vel[i]=0.;

		xyz=new double[3*n];
		for(int i=0;i<n;i++) {
			xyz[3*i]=tmp[8*i]; xyz[3*i+1]=tmp[8*i+1]; xyz[3*i+2]=tmp[8*i+2];
			for(int k=0;k<6;k++) pts[6*i+k]=tmp[8*i+2+k];
		}
		for(int i=0;i<ns;i++) pts[6*n+i]=tmp[8*n+i];
		delete[] tmp;
	}
	else {
		xyz=new double[3*n];
		safe_fread(xyz,sizeof(double),3*n,fp,"vertex positions");
		for(int i=0;i<Adof2;i++) vel[i]=0.;
		for(int i=0;i<n;i++) pts[6*i]=xyz[3*i+2];
	}
}

/** Saves the node positions as text to a file.
 * \param[in] fp a file handle to write to. */
void mesh::draw_nodes(FILE *fp) {
    for(int i=0;i<n;i++)
        fprintf(fp,"%g %g %g\n",xyz[3*i],xyz[3*i+1],pts[6*i]);
}

/** Saves the node positions as POV-Ray spheres to a file.
 * \param[in] fp a file handle to write to. */
void mesh::draw_nodes_pov(FILE *fp) {
    for(int i=0;i<n;i++) fprintf(fp,"sphere{<%g,%g,%g>,r}\n",xyz[3*i],xyz[3*i+1],pts[6*i]);
}

/** Draws the mesh of edges to a text file for plotting with Gnuplot. Since
 * Gnuplot is more efficient at drawing contiguous lines, the routine outputs
 * the mesh by tracing out contiguous paths of unprinted edges until it get
 * stuck.
 * \param[in] fp a file handle to write to. */
void mesh::draw_mesh_gnuplot(FILE *fp) {
    int i,j,*edp=edm,*edp2;
    for(i=0;i<n;i++) {
        while(edp<ed[i+1]) {

            // If this edge hasn't been marked, then trace a path starting from
            // this edge
            if((*edp&bflag)==0) {

                // Mark and print this edge
                j=edge_mark(i,edp);
                fprintf(fp,"%g %g %g\n%g %g %g\n",
                           xyz[3*i],xyz[3*i+1],pts[6*i],
                           xyz[3*j],xyz[3*j+1],pts[6*j]);

                // Follow and print as many unmarked edges as possible
                while(find_unmarked(j,edp2)) {
                    j=edge_mark(j,edp2);
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

/** Marks an edge by setting the last bit to 1. It also find and marks the
 * reverse edge.
 * \param[in] i the current vertex.
 * \param[in] edp a pointer to the edge to consider.
 * \return The vertex that the pointer pointed to. */
int mesh::edge_mark(int i,int *edp) {
    int j=*edp;
    *edp|=bflag;
    for(int *p=ed[j];p<ed[j+1];p++) if(*p==i) {
        *p|=bflag;return j;
    }
    fputs("Couldn't find reverse edge to mark\n",stderr);
    exit(1);
}

/** Draws the mesh in the POV-Ray mesh2 format.
 * \param[in] fp a file handle to write to.
 * \param[in] normals whether to include the normal vectors or not. */
void mesh::draw_mesh_pov(FILE *fp,bool normals) {
    int i,j,k,*edp=*ed,su=0;
    double *xy,*pp,x,y,z;

    // Output the vertex vectors
    fprintf(fp,"mesh2{\n\tvertex_vectors{\n\t\t%d",n);
    for(xy=xyz,pp=pts;pp<pts+6*n;xy+=3,pp+=6) fprintf(fp,",\n\t\t<%g,%g,%g>",*xy,xy[1],*pp);

    // Output the normal vectors
    if(normals) {
        fprintf(fp,"\n\t}\n\tnormal_vectors{\n\t\t%d",n);
        for(i=0;i<n;i++) {
            normal_vector(i,x,y,z);
            fprintf(fp,",\n\t\t<%g,%g,%g>",x,y,z);
        }
    }

    // Count the face indices
    for(i=0;i<n;i++) su+=(ncn[i]&~bflag)-(ncn[i]&bflag?1:0);
    fprintf(fp,"\n\t}\n\tface_indices{\n\t\t%d",su/3);

    // Output the face indices
    for(i=0;i<n;i++) if(edp!=ed[i+1]) {
        j=k=*(edp++);
        while(edp<ed[i+1]) {
            if(j>i&&*edp>i) fprintf(fp,",\n\t\t<%d,%d,%d>",i,j,*edp);
            j=*(edp++);
        }
        if((ncn[i]&bflag)==0&&j>i&&k>i) fprintf(fp,",\n\t\t<%d,%d,%d>",i,j,k);
    }
    fputs("\n\t}\n\ttexture{t_mesh}\n}\n",fp);
}

/** Draw the mesh using cylinders and spheres in POV-Ray format.
 * \param[in] fp a file handle to write to. */
void mesh::draw_cylinders_pov(FILE *fp) {
    int i,*edp=edm;
    double *p=pts, *xy=xyz;
    for(i=0;i<n;i++,p+=6,xy+=3) {
        fprintf(fp,"\tsphere{<%g,%g,%g>,s}\n",*xy,xy[1],*p);
        for(;edp<ed[i+1];edp++) if(i<*edp) cyl_print(fp,p,pts+*edp*6);
    }
}

/** Prints a cylinder connecting two vertices to a file, skipping any with
 * duplicate start and end points, since they cause an error in POV-Ray.
 * \param[in] fp the file handle to write to.
 * \param[in] (p,q) pointers to the vertex pair to consider. */
void mesh::cyl_print(FILE *fp,double *p,double *q) {
    char vbuf1[80],vbuf2[80];
    sprintf(vbuf1,"%g,%g,%g",*p,p[1],p[2]);
    sprintf(vbuf2,"%g,%g,%g",*q,q[1],q[2]);
    if(strcmp(vbuf1,vbuf2)!=0)
        fprintf(fp,"\tcylinder{<%s>,<%s>,s}\n",vbuf1,vbuf2);
}

/** Computes an approximate normal vector at a mesh point.
 * \param[in] i the index of the mesh point.
 * \param[out] (x,y,z) the normal vector. */
void mesh::normal_vector(int i,double &x,double &y,double &z) {
    const double eps=std::numeric_limits<double>::epsilon();

    // Loop over all of the triangles and sum the normal vectors of each
    x=y=z=0;
    if(ed[i]==ed[i+1]) return;
    int *edp=ed[i],j=*edp,k=j;
    while(edp<ed[i+1]) {
        n_vec_contrib(x,y,z,i,j,*edp);
        j=*(edp++);
    }
    if((ncn[i]&bflag)==0) n_vec_contrib(x,y,z,i,j,k);

    // Normalize the vector
    double rsq=x*x+y*y+z*z;
    if(rsq>eps*eps) {
        rsq=1./sqrt(rsq);
        x*=rsq;y*=rsq;z*=rsq;
    }
}

/** Computes a contribution to the normal vector at a mesh point, by evaluating
 * a normal vector to a triangle.
 * \param[in,out] (x,y,z) the vector on which to add the normal vector contribution.
 * \param[in] (i,j,k) the triplets of mesh points that make up a triangle. */
inline void mesh::n_vec_contrib(double &x,double &y,double &z,int i,int j,int k) {
    double ux=xyz[3*j]-xyz[3*i],uy=xyz[3*j+1]-xyz[3*i+1],uz=pts[6*j]-pts[6*i];
    double vx=xyz[3*k]-xyz[3*i],vy=xyz[3*k+1]-xyz[3*i+1],vz=pts[6*k]-pts[6*i];
    x+=uy*vz-uz*vy;
    y+=uz*vx-ux*vz;
    z+=ux*vy-uy*vx;
}

/** Sets up a directory for saving simulation output.
 * \param[in] odir_ the name of the directory. */
void mesh::setup_output_dir(const char *odir_) {
    size_t l=strlen(odir_)+1;

    // Allocate space for the directory filename and copy it. Allocate additional
    // space for assembling the output filenames.
    odir=new char[2*l+32];
    memcpy(odir,odir_,sizeof(char)*l);
    obuf=odir+l;

    // Make the output directory if it doesn't already exist
    mkdir(odir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Save the mesh topology
    sprintf(obuf,"%s/topo",odir);
    output_topology(obuf);
}

/** Outputs the mesh topology.
 * \param[in] fp a file handle to write to. */
void mesh::output_topology(FILE *fp) {

	if(output_refined) {
		int nx=static_cast<int>(sqrt(n)); // TODO - nx and ny may be different.
		const int nn = nx + (np-1)*(nx-1);
		char buf[50], buf1[50];
		float s=sed/np;
		printf("The refined side length is set to %f.\n",s);
		sprintf(buf,"./sheet_gen rec48 %f %d %d",s,nn,nn);
		mesh_param par(K,drag,false,false,false);
		std::system(buf);
		sprintf(buf1,"sh48_%dx%d.bin",nn,nn);
		mesh *mp_deluxe=new mesh(par,buf1);

		fwrite(&(mp_deluxe->n),sizeof(int),2,fp);
		fwrite(mp_deluxe->ncn,sizeof(int),mp_deluxe->n,fp);
		fwrite(mp_deluxe->edm,sizeof(int),mp_deluxe->nc,fp);

		printf("Successfully saved the topology info for a refined REC48 mesh.\n");
	}
	else {
		// Save the total vertices and edges
		fwrite(&n,sizeof(int),2,fp);

		// Save the vertex connection numbers
		fwrite(ncn,sizeof(int),n,fp);

		// Save the edge table
		fwrite(edm,sizeof(int),nc,fp);
	}
}

/** Currently this works for a rec48 mesh. */
void mesh::mesh_print_dense_del(int fr,double t_,double *in) {
    printf("# Output frame %d (t=%g)\n",fr,t_);

    sprintf(obuf,"%s/pts.%d",odir,fr);
    FILE *fp=safe_fopen(obuf,"wb");

	const int n_del=(np-2)*(np-1)*ntri/2 + (np-1)*ns + n;
	double *dvals=new double[3*n_del];
	int nx=static_cast<int>(sqrt(n)); // TODO - nx and ny may be different.
	const int nn = nx + (np-1)*(nx-1);
	if(nn*nn!=n_del) {
		printf("Error: total computed deluxe points is incorrect for the rec48 mesh.\n");
		exit(1);
	}
	
	// Compute solution at deluxe points
	get_dvals(dvals,nn,n_del,nx);
 
	fwrite(dvals,sizeof(double),3*n_del,fp);
	delete[] dvals;
    fclose(fp);
}

void mesh::mesh_print_dense(int fr,double t_,double *in) {
    printf("# Output frame %d (t=%g)\n",fr,t_);
    sprintf(obuf,"%s/pts.%d",odir,fr);
    FILE *fp=safe_fopen(obuf,"wb");
	
	double *Apts = new double[3*n]; 
	for (int i=0,j=0;j<6*n;i+=3,j+=6) {
		Apts[i] = xyz[i];
		Apts[i+1] = xyz[i+1];
		Apts[i+2] = in[j];
	}
	fwrite(Apts,sizeof(double),3*n,fp);
	delete[] Apts;
    fclose(fp);
}

void mesh::mesh_print_last_step() {
	if(wr_all_dofs) {
		sprintf(obuf,"%s/pts_Argyris.%d",odir,1);
		FILE *fp=safe_fopen(obuf,"wb");
		double *Apts=new double[8*n+ns];
		for(int i=0;i<n;i++) {
			Apts[8*i]=xyz[3*i]; Apts[8*i+1]=xyz[3*i+1];
			for(int k=0;k<6;k++) Apts[8*i+2+k]=pts[6*i+k];
		}
		for(int i=0;i<ns;i++) Apts[8*n+i]=pts[6*n+i];
		fwrite(Apts,sizeof(double),8*n+ns,fp);
		delete[] Apts;
		fclose(fp);
	}
	else {
		sprintf(obuf,"%s/pts.%d",odir,1);
		FILE *fp=safe_fopen(obuf,"wb");
		double *Apts = new double[3*n]; 
		for (int i=0,j=0;j<6*n;i+=3,j+=6) {
			Apts[i] = xyz[i];
			Apts[i+1] = xyz[i+1];
			Apts[i+2] = pts[j];
		}
		fwrite(Apts,sizeof(double),3*n,fp);
		delete[] Apts;
		fclose(fp);
	}
}

/** Outputs the mesh vertex positions.
 * \param[in] fp a file handle to write to. */
void mesh::output_positions(FILE *fp) {
    fwrite(pts,sizeof(double),Adof2,fp);
}

/** Outputs the mesh vertex positions into the output directory with a given
 * frame index.
 * \param[in] l the frame index. */
void mesh::output_positions(int l) {

    // Check that the output directory has been set up
    if(odir==NULL) fatal_error("output directory not set up",1);

    // Assemble the output filename and save the data
    sprintf(obuf,"%s/pts.%d",odir,l);
    output_positions(obuf);
}

/** Outputs the mesh vertex velocities.
 * \param[in] fp a file handle to write to. */
void mesh::output_velocities(FILE *fp) {
    fwrite(vel,sizeof(double),Adof2,fp);
}

/** Outputs the mesh vertex velocities into the output directory with a given
 * frame index.
 * \param[in] l the frame index. */
void mesh::output_velocities(int l) {

    // Check that the output directory has been set up
    if(odir==NULL) fatal_error("output directory not set up",1);

    // Assemble the output filename and save the data
    sprintf(obuf,"%s/vel.%d",odir,l);
    output_velocities(obuf);
}

/** Prints out the edge table for diagnostic purposes.
 * \param[in] fp a file handle to write to. */
void mesh::edge_diagnostic(FILE *fp) {
    int *edp=edm;
    for(int i=0;i<n;i++) {
        fprintf(fp,"%d -",i);
        while(edp<ed[i+1]) {
            fprintf(fp," %d",*edp);edp++;
        }
        fputc('\n',fp);
    }
}
