#include "mesh.hh"
#include "vec3.hh"

/** The constructor reads in a mesh from a file, and sets up the vertices and
 * edge tables.
 * \param[in] mp a mesh_param structure containing simulation constants.
 * \param[in] filename the file to read from. */
mesh::mesh(mesh_param &mp,const char* filename) : mesh_param(mp),
    n_ep(0), reg(NULL), odir(NULL) {
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
    mesh_param(mp), n_ep(0), reg(NULL), odir(NULL) {

    // Read in the mesh topology
    FILE *fp=safe_fopen(f_topo,"rb");
    read_topology(fp);
    fclose(fp);

    // Read in the vertex positions
    fp=safe_fopen(f_pts,"rb");
    read_positions(fp);
    fclose(fp);
}

/** The class destructor frees the dynamically allocated memory. */
mesh::~mesh() {
    if(odir!=NULL) delete [] odir;
    if(reg!=NULL) {
    if(bsheet_model) {delete [] ref; delete [] tom; delete [] to;}
    delete [] reg; delete [] eom; delete [] eo;
    }
    delete [] edm;delete [] ed;
    delete [] ncn;
    delete [] pts;
}

/** Sets up the spring network table and initializes the spring rest lengths to
 * be fully relaxed. */
void mesh::setup_springs() {

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

    // Set up relaxed connection lengths
    reg=new double[ns];
    double *regp=reg,dx,dy,dz,emax=0;
    for(i=0,edp=eom;i<n;i++) {
        while(edp<eo[i+1]) {
            dx=pts[3*i]-pts[3*(*edp)];
            dy=pts[3*i+1]-pts[3*(*edp)+1];
            dz=pts[3*i+2]-pts[3*(*(edp++))+2];
            *regp=sqrt(dx*dx+dy*dy+dz*dz);
        if(*regp>emax) emax=*regp;
            regp++;
        }
    }
    sigma=emax;

    if(bsheet_model) {

        // Set up relaxed areas of triangles joined by a "hinge" edge.
        // Determine number of "hinges" by subtracting boundary edges from total.
        int nb=0,j,j2,k,ct=0;
        for(i=0;i<n;i++) {
            if(ncn[i]&bflag) nb++;
        }
        nh=ns-nb;
        to=new int*[n+1];
        tom=new int[3*nh];
        ref=new double[nh];

        // Set up relaxed triangle areas as sums of triangle pairs
        int *top=tom;
        regp=ref; edp=*ed;
        for(i=0;i<n;i++) {
        to[i]=top;
        if(edp!=ed[i+1]) {

                // Cycle around the edges. Remember the first connected vertex.
                j=*(edp++);
                if(edp<ed[i+1]) {
                    j2=*edp;

                    // Loop over the other vertices
                    while(edp+1<ed[i+1]) {
                        k=*(edp++);
                        if(i<k) {
                            *(regp++)=edge_factor(pts,i,edp[-2],k,*edp);
                            *(top++)=edp[-2];
                            *(top++)=k;
                            *(top++)=*edp;
                            ct++;
                        }
                    }
                    k=*edp;

                    // Add triangles to make a complete loop, if this is not a boundary
                    // vertex
                    if((ncn[i]&bflag)==0) {
                        if(i<k) {
                            *(regp++)=edge_factor(pts,i,edp[-1],k,j);
                            *(top++)=edp[-1];
                            *(top++)=k;
                            *(top++)=j;
                            ct++;
                        }
                        if(i<j) {
                            *(regp++)=edge_factor(pts,i,k,j,j2);
                            *(top++)=k;
                            *(top++)=j;
                            *(top++)=j2;
                            ct++;
                        }
                    }
                    edp++;
                }
        }
        }
        to[n]=top;
        printf("%d %d\n",nh,ct);
    }
}

/** Perturbs the rest lengths of the springs by uniform random numbers
 * \param[in] (min_fac,max_fac) the range that the uniform numbers are sampled
 *                              from. */
void mesh::perturb_springs(double min_fac,double max_fac) {
    double *rp=reg,rfac=(max_fac-min_fac)/RAND_MAX,emax=0;
    int i,*ep=eom;
    unsigned int *np=ncn,*np2;

    // Loop over all of the springs
    for(i=0;i<n;i++,np++) while(ep<eo[i+1]) {
        np2=ncn+(*(ep++));

        // Randomly perturb all edges except those between boundary nodes, if boundary fixed
        if(fix_boundary&&(*np&bflag)&&(*np2&bflag)) rp++;
        else {
            *rp*=min_fac+static_cast<double>(rand())*rfac;
            if(*rp>emax) emax=*rp;
            rp++;
        }
    }
    sigma=emax;
}

/** Resets the current mesh configuration to be the new relaxed configuration. */
void mesh::reset_relaxed() {
    double *rp=reg,dx,dy,dz,emax=0; int i,*ep=eom,*tp=tom;
    for(i=0;i<n;i++) while(ep<eo[i+1]) {
        dx=pts[3*i]-pts[3*(*ep)];
        dy=pts[3*i+1]-pts[3*(*ep)+1];
        dz=pts[3*i+2]-pts[3*(*(ep++))+2];
        *rp=sqrt(dx*dx+dy*dy+dz*dz);
        if(*rp>emax) emax=*rp;
        rp++;
    }
    sigma=emax;

    if(bsheet_model) {
        rp=ref;
        for(i=0;i<n;i++) while(tp<to[i+1]) {
            *(rp++)=edge_factor(pts,i,*tp,tp[1],tp[2]);
            tp+=3;
        }
    }
}

void mesh::mesh_ff(double t_,double *in,double *out) {
    double *acc=out+3*n;
    int i;

    // Compute drag
    for(double *ap=acc,*vp=in+3*n;ap<acc+3*n;) *(ap++)=-*(vp++)*drag;

    // Add accelerations due to springs and external potentials
    acceleration(t_,in,acc);

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

/** Computes the acceleration due to springs and external potentials.
 * \param[in] t_ the time at which to evaluate the acceleration.
 * \param[in] in the mesh point positions.
 * \param[in] out the mesh point accelerations (cumulative). */
void mesh::acceleration(double t_,double *in,double *acc) {

    // Add repulsive forces if present
    //if(repulsion) accel_repulsive(in,acc);

    // Compute accelerations from the sheet
    bsheet_model?accel_rbsheet(in,acc):accel_springs(in,acc);

    // Add accelerations from the external potentials
    for(int i=0;i<n_ep;i++) ex_pot[i]->accel(t_,n,in,acc);
}

/** Computes the energy.
 * \param[in] in the mesh point positions. */
double mesh::energy(double t_,double *in) {

    // Compute energies due to edge springs
    double en=bsheet_model?energy_bsheet(in):energy_springs(in);

    // Add energies due to external potentials
    for(int i=0;i<n_ep;i++) en+=ex_pot[i]->energy(t_,n,in);
    return en;
}

/** Adds the edge spring forces for the shapeable sheet model.
 * \param[in] in the mesh point positions.
 * \param[in] acc the mesh point accelerations (cumulative). */
void mesh::accel_springs(double *in,double *acc) {

    int i,*ep=eo[0],of=ep-eom;
    double *rp=reg+of;

    for(i=0;i<n;i++) while(ep<eo[i+1]) {
        if(dashpot) damp_force(in,acc,i,*ep);
        stretch_force(in,acc,i,*(ep++),*(rp++));
    }
}

/** Computes the accelerations based on the bending sheet model with random mesh.
 * \param[in] in the mesh point positions.
 * \param[in] acc the mesh point accelerations (cumulative). */
void mesh::accel_rbsheet(double *in,double *acc) {
    int i,*ep=eo[0],*tp=to[0],of=ep-eom;
    double *rp=reg+of;

    // first compute the edge forces.
    for(i=0;i<n;i++) while(ep<eo[i+1]) {
        if(dashpot) damp_force(in,acc,i,*ep);
        stretch_force(in,acc,i,*(ep++),*(rp++));
    }

    // next add in bending forces.
    of=(tp-tom)/3; rp=ref+of;
    for(i=0;i<n;i++) while(tp<to[i+1]) {
        bend_force(in,acc,i,*tp,tp[1],tp[2],*(rp++));
        tp+=3;
    }
}

/** Computes the energy from the spring forces in the shapeable sheet
 * model.
 * \param[in] in the mesh point positions.
 * \return The energy. */
double mesh::energy_springs(double *in) {
    double en=0,*pp,*pp2,*rp=reg,rs,dx,dy,dz;int i,*ep=eom;
    for(pp=in,i=0;i<n;i++,pp+=3) while(ep<eo[i+1]) {

        // Find the vector to the neighboring vertex
        pp2=in+3*(*ep);
        dx=*pp-*(pp2++);
        dy=pp[1]-*(pp2++);
        dz=pp[2]-*pp2;

        // Compute the force
        rs=*(rp++)-sqrt(dx*dx+dy*dy+dz*dz);
        en+=rs*rs;
    }

    // Return the value, scaled by the spring constant
    return 0.5*K*en;
}

/** Computes the bandwidth of the banded hessian matrix. */
int mesh::bandwidth() {
    int i,e,*ep=eom,band=0;
    for(i=0;i<n;i++) {
        while(ep<eo[i+1]) {
            e=*ep;
            if((e-i)>band) band=e-i;
            ep++;
        }
    }
    band=3*band+2;
    return band;
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

/** Computes the accelerations based on the bending sheet model.
 * \param[in] in the mesh point positions.
 * \param[in] acc the mesh point accelerations (cumulative). */
void mesh::accel_bsheet(double *in,double *acc) {
    int i,*edp=*ed,j,j2,k;

    for(i=0;i<n;i++) if(edp!=ed[i+1]) {

        // Cycle around the edges. Remember the first connected vertex.
        j=*(edp++);
        if(i<j) edge_force(in,acc,i,j);

        if(edp<ed[i+1]) {
            j2=*edp;

            // Loop over the other vertices
            while(edp+1<ed[i+1]) {
                k=*(edp++);
                if(i<k) {
                    edge_force(in,acc,i,k);
                    triangle_force(in,acc,i,edp[-2],k,*edp);
                }
            }
            k=*edp;
            if(i<k) edge_force(in,acc,i,k);

            // Add triangles to make a complete loop, if this is not a boundary
            // vertex
            if((ncn[i]&bflag)==0) {
                if(i<k) triangle_force(in,acc,i,edp[-1],k,j);
                if(i<j) triangle_force(in,acc,i,k,j,j2);
            }
            edp++;
        }
    }
}

/** Computes the accelerations based on the bending sheet model.
 * \param[in] in the mesh point positions.
 * \param[in] acc the mesh point accelerations (cumulative). */
double mesh::energy_bsheet(double *in) {
    int i,*edp=*ed,j,j2,k;double ee=0.,et=0.;

    for(i=0;i<n;i++) if(edp!=ed[i+1]) {

        // Cycle around the edges. Remember the first connected vertex.
        j=*(edp++);
        if(i<j) ee+=edge_energy(in,i,j);

        if(edp<ed[i+1]) {
            j2=*edp;

            // Loop over the other vertices
            while(edp+1<ed[i+1]) {
                k=*(edp++);
                if(i<k) {
                    ee+=edge_energy(in,i,k);
                    et+=triangle_energy(in,i,edp[-2],k,*edp);
                }
            }
            k=*edp;
            if(i<k) ee+=edge_energy(in,i,k);

            // Add triangles to make a complete loop, if this is not a boundary
            // vertex
            if((ncn[i]&bflag)==0) {
                if(i<k) et+=triangle_energy(in,i,edp[-1],k,j);
                if(i<j) et+=triangle_energy(in,i,k,j,j2);
            }
            edp++;
        }
    }

    // Scale the accumulators by the physical constants
    return kappa*et+0.5*K*ee;
}

/** Adds the force due to an edge in the bendable sheet model.
 * \param[in] in the mesh point positions.
 * \param[in] acc the mesh point accelerations.
 * \param[in] (i,k) the two vertices that the edge connects. */
void mesh::edge_force(double *in,double *acc,int i,int k) {
    double *ii=in+3*i,*ik=in+3*k,*ai=acc+3*i,*ak=acc+3*k,
           dx=*ii-*ik,dy=ii[1]-ik[1],dz=ii[2]-ik[2],
           rs=1./sqrt(dx*dx+dy*dy+dz*dz)-1;

    // Add the force contributions to the two vertices
    dx*=rs*K;dy*=rs*K;dz*=rs*K;
    *ai+=dx;ai[1]+=dy;ai[2]+=dz;
    *ak-=dx;ak[1]-=dy;ak[2]-=dz;
}

/** Calculates the acceleration due to the angle between two adjacent triangles
 * (i,j,k) and (i,j,l) in the bendable sheet model. \param[in] in the mesh
 * point positions.
 * \param[in] (i,j,k,l) the four vertices defining the two triangles.
 * \return The energy. */
void mesh::triangle_force(double *in,double *acc,int i,int j,int k,int l) {
    double *ii=in+3*i,*ij=in+3*j,*ik=in+3*k,*il=in+3*l,ne,nf,sp,fac;
    vec3 b(*ij-*ii,ij[1]-ii[1],ij[2]-ii[2]),
         c(*ik-*ii,ik[1]-ii[1],ik[2]-ii[2]),
         d(*il-*ii,il[1]-ii[1],il[2]-ii[2]),e,f,u,v,a0,a1,a2,a3;

    // Compute normals to the triangles
    e=b*c;
    f=c*d;
    ne=1./mod_sq(e);
    nf=1./mod_sq(f);

    // Compute forces on the four vertices
    fac=kappa*sqrt(ne*nf);sp=dot(e,f);
    u=fac*(f-ne*sp*e);
    v=fac*(nf*sp*f-e);
    a1=c*u;
    a2=-b*u-d*v;
    a3=c*v;

    // Add forces to the acc array. Rather than compute the force on vertex i,
    // obtain it by subtracting the other three.
    a0=-a1-a2-a3;
    a0.add(acc+3*i);
    a1.add(acc+3*j);
    a2.add(acc+3*k);
    a3.add(acc+3*l);
}

/** Calculates the (unnormalized) energy due to an edge in the bendable sheet
 * model.
 * \param[in] in the mesh point positions.
 * \param[in] (i,k) the two vertices that the edge connects.
 * \return The energy. */
double mesh::edge_energy(double *in,int i,int k) {
    double *ii=in+3*i,*ik=in+3*k,
           dx=*ii-*ik,dy=ii[1]-ik[1],dz=ii[2]-ik[2],
           val=sqrt(dx*dx+dy*dy+dz*dz)-1;
    return val*val;
}

/** Calculates the (unnormalized) energy due to the angle between two adjacent
 * triangles (i,j,k) and (i,j,l) in the bendable sheet model.
 * \param[in] in the mesh point positions.
 * \param[in] (i,j,k,l) the four vertices defining the two triangles.
 * \return The energy. */
double mesh::triangle_energy(double *in,int i,int j,int k,int l) {
    double *ii=in+3*i,*ij=in+3*j,*ik=in+3*k,*il=in+3*l;
    vec3 b(*ij-*ii,ij[1]-ii[1],ij[2]-ii[2]),
         c(*ik-*ii,ik[1]-ii[1],ik[2]-ii[2]),
         d(*il-*ii,il[1]-ii[1],il[2]-ii[2]),e,f;

    // Compute normals to the triangles
    e=b*c;
    f=c*d;

    // Compute the energy contribution
    return 1-dot(e,f)/sqrt(mod_sq(e)*mod_sq(f));
}

double mesh::edge_factor(double *in,int i,int j,int k,int l) {
    double *ii=in+3*i,*ij=in+3*j,*ik=in+3*k,*il=in+3*l;
    vec3 e1(*ij-*ii,ij[1]-ii[1],ij[2]-ii[2]),
         e0(*ik-*ii,ik[1]-ii[1],ik[2]-ii[2]),
         e2(*il-*ii,il[1]-ii[1],il[2]-ii[2]),n1,n2;

    // Compute normals to the triangles
    n1=e0*e1;
    n2=-e0*e2;

    // Return the ratio of the hinge length squared to sum of adjacent triangle areas.
    return 2*mod_sq(e0)/(sqrt(mod_sq(n1))+sqrt(mod_sq(n2)));
}

void mesh::stretch_force(double *in,double *acc,int i,int k,double sf) {
    double *ii=in+3*i,*ik=in+3*k,*ai=acc+3*i,*ak=acc+3*k,
           dx=*ii-*ik,dy=ii[1]-ik[1],dz=ii[2]-ik[2],
           rs=sf/sqrt(dx*dx+dy*dy+dz*dz)-1;

    // Add the force contributions to the two vertices
    dx*=rs*K;dy*=rs*K;dz*=rs*K;
    *ai+=dx;ai[1]+=dy;ai[2]+=dz;
    *ak-=dx;ak[1]-=dy;ak[2]-=dz;
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

/** Calculates the acceleration due to the angle between two adjacent triangles
 * (i,j,k) and (i,j,l) in the bendable sheet model with random mesh.
 * \param[in] in the mesh point positions.
 * \param[in] (i,j,k,l) the four vertices defining the two triangles.
 * \return The energy. */
void mesh::bend_force(double *in,double *acc,int i,int j,int k,int l,double ef) {
    double *ii=in+3*i,*ij=in+3*j,*ik=in+3*k,*il=in+3*l,ne,nf,sp,fac;
    vec3 b(*ij-*ii,ij[1]-ii[1],ij[2]-ii[2]),
         c(*ik-*ii,ik[1]-ii[1],ik[2]-ii[2]),
         d(*il-*ii,il[1]-ii[1],il[2]-ii[2]),e,f,u,v,a0,a1,a2,a3;

    // Compute normals to the triangles
    e=b*c;
    f=c*d;
    ne=1./mod_sq(e);
    nf=1./mod_sq(f);

    // Compute forces on the four vertices
    fac=sqrt(3)/2*kappa*ef*sqrt(ne*nf);sp=dot(e,f);
    u=fac*(f-ne*sp*e);
    v=fac*(nf*sp*f-e);
    a1=c*u;
    a2=-b*u-d*v;
    a3=c*v;

    // Add forces to the acc array. Rather than compute the force on vertex i,
    // obtain it by subtracting the other three.
    a0=-a1-a2-a3;
    a0.add(acc+3*i);
    a1.add(acc+3*j);
    a2.add(acc+3*k);
    a3.add(acc+3*l);
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
