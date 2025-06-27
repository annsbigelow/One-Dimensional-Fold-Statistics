#include "ext_potential.hh"

#include <cstdio>
#include <cstdlib>
#include <cmath>

/** Adds accelerations to the mesh points due to a centering potential.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] is the starting mesh point.
 * \param[in] ie the ending mesh point (excluded).
 * \param[in] in the mesh point positions.
 * \param[in] acc the mesh point accelerations (cumulative). */
void ep_centering::accel(double t,int n,int is,int ie,double *in,double *acc) {
    double fac,dx=0,dy=0,dz=0,*ap,*pp;

    // Compute the centroid
    dx=dy=dz=0;
    for(pp=in;pp<in+3*n;) {
        dx+=*(pp++);
        dy+=*(pp++);
        dz+=*(pp++);
    }
    fac=cforce/static_cast<double>(n);dx*=fac;dy*=fac;dz*=fac;

    // Use the centroid to create a drift toward the origin
    for(ap=acc+3*is;ap<acc+3*ie;ap+=3) {
        *ap-=dx;
        ap[1]-=dy;
        ap[2]-=dz;
    }
}

/** Computes the potential energy due to the centering force.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] in the mesh point positions. */
double ep_centering::energy(double t,int n,double *in) {
    double dx=0,dy=0,dz=0,*pp;

    // Compute the centroid
    dx=dy=dz=0;
    for(pp=in;pp<in+3*n;) {
        dx+=*(pp++);
        dy+=*(pp++);
        dz+=*(pp++);
    }

    // Use the centroid to create a drift toward the origin
    return 0.5*cforce/static_cast<double>(n)*(dx*dx+dy*dy+dz*dz);
}

/** Adds accelerations to the mesh points due to the quadratic potential.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] is the starting mesh point.
 * \param[in] ie the ending mesh point (excluded).
 * \param[in] in the mesh point positions.
 * \param[in] acc the mesh point accelerations (cumulative). */
void ep_quadratic::accel(double t,int n,int is,int ie,double *in,double *acc) {
    t*=-lambda;
    double *pp,*ap;
    for(pp=in+3*is,ap=acc+3*is;ap<acc+3*ie;pp+=3,ap+=3) {
        *ap+=*(pp)*t;
        ap[1]+=pp[1]*t;
        ap[2]+=pp[2]*t;
    }
}

/** Computes the energy due to the piston potential.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] in the mesh point positions. */
double ep_quadratic::energy(double t,int n,double *in) {
    double en=0.;
    for(double *pp=in;pp<in+3*n;pp++) en+=(*pp)*(*pp);
    return 0.5*lambda*t*en;
}

/** Initialize potential parameters.
 * \param[in] r_cut_ cutoff radius for spherical potential.
 * \param[in] C_ strength of potential. */
void ep_spherical::set_params(double r_cut_,double C_) {
    r_cut=r_cut_;
    r_cut2=r_cut*r_cut;
    C=C_;
}

/** Adds accelerations to spherical potential.
 * \param[in] t the current time.
 * \param[in] is the starting mesh point.
 * \param[in] ie the ending mesh point (excluded).
 * \param[in] in the mesh point positions.
 * \param[in] acc the mesh point accelerations (cumulative). */
void ep_spherical::accel(double t,int n,int is,int ie,double *in,double *acc) { // In parallel
    double rsq,fac,*pp,*ap;
    for(pp=in+3*is,ap=acc+3*is;ap<acc+3*ie;pp+=3,ap+=3) {
        rsq=(*pp)*(*pp)+pp[1]*pp[1]+pp[2]*pp[2];
        if(rsq>r_cut2) {
            fac=C*(1-r_cut/sqrt(rsq));
            *ap-=fac*(*pp); ap[1]-=fac*pp[1]; ap[2]-=fac*pp[2];
        }
    }
}

/** Computes the energy due to the spherical potential.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] in the mesh point positions. */
double ep_spherical::energy(double t,int n,double *in) {
    double en=0.,fac,rsq;
    for(double *pp=in;pp<in+3*n;pp+=3) {
        rsq=(*pp)*(*pp)+pp[1]*pp[1]+pp[2]*pp[2];
        if(rsq>r_cut2) {
            fac=sqrt(rsq)-r_cut;
            en+=fac*fac;
        }
    }
    return 0.5*C*en;
}

/** Adds accelerations to the mesh points due to the piston potential.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] is the starting mesh point.
 * \param[in] ie the ending mesh point (excluded).
 * \param[in] in the mesh point positions.
 * \param[in] acc the mesh point accelerations (cumulative). */
void ep_piston::accel(double t,int n,int is,int ie,double *in,double *acc) {
    double r,rsq,expr,*pp,*ap;

    for(pp=in+3*is,ap=acc+3*is;ap<acc+3*ie;pp+=3,ap+=3) {
        rsq=(*pp)*(*pp)+pp[1]*pp[1];

        // Constrain at the walls
        if(rsq>r_piston2) {
            r=sqrt(rsq);
            expr=exp(r-r_piston);
            *ap-=expr*(*pp)/r;
            ap[1]-=expr*pp[1]/r;
        }

        // Constrain at the bottom boundary
        if(pp[2]<z_bot) ap[2]+=exp(z_bot-pp[2]);

        // Constrain at the top boundary
        if(pp[2]>(z_top-delta)) ap[2]-=(((sin(freq*t+20.3)-0.998)<0)?0.05:0);
    }
}

/** Computes the energy due to the piston potential.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] in the mesh point positions. */
double ep_piston::energy(double t,int n,double *in) {
    fputs("Piston energy function not implemented yet\n",stderr);
    exit(1);
}

/** Initializes the piston parameters. */
void ep_piston::set_params(double r_piston_,double z_piston_,double delta_,double freq_) {
    r_piston=r_piston_;
    r_piston2=r_piston*r_piston;
    z_top=z_piston_;
    z_bot=-z_piston_;
    delta=delta_;
    freq=freq_;
}

/** Adds accelerations to the mesh points due to axial compression.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] is the starting mesh point.
 * \param[in] ie the ending mesh point (excluded).
 * \param[in] in the mesh point positions.
 * \param[in] acc the mesh point accelerations (cumulative). */
void ep_axial::accel(double t,int n,int is,int ie,double *in,double *acc) {
    double r,rsq,expr,pz,*pp,*ap;

    for(pp=in+3*is,ap=acc+3*is;ap<acc+3*ie;pp+=3,ap+=3) {
        rsq=(*pp)*(*pp)+pp[1]*pp[1];

        // Constrain at the walls
        if(rsq>r_piston2) {
            r=sqrt(rsq);
            expr=exp(r-r_piston);
            *ap-=expr*(*pp)/r;
            ap[1]-=expr*pp[1]/r;
        }

        // Constrain at the bottom boundary
        if(pp[2]<z_bot) ap[2]+=az;//ap[2]+=exp(z_bot-pp[2]);

        // Constrain at the top boundary
        pz=z_top-vz*t+0.5*az*t*t;
        if(pp[2]>pz) ap[2]-=az;
    }
}

/** Computes the energy due to axial compression.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] in the mesh point positions. */
double ep_axial::energy(double t,int n,double *in) {
    fputs("Axial compression energy function not implemented yet\n",stderr);
    exit(1);
}

/** Initializes the axial compression parameters. */
void ep_axial::set_params(double r_piston_,double z_piston_,double tf_) {
    r_piston=r_piston_;
    r_piston2=r_piston*r_piston;
    z_top=z_piston_;
    z_bot=-z_piston_;
    tf=tf_;
    vz=2*z_top/tf;
    az=vz/tf;
}

/** Initialize potential parameters.
 * \param[in] r_cut_ cutoff radius for spherical potential with drift in -z.
 * \param[in] C_ strength of potential.
 * \param[in] vz_ magnitude of initial velocity of sphere center, moving in -z direction.
 * \param[in] tf_ time for sphere to displace distance dz, denting the sheet. */
void ep_spherical_drift::set_params(double r_cut_,double C_,double dz_,double tf_) {
    r_cut=r_cut_;
    r_cut2=r_cut*r_cut;
    C=C_;
    dz=dz_;
    tf=tf_;
    vz=2*dz/tf;
    az=vz/tf;
}

/** Adds accelerations to spherical potential with drift in -z.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] is the starting mesh point.
 * \param[in] ie the ending mesh point (excluded).
 * \param[in] in the mesh point positions.
 * \param[in] acc the mesh point accelerations (cumulative). */
void ep_spherical_drift::accel(double t,int n,int is,int ie,double *in,double *acc) {
    double rsq,fac,*pp,*ap;
    for(pp=in+3*is,ap=acc+3*is;ap<acc+3*ie;pp+=3,ap+=3) {
        rsq=(*pp)*(*pp)+pp[1]*pp[1]+(pp[2]-(r_cut-vz*t+0.5*az*t*t))*(pp[2]-(r_cut-vz*t+0.5*az*t*t));
        if(rsq<r_cut2) {
            fac=-C*(1-r_cut/sqrt(rsq));
            *ap+=fac*(*pp); ap[1]+=fac*pp[1]; ap[2]+=(fac*pp[2]-az);
        }
    }
}

/** Computes the energy due to the spherical potential with drift in -z.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] in the mesh point positions. */
double ep_spherical_drift::energy(double t,int n,double *in) {
    fputs("Energy function not implemented yet\n",stderr);
    exit(1);
}

/** Adds tensile force in x direction.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] is the starting mesh point.
 * \param[in] ie the ending mesh point (excluded).
 * \param[in] in the mesh point positions.
 * \param[in] acc the mesh point accelerations (cumulative). */
void ep_tensilex::accel(double t,int n,int is,int ie,double *in,double *acc) {
    double *vp,*ap;
    int *gp;
    for(vp=in+3*n+3*is,ap=acc+3*is,gp=grip+is;ap<acc+3*ie;vp+=3,ap+=3,gp++) if(*gp!=0) {
    *ap=ap[1]=ap[2]=0;
    *vp=gamma*(*gp);
    vp[1]=vp[2]=0;
    }
}

/** Computes the energy due to the tensile force in x direction.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] in the mesh point positions. */
double ep_tensilex::energy(double t,int n,double *in) {
    fputs("Energy function not implemented yet\n",stderr);
    exit(1);
}

/** Adds tensile force in y direction.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] is the starting mesh point.
 * \param[in] ie the ending mesh point (excluded).
 * \param[in] in the mesh point positions.
 * \param[in] acc the mesh point accelerations (cumulative). */
void ep_tensiley::accel(double t,int n,int is,int ie,double *in,double *acc) {
    double *vp,*ap;
    int *gp;
    for(vp=in+3*n+3*is,ap=acc+3*is,gp=grip+is;ap<acc+3*ie;vp+=3,ap+=3,gp++) if(*gp!=0) {
        *ap=ap[1]=ap[2]=0;
        vp[1]=gamma*(*gp);
        *vp=vp[2]=0;
    }
}

/** Computes the energy due to the tensile force in y direction.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] in the mesh point positions. */
double ep_tensiley::energy(double t,int n,double *in) {
    fputs("Energy function not implemented yet\n",stderr);
    exit(1);
}

/** Initialize potential parameters.
 * \param[in] R_ the initial shell radius.
 * \param[in] t_ initial time. */
void ep_shell::set_params(double r_,double t_) {
    r0=r_;
    t0=t_;
}

/** Adds accelerations to shell potential.
 * \param[in] t the current time.
 * \param[in] is the starting mesh point.
 * \param[in] ie the ending mesh point (excluded).
 * \param[in] in the mesh point positions.
 * \param[in] acc the mesh point accelerations (cumulative). */
void ep_shell::accel(double t,int n,int is,int ie,double *in,double *acc) { // In parallel
    double rsq,fac,*pp,*ap,r=r0-gamma*(t-t0);
    for(pp=in+3*is,ap=acc+3*is;ap<acc+3*ie;pp+=3,ap+=3) {
        rsq=(*pp)*(*pp)+pp[1]*pp[1]+pp[2]*pp[2];
        if(rsq>r*r) {
            fac=C*(1-r/sqrt(rsq));
            *ap-=fac*(*pp); ap[1]-=fac*pp[1]; ap[2]-=fac*pp[2];
        }
    }
}

/** Computes the energy due to the shell potential.
 * \param[in] t the current time.
 * \param[in] n the number of mesh points.
 * \param[in] in the mesh point positions. */
double ep_shell::energy(double t,int n,double *in) {
    double en=0.,fac,rsq,r=r0-gamma*(t-t0);
    for(double *pp=in;pp<in+3*n;pp+=3) {
        rsq=(*pp)*(*pp)+pp[1]*pp[1]+pp[2]*pp[2];
        if(rsq>r*r) {
            fac=sqrt(rsq)-r;
            en+=fac*fac;
        }
    }
    return 0.5*C*en;
}
