#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "p_grid.hh"

/** Initializes the proximity grid class for rapidly detecting points that are
 * spatially near each other. */
proximity_grid::proximity_grid() : mno(0), count(0),
    mem_ov(pgrid_init_overflow_buffer), ov(new oflow_info[mem_ov]) {}

/** The class destructor frees any dynamically allocated memory. */
proximity_grid::~proximity_grid() {
    if(mno>0) free_grid();
    delete [] ov;
}

/** Sets up the grid data structure. It first periodically rebuilds the grid
 * dimensions. Then is populates the grid with the point positions.
 * \param[in] q a pointer to the point (x,y,z) positions.
 * \param[in] c the number of points to consider. */
void proximity_grid::setup(double *q,int c) {
    if(mno>0) {
        if(++count==pgrid_build_interval) {allocate(q,c);count=0;}
    } else allocate(q,c);
}

/** Populates the grid data structure with point positions.
 * \param[in] q a pointer to the point (x,y,z) positions.
 * \param[in] c the number of points to consider. */
void proximity_grid::populate(double *q,int c) {
    clear();
    co_ov=0;
#pragma omp parallel for
    for(int i=0;i<c;i++) put(i,q+3*i);
    resolve_overflows();
}

/** Puts a point into the grid structure.
 * \param[in] id the numerical ID of this point.
 * \param[in] qp a pointer to the point's (x,y,z) position. */
void proximity_grid::put(int id,double *qp) {
    // Find the block that this point is inside
    int i=(*qp-ax)*xsp,j=(qp[1]-ay)*ysp,k=(qp[2]-az)*zsp;
    map_inside(i,j,k);
    int ijk=i+m*(j+n*k),s;
#pragma omp atomic capture
    s=co[ijk]++;

    // If there is an available memory slot for this point, then just add it
    if(s<mem[ijk]) p[ijk][s].set(id,qp);
    else {

        // If there isn't an available memory slot, then add it to the overflow
        // buffer. This will happen infrequently, but needs to be done in a
        // critical block to prevent race conditions between threads.
#pragma omp critical
        {
            if(co_ov==mem_ov) add_overflow_memory();
            ov[co_ov++].set(ijk,s,id,qp);
	}
    }
}

/** Thread-specific version of clear.
 * \param[in] nt the total number of threads.
 * \param[in] j the current thread index. */
void proximity_grid::clear(int nt,int j) {
    if(j==0) co_ov=0;
    int ld=mno/nt,is=j*ld,ie=(j==nt-1)?mno:(j+1)*ld;
    for(int *cop=co+is;cop<co+ie;cop++) *cop=0;
}

/** Moves any points in the overflow buffer into the main grid structure. */
void proximity_grid::resolve_overflows() {
    for(oflow_info *op=ov,*oe=ov+co_ov;op<oe;op++) {
        while(co[op->ijk]>=mem[op->ijk]) add_memory(op->ijk);
        op->put(p[op->ijk]+op->s);
    }
}

/** Allocates memory for the grid structure.
 * \param[in] q a pointer to the point (x,y,z) positions that will be added to
 *              the structure.
 * \param[in] c the number of points to consider. */
void proximity_grid::allocate(double *q,int c,int imem,double opt_points) {
    const double apad=1e-3,rpad=0.05;

    // Find the extremal (x,y,z) values of the points
    double lx=*q,ux=*q,ly=q[1],uy=q[1],lz=q[2],uz=q[2];
    for(double *qp=q+3,*qe=q+3*c;qp<qe;qp+=3) {
        if(*qp<lx) lx=*qp;else if(*qp>ux) ux=*qp;
        if(qp[1]<ly) ly=qp[1];else if(qp[1]>uy) uy=qp[1];
        if(qp[2]<lz) lz=qp[2];else if(qp[2]>uz) uz=qp[2];
    }

    // Compute the dimensions of the grid data structure by padding
    // the extremal (x,y,z) ranges of points.
    double wx=ux-lx,wy=uy-ly,wz=uz-lz,
           cx=0.5*(ux+lx),cy=0.5*(uy+ly),cz=0.5*(uz+lz);
    wx+=wx*rpad;wx+=apad;
    wy+=wy*rpad;wy+=apad;
    wz+=wz*rpad;wz+=apad;
    ax=cx-0.5*wx;ay=cy-0.5*wy;az=cz-0.5*wz;

    // Determine the number of blocks to use in the grid data structure
    double ilscale=pow(c/(opt_points*wx*wy*wz),1/3.0);
    m=int(wx*ilscale+1);dx=wx/m;xsp=m/wx;
    n=int(wy*ilscale+1);dy=wy/n;ysp=n/wy;
    o=int(wz*ilscale+1);dz=wz/o;zsp=o/wz;
    printf("# Grid %d %d %d\n",m,n,o);

    // Allocate the memory for the grid blocks
    int nmno=m*n*o;
    if(mno!=nmno) {
        if(mno>0) free_grid();
        co=new int[mno=nmno];
        mem=new int[mno];
        p=new point_info*[mno];
        for(int k=0;k<mno;k++)
            p[k]=new point_info[(mem[k]=pgrid_default_memory)];
    }
}

/** Increases the memory in a grid block.
 * \param[in] ijk the block to consider. */
void proximity_grid::add_memory(int ijk) {
    if(mem[ijk]>=pgrid_max_memory) {
        fprintf(stderr,"Maximum memory allocation exceeded in block %d\n",ijk);
        exit(1);
    }
    int nmem=mem[ijk]<<1;
    point_info *np=new point_info[nmem];
    memcpy(np,p[ijk],sizeof(point_info)*mem[ijk]);
    delete [] p[ijk];
    p[ijk]=np;
    mem[ijk]=nmem;
}

/** Increases memory in the overflow buffer. */
void proximity_grid::add_overflow_memory() {
    if(mem_ov>=pgrid_max_memory) {
        fputs("Maximum memory allocation exceeded in overflow buffer\n",stderr);
        exit(1);
    }
    int nmem=mem_ov<<1;
    oflow_info *nov=new oflow_info[nmem];
    memcpy(nov,ov,sizeof(oflow_info)*mem_ov);
    delete [] ov;
    ov=nov;
    mem_ov=nmem;
}

/** Frees the dynamically allocated memory for the grid. */
void proximity_grid::free_grid() {
    for(int k=mno-1;k>=0;k--) delete [] p[k];
    delete [] p;
    delete [] mem;
    delete [] co;
}
