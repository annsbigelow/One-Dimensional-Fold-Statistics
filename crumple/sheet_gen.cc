#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>
#include "voro++.hh"
using namespace voro;

// Bit flag to signify a mesh point without a full loop of triangles
const unsigned int bflag=1<<31;

// This function returns a random floating point number between 0 and 1
inline double rnd() {return static_cast<double>(rand())*(1./RAND_MAX);}

/** Takes an array of edge connections at a vertex, with zero signifying an
 * unused connection, and extracts out a contiguous list of connections,
 * possibly wrapping around in the given array.
 * \param[in] e an array of connections.
 * \param[in] n the number of elements in e.
 * \param[in] l the index of this vertex.
 * \param[in,out] ncp a pointer to an array in which to store the number of
 *                    connections.
 * \param[in,out] edp a pointer to an array in which to store the contiguous
 *                    list. */
void store_contiguous(int *e,unsigned int n,unsigned int l,unsigned int *&ncp,unsigned int *&edp) {
    unsigned int *edc=edp,j,k;
    if(*e!=0) {
        k=n;

        // The list of connections may wrap around in the e array
        while(k>1&&e[k-1]!=0) k--;
        for(j=k;j<n;j++) *(edp++)=e[j]+l;
        for(j=0;j<k&&e[j]!=0;j++) *(edp++)=e[j]+l;
    } else {

        // The list of connections does not wrap around
        k=1;
        while(k<n&&e[k]==0) k++;
        while(k<n&&e[k]!=0) *(edp++)=e[k++]+l;
    }

    // Determine the total number of connections. If it's fewer than 6, then
    // apply the flag signifying this vertex doesn't have a complete loop
    // of triangles.
    *ncp=static_cast<unsigned int>(edp-edc);
    if(*ncp<n) *ncp|=bflag;
    ncp++;
}

/** Outputs a rectangular sheet mesh to a file, using a regular hexagonal
 * mesh topology.
 * \param[in] fb a file handle to write to.
 * \param[in] s the edge length.
 * \param[in] (nx,ny) the dimensions of the mesh. */
void output_rect_mesh(FILE *fb,double s,int nx,int ny) {
    double h=0.5*sqrt(3)*s;
    int i,j,e[6],nt,l=0;

    // Compute total nodes and triangles. Allocate memory.
    unsigned int n=nx*ny+(ny>>1),ned=nx*(6*ny-4)-ny-(ny&1),
                 *ibuf=new unsigned int[2+n+ned],
                 *ncp=ibuf+2,*edp=ncp+n;
    *ibuf=n;ibuf[1]=ned;

    // Compute and output the edge table
    for(j=0;j<ny;j++) {
        nt=nx+(j&1);
        for(i=0;i<nt;i++,l++) {
            if(i>0&&i<nt-1&&j>0&&j<ny-1) {

                // Deal with the interior case
                *(ncp++)=6;
                *(edp++)=-nx-1+l;*(edp++)=-nx+l;
                *(edp++)=1+l;*(edp++)=nx+1+l;
                *(edp++)=nx+l;*(edp++)=-1+l;
            } else {

                // Assemble top and bottom row connections (if present)
                if(j==0) *e=e[1]=0;
                else {*e=-nx-1;e[1]=-nx;}
                if(j==ny-1) e[3]=e[4]=0;
                else {e[3]=nx+1;e[4]=nx;}

                // Assemble left and right connections (if present), accounting
                // for possible additional missed neighbors on odd rows
                if(i==0) {
                    e[5]=0;
                    if(j&1) *e=e[4]=0;
                } else e[5]=-1;
                if(i==nt-1) {
                    e[2]=0;
                    if(j&1) e[1]=e[3]=0;
                } else e[2]=1;

                // Store a contiguous loop of edges into the table
                store_contiguous(e,6,l,ncp,edp);
            }
        }
    }

    // Write all of the data
    fwrite(ibuf,sizeof(int),2+n+ned,fb);
    delete [] ibuf;

    // Compute node positions
    double *dbuf=new double[3*n],*dp=dbuf;
    for(j=0;j<ny;j++) {
        if(j&1) {

            // For odd rows, adjust the end points
            *(dp++)=0.;*(dp++)=j*h;*(dp++)=0.;
            for(i=0;i<nx-1;i++) {
                *(dp++)=(0.5+i)*s;*(dp++)=j*h;*(dp++)=0.;
            }
            *(dp++)=(nx-1)*s;*(dp++)=j*h;*(dp++)=0.;
        } else {

            // For even rows, no adjustments are needed
            for(i=0;i<nx;i++) {
                *(dp++)=i*s;*(dp++)=j*h;*(dp++)=0.;
            }
        }
    }

    // Output the node positions and free memory
    fwrite(dbuf,sizeof(double),3*n,fb);
    delete [] dbuf;
}

/** Outputs a triangular sheet mesh to a file, using a regular hexagonal
 * mesh topology.
 * \param[in] fb a file handle to write to.
 * \param[in] s the edge length.
 * \param[in] nx the number of nodes along a side of the triangle. */
void output_tri_mesh(FILE *fb,double s,int nx) {
    double h=0.5*sqrt(3)*s;
    int i,j,e[6],l=0;

    // Compute total nodes and triangles. Allocate memory.
    unsigned int n=(nx*(nx+1))>>1,ned=3*nx*(nx-1),
                 *ibuf=new unsigned int[2+n+ned],
                 *ncp=ibuf+2,*edp=ncp+n;
    *ibuf=n;ibuf[1]=ned;

    for(j=nx-1;j>=0;j--) {
        for(i=0;i<=j;i++,l++) {
            if(false) {//i==0||i==j||j==0) {

                // Deal with the interior case
                *(ncp++)=6;
                *(edp++)=-j-1+l;*(edp++)=-j+l;
                *(edp++)=1+l;*(edp++)=j+l;
                *(edp++)=j-1+l;*(edp++)=-1+l;
            } else {

                // Assemble bottom row connections (if present)
                if(j==nx-1) *e=e[1]=0;
                else {*e=-j-2;e[1]=-j-1;}

                // Assemble up-left and up-right connections (if present)
                if(i==0) {e[5]=e[4]=0;}
                else {e[5]=-1;e[4]=j;}
                if(i==j) {e[2]=e[3]=0;}
                else {e[2]=1;e[3]=j+1;}

                // Store a contiguous loop of edges into the table
                store_contiguous(e,6,l,ncp,edp);
            }
        }
    }

    // Write all of the data
    fwrite(ibuf,sizeof(int),2+n+ned,fb);
    delete [] ibuf;

    // Compute node positions
    double *dbuf=new double[3*n],*dp=dbuf;
    for(j=nx-1;j>=0;j--) {
        for(i=0;i<=j;i++) {
            *(dp++)=s*(i+(nx-1-j)*0.5);
            *(dp++)=h*(nx-1-j);
            *(dp++)=0.;
        }
    }

    // Output the node positions and free memory
    fwrite(dbuf,sizeof(double),3*n,fb);
    delete [] dbuf;
}

/** Outputs a rectangular sheet mesh to a file, using alternating 4 and 8 point
 * vertices.
 * \param[in] fb a file handle to write to.
 * \param[in] s the edge length.
 * \param[in] (nx,ny) the dimensions of the mesh. */
void output_rect48_mesh(FILE *fb,double s,int nx,int ny) {
    int i,j,e[8],l=0;

    // Compute total nodes and triangles. Allocate memory.
    unsigned int n=nx*ny,ned=6*nx*ny-4*(nx+ny)+2,
                 *ibuf=new unsigned int[2+n+ned],
                 *ncp=ibuf+2,*edp=ncp+n;
    *ibuf=n;ibuf[1]=ned;

    for(j=0;j<ny;j++) {
        for(i=0;i<nx;i++,l++) {
            if((i+j)&1) {

                // Deal with the 4-edge vertex
                *e=i==nx-1?0:1;
                e[1]=j==ny-1?0:nx;
                e[2]=i==0?0:-1;
                e[3]=j==0?0:-nx;
                store_contiguous(e,4,l,ncp,edp);

            } else {

                // Deal with the 8-edge vertex. Set diagonal connections.
                e[1]=nx+1;e[3]=nx-1;e[5]=-nx-1;e[7]=-nx+1;

                // Set orthogonal connections and remove diagonals if needed
                if(i==nx-1) *e=e[1]=e[7]=0;else *e=1;
                if(j==ny-1) e[1]=e[2]=e[3]=0;else e[2]=nx;
                if(i==0) e[3]=e[4]=e[5]=0;else e[4]=-1;
                if(j==0) e[5]=e[6]=e[7]=0;else e[6]=-nx;

                // Store a contiguous loop of edges into the table
                store_contiguous(e,8,l,ncp,edp);
            }
        }
    }

    // Write all of the data
    fwrite(ibuf,sizeof(int),2+n+ned,fb);
    delete [] ibuf;

    // Compute node positions
    double *dbuf=new double[3*n],*dp=dbuf;
    for(j=0;j<ny;j++) {
        for(i=0;i<nx;i++) {
            *(dp++)=s*i;
            *(dp++)=s*j;
            *(dp++)=0.;
        }
    }

    // Output the node positions and free memory
    fwrite(dbuf,sizeof(double),3*n,fb);
    delete [] dbuf;
}

/** Outputs a hexagonal sheet with a dislocation in the center.
 * \param[in] fb a file handle to write to.
 * \param[in] s the edge length.
 * \param[in] nx the dimension of the mesh. */
void output_disloc_mesh(FILE *fb,double s,int nx) {
    int i,j,e[7],l=0,w,q=0;
    double h=0.5*sqrt(3)*s;

    // Compute total nodes and triangles. Allocate memory.
    unsigned int n=1+nx*(2+3*nx),ned=18*nx*nx+2,
                 *ibuf=new unsigned int[2+n+ned],
                 *ncp=ibuf+2,*edp=ncp+n;
    *ibuf=n;ibuf[1]=ned;

    // Part I: the lower half of the hexagon
    for(j=0;j<nx;j++) {
        w=j+nx;
        for(i=0;i<w;i++,l++) {

            // Set connections
            *e=1;e[1]=-w+1;e[2]=-w;
            e[3]=-1;e[4]=w+q;e[5]=w+1+q;

            // Remove unused connections
            if(j==0) e[1]=e[2]=0;
            if(i==0) e[2]=e[3]=0;
            if(i==w-1) *e=e[1]=0;

            // Store a contiguous loop of edges into the table
            if(j==nx-1&&i==nx-1) {
                e[6]=w+2;q=1;
                store_contiguous(e,7,l,ncp,edp);
            } else store_contiguous(e,6,l,ncp,edp);
        }
    }

    // Part II: the upper half of the hexagon
    q=2;
    for(j=nx;j<2*nx+1;j++) {
        w=3*nx+1-j;
        for(i=0;i<w;i++,l++) {

            // Set connections
            *e=1;e[1]=-w+q;e[2]=-w-1+q;
            e[3]=-1;e[4]=w-1;e[5]=w;

            // Remove unused connections
            if(j==2*nx) e[4]=e[5]=0;
            if(i==0) {e[3]=e[4]=0;if(j==nx) e[2]=0;}
            if(i==w-1) {*e=e[5]=0;if(j==nx) {e[1]=0;q=0;}}

            // Store a contiguous loop of edges into the table
            if(j==nx&&i==nx) {
                e[1]=*e;q=1;
                store_contiguous(e+1,5,l,ncp,edp);
            } else store_contiguous(e,6,l,ncp,edp);
        }
    }

    // Write all of the data
    fwrite(ibuf,sizeof(int),2+n+ned,fb);
    delete [] ibuf;

    // Compute node positions
    double *dbuf=new double[3*n],*dp=dbuf,ax,dx;
    for(j=0;j<2*nx+1;j++) {
        w=j<nx?nx+j:3*nx+1-j;
        ax=-0.5*s*(j<nx?nx+j:3*nx-j);
        dx=-2*ax/(w-1);
        for(i=0;i<w;i++) {
            *(dp++)=ax+dx*i;
            *(dp++)=h*(j-nx);
            *(dp++)=0.;
        }
    }

    // Output the node positions and free memory
    fwrite(dbuf,sizeof(double),3*n,fb);
    delete [] dbuf;
}

/** Outputs a rectangular sheet mesh to a file, using a random mesh topology.
 * \param[in] fb a file handle to write to.
 * \param[in] seed a seed for the random node position generation.
 * \param[in] num_nodes the number of nodes in the mesh. */
void output_rand_mesh(FILE *fb,unsigned int seed,unsigned int num_nodes) {
    const unsigned int n=num_nodes;
    int i,j,l; double x,y;
    int npg=16*(1+(n/100/16)); // initial memory allocation per grid square; round up to nearest multiple of 16
    // Initialize the container class to be the unit square, with
    // non-periodic boundary conditions. Divide it into a 10 by 10 grid,
    // with an initial memory allocation of npg particles per grid square.
    container_2d con(0,1,0,1,10,10,false,false,npg);

    // Create y positions and sort them
    double *ys=new double[n];
    srand(seed);
    for(i=0;i<(signed int) n;i++) ys[i]=rnd();
    std::sort(ys,ys+n);

    // Add n random points to the container
    particle_order po;
    for(i=0;i<(signed int) n;i++) con.put(po,i,rnd(),ys[i]);

    // Delete allocated memory
    delete [] ys;

    // Voro++ loop object
    container_2d::iterator cli;

    // Voronoi cell with neighbor information
    voronoicell_neighbor_2d c(con);
    std::vector<int> nei;

    // Loop through once to get the total number of edge connections.
    unsigned int ned=0;
    for(cli=con.begin();cli<con.end();cli++) if(con.compute_cell(c,cli)) {
        c.neighbors_sorted(nei);
        j=0;
        for(l=0;l<(signed int) nei.size();l++) {
            if(nei[l]>=0) {ned++; j++;}
        }
        if(j<2) printf("Warning: found a node with only one neighbor.\n");
    }

    // Allocate memory for edge connections and node positions.
    unsigned int *ibuf=new unsigned int[2+n+ned],
             *ncp=ibuf+2,*edp=ncp+n;
    *ibuf=n;ibuf[1]=ned;
    double *dbuf=new double[3*n],*dp=dbuf;
    bool edge;
    std::vector<int> nei2;
    std::vector<int>::iterator it;

    // Loop over all particles and print the neighbors of each (defined in
    // terms of particles that share a Voronoi cell edge)
    for(cli=con.begin();cli<con.end();cli++) if(con.compute_cell(c,cli)) {

        // This function should ensure that the neighbors are printed in
        // an anticlockwise order. Negative values from -1 to -4 correspond
        // to edges that touch the walls
        c.neighbors_sorted(nei2);
        nei.erase(nei.begin(),nei.end()); // clear array
        it = nei.begin(); // get iterator at the start

        // Print particle ID and position, plus the neighbor list
        con.pos(cli,x,y);
        *(dp++)=x;*(dp++)=y;*(dp++)=0;
        edge=false;
        *ncp=0;
        for(l=0;l<(signed int) nei2.size();l++) {
        if(nei2[l]>=0) {
                it=nei.insert(it,nei2[l]);
                it++; // increment iterator
            }
            else {
                edge=true;
                it=nei.begin();
            }
        }
        // Now loop through nei to store in correct order
        for(l=0;l<(signed int) nei.size();l++) {
            (*ncp)++;
            *(edp++)=(unsigned int) nei[l];
        }
        if(edge) *ncp|=bflag;
        ncp++;

    }

    // Write all of the data
    fwrite(ibuf,sizeof(int),2+n+ned,fb);
    delete [] ibuf;

    fwrite(dbuf,sizeof(double),3*n,fb);
    delete [] dbuf;
}

int main(int argc,char **argv) {

    // Check for the correct number of command-line arguments
    if(argc<4||argc>5) {
        fputs("Syntax: ./sheet_gen rec <side_length> <N_long> <N_short>\n"
              "        ./sheet_gen rec48 <side_length> <N_long> <N_short>\n"
              "        ./sheet_gen tri <side_length> <N>\n"
              "        ./sheet_gen disloc <side_length> <N>\n"
              "        ./sheet_gen rand <seed> <n>\n",stderr);
        return 1;
    }

    // Read the mesh type
    int mode,eargs;
    if(strcmp(argv[1],"rec")==0) {mode=0;eargs=5;}
    else if(strcmp(argv[1],"rec48")==0) {mode=1;eargs=5;}
    else if(strcmp(argv[1],"tri")==0) {mode=2;eargs=4;}
    else if(strcmp(argv[1],"disloc")==0) {mode=3;eargs=4;}
    else if(strcmp(argv[1],"rand")==0) {mode=4;eargs=4;}
    else {
        fputs("Unknown mesh type\n",stderr);
        return 1;
    }

    // Read the command-line arguments and do some basic checks
    if(argc!=eargs) {
        fprintf(stderr,"Expected %d arguments, but got %d\n",eargs,argc);
        return 1;
    }
    double s=atof(argv[2]);
    unsigned int seed=static_cast<unsigned int>(s); // only for random sheet
    int nx=atoi(argv[3]),ny=argc==4?2:atoi(argv[4]);
    if(nx<=1||ny<=1) {
        fputs("Dimensions must both be greater than 1\n",stderr);
        return 1;
    }

    // For the rectangular meshes, reorder if short side > long side
    if(mode!=1&&ny>nx) {int tmp=nx;nx=ny;ny=tmp;}

    // Open the file
    char buf[128];
    if(mode==0) sprintf(buf,"sheet_%g_%dx%d.bin",s,nx,ny);
    else if(mode==1) sprintf(buf,"sh48_%dx%d.bin",nx,ny);
    else if(mode==2) sprintf(buf,"tri_%d.bin",nx);
    else if(mode==3) sprintf(buf,"disloc_%d.bin",nx);
    else sprintf(buf,"rsheet_%d_%d.bin",nx,seed);
    FILE *fb=fopen(buf,"wb");
    if(fb==NULL) {
        fprintf(stderr,"Can't open file\n");
        exit(1);
    }

    // Output the mesh and close the file
    if(mode==0) output_rect_mesh(fb,s,nx,ny);
    else if(mode==1) output_rect48_mesh(fb,s,nx,ny);
    else if(mode==2) output_tri_mesh(fb,s,nx);
    else if(mode==3) output_disloc_mesh(fb,s,nx);
    else output_rand_mesh(fb,seed,nx);
    fclose(fb);
}
