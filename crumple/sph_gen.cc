#include "wall_shell.hh"

#include "voro++.hh"
using namespace voro;

#include <vector>
using namespace std;

// Set up constants for the container geometry
const double boxl=1.2;

// Set up the number of blocks that the container is divided into
const int bl=14;

// Set the number of particles that are going to be randomly introduced
const int particles=2000;

// Maximum number of edges to collect face information for
const int nface=11;

// This function returns a random double between -1 and 1
inline double rnds() {return static_cast<double>(rand())*(2./RAND_MAX)-1;}

int main(int argc,char **argv) {

    // Check command line arguments
    if(argc!=3) {
        fputs("Syntax: ./sph_gen <num_particles> <lloyd_iterations>\n",stderr);
        return 1;
    }
    int particles=atoi(argv[1]);
    if(particles<=1||particles>16777216) {
        fputs("Particle number out of range\n",stderr);
        return 1;
    }
    int lloyd=atoi(argv[2]);
    if(lloyd<0||lloyd>16777216) {
        fputs("Lloyd iterations out of range\n",stderr);
        return 1;
    }

    // Initialize memory
    char buf[128];
    int i=0,j,k,l,ll,o;
    double x,y,z,r,dx,dy,dz;
    int faces[nface],*fp;
    double *p=new double[3*particles];

    // Create a container with the geometry given above, and make it
    // non-periodic in each of the three coordinates. Allocate space for eight
    // particles within each computational block
    container con(-boxl,boxl,-boxl,boxl,-boxl,boxl,bl,bl,bl,false,false,false,8);

    // Add the custom shell wall object
    wall_shell ws(0,0,0,1,0.00001);
    con.add_wall(ws);

    // Randomly add particles into the container
    while(i<particles) {
        x=rnds();y=rnds();z=rnds();
        r=x*x+y*y+z*z;
        if(r>1e-5&&r<1) {
            r=1/sqrt(r);x*=r;y*=r;z*=r;
            con.put(i,x,y,z);
            i++;
        }
    }

    // Perform the Lloyd iterations
    for(l=0;l<lloyd;l++) {
        c_loop_all vl(con);
        voronoicell c;
        for(fp=faces;fp<faces+nface;fp++) *fp=0;

        // Move each Voronoi cell to its centroid position
        if(vl.start()) do if(con.compute_cell(c,vl)) {
            vl.pos(i,x,y,z,r);
            c.centroid(dx,dy,dz);
            p[3*i]=x+dx;
            p[3*i+1]=y+dy;
            p[3*i+2]=z+dz;

            // Tally the number of faces
            i=c.number_of_faces()-4;
            if(i<0) i=0;else if(i>=nface) i=nface-1;
            faces[i]++;
        } while (vl.inc());
        con.clear();

        // Apply a random perturbation to the points
        double fac=l<10*lloyd/9?0.1/sqrt(double(l)):0.;
        for(i=0;i<particles;i++) con.put(i,p[3*i]+fac*rnds(),
                p[3*i+1]+fac*rnds(),p[3*i+2]+fac*rnds());
        printf("%d",l);

        // Print the frequencies of each type of face
        for(fp=faces;fp<faces+nface;fp++) printf(" %d",*fp);
        puts("");
    }

    // Output the particle positions in gnuplot format
    sprintf(buf,"sph_%d_p.gnu",particles);
    con.draw_particles(buf);

    // Output the Voronoi cells in gnuplot format
    sprintf(buf,"sph_%d_v.gnu",particles);
    con.draw_cells_gnuplot(buf);

    // Allocate memory for neighbor relations
    int *q=new int[particles*nface],*qn=new int[particles],*qp;
    for(l=0;l<particles;l++) qn[l]=0;

    // Create a table of all neighbor relations
    vector<int> vi;
    voronoicell_neighbor c;
    c_loop_all vl(con);
    if(vl.start()) do if(con.compute_cell(c,vl)) {
        i=vl.pid();qp=q+i*nface;
        c.neighbors(vi);
        if(vi.size()>nface+2) voro_fatal_error("Too many faces; boost nface",5);
        for(l=0;l<(signed int) vi.size();l++) if(vi[l]>=0) qp[qn[i]++]=vi[l];
    } while (vl.inc());

    // Sort the connections in anti-clockwise order
    bool connect;
    int tote=0;
    for(l=0;l<particles;l++) {
        tote+=qn[l];
        for(i=0;i<qn[l]-2;i++) {
            o=q[l*nface+i];
            j=i+1;
            while(j<qn[l]-1) {
                ll=q[l*nface+j];
                connect=false;
                for(k=0;k<qn[ll];k++) {
                    if(q[ll*nface+k]==o) {connect=true;break;}
                }
                if(connect) break;
                j++;
            }

            // Swap the connected vertex into this location
            o=q[l*nface+i+1];
            q[l*nface+i+1]=q[l*nface+j];
            q[l*nface+j]=o;
        }

        // Reverse all connections if the have the wrong handedness
        j=3*l;k=3*q[l*nface];o=3*q[l*nface+1];
        x=p[j]-p[k];dx=p[j]-p[o];
        y=p[j+1]-p[k+1];dy=p[j+1]-p[o+1];
        z=p[j+2]-p[k+2];dz=p[j+2]-p[o+2];
        if(p[j]*(y*dz-z*dy)+p[j+1]*(z*dx-x*dz)+p[j+2]*(x*dy-y*dx)<0) {
            for(i=0;i<qn[l]/2;i++) {
                o=q[l*nface+i];
                q[l*nface+i]=q[l*nface+qn[l]-1-i];
                q[l*nface+qn[l]-1-i]=o;
            }
        }
    }

    // Output the mesh network based on Voronoi neighbors
    sprintf(buf,"sph_%d.net",particles);
    FILE *ff=safe_fopen(buf,"w");
    int *mp=new int[particles],*mpi=new int[particles];
    for(i=0;i<particles;i++) mp[i]=-1;
    *mpi=0;*mp=0;l=1;o=0;
    while(o<l) {
        i=mpi[o];
        for(j=0;j<qn[i];j++) {
            k=q[i*nface+j];
            if(mp[k]==-1) {
                mpi[l]=k;
                mp[k]=l++;
            }
        }
        o++;
    }
    fclose(ff);

    // Save binary representation of the mesh
    sprintf(buf,"sph_%d.bin",particles);
    FILE *fb=safe_fopen(buf,"wb");

    // Write header
    int *ibuf=new int[2+particles+tote],*ip=ibuf;
    *(ip++)=particles;
    *(ip++)=tote;

    // Assemble the connections and write them
    for(l=0;l<particles;l++) *(ip++)=qn[mpi[l]];
    for(l=0;l<particles;l++) {
        i=mpi[l];
        for(j=0;j<qn[i];j++) *(ip++)=mp[q[i*nface+j]];
    }
    fwrite(ibuf,sizeof(int),2+particles+tote,fb);

    // Write particle positions
    double *pm=new double[3*particles],*a=pm,*b;
    for(i=0;i<particles;i++) {
        b=p+3*mpi[i];
        *(a++)=*(b++);*(a++)=*(b++);*(a++)=*b;
    }
    fwrite(pm,sizeof(double),3*particles,fb);
    delete [] pm;

    // Free dynamically allocated arrays
    delete [] ibuf;
    delete [] mpi;
    delete [] mp;
    delete [] qn;
    delete [] q;
    delete [] p;
}
