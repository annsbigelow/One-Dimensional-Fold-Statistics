#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "mesh.hh"

int main(int argc,char **argv) {

    // Check for the correct number of command-line arguments
    if(argc<4||argc>5) {
        fputs("Syntax: ./unpack <mode> <directory> <frame> <output_filename>\n"
              "        ./unpack <mode> <input_filename> <output_filename>\n\n"
              "Mode: \"msh\" for POV-Ray mesh (with normals)\n"
              "      \"mtr\" for POV-Ray mesh (with flat triangles, no normals)\n"
              "      \"gnu\" for Gnuplot mesh\n"
              "      \"cyl\" for POV-Ray cylinders\n"
              "      \"sph\" for POV-Ray spheres\n"
              "      \"txt\" for plain text vertices\n"
              "      \"edg\" for plain text edge table\n\n"
              "If the output filename is \"-\" then the data will written to standard\n"
              "output\n",stderr);
        return 1;
    }

    // Find output mode
    int mode,fnum;
    if(strcmp(argv[1],"msh")==0) mode=0;
    else if(strcmp(argv[1],"mtr")==0) mode=1;
    else if(strcmp(argv[1],"gnu")==0) mode=2;
    else if(strcmp(argv[1],"cyl")==0) mode=3;
    else if(strcmp(argv[1],"sph")==0) mode=4;
    else if(strcmp(argv[1],"txt")==0) mode=5;
    else if(strcmp(argv[1],"edg")==0) mode=6;
    else {
        fprintf(stderr,"Mode type \"%s\" not known\n",argv[1]);
        return 1;
    }

    // Read in the mesh
    mesh_param par(0.05,0.02,false);
    mesh *mp;
    if(argc==5) {

        // If the are five command line arguments, then look for mesh
        // information and topology separately in an output directory. First,
        // check the frame number is sensible.
        fnum=atoi(argv[3]);
        if(fnum<0||fnum>16777216) {
            fputs("Frame number out of range\n",stderr);
            return 1;
        }

        // Assemble input filenames and read in the data
        size_t l=strlen(argv[2])+16;
        char *f_topo=new char[2*l],*f_pts=f_topo+l;
        sprintf(f_topo,"%s/topo",argv[2]);
        sprintf(f_pts,"%s/pts.%d",argv[2],fnum);
        mp=new mesh(par,f_topo,f_pts);

        // Free the dynamically allocated memory
        delete [] f_topo;
    } else mp=new mesh(par,argv[2]);

    // Output the relevant data, checking for when the output filename is "-"
    // to signal that the data should be sent to standard output
    bool std=strcmp(argv[argc-1],"-")==0;
    FILE *fp=std?stdout:safe_fopen(argv[argc-1],"w");
    switch(mode) {
        case 0: mp->draw_mesh_pov(fp);break;
        case 1: mp->draw_mesh_pov(fp,false);break;
        case 2: mp->draw_mesh_gnuplot(fp);break;
        case 3: mp->draw_cylinders_pov(fp);break;
        case 4: mp->draw_nodes_pov(fp);break;
        case 5: mp->draw_nodes(fp);break;
        case 6: mp->edge_diagnostic(fp);
    }
    if(!std) fclose(fp);
    else fflush(stdout);

    // Free the dynamically allocated mesh memory
    delete mp;
}
