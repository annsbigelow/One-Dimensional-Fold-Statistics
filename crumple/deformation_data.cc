#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "mesh.hh"
#include "common.hh"

int main(int argc,char **argv) {

	// Check for the correct command-line arguments
	if (argc!=8) {
		fputs("Syntax: ./deformation_data <mode> <input_directory> <frame> <nx> <ny> <s> <topo>\n"
			"Mode: \"0\" for single frame\n"
			"		\"1\" for many frames\n"
			"Frame: frame number for single; num_frames for many\n"
			"nx, ny: dimensions of the sheet\n"
			"s: the side edge length\n"
			"topo: \"0\" for regular hexagonal topology\n", stderr);
		return 1;
	}
	if(strcmp(argv[7],"0")!=0) {
		fputs("This mesh topology is not supported yet.\n",stderr);
		return 1;
	}

	// Set dimensions of the mesh
	int nx=atoi(argv[4]),ny=atoi(argv[5]);
	double sed=atof(argv[6]);

	// Find output mode
	int mode;
	if(strcmp(argv[1],"0")==0) mode=0;
	else mode=1;

	// Read in the mesh
	mesh_param par(0.5,0.01,0,0.2,false,true,1.3,sed);
	mesh *mp;

	int fnum = atoi(argv[3]);
	if (fnum < 0 || fnum>16777216) {
		fputs("Frame number out of range\n", stderr);
		return 1;
	}

	// Assemble input filenames
	size_t l = strlen(argv[2])+16;
	char* f_topo=new char[2*l], *f_pts=f_topo+l;

	if (mode == 0) {
		// Read in the data
		sprintf(f_topo, "%s/topo", argv[2]);
		sprintf(f_pts, "%s/pts.%d", argv[2], fnum);
		mp = new mesh(par, f_topo, f_pts);

		// Setup triangle info and subsheet boundaries
		mp->setup_springs();
		// Print the roughness metrics
		double Sq,Sa; mp->Sq_Sa(Sq,Sa,nx,ny);
		printf("Root mean square height: %g\nArithmetical mean height: %g\n",Sq,Sa);
		printf("Sheet area: %g\n", mp->tot_area_rec(nx,ny));
		delete mp;
	}
	else {
		double* sdevs = new double[fnum];
		double* Sas=new double[fnum];
		double* area_arr = new double[fnum];
		// Read in the data for each frame
		for(int j=0;j<fnum;j++) {
			sprintf(f_topo, "%s/topo", argv[2]);
			sprintf(f_pts, "%s/pts.%d", argv[2], j);
			mp=new mesh(par, f_topo, f_pts);
			mp->setup_springs();

			double Sq,Sa;
			mp->Sq_Sa(Sq,Sa,nx,ny);
			sdevs[j]=Sq;
			Sas[j]=Sa;
			area_arr[j]=mp->tot_area_rec(nx,ny);

			delete mp;
		}

		// Write the data to a files
		FILE* fp=safe_fopen("sdevs.bin", "wb");
		fwrite(sdevs,sizeof(double),fnum,fp);
		fclose(fp);
		FILE* fp1 = safe_fopen("areas.bin", "wb");
		fwrite(area_arr, sizeof(double), fnum, fp1);
		fclose(fp1);
		FILE* fp2=safe_fopen("mean_heights.bin","wb");
		fwrite(Sas,sizeof(double),fnum,fp2);
		fclose(fp2);

		delete[] area_arr;
		delete[] sdevs;
		delete[] Sas;
	}
	delete[] f_topo;
}
