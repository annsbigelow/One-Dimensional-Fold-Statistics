#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "mesh.hh"
#include "common.hh"

int main(int argc,char **argv) {

	// Check for the correct number of command-line arguments
	if (argc!=4) {
		fputs("Syntax: ./deformation_data <mode> <input_directory> <frame>\n"
			"Mode: \"0\" for single frame\n"
			"		\1\" for many frames\n"
			"Frame: frame number for single; num_frames for many\n", stderr);
		return 1;
	}

	// Set dimensions of the mesh
	int nx=100,ny=100; 

	// Find output mode
	int mode;
	if(strcmp(argv[1],"0")==0) mode=0;
	else mode=1;

	bool area=true;

	// Read in the mesh
	mesh_param par(0.5,0.01,0,0.2,false,true,0.001,1.3);
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

		double frac=.2;
		// Choose a percentage of the nodes 
		//mp->setup_springs();
		//mp->select_subsheet(frac,nx,ny);

		// Print the standard deviation of z-coordinates of the nodes
		printf("Roughness measure: %g\n", mp->sdev(frac,nx,ny));

		// Setup triangle info and print the sheet area
		if(area) {
			mp->setup_springs();
			printf("Sheet area: %g\n", mp->tot_area(frac,nx,ny));
		}

		delete mp;
	}
	else {
		double* sdevs = new double[fnum];
		double* area_arr = new double[fnum];
		// Read in the data for each frame
		for(int j=0;j<fnum;j++) {
			sprintf(f_topo, "%s/topo", argv[2]);
			sprintf(f_pts, "%s/pts.%d", argv[2], j);
			mp=new mesh(par, f_topo, f_pts);

			sdevs[j]=mp->sdev(.2,nx,ny);

			if (area) {
				mp->setup_springs();
				area_arr[j]=mp->tot_area(.2,nx,ny);
			}

			delete mp;
		}

		// Write the data to a file
		FILE* fp=safe_fopen("sdevs.bin", "wb");
		fwrite(sdevs,sizeof(double),fnum,fp);
		fclose(fp);
		FILE* fp1 = safe_fopen("areas.bin", "wb");
		fwrite(area_arr, sizeof(double), fnum, fp1);
		fclose(fp1);

		delete[] area_arr;
		delete[] sdevs;
	}
	delete[] f_topo;
}
