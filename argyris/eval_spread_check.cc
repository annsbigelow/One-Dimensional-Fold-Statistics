#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>

#include "mesh.hh"

int main() {
	// Set parameters
	bool fix_boundary=true;
	char buf[50], buf1[50];
	bool r_all_dofs=false, wr_all_dofs=false;
	const double K=.05,drag=.02,kb=.001;

	std::ofstream File("Kevals.csv");
	
	int nx_array[9]={7,8,10,11,13,15,18,22,25};
	for(int k=0;k<9;k++) {
		int nx=nx_array[k];
		float s=5.0/nx;

		// Create mesh
		mesh_param par(K,drag,kb,false,fix_boundary,r_all_dofs,wr_all_dofs,s);
		sprintf(buf,"./sheet_gen rec48 %f %d %d",s,nx,nx);
		std::system(buf);
		sprintf(buf1,"sh48_%dx%d.bin",nx,nx);
		mesh_rk4 mp(par,buf1);

		// Centralize and scale the mesh
		double wx, wy, wz;
		mp.centralize(wx, wy, wz);

		// Build FEM matrices.
		mp.output_refined=false;
		mp.setup_springs(); mp.setup_fem(); mp.assemble_K();
		mp.checkK(File);

		printf("Eigenvalue spectrum check finished for nx=%d (n=%d).\n",nx,mp.n);
	}
	File.close();

	//// Check the kernel of the stiffness matrix (the biharmonic operator). ////

	// Create mesh
	/*int nx=7; float s=5.0/nx;
	mesh_param par(K,drag,kb,false,fix_boundary,r_all_dofs,wr_all_dofs,s);
	sprintf(buf,"./sheet_gen rec48 %f %d %d",s,nx,nx);
	std::system(buf);
	sprintf(buf1,"sh48_%dx%d.bin",nx,nx);
	mesh_rk4 mp(par,buf1);

	// Centralize and scale the mesh
	double wx, wy, wz;
	mp.centralize(wx, wy, wz);

	// Build FEM matrices.
	mp.output_refined=false;
	mp.setup_springs(); mp.setup_fem(); mp.assemble_K();

	mp.check_kernel_K();*/
}