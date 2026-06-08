#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "mesh.hh"

int main() {

	// Create mesh and initialize acceleration
	mesh_param par(0.05, 0.03, 0.001, false, false);
	mesh_rk4 mp(par, "../crumple/sh48_43x43.bin");

	// Centralize and scale the mesh
	double wx, wy, wz;
	mp.centralize(wx, wy, wz);

	/* Set up springs. Assign lengths manually, since initializing out of
	equilibrium configuration. */
	mp.lump = false;
	mp.PCG = false; 
	mp.CG = false;
	mp.setup_springs();
}
