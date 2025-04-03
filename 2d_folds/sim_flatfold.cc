#include "sim_flatfold.hh"
#include <cstdlib>
#include <cmath>
#include <vector>

//const double M_PI=3.1415926535897932384626433832795;

/** Initializes the flatfold class to be a square sheet covering [-1,1]^2. */
sim_flatfold::sim_flatfold() : rng(gsl_rng_alloc(gsl_rng_taus)) {
	f.push_back(new facet(-1,1,-1,1));
	cx=0;cy=0;crsq=2.;cr=sqrt(2.);
//	v={-1.,1., 1.,1., 1.,-1., -1.,-1.};
}

/** Initializes the flatfold class to be a rectangular sheet.
 * \param[in] (ax,bx) the x-dimensions of the sheet.
 * \param[in] (ay,by) the y-dimensions of the sheet. */
sim_flatfold::sim_flatfold(double ax,double bx,double ay,double by) :
	rng(gsl_rng_alloc(gsl_rng_taus)) {
	f.push_back(new facet(ax,bx,ay,by));
	cx=0.5*(bx+ax);cy=0.5*(by+ay);
	double lx=bx-ax,ly=by-ay;
	crsq=0.25*(lx*ly+ly*ly);
	cr=sqrt(crsq);
}

/** The class destructor frees the dynamically allocated facets. */
sim_flatfold::~sim_flatfold() {
	for(unsigned int i=0;i<f.size();i++) delete f[i];
	gsl_rng_free(rng);
}

/** Computes a bounding circle to the current sheet, by first finding the
 * center of mass, and then finding the maximum distance to a vertex. */
void sim_flatfold::compute_bounds() {
	cx=0;cy=0;
	double sa=0,rsq;

	// Find the center of mass
	for(unsigned int i=0;i<f.size();i++)
		f[i]->weight_contrib(cx,cy,sa);
	cx/=sa;cy/=sa;

	// Compute the maximum distance of a vertex to the center of mass. For
	// quick lookups, both the distance and the square distance are stored.
	crsq=0.;
	for(unsigned int i=0;i<f.size();i++)
		if((rsq=f[i]->max_rad_sq(cx,cy))>crsq) crsq=rsq;
	cr=sqrt(crsq);
}

/** Finds a random point on a random edge.
* \param[out] (epx,epy) The coordinates of the point. */
void sim_flatfold::ed_pts(double &epx, double &epy) {
	std::vector<int> ct; // Use a vector because the size isn't initially obvious
	std::vector<double> v; // use a 2D array?

	// Ignore stacked edges
	int j=0; 
	for (int l=0; l<f.size(); l++) {
		int k=0; 
		do {
			v.push_back(f[l]->c.pts[2*k]); v.push_back(f[l]->c.pts[2*k+1]);
			int q=f[l]->c.ed[2*k];
			v.push_back(f[l]->c.pts[2*q]); v.push_back(f[l]->c.pts[2*q+1]);
			ct.push_back(1);
			// Check for recurring edges from any other facet
			for (int m=0; m<v.size();) {
				if (j!=m && v[j]==v[m]&&v[j+1]==v[m+1]&&v[j+2]==v[m+2]&&v[j+3]==v[m+3]) ct[j/4]=0;
				m+=4;
			}
			k=f[l]->c.ed[2*k];
			j+=4;
		} while (k!=0);
	}

	// Pick a random edge using nonzero entries of ct 
	std::vector<int> nz;
	for (unsigned int k=0; k<ct.size(); k++) { if (ct[k]!=0) nz.push_back(k); }
	if (nz.empty()) { fputs("Unable to find edges\n",stderr); exit(1); }
	unsigned int idx = gsl_rng_uniform_int(rng, nz.size());
	idx=nz[idx]*4; // v has four times the size of ct.
	double v1x=v[idx], v1y=v[idx+1], v2x=v[idx+2], v2y=v[idx+3];

	epx = v1x + (v2x-v1x)*gsl_rng_uniform(rng); epy = v1y + (v2y-v1y)*gsl_rng_uniform(rng);
}

/** Uses two points along the boundary to fold the sheet.
* \param[in] rand_sign whether to choose a random sign for the fold or not. */
void sim_flatfold::random_fold5(bool rand_sign) {
	/* Loop a few times in case two points are chosen 
	along the same edge and the fold is trivial. */
	for (int k=0; k<10; k++) {
		bool p1=true, p2=true;
		int j = 0;
		while (j<sim_flatfold_max_attempts && (p1||p2)) {
			if (p1) { find_ed_pts(p1x,p1y); p1 = point_inside(p1x,p1y); }
			if (p2) { find_ed_pts(p2x,p2y); p2 = point_inside(p2x,p2y); }
			j++;
		}
		if (p1||p2) { fputs("Too many attempts to find two boundary points.\n", stderr); exit(1); }
		if (random_flatfold5(rand_sign)) return;
	}
	fputs("Too many trivial flat fold attempts.\n", stderr); exit(1);
}

/** Applies a random flat fold to the sheet using two edge points.
* \param[in] rand_sign whether to choose a random sign for the fold or not.
* \return Whether the fold intersected the sheet or not.*/
bool sim_flatfold::random_flatfold5(bool rand_sign) {
	double magn = 1/sqrt((p1y-p2y)*(p1y-p2y) + (p2x-p1x)*(p2x-p1x));
	double nx = magn*(p1y-p2y), ny = magn * (p2x-p1x), di = p2x*nx+p2y*ny;
	return flatfold(nx, ny, di, rand_sign?random_sign():1);
}

/** Applies a fold by picking a random angle followed by a uniformly chosen displacement
* using the bounding circle.
* \param[in] rand_sign whether to choose a random sign for the fold or not. */
void sim_flatfold::random_fold4(bool rand_sign) {
	double th = 2*M_PI*gsl_rng_uniform(rng);
	double nx = cos(th);
	double ny = sin(th);

	// Choose an angle uniformly using the radius of the bounding circle.
	// Now we must check that the fold is non-trivial.
	for (int k=0; k<sim_flatfold_max_attempts; k++) {
		if (flatfold(nx,ny,cx*nx+cy*ny+cr*(-1+2*gsl_rng_uniform(rng)),rand_sign?random_sign():1)) return;
	}
	fputs("Too many fold attempts in fold4.\n",stderr);
	exit(1);
}

/** Applies a fold by choosing a random angle followed by 
* a random displacement.
* \param[in] rand_sign whether to choose a random sign for the fold or not. */
void sim_flatfold::random_fold3(bool rand_sign) {
	double th = 2*M_PI*gsl_rng_uniform(rng);
	double nx = cos(th);
	double ny = sin(th);
	// Find extremal vertices with respect to the normal vector 
	// Initializing values may need to be adjusted
	double min=10., max=-10.;
	double x,y,d;
	for (unsigned int i=0; i<f.size(); i++) {
		int k=0;
		do {
			x=f[i]->c.pts[2*k]; y=f[i]->c.pts[2*k+1];
			d=x*nx+y*ny;
			if (d<=min) min=d;
			else if (d>=max) max=d;
			k = f[i]->c.ed[2*k];
		} while (k!=0);
	}

	//for (int k=0; k<10; k++) {
		double di = min+(max-min)*gsl_rng_uniform(rng); // this displacement should make a valid fold. 
		if (flatfold(nx, ny, di, rand_sign?random_sign():1)) return; 
	//}
	fputs("random_fold3 unsuccessful\n", stderr);
	exit(1);
}

/** Computes the total area of the facets, which should remain constant
 * during folding.
 * \return The total area. */
double sim_flatfold::tot_area() {
    double t_area=f[0]->area();
	for(unsigned int k=1; k<f.size(); k++) t_area+=f[k]->area();
    return t_area;
}

/** Chooses a random point in a facet weighted via area and applies a fold there.
 * \param[in] t_area the precomputed total area of the sheet.
 * \param[in] rand_sign whether to choose a random sign for the fold or not. */
void sim_flatfold::random_fold2(double t_area,bool rand_sign) {

    // Select a random facet, weighted according to their area
	double x=t_area*gsl_rng_uniform(rng),cu=0;
    unsigned int k=0;
    while(k<f.size()-1) {
        cu+=f[k]->area();
        if(cu>x) break;
        k++;
    }
    facet *rf=f[k];

	// Find a random point in the chosen facet
	double lx,ux,ly,uy;
    rf->rect_bounds(lx,ux,ly,uy);
	for(int k=0; k<sim_flatfold_max_attempts; k++) {
        px=lx+(ux-lx)*gsl_rng_uniform(rng);
        py=ly+(uy-ly)*gsl_rng_uniform(rng);
		if(rf->point_inside(px, py)) {
			random_flatfold_point(rand_sign);
			return;
		}
	}
	fputs("Too many attempts to find a point in the facet in random_fold2\n", stderr);
	exit(1);
}


/** Applies a random flat fold to the sheet using the random point (px,py)
* and a random angle.
* \param[in] rand_sign whether to choose a random sign for the fold or not.
*/
void sim_flatfold::random_flatfold_point(bool rand_sign) {
	double th=2*M_PI*gsl_rng_uniform(rng);
	double nx=-sin(th);
	double ny=cos(th);
	if (ny<0) {nx=-nx; ny=-ny;}
	double di=px*nx+py*ny;
	flatfold(nx,ny,di,rand_sign?random_sign():1);
}

/** Chooses a random point, checks that it is on the sheet, and applies a fold.
* The option for radial folds has not yet been added.
* \param[in] rand_sign whether to choose a random sign for the fold or not. */
void sim_flatfold::random_fold1(bool rand_sign) {
	for (int k=0; k < sim_flatfold_max_attempts; k++) {
		double th_p = 2*M_PI*gsl_rng_uniform(rng);
		double r_p = cr*sqrt(gsl_rng_uniform(rng));
		px = r_p * cos(th_p) + cx;
		py = r_p * sin(th_p) + cy;
		if (point_inside(px, py)) break;
		if (k == sim_flatfold_max_attempts - 1) {
			fputs("Too many attempts to find a point in the sheet in random_fold1\n", stderr);
			exit(1);
		}
	}
	random_flatfold_point(rand_sign);
}

/** Applies a random flat fold to the sheet, choosing a random angle and
 * displacement.
 * \param[in] rand_sign whether to choose a random sign for the fold or not.
 * \return Whether the fold intersected the sheet or not. */
bool sim_flatfold::random_flatfold(bool rand_sign) {
	double th=2*M_PI*gsl_rng_uniform(rng);
	double nx=cos(th),ny=sin(th);
	return flatfold(nx,ny,cx*nx+cy*ny+cr*(-1+2*gsl_rng_uniform(rng)),
			rand_sign?random_sign():1);
}

/** Applies a random radial fold to the sheet.
 * \param[in] rand_sign whether to choose a random sign for the fold or not. */
void sim_flatfold::random_radial_fold(bool rand_sign) {
	int k=0;
	double x,y;
	while(k<sim_flatfold_max_attempts) {

		// Find a point inside the bounding circle
		x=-1+2*gsl_rng_uniform(rng);
		y=-1+2*gsl_rng_uniform(rng);
		if(x*x+y*y>crsq) continue;

		// Scale the point and check if it is inside the folded sheet.
		// If so, apply a random radial fold.
		x=cx+cr*x;
		y=cy+cr*y;
		if(point_inside(x,y)) {
			radial_fold(x,y,2*M_PI*gsl_rng_uniform(rng),
				    M_PI*gsl_rng_uniform_pos(rng),
				    M_PI*gsl_rng_uniform_pos(rng),
				    rand_sign?random_sign():1);
			return;
		}
		k++;
	}
	fputs("Too many attempts to find a point inside the sheet\n",stderr);
	exit(1);
}

/** Applies a flat fold to the sheet.
 * \param[in] (nx,ny) the normal vector of the fold to apply.
 * \param[in] di the displacement of the fold to apply.
 * \param[in] fsign the sign of the fold.
 * \return Whether the line given actually folded the sheet or not. */
bool sim_flatfold::flatfold(double nx,double ny,double di,int fsign) {
	unsigned int i=0,n=f.size();

	// Check to see if the fold intersects any facet. If not, then skip it
	for(i=0;i<n;i++) if(f[i]->intersects(nx,ny,di)) break;
	if(i==n) return false;

	// If the code reaches here, then we know that the fold intersects a
	// facet. Perform the folding by splitting the facet into two, and
	// mirroring one about the line of the fold.
	for(i=0;i<n;i++) {
		facet *nf=new facet(*f[i]);
		bool p1=f[i]->c.nplane(nx,ny,di,fsign);
		if(nf->c.nplane(-nx,-ny,-di,0)) {
			nf->mirror(nx,ny,di);

			// Update the bounding circle, since the mirrored facet
			// might extend further than before
			update_bounding_circle(nf);
			if(p1) f.push_back(nf);
			else {
				delete f[i];
				f[i]=nf;
			}
		} else {
			delete nf;
			if(!p1) {
				fputs("Error: both facets removed by fold\n",stderr);
				exit(1);
			}
		}
	}
	return true;
}

/** Applies a radial fold to the sheet.
 * \param[in] (x,y) the focus point of the radial fold.
 * \param[in] ro the overall rotation of the fold position.
 * \param[in] (al,be) the first two angles making up the radial fold.
 * \param[in] fsign the sign of the three folds that share the same sign. */
void sim_flatfold::radial_fold(double x,double y,double ro,double al,double be,int fsign) {
	unsigned int i=0,j,n=f.size();
	facet *qf[4];

	// Compute the normals to the four planes making up this radial fold
	double nxa=cos(ro),nya=sin(ro),nxb=cos(ro+al),nyb=sin(ro+al),
	       nxc=cos(ro+al+be),nyc=sin(ro+al+be),nxd=cos(ro+be+M_PI),nyd=sin(ro+be+M_PI),
	       da=nxa*x+nya*y,db=nxb*x+nyb*y,dc=nxc*x+nyc*y,dd=nxd*x+nyd*y;

	for(;i<n;i++) {
		j=0;
		facet *nbe=new facet(*f[i]),*nga=new facet(*f[i]),*nde=new facet(*f[i]);

		// Compute the first facet
		if(f[i]->c.nplane(-nxa,-nya,-da,fsign)&&f[i]->c.nplane(nxb,nyb,db,fsign)) qf[j++]=f[i];
		else delete f[i];

		// Compute the second facet
		if(nbe->c.nplane(-nxb,-nyb,-db,0)&&nbe->c.nplane(nxc,nyc,dc,0)) {
			nbe->mirror(nxb,nyb,db);
			update_bounding_circle(nbe);
			qf[j++]=nbe;
		} else delete nbe;

		// Compute the third facet
		if(nga->c.nplane(-nxc,-nyc,-dc,fsign)&&nga->c.nplane(nxd,nyd,dd,0)) {
			nga->mirror(nxc,nyc,dc);
			nga->mirror(nxb,nyb,db);
			update_bounding_circle(nga);
			qf[j++]=nga;
		} else delete nga;

		// Compute the fourth facet
		if(nde->c.nplane(-nxd,-nyd,-dd,-fsign)&&nde->c.nplane(nxa,nya,da,0)) {
			nde->mirror(nxa,nya,da);
			update_bounding_circle(nde);
			qf[j++]=nde;
		} else delete nde;

		// At least once facet should be present. If not, print an
		// error message
		if(j==0) {
			fputs("Error: all four facets removed by radial fold\n",stderr);
			exit(1);
		}

		// Store the available facets. One can be placed at the
		// location of the original facet. Others need to be added.
		f[i]=*qf;
		while(j>1) f.push_back(qf[--j]);
	}
}


/** Applies a random fold to the sheet.
 * \param[in] frac the probability of the fold being a radial fold.
 * \param[in] rand_sign whether to choose a random sign for the fold or not. */
void sim_flatfold::random_fold(double frac,bool rand_sign) {
	if(gsl_rng_uniform(rng)>frac) {
		for(int k=0;k<sim_flatfold_max_attempts;k++)
			if(random_flatfold(rand_sign)) return;
		fputs("Too many flatfold attempts in random_fold\n",stderr);
		exit(1);
	} else random_radial_fold(rand_sign);
}

/** Updates the bounding circle to ensure in encompasses a given facet.
 * \param[in] nf a pointer to the facet. */
void sim_flatfold::update_bounding_circle(facet *nf) {
	double nrsq=nf->max_rad_sq(cx,cy);
	if(nrsq>crsq) {crsq=nrsq;cr=sqrt(crsq);}
}

/** Checks to see if a given position is inside any facet.
 * \param[in] (x,y) the position to consider.
 * \return Whether the position is inside or not. */
bool sim_flatfold::point_inside(double x,double y) {
	for(unsigned int i=0;i<f.size();i++)
		if(f[i]->point_inside(x,y)) return true;
	return false;
}

/** Outputs the positions of creases in the original coordinate system.
 * \param[in] fp the file handle to write to.
 * \param[in] positive whether to write positive or negative creases. */
void sim_flatfold::crease_map(FILE *fp,bool positive) {
	for(unsigned int i=0;i<f.size();i++)
		f[i]->unfolded_crease(fp,positive);
}

/** Outputs the total crease mileage.
 * \param[out] pos the positive crease mileage.
 * \param[out] neg the negative crease mileage. */
void sim_flatfold::crease_mileage(double &pos,double &neg) {
	pos=neg=0;
	for(unsigned int i=0;i<f.size();i++)
		f[i]->crease_mileage(pos,neg);
}

/** Outputs the facet positions in the current coordinate system.
 * \param[in] fp the file handle to write to. */
void sim_flatfold::output(FILE *fp) {
	for(unsigned int i=0;i<f.size();i++) f[i]->output(fp);
}

/** Updates the vertices of the sheet after a fold. */
//void sim_flatfold::update_boundary(){
	// Define edges of sheet pre-fold outside of here?
	// Define the fold as a line segment
	// Add vertices using p1x,p1y,p2x,p2y and where an old edge intersects the fold
	// Loop through vertices and check if they are in the sheet
//}
