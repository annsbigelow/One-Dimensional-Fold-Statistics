#include <iostream>
using namespace std;
#include <vector>

#ifndef folds_stats_h
#define folds_stats_h

#include <gsl/gsl_rng.h>

class folds_stats {
	public: 
		// Initializing variables 
		int direct;
		long size, sizeo;
		double LO, HI, seg_l, seg_r, x, min, max;
		vector<double> segs_i, segs_o, c, logc, cavg;
		//vector < double> segs_out; // Initializing value?? 
		vector<int> f;
		// instance = number of iterations to take avg over (TEST GSL)
		const int instance = 4;
		// n = number of folds (TEST GSL)
		const int n = 5;

		// Class constructor "dynamically allocates memory for initial number of slots." Written as an initializer list: same as writing rng  = gsl_rng_alloc(...)
		folds_stats() : direct(0), HI(1), LO(0), seg_l(0), seg_r(1), x(0), min(0), max(1), segs_i({ 0 }), segs_o({ 0 }), c({ 0 }), logc({ 0 }), cavg({ 0 }),
			size(1), sizeo(1), rng(gsl_rng_alloc(gsl_rng_taus2)) {} // Creates instance of the Tausworthe generator
		// initializing variables without values is sketchy. If you do, the compiler gives these some val anyway or gives error
		
		// Class destructor " frees dynamically allocated memory " 
		~folds_stats() {
			gsl_rng_free(rng);
		}

		inline void seed(long int k) {
			gsl_rng_set(rng, k);
		}

		void logavg();
		vector<double> fold(vector<double> &segs_in);
		double getmin(vector<double> &segs_in);
		double getmax(vector<double> &segs_in);
		void display(vector<double> arr);

private: // Only accessible in class
	gsl_rng* rng;
};

#endif