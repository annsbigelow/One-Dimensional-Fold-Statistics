#include <iostream>
using namespace std;
#include <vector>

#ifndef folds_stats_h
#define folds_stats_h

#include "gsl/gsl_rng.h"

class folds_stats {
	public: 
		// Initializing variables 
		int direct;
		long size, sizeo;
		double LO, HI, seg_l, seg_r, x, min, max;
		vector<double> segs_i, segs_o, c, logc, cavg; 
		vector<int> f;
		// instance = number of iterations to take avg over (TEST GSL)
		const int instance = 4;
		// n = number of folds (TEST GSL)
		const int n = 5;

		// Class constructor " dynamically allocates memory for initial number of slots " ( NOTATION ?? )
		folds_stats() : rng(gsl_rng_alloc(gsl_rng_taus2)) {} 
		
		// Class destructor " frees dynamically allocated memory " 
		~folds_stats() {
			gsl_rng_free(rng);
		}

		inline void seed(int k) {
			gsl_rng_set(rng, k);
		}

		void logavg();
		vector<double> fold(vector<double> &segs_in);
		double getmin(vector<double> &segs_in);
		double getmax(vector<double> &segs_in);
		void display(vector<double> arr);

private:
	gsl_rng* rng;
};

#endif