#include <iostream>
using namespace std;
#include <vector>
#include <fstream>

#ifndef folds_stats_h
#define folds_stats_h

#include <gsl/gsl_rng.h>

class folds_stats {
	public: 
		// Initializing variables
        const int th;
		int direct;
		long size, sizeo;
		double LO, HI, seg_l, seg_r, x, min, max, step, lo, hi, sum;
		vector<double> segs_i, segs_o, c, logc, cavg, linpoints, domain, range,pavg,normx, normy;
		//vector < double> segs_out; // Initializing value?? 
		vector<int> f;
		// instance = number of iterations to take avg over
		const int instance = 10000;
		// n = number of folds 
		const int n = 47;

		// Class constructor 
		// seeds rng upon class instance creation
		folds_stats(int th_=1) : th(th_), direct(0), HI(1), LO(0), seg_l(0), seg_r(1), x(0), min(0), max(1), step(0), lo(0), hi(0), sum(0), segs_i({ 0 }), segs_o({ 0 }), c({ 0 }), logc({ 0 }), cavg({ 0 }),
			linpoints({ 0 }), domain({ 0 }), range({ 0 }), pavg({ 0 }), normx({ 0 }), normy({ 0 }), size(1), sizeo(1), rng(gsl_rng_alloc(gsl_rng_taus2))
		{seed(th);}
		 
		// Class destructor
		~folds_stats() {
			gsl_rng_free(rng);
		}

		inline void seed(long int k) {
			gsl_rng_set(rng, k);
		}

		void segdens(int &numplaces);
		void log_fixedn();
		void logavg();
		vector<double> fold(vector<double> &segs_in);
		double getmin(vector<double> &segs_in);
		double getmax(vector<double> &segs_in);
		void display(vector<double> arr);
		vector<double> linspace(const double& start, const double& end, int& num);

		ofstream data;
		ofstream densData;
		ofstream densData1;

//private: // Only accessible in class
	gsl_rng* rng;
};

#endif
