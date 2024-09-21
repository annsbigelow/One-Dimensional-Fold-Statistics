#include <iostream>
using namespace std;
#include <vector>

#ifndef folds_stats_h
#define folds_stats_h

class folds_stats {
	public: 
		// Initializing variables 
		int direct;
		long size, sizeo;
		double LO, HI, seg_l, seg_r, x, min, max;
		vector<double> segs_i, segs_o, c, logc, cavg; 
		vector<int> f;
		// instance = number of iterations to take avg over
		const int instance = 500;
		// n = number of folds
		const int n = 40;


		void logavg();
		vector<double> fold(vector<double> &segs_in);
		double getmin(vector<double> &segs_in);
		double getmax(vector<double> &segs_in);
		void display(vector<double> arr);
};

#endif