#include <iostream>
using namespace std;
#include <vector>

#ifndef folds_stats_h
#define folds_stats_h

class folds_stats {
	public: 
		// Initializing variables 
		int direct, size;
		float LO, HI, seg_l, seg_r, x, min, max, sizeo;
		vector<float> segs_i, segs_o, c, logc, cavg, f;
		const int instance = 350;
		// n = number of folds
		const int n = 26;


		void logavg();
		vector<float> fold(vector<float> &segs_in);
		float getmin(vector<float> &segs_in);
		float getmax(vector<float> &segs_in);
		void display(vector<float> arr);
};

#endif