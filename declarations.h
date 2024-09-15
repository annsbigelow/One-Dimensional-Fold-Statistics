#include <iostream>
using namespace std;

#ifndef folds_stats_h
#define folds_stats_h

class folds_stats {
	public: 
		// Initializing variables that shouldn't be re-initialized in loops.
		vector<float> segs_in, segs_out;
		int direct, size;
		float LO, HI, seg_l, seg_r, x, min, max;

		vector<float> logavg();
		vector<float> fold();
		float getmin();
		float getmax();
		void display(vector<float> arr);
};

#endif