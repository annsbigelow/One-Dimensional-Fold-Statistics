#include <iostream>
using namespace std;
#include <vector>

#ifndef folds_stats_h
#define folds_stats_h

class folds_stats {
	public: 
		// Initializing variables that shouldn't be re-initialized in loops.
		int direct, size;
		float LO, HI, seg_l, seg_r, x, min, max;

		vector<float> logavg();
		vector<float> fold(vector<float> segs_in);
		float getmin(vector<float> segs_in);
		float getmax(vector<float> segs_in);
		void display(vector<float> arr);
};

#endif