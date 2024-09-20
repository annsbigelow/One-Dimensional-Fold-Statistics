#include <iostream>
using namespace std;
#include <algorithm>
#include <fstream>
#include <vector>
#include <ctime>
#include <cmath>
#include "declarations.h"


void folds_stats::logavg() {
    // Seed the RNG
    srand(static_cast <unsigned> (time(0)));

    // fold count is just a vector of integers increasing by 1.
    for (int i = 0; i < n; i++) {
        f.insert(f.end(), {i + 1});
    }

    for (int i = 0; i < n; i++) {
        cavg.insert(cavg.end(), { 0 });
        c.insert(c.end(), { 0 });
        logc.insert(logc.end(), { 0 });
    } 

    // Loop (instance) times to average log(crease values) for accurate statistics. 
    for (int j = 1; j <= instance; j++) {
        segs_i = { 0,1 };
        // Fold n times. Count creases and take their logs
        for (int i = 0; i < n; i++) {
            segs_o = fold(segs_i);
            sizeo = segs_o.size();
            // Number of creases is number of segments - 1
            c[i] = (sizeo / 2) - 1;
            logc[i] = log(c[i]);
            cavg[i] += logc[i];
            // Avoid the case when no fold happened the first time. Then log(0) = -inf
            if (c[i] == 0.0) {
                c[i] = 1;
                logc[i] = 0;
                cavg[i] = 0; 
            }
            segs_i = segs_o;
        }
    }
    // Print to file to be used in gnuplot
    ofstream data;
    data.open("data.txt");
    // Average the logs of crease values 
    for (int i = 0; i < n; i++) {
        cavg[i] /= instance;
        data << f[i] << ' ' << cavg[i] << '\n';
    }
    data.close();
}
 
vector<double> folds_stats::fold(vector<double> &segs_in) {
    vector<double> segs_out;
    // 1 for fold direction left; 0 for right.
    direct = rand() % 2;
   
    // generates random fold pos.
    LO = getmin(segs_in);
    HI = getmax(segs_in);
    x = LO + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (HI - LO)));

    // check if x in outermost segment
    if (LO < x && x < HI) {

        // loop through each segment
        size = segs_in.size();
        for (int i = 0; i < size; i += 2) {

            // identify i'th segment
            seg_l = segs_in[i];
            seg_r = segs_in[i + 1];

            // right fold
            if (direct == 0) {
                if (seg_l < x && x < seg_r) {
                    segs_out.insert(segs_out.end(), { x, (2 * x) - seg_l });
                    segs_out.insert(segs_out.end(), { x, seg_r });
                }
                else if (x > seg_r) {
                    segs_out.insert(segs_out.end(), { (2 * x) - seg_r, (2 * x) - seg_l });
                }
                else if (x <= seg_l) {
                    segs_out.insert(segs_out.end(), { seg_l, seg_r });
                }
                else if (x == seg_r) {
                    segs_out.insert(segs_out.end(), { x, (2 * x) - seg_l });
                }
            }
            // left fold
            else if (direct == 1) {
                if (seg_l < x && x < seg_r) {
                    segs_out.insert(segs_out.end(), { seg_l, x });
                    segs_out.insert(segs_out.end(), { (2 * x) - seg_r, x });
                }
                else if (x < seg_l) {
                    segs_out.insert(segs_out.end(), { (2 * x) - seg_r, (2 * x) - seg_l });
                }
                else if (x >= seg_r) {
                    segs_out.insert(segs_out.end(), { seg_l, seg_r });
                }
                else if (x == seg_l) {
                    segs_out.insert(segs_out.end(), { (2 * x) - seg_r, x });
                }
            }
        }
    }
    // else: don't fold at outermost endpoints. Then crease # does not change
    else {
        segs_out = segs_in;
    }
    return segs_out;
}

// Find min element at beginning of fold() 
double folds_stats::getmin(vector<double> &segs_in) {
    size = segs_in.size();
    min = segs_in[0];
    for (int i = 1; i < size; i++) {
        if (segs_in[i] < min) {
            min = segs_in[i];
        }
    }
    return min;
}

// Find max element at beginning of fold()
double folds_stats::getmax(vector<double> &segs_in) {
    size = segs_in.size();
    max = segs_in[0];
    for (int i = 1; i < size; i++) {
        if (segs_in[i] > max) {
            max = segs_in[i];
        }
    }
    return max;
}


// Display vec function. For debug purposes only
void folds_stats::display(vector<double> arr) {
    for (int i = 0; i < arr.size(); i++) {
        cout << arr[i] << '\n';
    }
}


