#include <iostream>
using namespace std;
#include <algorithm>
#include <fstream>
#include <vector>
#include <ctime>
#include <cmath>
#include <omp.h>
#include "declarations.h"
#include <gsl/gsl_rng.h>


void folds_stats::segdens(int &numplaces) {
    pavg.clear();
    normx.clear();
    normy.clear();
    for (int i = 0; i < numplaces; i++) {
        pavg.push_back(0);
        normx.push_back(0);
        normy.push_back(0);
    }

    for (int w = 0; w < instance; w++) {
        //cout << omp_get_thread_num() << '\n';
        segs_i = { 0,1 };
        for (int i = 0; i < n; i++) {
            segs_o = fold(segs_i);
            segs_i = segs_o;
        }

        lo = getmin(segs_o);
        hi = getmax(segs_o);
        domain = linspace(lo, hi, numplaces); // each step is 1/numplaces
        range.clear();
        for (int i = 0; i < numplaces; i++) {
            range.push_back(0);
        }

        // scan across final segment
        for (int i = 0; i < segs_o.size(); i += 2) {
            for (int m = 0; m < numplaces; m++) {
                if (segs_o[i] <= domain[m] && domain[m] <= segs_o[i + 1]) {
                    range[m] += 1;
                }
            }
        }
        for (int j = 0; j < numplaces; j++) {
            sum += range[j]; // the un-normalized integral
        }
        sum /= numplaces;
        for (int i = 0; i < numplaces; i++) {
            //domavg[i] += domain[i];
            normy[i] = range[i] / sum; // normalize y-axis to integrate to 1.
            pavg[i] += normy[i];
        }
        //if (w == 0) { cout << gsl_rng_uniform(rng) << '\n'; }
        sum = 0;
    }

    /*char buf[128];
    sprintf(buf,"NormSegDens_testomp%d.txt",th);
    densData.open(buf);*/

    densData.open("NormSegDens_omptest.txt");
    for (int j = 0; j < numplaces; j++) {
        //domavg[j] /= instance;
        normx[j] = (domain[j] - lo) / (hi - lo); // normalize x-axis
        pavg[j] /= instance;
        densData << normx[j] << ' ' << pavg[j] << '\n';
    } 
    densData.close();
}

void folds_stats:: log_fixedn(){
    // Essentially the same as logavg(), but now we print logged crease vals at fixed n 
    c.clear();
    logc.clear();

    for (int i = 0; i < n; i++) {
        c.push_back(0);
        logc.push_back(0);
    }
 
    ofstream data;
    data.open("logc.txt");


    for (int j = 1; j <= instance; j++) {
        segs_i = { 0,1 };
        for (int i = 0; i < n; i++) {
            segs_o = fold(segs_i);
            sizeo = segs_o.size();
            c[i] = (sizeo / 2) - 1;
            logc[i] = log(c[i]);
            if (c[i] == 0.0) {
                c[i] = 1;
                logc[i] = 0;
            }
            segs_i = segs_o;
        }
        for (int i = 0; i < n; i++) {
            data << logc[n-1] << '\n';
        }
    }
    data.close();
}

void folds_stats::logavg() {
    // Make sure vecs are sized n 
    f.clear();
    c.clear(); 
    logc.clear();
    cavg.clear();

    // fold count is just a vector of integers increasing by 1.
    for (int i = 0; i < n ; i++) {
        f.push_back(i + 1);
    }

    for (int i = 0; i < n ; i++) {
        cavg.push_back(0);
        c.push_back(0);
        logc.push_back(0);
    } 
    long int k = static_cast<long int> (time(NULL)); 
    // Loop (instance) times to average log(crease values) for accurate statistics. 
    for (int j = 1; j <= instance; j++) {
        segs_i = { 0,1 };
        // Fold n times. Count creases and take their logs
        for (int i = 0; i < n; i++) {
            segs_o = fold(segs_i);
            //display(segs_o);
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
    // Print logc file 
    ofstream data;
    data.open("test.txt");
    // Average the logs of crease values 
    for (int i = 0; i < n ; i++) {
        cavg[i] /= instance;
        data << f[i] << ' ' << cavg[i] << '\n';
    }
    data.close();
    display(c); 
}

vector<double> folds_stats::fold(vector<double> &segs_in) {
    vector<double> segs_out;  // Include in class construction ?? If you do, you get very large crease values
    // 1 for fold direction left; 0 for right.
    direct = rand() % 2; // convert to GSL RNG?
   
    // generates random fold pos.
    LO = getmin(segs_in);
    HI = getmax(segs_in);

    x = (HI - LO)*gsl_rng_uniform(rng) + LO;
    // NOTE: UNIFORM RANGE INCLUDES 0.0 BUT EXCLUDES 1.0

    // check if x in outermost segment
    if (LO < x && x < HI) {

        // loop through each segment
        size = segs_in.size();
        for (long i = 0; i < size; i += 2) {

            // identify i'th segment
            seg_l = segs_in[i];
            seg_r = segs_in[i + 1];

            // right fold
            if (direct == 0) {
                if (seg_l < x && x < seg_r) {
                    segs_out.push_back(x); // IS segs_out AN EMPTY VEC AT FIRST?
                    segs_out.push_back((2 * x) - seg_l);
                    segs_out.push_back(x);
                    segs_out.push_back(seg_r);
                }
                else if (x > seg_r) {
                    segs_out.push_back((2 * x) - seg_r);
                    segs_out.push_back((2 * x) - seg_l);
                }
                else if (x <= seg_l) {
                    segs_out.push_back(seg_l);
                    segs_out.push_back(seg_r);
                }
                else if (x == seg_r) {
                    segs_out.push_back(x);
                    segs_out.push_back((2 * x) - seg_l);
                }
            }
            // left fold
            else if (direct == 1) {
                if (seg_l < x && x < seg_r) {
                    segs_out.push_back(seg_l);
                    segs_out.push_back(x);
                    segs_out.push_back((2 * x) - seg_r);
                    segs_out.push_back(x);
                }
                else if (x < seg_l) {
                    segs_out.push_back((2 * x) - seg_r);
                    segs_out.push_back((2 * x) - seg_l);
                }
                else if (x >= seg_r) {
                    segs_out.push_back(seg_l);
                    segs_out.push_back(seg_r);
                }
                else if (x == seg_l) {
                    segs_out.push_back((2 * x) - seg_r);
                    segs_out.push_back(x);
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
    for (long i = 1; i < size; i++) {
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
    for (long i = 1; i < size; i++) {
        if (segs_in[i] > max) {
            max = segs_in[i];
        }
    }
    return max;
}


// Display vec function. For debug purposes only
void folds_stats::display(vector<double> arr) {
    for (long i = 0; i < arr.size(); i++) {
        cout << arr[i] << '\n';
    }
}

// Similar to numpy.linspace()
vector<double> folds_stats:: linspace(const double &start, const double &end, int &num) {
    linpoints.clear();
    if (num == 0) return linpoints;
    if (num == 1) {
        linpoints.push_back(start);
        return linpoints;
    }
    step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        linpoints.push_back(start + i * step);
    }
    return linpoints;
}


