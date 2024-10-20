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

void folds_stats::altFold_segdens(int& numplaces){

    double* pavg = new double[numplaces];
    double* range = new double[numplaces];
    for (int i = 0; i < numplaces; i++) pavg[i] = 0;

    for (int w = 0; w < instance; w++) {
        segs_o = altFold();
        lo = getmin(segs_o);
        hi = getmax(segs_o);

        // Clear range array; this needs to be done each time
        for(int i=0;i<numplaces;i++) range[i]=0;

        // scan across final segment
        double fac = numplaces / (hi - lo);
        for (unsigned int i = 0; i < segs_o.size(); i += 2) {

            double x = (segs_o[i] - lo) * fac, y = (segs_o[i + 1] - lo) * fac;
            if (x > y) { double z = x; x = y; y = z; }
            int mlo = int(x), mhi = int(y);
            if (mlo < 0) mlo = 0; else if (mlo >= numplaces) mlo = numplaces - 1;
            if (mhi < 0) mhi = 0; else if (mhi >= numplaces) mhi = numplaces - 1;
            range[mlo] -= x - mlo;   // Account for first bin only being partially covered
            for (int m = mlo; m < mhi; m++) range[m] += 1;
            range[mhi] -= mhi - y;  // Account for last bin only being partially covered
        }

        // Compute the normalizing factor
        double sum=0;
        for (int j = 0; j < numplaces; j++) sum += range[j];

        // Add the normalized distribution to the running total
        double nor=numplaces/sum;
        for (int i = 0; i < numplaces; i++) pavg[i] += range[i]*nor;
    }

    densData.open("Altfold_segdens_n10.txt");
    double h = 1.0 / (double)numplaces;
    for (int j = 0; j < numplaces; j++) {
        pavg[j] /= instance;
        densData << (j + 0.5) * h << ' ' << pavg[j] << '\n';
    }
    densData.close();
    delete[] pavg;
    delete[] range;
}

void folds_stats::segdens(int &numplaces) {
    double *pavg = new double[numplaces],
           *range = new double[numplaces];
    for (int i = 0; i < numplaces; i++) pavg[i] = 0;

    for (int w = 0; w < instance; w++) {
        //cout << omp_get_thread_num() << '\n';
        segs_i = { 0,1 };
        for (int i = 0; i < n; i++) {
            segs_o = fold(segs_i);
            segs_i = segs_o;
        }

        lo = getmin(segs_o);
        hi = getmax(segs_o);

        // Clear range array; this needs to be done each time
        for(int i=0;i<numplaces;i++) range[i]=0;

        // scan across final segment
        double fac=numplaces/(hi-lo);
        for (unsigned int i = 0; i < segs_o.size(); i += 2) {

            double x=(segs_o[i]-lo)*fac,y=(segs_o[i+1]-lo)*fac;
            if(x>y) {double z=x;x=y;y=z;}
            int mlo=int(x),mhi=int(y);
            if(mlo<0) mlo=0;else if(mlo>=numplaces) mlo=numplaces-1;
            if(mhi<0) mhi=0;else if(mhi>=numplaces) mhi=numplaces-1; // Ann: aren't mlo/mhi never less than 0?
            range[mlo]-=x-mlo;   // Account for first bin only being partially covered
            for(int m=mlo;m<mhi;m++) range[m]+=1;
            range[mhi]-=mhi-y;  // Account for last bin only being partially covered

        }

        // Compute the normalizing factor
        double sum=0;
        for (int j = 0; j < numplaces; j++) sum += range[j];

        // Add the normalized distribution to the running total
        double nor=numplaces/sum;
        for (int i = 0; i < numplaces; i++) pavg[i] += range[i]*nor;
    }

    densData.open("NormSegDens_12new.txt");
    double h = 1.0 / (double)numplaces;
    for (int j = 0; j < numplaces; j++)
        densData << (j+0.5)*h << ' ' << pavg[j]/instance << '\n';
    densData.close();
    delete[] pavg; // free memory for arrays
    delete[] range;
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
    //long int k = static_cast<long int> (time(NULL));
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

vector<double> folds_stats::altFold() {
    vector<double> segs_out;
    vector<double> segsin = { 0,1 };
    vector<double> q;
    vector<double> qmap;
    bool* dirvec = new bool[n],dir;

    for (int k = 0; k < n; k++) {
        q.push_back(gsl_rng_uniform(rng));
        dir = rand_bit();
        dirvec[k] = dir; // Store the directions

        // Map q[k]
        // q.size() is the number of folds that have already happened
        qmap.push_back(q[k]);
        for (unsigned int m = 0; m < q.size()-1; m++) {
            if (dirvec[m]) { // right fold
                if (qmap[q.size() - 1] < qmap[m]) { qmap[q.size() - 1] = 2 * qmap[m] - qmap[q.size() - 1]; }
            } else {
                if (qmap[m] < qmap[q.size() - 1]) {qmap[q.size() - 1] = 2*qmap[m] - qmap[q.size() - 1]; }
            }
        }
        LO = getmin(segsin);
        HI = getmax(segsin);
        x = qmap[k];
        // check if x in outermost segment
        if (LO < x && x < HI) {
            // loop through each segment
            size = segsin.size();
            for (long i = 0; i < size; i += 2) {
                // identify i'th segment
                double seg_l = segsin[i],seg_r = segsin[i + 1];
                // right fold
                if (dir) {
                    if (seg_l < x && x < seg_r) {
                        segs_out.push_back(x);
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
                else {
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
            segs_out = segsin;
        }

        segsin = segs_out;
    }
    return segs_out;
    delete[] dirvec;
}

vector<double> folds_stats::fold(vector<double> &segs_in) {
    vector<double> segs_out;

    // generates random fold pos.
    LO = getmin(segs_in);
    HI = getmax(segs_in);
    bool dir=rand_bit();

    x = (HI - LO)*gsl_rng_uniform(rng) + LO;
    // NOTE: UNIFORM RANGE INCLUDES 0.0 BUT EXCLUDES 1.0

    // check if x in outermost segment
    // CHR: I don't see the need for this check, since x is chosen to be
    // between LO and HI
    if (LO < x && x < HI) {

        // loop through each segment
        size = segs_in.size();
        for (long i = 0; i < size; i += 2) {

            // identify i'th segment
            double seg_l = segs_in[i],seg_r = segs_in[i + 1];

            // right fold
            // CHR: I cleaned up this section to require fewer comparisons
            if (dir) {
                if (x > seg_r) {
                    segs_out.push_back((2 * x) - seg_r);
                    segs_out.push_back((2 * x) - seg_l);
                } else if (x > seg_l) {
                    segs_out.push_back(x);
                    segs_out.push_back((2 * x) - seg_l);
                    segs_out.push_back(x);
                    segs_out.push_back(seg_r);
                } else {
                    segs_out.push_back(seg_l);
                    segs_out.push_back(seg_r);
                }
            }
            // left fold
            else {
                if (x > seg_r) {
                    segs_out.push_back(seg_l);
                    segs_out.push_back(seg_r);
                } else if(x > seg_l) {
                    segs_out.push_back(seg_l);
                    segs_out.push_back(x);
                    segs_out.push_back((2 * x) - seg_r);
                    segs_out.push_back(x);
                } else {
                    segs_out.push_back((2 * x) - seg_r);
                    segs_out.push_back((2 * x) - seg_l);
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
    for (unsigned long i = 0; i < arr.size(); i++) {
        cout << arr[i] << '\n';
    }
}

// Similar to numpy.linspace()
vector<double> folds_stats:: linspace(const double &start, const double &end, int &num) {
    vector<double> linpoints;
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
