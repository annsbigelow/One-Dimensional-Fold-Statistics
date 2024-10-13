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

    // CHR: pavg, normx, and normy only exist within this function so you
    // declare them locally instead of as class members. In addition, you know
    // their size exactly, so this is a case where the extensibility
    // of vectors isn't needed. 

    double* pavg = new double[numplaces]; // allocates memory for numplaces doubles
    double* normx = new double[numplaces];
    double* normy = new double[numplaces];
    double* range = new double[numplaces];
    vector<double> domain;
    for (int i = 0; i < numplaces; i++) { pavg[i] = 0; normx[i] = 0; normy[i] = 0; range[i] = 0; }

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

        // scan across final segment
        double fac=numplaces/(hi-lo);
        for (int i = 0; i < segs_o.size(); i += 2) {

            // CHR: In general, I think that this binning procedure works. But
            // it seems to assume that all of your segments are properly aligned,
            // so segs_o[i]<segs_o[i+1]. If it's the other way round, then the
            // segment won't get binned anywhere.
            //
            // Another issue here is that you have to do a for-loop over all
            // bins. It can be done more efficiently than this.
           /* for (int m = 0; m < numplaces; m++) {
                if (segs_o[i] <= domain[m] && domain[m] <= segs_o[i + 1]) {
                    range[m] += 1;
                }
            }*/

            //// CHR: Here's a more efficient way to do it, without looping over all
            //// of the data. I'm also assuming that the segments might not have
            //if(false) {
            //    int mlo=int((segs_o[i]-lo)*fac),mhi=int((segs_o[i+1]-lo)*fac);
            //    if(mhi>mlo) {int o=mlo;mlo=mhi;mhi=o;}
            //    if(mlo<0) mlo=0;else if(mlo>numplaces) mlo=numplaces;
            //    if(mhi<0) mhi=0;else if(mhi>numplaces) mhi=numplaces;
            //    for(int m=mlo;m<mhi;m++) range[m]+=1;
            //}

            // CHR: Here's an even better way to do it. The previous approaches
            // only bin segments when they specific values (given by the
            // domain[m] values in the original version). But we might have
            // many segments that lie wholly between domain[m] and domain[m+1].
            // Currently they'd just be omitted. We could improve this by
            // binning them according to the fraction of the range between
            // domain[m] and domain[m+1] that they cover.
           
                double x=(segs_o[i]-lo)*fac,y=(segs_o[i+1]-lo)*fac,z;
                if(x>y) {double z=x;x=y;y=z;}
                int mlo=int(x),mhi=int(y);
                if(mlo<0) mlo=0;else if(mlo>=numplaces) mlo=numplaces-1;
                if(mhi<0) mhi=0;else if(mhi>=numplaces) mhi=numplaces-1; // Ann: aren't mlo/mhi never less than 0? 
                range[mlo]-=x-mlo;   // Account for first bin only being partially covered
                for(int m=mlo;m<=mhi;m++) range[m]+=1;
                range[mhi]-=mhi+1-y;  // Account for last bin only being partially covered
            
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

    densData.open("NormSegDens_40new.txt");
    double h = 1.0 / (double)numplaces;
    for (int j = 0; j < numplaces; j++) {
        //domavg[j] /= instance;
        //
        // CHR: It's unclear to me why hi and lo enter in here. They will be set by
        // the final instance. I don't see why those final values would be relevant
        // for the outputting this data.
        //normx[j] = (domain[j] - lo) / (hi - lo); // normalize x-axis
        pavg[j] /= instance;
       // densData << normx[j] << ' ' << pavg[j] << '\n';

        // CHR: For the third method, if you have spacing h=1/numplaces, each bin covers
        // [j*h,(j+1)*h]. So it would make sense to use an x coordinate halfway
        // between them, like
         densData << (j+0.5)*h << ' ' << pavg[j] << '\n';
    } 
    densData.close();
    delete[] pavg; // free memory for arrays
    delete[] normx;
    delete[] normy;
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


