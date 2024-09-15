#include <iostream>
using namespace std;
#include <algorithm>
#include <vector>
#include <ctime>
#include "declarations.h"


vector<float> folds_stats::logavg() {

    // Initialize values for later to prevent re-initializing
    // n = number of folds
    int n = 25;
    float sizeout; 

    // fold count is just a vector of integers increasing by 1
    vector<int> f;
    for (int i = 0; i < n; i++) {
        f.insert(f.end(), {i + 1});
    }

    // Initialize vecs: creases, log(creases), average(log(creases))
    vector<float> c;
    vector<float> logc;
    vector<float> cavg;
    for (int i = 0; i < n; i++) {
        cavg.insert(cavg.end(), { 0 });
        c.insert(c.end(), { 0 });
        logc.insert(logc.end(), { 0 });
    } 

    // loop 4 times to average log(crease values) for accurate statistics.
    int instance = 4; 
    for (int j = 1; j <= instance; j++) {
        segs_in = { 0,1 };
        // Fold 25 times. Count creases and take their logs
        for (int i = 0; i < n; i++) {
            segs_out = fold();
            sizeout = segs_out.size();
            // number of creases is number of segments - 1
            c[i] = (sizeout / 2) - 1;
            logc[i] = log(c[i]);
            cavg[i] += logc[i];
            segs_in = segs_out;
        }
    }
    // Average the logs of crease values
    for (int i = 0; i < n; i++) {
        cavg[i] /= instance;
    }
    return cavg;
}

// fold function.  
vector<float> folds_stats::fold() {
    
    // 1 for fold direction left; 0 for right. rand() is seeded in main()
    direct = rand() % 2;
   
    // generates random fold pos.
    LO = getmin();
    HI = getmax();
    x = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));

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
float folds_stats::getmin() {
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
float folds_stats::getmax() {
    size = segs_in.size();
    max = segs_in[0];
    for (int i = 1; i < size; i++) {
        if (segs_in[i] > max) {
            max = segs_in[i];
        }
    }
    return max;
}


// display vec function. For debug purposes only
void folds_stats::display(vector<float> arr) {
    for (int i = 0; i < arr.size(); i++) {
        cout << arr[i] << '\n';
    }
}

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

