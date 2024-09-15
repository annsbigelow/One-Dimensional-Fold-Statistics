#include <iostream>
using namespace std;
#include <algorithm>
#include <vector>
#include "declarations.h"
#include <cstdlib>


int main() {
    folds_stats fs; 
    // for random num generating. Only needs to be called once.
    srand(static_cast <unsigned> (time(0)));
    
    vector<float> cavg = fs.logavg();
    fs.display(cavg); 
    
    return 0;
}


