#include <iostream>
using namespace std;
#include <vector>
#include <ctime>
#include "declarations.h"


int main() {
    folds_stats fs; 
    // for random num generating. Only needs to be called once.
    srand(static_cast <unsigned> (time(0)));
    
    vector<float> cavg = fs.logavg();
    fs.display(cavg); 
    
    return 0;
}


