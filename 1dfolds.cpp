#include <iostream>
using namespace std;
#include "declarations.h"
#include <omp.h>

int main() {
    // User input
    int seed, inst, n;
    string filename;
    std::cout << "Type any number to seed the RNG:";
    std::cin >> seed;
    std::cout << "Type the number of instances for averaging:";
    std::cin >> inst;
    std::cout << "Type the number of folds:";
    std::cin >> n;
    std::cout << "Type the filename to write data to:";
    std::cin >> filename;
        folds_stats fs(seed,inst,n); // if you want new sequence of random nums, change fs() seed input.

    //fs.segdens(s_divide,"SegDens_5new.txt");
        double t0 = omp_get_wtime();
        fs.altFold(filename);
        printf("Total time elapsed is %g s\n", omp_get_wtime() - t0);

    return 0;
}
