#include <iostream>
using namespace std;
#include "declarations.h"
#include <omp.h>

int main() {
    // User input
    int seed, inst, n;
    string filename;
    std::cout << "Input any number to seed the RNG:";
    std::cin >> seed;
    std::cout << "Input the number of instances for averaging:";
    std::cin >> inst;
    std::cout << "Input the number of folds:";
    std::cin >> n;
    std::cout << "Input the filename to write data to:";
    std::cin >> filename;
        
	folds_stats fs(seed,inst,n); // if you want new sequence of random nums, change fs() seed input.
	double t0 = omp_get_wtime();

	// Segment Density
	// Specify the bins spacing to find the density of facets. Usually we use ~200.
	int s_divide = 200;
    fs.segdens(s_divide,filename);
         
	// Test fold protocols
      //  fs.altFold(filename);
        printf("Total time elapsed is %g s\n", omp_get_wtime() - t0);

    return 0;
}
