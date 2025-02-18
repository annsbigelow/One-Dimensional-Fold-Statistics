#include "sim_flatfold.hh"

#include <cstdlib>
#include <cmath>

int main(int argc,char **argv) {

	// Check for the correct number of command-line arguments
	if(argc<2||argc>6) {
		fputs("Syntax: ./flatfold_gen <num_folds> [<seed>] [<percent_radial>] [<sign>] [<fold_option>]\n\n"
		      "sign=0 for positive folds, sign=1 for folds with random sign\n"
			"fold_option=0 for standard random fold, fold_option=1 for an updated protocol\n",stderr);
		return 1;
	}

	// Check for a sensible number of folds
	int folds=atoi(argv[1]),i=0;
	if(folds<0||folds>100) {
		fputs("Fold number out of bounds\n",stderr);
		return 1;
	}

	// Check for optional arguments
	unsigned long seed=1;
	double frac=0.;
	bool rand_sign=false;
	int fold_option=1;
	if(argc>2) {

		// Read the random seed
		seed=atol(argv[2]);
		if(argc>3) {

			// Read the fraction of radial folds
			frac=0.01*atof(argv[3]);
			if(frac<0) frac=0;
			else if(frac>1) frac=1;

			// Read whether to use random signs on the folds
			if(argc>4) {
				rand_sign=atoi(argv[4])==1;
				
				// Read the random fold protocol choice
				if(argc>5) fold_option=atoi(argv[5]);
			}
		}
	}

	// Apply random folds, periodically reevaluating the bounding circle
	// for better efficiency in sampling folds
	sim_flatfold ff;
	ff.seed(seed);
	while(i<folds) {
		std::printf("fold: %d \n", i);
		if(fold_option==0) ff.random_fold(frac,rand_sign);
		else if(fold_option==1) ff.random_fold1(rand_sign);
		if(++i%2==0) ff.compute_bounds();
	}

	// Output positive and negative creases, plus the shape of the folded
	// sheet
	ff.crease_map("pos.dat",true);
	ff.crease_map("neg.dat",false);
	ff.output("ff.dat");
}
