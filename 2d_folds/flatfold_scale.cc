#include "sim_flatfold.hh"

#include <cstdlib>
#include <cmath>

const int max_fold=30;
const int trials=400;

int main(int argc,char **argv) {

	// Obtain seed
	if(argc!=3) {
		fputs("Syntax: ./flatfold_scale <seed> <fold_option>\n"
			"fold_option=0 for standard random fold; 1 to choose a random point then angle; " 
				"2 to choose a random point uniformly on [-1,1]^2, 3 to choose a random angle then displacement; "
				"4 to pick a random angle followed by a uniformly chosen displacement\n",stderr);
		return 1;
	}
	int seed=atoi(argv[1]);
	int fold_option=atoi(argv[2]);

	// Initialize empty counter arrays
	double sl[6*max_fold],*sll=sl+max_fold,
	       *sm=sll+max_fold,*smm=sm+max_fold,
	       *sd=smm+max_fold,*sdd=sd+max_fold;
	for(int k=0;k<6*max_fold;k++) sl[k]=0;
	unsigned int mxf = 0, mnf = 3000000, fsz;
	// Create many random folding instances, and count the number of facets
	// as a function of the folds
#pragma omp parallel
	{
		long fo[max_fold];*fo=1;
		double pos[max_fold],neg[max_fold];*pos=*neg=0;
#pragma omp for schedule(dynamic)
		for(int j=0;j<trials;j++) {
			sim_flatfold ff;
			ff.seed(1024*seed+1+j);

			// Perform random folds, according to the chosen random fold protocol,
			// and store the number of facets after each
			if(fold_option==1){
				for (int i=0;i<max_fold-1;){ 
					ff.random_fold1();
					ff.compute_bounds();
					fo[++i] = ff.f.size();    // CHR: the ++i needs to happen here because the stats need to be binned into the next array element
					ff.crease_mileage(pos[i], neg[i]);
				}
			}
			else if (fold_option==0){
				for(int i=0;i<max_fold-1;) if(ff.random_flatfold()) {
					if(++i%3==0) ff.compute_bounds();
					fo[i]=ff.f.size();
					ff.crease_mileage(pos[i],neg[i]);
				}
			}
			else if (fold_option==2){
                double t_area=ff.tot_area();
				for (int i=0; i<max_fold-1;){ 
					ff.random_fold2(t_area);
                    ++i;
					fo[i] = ff.f.size();
					ff.crease_mileage(pos[i], neg[i]);
				}
			}
			else if (fold_option==3){
				for (int i=0; i<max_fold-1;){
					ff.random_fold3();
					++i;
					fo[i] = ff.f.size();
					ff.crease_mileage(pos[i], neg[i]);
				}
			}
			else if (fold_option==4) {
				for (int i=0; i<max_fold-1;) {
					ff.random_fold4();
					if (++i%2==0) ff.compute_bounds();
					fo[i]=ff.f.size();
					ff.crease_mileage(pos[i], neg[i]);
				}
			}

#pragma omp critical
			{	
				// Store the facet numbers and print a diagnostic message
				printf("%d",j);
				for(int i=0;i<max_fold;i++) {
					printf(" %ld",fo[i]);
					double fd=log(fo[i]);
					sl[i]+=fd;sll[i]+=fd*fd;
					fd=log(pos[i]+neg[i]);
					sm[i]+=fd;smm[i]+=fd*fd;
					fd=pos[i]-neg[i];
					sd[i]+=fd;sdd[i]+=fd*fd;
				}
				putchar('\n');
			}
			fsz=fo[max_fold-2];
			if (fsz<mnf) mnf=fsz; else if (fsz>mxf) mxf=fsz;
		}
	}
	printf("Minimum number of facets: %d\nMaximum number of facets: %d\n",mnf,mxf); // may be out of scope
	// Output the mean and standard deviation of the number of facets
	char buf[256];
	sprintf(buf,"ff%d_scale_%d.dat",fold_option,seed);
	FILE *fp=fopen(buf,"wb");
	if(fp==NULL) {
		fputs("Can't open output file\n",stderr);
		return 1;
	}
	fwrite(&trials,sizeof(int),1,fp);
	fwrite(sl,sizeof(double),6*max_fold,fp);
	fclose(fp);
}
