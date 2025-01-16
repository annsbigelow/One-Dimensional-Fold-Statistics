#ifndef FOLDS_STATS_H
#define FOLDS_STATS_H
// ^^^ CHR: The "include guard" lines should occur first, before any library includes

//using namespace std;
// ^^^ CHR: It's generally a good idea not to put any "using namespace" command
// into header files, since it means that any program that includes this file
// is forced to have all those functions in the namespace (e.g. you write a
// program with a local variable called "ofstream", and this command now means
// that std::ofstream will clash with your local one.) Namespace commands are
// okay in the cpp files.
#include <iostream>
#include <vector>
#include <fstream>

#include <gsl/gsl_rng.h>

class folds_stats {
    public:
        const int th,instance,n;
        // CHR: a number of the variables in the lists below are mainly used
        // for local manipulations within functions. They would be better declared
        // within those functions themselves for several reasons: (1) having
        // them declared locally reduces the chances of bugs due to variable
        // crosstalk between functions, (2) the compiler can better optimize
        // when it knows a variable is local. You want to include variables
        // here that have a global purpose (e.g. rng, rbits, etc.).

        // Class constructor
        // seeds rng upon class instance creation
        folds_stats(int th_=1,const int inst=1000, const int n_=28) : th(th_), rng(gsl_rng_alloc(gsl_rng_taus2)), n_rbits(0), instance(inst),n(n_)
        {seed(th);}

        // Class destructor
        ~folds_stats() {
            gsl_rng_free(rng);
        }

        inline void seed(long int k) { 
            gsl_rng_set(rng, k);
        }

        void altFold_segdens(int &numplaces,string txt);
        void segdens(int &numplaces,string txt);
        void log_fixedn(string txt);
        void logavg(string txt);
        vector<double> altFold(string txt);
        vector<double> fold(vector<double> &segs_in);
        double getmin(vector<double> &segs_in);
        double getmax(vector<double> &segs_in);
        void display(vector<double> arr);

        std::ofstream data;
        std::ofstream densData;
        std::ofstream densData1;
    private:
        gsl_rng* rng;
        /** Generates a random bit of information.
         * \return The random bit, with true and false having equal
         * probability. */
        inline bool rand_bit() {

            // If there are no random bits available in our previous random
            // number then generate a new one.
            if(n_rbits==0) {
                rbits=gsl_rng_get(rng);
                n_rbits=8*sizeof(int)-1;
            } else n_rbits--;

            // Pick off a random bit from the number, and perform a binary
            // shift to remove it from the number
            bool res=rbits&1;
            rbits>>=1;
            return res;
        }
        /** The number of random bits available. */
        int n_rbits;
        /** The random number to use for generating random bits. */
        unsigned int rbits;
};

#endif
