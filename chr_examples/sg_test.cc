#include <cstdio>
#include <cstring>

#include "seg_collect.hh"

int main() {

    // Add 18 entries to the seg_collect, which will trigger some memory
    // reallocations
    seg_collect sc;
    for(int i=0;i<18;i++) sc.add(i*0.01,1-i*0.01);

    // Print out a list of all of the entries
    sc.print_segs();
}
