// Copyright (c) 2023, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>
#include <vector>

#include <iostream>
#include <sdsl/lcp.hpp>
#include <sdsl/suffix_arrays.hpp>

using namespace std;
using namespace sdsl;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " file" << endl;
        return 1;
    }
    string file = argv[1];

	cache_config cc(true); // delete temp files after lcp construction
    lcp_wt<> lcp;
    construct(lcp, file, 1);

    vector<int64_t> dk(lcp.size(),0);

    //cout << "sa = "; for(auto x:csa) cout << x << " "; cout << endl;
    //cout << "lcp = "; for(auto x:lcp) cout << x << " "; cout << endl;

    dk[1]++; //this is dk[lcp[0]+1]++

    for(uint64_t i=1;i<lcp.size();++i){
    	dk[lcp[i]+1]++; // let lcp[i]=k. Here we see a new factor of length k+1, k+2, k+3, ...
    	dk[i]--; // all suffixes of length <= i (ended by $) create a spurious distinct factor of length i
    }

    for(uint64_t i=1;i<dk.size();++i) dk[i] += dk[i-1];
    
    double delta = 0;
	double argmax = 0;

    for(uint64_t i=1;i<dk.size();++i){

    	if(double(dk[i])/i>delta) argmax = i;

    	delta = max(delta, double(dk[i])/i);

    }


    cout << "delta = " << delta << endl; 
    cout << "argmax_k(d_k/k) = " << argmax << endl; 

    cout << endl << "k\td_k\td_k/k" << endl;

    for(int i=1;i<min(100,int(dk.size()));++i) 
    	cout << i << "\t" << dk[i] << "\t" << double(dk[i])/i << endl;

}