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

    if (argc < 2 or argc > 4) {
        cout << "Usage: " << argv[0] << " <file> [N] [K]" << endl;
        cout << "where <file> contains the input string, N is the number of delta_k "
        	 << "values to print (default: 0), and K is the maximum k-mer length to" 
        	 << "consider when estimating delta (default: 1000)." << endl;
        return 1;
    }
    string file = argv[1];
    uint64_t N = 0;
    uint64_t K = 1000;

    if(argc>=3) N = atoi(argv[2]);
    if(argc==4) K = atoi(argv[3]);

	cache_config cc(true); // delete temp files after lcp construction
    lcp_wt<> lcp;
    construct(lcp, file, 1);

    //NOTE: use wider integers on very large files
    vector<uint32_t> dk(lcp.size(),0);

    //cout << "lcp = "; for(auto x:lcp) cout << x << " "; cout << endl;

    dk[1]++; //this is dk[lcp[0]+1]++

    for(uint64_t i=1;i<lcp.size();++i)
    	dk[lcp[i]+1]++; // let lcp[i]=k. Here we see a new factor of length k+1, k+2, k+3, ...
    					// note: "spurious" factors of length <=i contining a $ will be removed in the next loop
        
    double delta = 0;
    double delta_leqK = 0;
	double delta_i = 0;
	double dk_value = 0;
	uint64_t argmax = 0;
	uint64_t lrs = 0;

	if(N>0) cout << endl << "k\td_k\td_k/k" << endl;

    for(uint64_t i=1;i<dk.size();++i){

    	dk_value = dk_value + dk[i] - 1;
		delta_i = dk_value/i;

		if(dk_value < lcp.size() - i) 
			lrs = i; // dk < n - i + 1: there exists repeated string of length i

    	if(delta_i>delta){
    		argmax = i;
			delta = delta_i;
		}

    	if(i <= K and delta_i>delta_leqK)
			delta_leqK = delta_i;

		if(i<=N) cout << i << "\t" << dk_value << "\t" << delta_i << endl;

    }

	if(N>0) cout << endl;

	cout << "string length n = " << lcp.size() - 1 << endl;
    cout << "delta = " << int(delta) << endl;
    cout << "delta_{<=K} (with K = " << K << ") = " << delta_leqK << endl;
    cout << "ratio n/delta = " << double(lcp.size() - 1)/delta << endl; 
    cout << "argmax_k(d_k/k) = " << argmax << endl; 
    cout << "length of longest repeated substring = " << lrs << endl; 
    /*
    ofstream output("exact_delta.txt");
    output << int(delta);
    output.close();
    */
}