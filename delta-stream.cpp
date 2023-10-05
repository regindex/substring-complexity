// Copyright (c) 2023, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>
#include <vector>
#include <iostream>
#include <sdsl/lcp.hpp>
#include <sdsl/suffix_arrays.hpp>
#include "hyperloglog.hpp"

using namespace std;
using namespace sdsl;
using namespace hll;

//template on the underlying count-distinct sketch
template<class hll_t = HyperLogLogHIP>
class sketch{

	/*vector<uint64_t> v {1,2,2,2,3,4,5,5,5,6,6,7,7,7,7,8,9,9,9,10,10};
	//HyperLogLog hll(10);
	HyperLogLogHIP hll(10);
	for(auto x:v) hll.add(reinterpret_cast<char*>(&x),sizeof(x));
	cout << "cardinality = " << hll.estimate() << endl;*/

public:

	//total sketch size is <= (_e + (_u/log2(a))) * 2^_r HLL registers 
	//(_u/log2(a) = log_a(U), where U=2^_u is the upper bound to the stream's length)

	static constexpr uint8_t default_u = 32; 
	static constexpr double default_a = 1.001; //logarithm base for the sampling of factor lengths
	static constexpr uint8_t default_r = 8;
	static constexpr uint64_t default_e = 100; // the first default_e k-mer lengths are 1..default_e

	sketch(	uint8_t u = default_u, 
			double a = default_a,
			uint8_t r = default_r, 
			uint64_t e = default_e) : _u(u), _a(a), _r(r), _e(e) {

		//total number of HLL sketches: at least default_exact
		//uint64_t number_hll_sketches = max(_u*_p,uint64_t(default_exact));

		//hll_sketches = vector<hll_t>(number_hll_sketches,{_r});

	}

	//right-extend text by one character (integer alphabet, we use uint64_t)
	void extend(uint64_t c){

	}

	//merge this sketch with s
	void merge(const sketch& s){

	}

	//get estimate of measure delta
	double estimate_delta() const {
		return 0;
	}

	//store sketch to output stream
	void store(ostream& os) const {

	}

	//load sketch from input stream
	void load(const istream& is){

	}

	//return current stream length
	uint64_t get_stream_length() const {
		return stream_length;
	}

	//return log2(upper bound to stream length)
	uint8_t get_upper_bound() const {
		return _u;
	}

	//return log2(number of registers used by HLL)
	uint8_t get_hll_registers() const {
		return _r;
	}


private:

	vector<hll_t> hll_sketches;	// the count-distinct sketches
	vector<uint64_t> lengths;  // sampled k-mer lengths corresponding to the sketches

	//current stream length (or sum of streams lengths if the sketch is the result of union of streams)
	uint64_t stream_length = 0; 

	uint8_t _u;  // log2(upper bound to stream length)
	double _a;   //logarithm base for the sampling of factor lengths
	uint8_t _r;  // log2(number of registers used by each HLL sketch)
	uint64_t _e; // the first default_e k-mer lengths are 1..default_e

};


/*
	build sketch on the input stream and save it to outfile, if outfile is not empty
*/
void stream_delta(string outfile = {}){

	cout << "file: " << outfile<<endl;

	sketch<> s;


}

int main(int argc, char* argv[]){

    if (argc < 2) {
        cout	<< "Usage: " << argv[0] << " [options] [< input stream]" << endl << endl

        		<< "Tool to compute and compare compressibility sketches based on the delta measure." << endl << endl

        	 	<< "	-s, --stream" << endl 
        		<< "		Computes the sketch of the input stream (pipeline or redirect). Outputs delta(stream). Can be combined with -o to save the sketch." << endl << endl

        		<< "	-d s, --delta s" << endl
        		<< "		Given file s containing a sketch, outputs delta(s)." << endl << endl

        		<< "	-o s, --output s" << endl
        		<< "		Store the resulting sketch to file s." << endl << endl
        		
        		<< "	-m s1 s2, --merge s1 s2" << endl
        		<< "		Merge the sketches contained in files s1 and s2. Outputs delta(s1,s2). The two input sketches must have the same parameters. Can be combined with -o to save the merged sketch." << endl << endl

        		<< "	-u u, --upper-bound u" << endl
        		<< "		Logarithm in base 2 of the upper bound to the stream length (used with -s). The tool uses working space O(sqrt(2^u) * u). Default: u = " << uint64_t(sketch<>::default_u) << "." << endl << endl

        		<< "	-a a, --sample-rate a" << endl
        		<< "		Sample rate. Samples log_a(stream length) factor lengths. Must be a double a>1. Default: a = " << sketch<>::default_a << "." << endl << endl

        		<< "	-e E, --exact E" << endl
        		<< "		Factor lengths {1,2,...,E} are always sampled. Default: E = " << sketch<>::default_e << endl << endl


        		<< "	-r R, --registers R" << endl
        		<< "		Logarithm in base 2 of the number of registers used by each HLL sketch. Default: R = " << uint64_t(sketch<>::default_r) << ". Range of R: [4,30]" << endl << endl

        		<< endl;
        return 1;
    }
	

    if(strcmp(argv[1], "-s") == 0 or strcmp(argv[1], "--stream")){

    	if(argc > 2) stream_delta(argv[2]);

    }

}