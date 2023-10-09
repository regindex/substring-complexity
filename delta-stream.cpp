// Copyright (c) 2023, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>
#include <vector>
#include <iostream>
#include <set>
#include <sdsl/lcp.hpp>
#include <sdsl/suffix_arrays.hpp>
#include "hyperloglog.hpp"

using namespace std;
using namespace sdsl;
using namespace hll;


//computes z^e mod q using fast exponentiation
uint64_t fexp(uint64_t z, uint64_t e, uint64_t q){

	__uint128_t z_2_i = z % q; // z^(2^i) mod q (initially, i=0)
	__uint128_t res = 1;

	while(e>0){
		res = (res * (e&1 ? z_2_i : 1)) % q;
		e = e>>1;
		z_2_i = (z_2_i * z_2_i) % q; 
	}

	return res;

}

//template on the underlying count-distinct sketch and on the char type
template<class hll_t = HyperLogLogHIP, class char_t = uint8_t>
class sketch{

public:

	//total sketch size is <= (_e + (_u/log2(a))) * 2^_r HLL registers 
	//(_u/log2(a) = log_a(U), where U=2^_u is the upper bound to the stream's length)

	static constexpr uint8_t default_u = 32; 
	static constexpr double default_a = 1.2; //logarithm base for the sampling of factor _lengths
	static constexpr uint8_t default_r = 10;
	static constexpr uint64_t default_e = 20; // the first default_e k-mer _lengths are 1..default_e

	static constexpr uint64_t q = (uint64_t(1)<<61) - 1; // prime for Karp-Rabin fingerprinting (M61)

	sketch(	uint64_t z, //random base for KR hashing
			uint8_t u = default_u, 
			double a = default_a,
			uint8_t r = default_r, 
			uint64_t e = default_e) : _z(z), _u(u), _a(a), _r(r), _e(e) {

		uint64_t U = uint64_t(1)<<_u;
		double exp = 1;//exponential
		uint64_t k = 1;//last sampled length
		uint64_t sampled_lenghts = _e;//number of sampled kmer _lengths
		
		//TO DO adjust the window size once theory is ready!!
		//for now, sqrt(max stream length)*log(max stream length)
		_window_size = (uint64_t(1)<<(_u/2))*_u;

		while(exp < U){
			if(uint64_t(exp)>k and uint64_t(exp) > _e){
				sampled_lenghts++;
				k = exp;
			}
			exp *= a;
		}

		_hll_sketches = vector<hll_t>(sampled_lenghts,{_r});
		_lengths = vector<uint64_t>(sampled_lenghts,0);

		uint64_t i=0;
		while(i++<_e) _lengths[i-1] = i;
		i--;

		exp = 1;
		k = _e;

		while(exp < U){
			if(uint64_t(exp)>k and uint64_t(exp) > _e){
				_lengths[i++]=exp;
				k = exp;
			}
			exp *= a;
		}

	}

	//right-extend text by one character
	//WARNING: cannot be called on a sketch loaded with load() or after running merge(..).
	void extend(char_t c){

		if(_stream_length==0){ // first stream character: initialize arrays

			_window = vector<char_t>(_window_size,0);
			_window_bookmarks = vector<uint64_t>(_lengths.size(),0);
			_fingerprints =	vector<uint64_t>(_lengths.size(),0);
			_exponents = vector<uint64_t>(_lengths.size(),0);

			for(uint64_t i=0;i<_lengths.size();++i){
				assert(i<_exponents.size());
				_exponents[i] = fexp(_z,_lengths[i]-1,q); //fast exponentiation: z^(_lengths[i]-1) mod q
			}

		}

		for(uint64_t i=0;i<_lengths.size();++i){

			/*DEBUG*/ //cout << "checking length " << _lengths[i] << " (window size = " << _window_size << ")" << endl;

			assert(i<_lengths.size());
			if(_lengths[i] <= _window_size){ //ignore k-mers that do not fit in window

				// remove from the active fingerprints the oldest character
				// and move forward their iterators
				
				/*DEBUG*/ //cout << "  stream len = " << _stream_length << ", lengths[i] = " << _lengths[i] << endl;

				if(_stream_length >= _lengths[i]){
					
					/*DEBUG*/ //cout << "  enter IF (remove oldest char)" << endl;

					char_t b = _window[_window_bookmarks[i]]; //get oldest character

					__uint128_t remove = (__uint128_t(b)*__uint128_t(_exponents[i])) % q;
					_fingerprints[i] = ((__uint128_t(_fingerprints[i]) + q) - remove) % q;

					/*DEBUG if(_lengths[i]==2){


						fp_alt.insert(hash_debug);

						cout << "exponent = " << _exponents[i] << endl;
						cout << "stream len = " << _stream_length << endl;
 						cout << "fingerprints: " << _fingerprints[i] << " / " << hash_debug << endl;

					}*/

					// shift forward the bookmark
					_window_bookmarks[i] = (_window_bookmarks[i]+1) % _window_size;

				}

				//append new character to all fingerprints:
				//since we have already removed the oldest character (done in previous if),
				//we only need to left-shift fingerprint (multiply by _z) and add c
				_fingerprints[i] = (__uint128_t(_fingerprints[i])*__uint128_t(_z) + c) % q;

				//_stream_length increased by 1 in the fingerprints: if the updated stream length
				// is >= the current k-mer length, insert the k-mer fingerprint in the HLL sketch
				if(_stream_length+1 >= _lengths[i]){

					/*DEBUG*/ //cout << "  inserting fingerprint " << _fingerprints[i] << " in sketch for length " << _lengths[i] << endl;

					assert(i<_hll_sketches.size());
					_hll_sketches[i].add(reinterpret_cast<char*>(&_fingerprints[i]),sizeof(_fingerprints[i]));

				}

			}

		}

		//if empty stream, leave _window_head to 0 (_window_head always points to
		//most recent stream character). Otherwise, increment it.
		_window_head = (_window_head + (_stream_length>0))%_window_size;

		//store new character in updated _window_head
		_window[_window_head] = c;

		//update stream length
		_stream_length++;

	}

	//merge this sketch with s
	void merge(const sketch& s){

	}

	//get estimate of measure delta
	double estimate_delta() const {
		
		double delta = 0;

		for(uint64_t i=0;i<_lengths.size();++i){

			/*DEBUG*/ //cout << "d_" << _lengths[i] << " = " << _hll_sketches[i].estimate() << endl;
			delta = max(delta,_hll_sketches[i].estimate()/_lengths[i]);

		}

		return delta;
	}

	//store sketch to output stream
	void store(ostream& os) const {

	}

	//load sketch from input stream
	void load(const istream& is){

	}

	//return current stream length
	uint64_t stream_length() const {
		return _stream_length;
	}

	//return log2(upper bound to stream length)
	uint8_t get_upper_bound() const {
		return _u;
	}

	//return log2(number of registers used by HLL)
	uint8_t get_hll_registers() const {
		return _r;
	}

	//number of sampled factor _lengths (for each, we store a HLL sketch)
	uint64_t get_number_of_samples(){
		return _lengths.size();
	}


private:

	vector<uint64_t> _fingerprints; // Karp-Rabin fingerprints of the windows

	vector<uint64_t> _exponents; // z^_lengths[0], z^_lengths[1], z^_lengths[2], ...

	//store a window of the last 
	vector<char_t> _window;
	uint64_t _window_head = 0;
	vector<uint64_t> _window_bookmarks;

	vector<hll_t> _hll_sketches;	// the count-distinct sketches
	vector<uint64_t> _lengths;  // sampled k-mer _lengths corresponding to the sketches

	//current stream length (or sum of streams _lengths if the sketch is the result of union of streams)
	uint64_t _stream_length = 0; 

	uint64_t _window_size;

	uint64_t _z; // base of Karp-Rabin hashing (should be a random number in (0,q))
	uint8_t _u;  // log2(upper bound to stream length)
	double _a;   //logarithm base for the sampling of factor _lengths
	uint8_t _r;  // log2(number of registers used by each HLL sketch)
	uint64_t _e; // the first default_e k-mer _lengths are 1..default_e

};




/*
	build sketch on the input stream and save it to outfile, if outfile name is not empty
*/
void stream_delta(string outfile = {}){

	//TODO replace with random integer

	//sketch<HyperLogLog> s(324289893284831);	
	sketch<> s(324289893284831);

	uint64_t i = 0;

	while(cin){
		uint8_t c = cin.get();
		if(cin){
			s.extend(c);
			i++;
			if(i%100000==0) cout << "Processed " << i << " characters." << endl;
		}
	}

	cout << "number of sampled lengths : " << s.get_number_of_samples() << endl;
	cout << "stream length = " << s.stream_length() << endl;
	cout << "delta = " << s.estimate_delta() << endl;

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

        		<< "	-c s1 s2, --ncd s1 s2" << endl
        		<< "		Outputs the normalized compression distance (value in [0,1]) of sketches s1 and s2. The two input sketches must have the same parameters." << endl << endl

        		<< "	-u u, --upper-bound u" << endl
        		<< "		Logarithm in base 2 of the upper bound to the stream length (used with -s). The tool uses working space O(2^(u/2) * u) to build the sketch of the stream. Default: u = " << uint64_t(sketch<>::default_u) << "." << endl << endl

        		<< "	-a a, --sample-rate a" << endl
        		<< "		Sample rate. Samples log_a(stream length) factor _lengths. Must be a double a>1. Default: a = " << sketch<>::default_a << "." << endl << endl

        		<< "	-e E, --exact E" << endl
        		<< "		Factor _lengths {1,2,...,E} are always sampled. Default: E = " << sketch<>::default_e << endl << endl


        		<< "	-r R, --registers R" << endl
        		<< "		Logarithm in base 2 of the number of registers used by each HLL sketch. Default: R = " << uint64_t(sketch<>::default_r) << ". Range of R: [4,30]" << endl << endl

        		<< endl;
        return 1;
    }
	

    if(strcmp(argv[1], "-s") == 0 or strcmp(argv[1], "--stream")){

    	stream_delta();

    }

}