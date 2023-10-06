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

//template on the type of char
template<class char_t = char>
class window_stream{

public:

	struct iterator{

        using iterator_category = std::forward_iterator_tag;
        using value_type        = char_t;
        using pointer           = char_t*;
        using reference         = char_t&;

        // constructor from a position in the stream
        iterator(window_stream& ws, uint64_t pos) : _ws(ws), _pos(pos) {};

        reference operator*(){ return _ws.at(_pos); };
        
        // Prefix increment. Increments the iterator and returns it. 
        iterator& operator++(){
        	_pos = (_pos + 1) % _ws._w;
        	return *this;
        }
        
    private:

        window_stream& _ws;
        uint64_t _pos;
    
    };

	static const uint64_t default_w = 10e6; // 1M characters

	window_stream(uint64_t window_size = default_w) : _w(window_size), _head(0) {

		_window = vector<char_t>(_w,0);

	}

	//shifts window to the right by 1 position and appends a new character at the beginning of stream
	void push(char_t c){
		_head = (_head+1)%_w;
		_window[_head] = c;
	}

	//get most recent character on the stream. Do not shift window.
	char& head(){

		return _window[_head];
			
	}

	//random access (within the window)
	char& at(uint64_t i){ return _window[i]; }

	uint64_t window_size(){return _w;}

	iterator head_iterator(){
		return {*this,_head};
	}

private:

	uint64_t _w; //window size

	vector<char_t> _window;
	uint64_t _head; //most recent stream character

};

//template on the underlying count-distinct sketch and on the char type
template<class hll_t = HyperLogLogHIP, class char_t = char>
class sketch{

public:

	//total sketch size is <= (_e + (_u/log2(a))) * 2^_r HLL registers 
	//(_u/log2(a) = log_a(U), where U=2^_u is the upper bound to the stream's length)

	static constexpr uint8_t default_u = 32; 
	static constexpr double default_a = 1.001; //logarithm base for the sampling of factor lengths
	static constexpr uint8_t default_r = 8;
	static constexpr uint64_t default_e = 40; // the first default_e k-mer lengths are 1..default_e

	static constexpr uint64_t q = (uint64_t(1)<<61) - 1; // prime for Karp-Rabin fingerprinting (M61)

	sketch(	uint64_t z, //random base for KR hashing
			uint8_t u = default_u, 
			double a = default_a,
			uint8_t r = default_r, 
			uint64_t e = default_e) : _z(z), _u(u), _a(a), _r(r), _e(e) {

		uint64_t U = uint64_t(1)<<_u;
		double exp = 1;//exponential
		uint64_t k = 1;//last sampled length
		uint64_t sampled_lenghts = _e;//number of sampled kmer lengths
		
		while(exp < U){
			if(uint64_t(exp)>k and uint64_t(exp) > _e){
				sampled_lenghts++;
				k = exp;
			}
			exp *= a;
		}

		hll_sketches = vector<hll_t>(sampled_lenghts,{_r});
		lengths = vector<uint64_t>(sampled_lenghts,0);

		uint64_t i=0;
		while(i++<_e) lengths[i-1] = i;
		i--;

		exp = 1;
		k = _e;

		while(exp < U){
			if(uint64_t(exp)>k and uint64_t(exp) > _e){
				lengths[i++]=exp;
				k = exp;
			}
			exp *= a;
		}

	}

	//right-extend text by one character
	//WARNING: cannot be called on a sketch loaded with load() or after running merge(..).
	void extend(char_t c){

		if(stream_length==0){ // first stream character: initialize arrays

			_fingerprints =	{lengths.size(),0};
			_iterators = {lengths.size(),_ws.head_iterator()};
			_exponents = {lengths.size(),0};

			for(uint64_t i=0;i<lengths.size();++i)
				_exponents[i] = fexp(_z,lengths[i]-1,q); //fast exponentiation: z^(lengths[i]-1) mod q

		}

		for(uint64_t i=0;i<lengths.size();++i){

			if(lengths[i] <= _ws.window_size()){ //ignore k-mers that do not fit in window

				// remove from the active fingerprints the oldest character
				// and move forward their iterators
				if(stream_length >= lengths[i]){

					char_t b = *_iterators[i]; //get oldest character
					__uint128_t remove = (__uint128_t(b)*__uint128_t(_exponents[i])) % q;
					_fingerprints[i] = ((__uint128_t(_fingerprints[i]) + q) - remove) % q;
					++_iterators[i];

				}

				//append new character to all fingerprints
				_fingerprints[i] = (__uint128_t(_fingerprints[i])*__uint128_t(_z) + c) % q;

				if(stream_length >= lengths[i]){

					//if updated fingerprint is active, insert it in the HLL sketch
					hll_sketches[i].add(reinterpret_cast<char*>(&_fingerprints[i]),sizeof(_fingerprints[i]));

				}

			}

		}

		// append the new character to the internal window-stream

		_ws.push(c);
		stream_length++;

	}

	//merge this sketch with s
	void merge(const sketch& s){

	}

	//get estimate of measure delta
	double estimate_delta() const {
		
		double delta = 0;

		for(uint64_t i=0;i<lengths.size();++i){

			delta = max(delta,hll_sketches[i].estimate()/lengths[i]);

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

	//number of sampled factor lengths (for each, we store a HLL sketch)
	uint64_t get_number_of_samples(){
		return lengths.size();
	}


private:

	vector<uint64_t> _fingerprints; // Karp-Rabin fingerprints of the windows

	vector<uint64_t> _exponents; // z^lengths[0], z^lengths[1], z^lengths[2], ...

	window_stream<char_t> _ws;
	vector<typename window_stream<char_t>::iterator> _iterators;

	vector<hll_t> hll_sketches;	// the count-distinct sketches
	vector<uint64_t> lengths;  // sampled k-mer lengths corresponding to the sketches

	//current stream length (or sum of streams lengths if the sketch is the result of union of streams)
	uint64_t stream_length = 0; 

	uint64_t max_length = 0; //largest integer such that lengths[max_length] <= stream_length (except with empty sketch) 

	uint64_t _z; // base of Karp-Rabin hashing (should be a random number in (0,q))
	uint8_t _u;  // log2(upper bound to stream length)
	double _a;   //logarithm base for the sampling of factor lengths
	uint8_t _r;  // log2(number of registers used by each HLL sketch)
	uint64_t _e; // the first default_e k-mer lengths are 1..default_e

};




/*
	build sketch on the input stream and save it to outfile, if outfile name is not empty
*/
void stream_delta(string outfile = {}){

	//TODO replace with random integer
	
	sketch<> s(324231);

	s.extend('a');
	s.extend('b');
	s.extend('c');
	s.extend('a');
	s.extend('b');
	s.extend('c');
	s.extend('a');
	s.extend('b');

	cout << "number of sampled lengths : " << s.get_number_of_samples() << endl;
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
        		<< "		Sample rate. Samples log_a(stream length) factor lengths. Must be a double a>1. Default: a = " << sketch<>::default_a << "." << endl << endl

        		<< "	-e E, --exact E" << endl
        		<< "		Factor lengths {1,2,...,E} are always sampled. Default: E = " << sketch<>::default_e << endl << endl


        		<< "	-r R, --registers R" << endl
        		<< "		Logarithm in base 2 of the number of registers used by each HLL sketch. Default: R = " << uint64_t(sketch<>::default_r) << ". Range of R: [4,30]" << endl << endl

        		<< endl;
        return 1;
    }
	

    if(strcmp(argv[1], "-s") == 0 or strcmp(argv[1], "--stream")){

    	stream_delta();

    }

}