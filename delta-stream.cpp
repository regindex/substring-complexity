// Copyright (c) 2023, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include "internal/delta-stream.hpp"
#include "internal/delta-stream-mt.hpp"

// struct containing command line parameters and other globals
struct Args {
  //string outname;
  bool stream = false;
  bool delta = false;
  bool merge = false;
  uint8_t registers = 10;
  int threads = 1;
  string output = string();
  string sketch1 = string();
  string sketch2 = string();
  bool verbose = false; // verbosity level
};

// function that prints the instructions for using the tool
void print_help(char** argv) { 
	cout << 
		"Usage: " << argv[0] << " [options] [< input stream]" << endl << endl

		<< "Tool to compute and compare compressibility sketches based on the delta measure." << endl << endl

	 	<< "	-s, --stream" << endl 
		<< "		Computes the sketch of the input stream (pipeline or redirect). Outputs delta(stream). Can be combined with -o to save the sketch." << endl << endl

		<< "	-d s, --delta s" << endl
		<< "		Given file s containing a sketch, outputs delta(s)." << endl << endl

		<< "	-o s, --output s" << endl
		<< "		Store the resulting sketch to file s (only for single thread so far)." << endl << endl
		
		<< "	-m s1 s2, --merge s1 s2" << endl
		<< "		Merge the sketches contained in files s1 and s2. Outputs delta(s1,s2). The two input sketches must have the same parameters. Can be combined with -o to save the merged sketch." << endl << endl

		//<< "	-c s1 s2, --ncd s1 s2" << endl
		//<< "		Outputs the normalized compression distance (value in [0,1]) of sketches s1 and s2. The two input sketches must have the same parameters." << endl << endl

		//<< "	-u u, --upper-bound u" << endl
		//<< "		Logarithm in base 2 of the upper bound to the stream length (used with -s). The tool uses working space O(2^(u/2) * u) to build the sketch of the stream. Default: u = " << uint64_t(sketch<>::default_u) << "." << endl << endl

		//<< "	-a a, --sample-rate a" << endl
		//<< "		Sample rate. Samples log_a(stream length) factor _lengths. Must be a double a>1. Default: a = " << sketch<>::default_a << "." << endl << endl

		//<< "	-e E, --exact E" << endl
		//<< "		Factor _lengths {1,2,...,E} are always sampled. Default: E = " << sketch<>::default_e << endl << endl

		<< "	-r R, --registers R" << endl
		<< "		Logarithm in base 2 of the number of registers used by each HLL sketch. Default: R = " << uint64_t(sketch<>::default_r) << ". Range of R: [4,30]" << endl << endl

		<< "	-t T, --threads T" << endl
		<< "		Number of threads used to construct the sketches. Default: T = 1" << ". Range of R: [2...]" << endl << endl

		<< endl;
}

// function for parsing the input arguments
void parseArgs(int argc, char** argv, Args& arg) {

	if(argc < 2){ print_help(argv); }

	for(size_t i=1;i<argc;++i)
	{

		string param = argv[i];

		if( param == "-s" or param == "--stream" )
		{
			arg.stream = true;
			// read a stream
		}
		else if( param == "-r" or param == "--registers" )
		{
			i++;
		    arg.registers = atoi( argv[i] );
		    if( arg.registers < 4 or arg.registers > 30 )
		    {
			    cerr << "The number of register range is [4,30]." << endl;
			    exit(1);
		    }
	  	}
		else if( param == "-t" or param == "--threads" )
		{
			i++;
		    arg.threads = atoi( argv[i] );
		    if( arg.registers < 2 )
		    {
			    cerr << "Select at least 2 threads." << endl;
			    exit(1);
		     }
	   	}
		else if( param == "-d" or param == "--delta" )
		{
			i++;
		    arg.sketch1 = string( argv[i] );
		    arg.delta = true;
	   	}
		else if( param == "-m" or param == "--merge" )
		{
			i++;
		    arg.sketch1 = string( argv[i++] );
		    arg.sketch2 = string( argv[i] );
		    arg.merge = true;
	   	}
		else if( param == "-o" or param == "--output" )
		{
			i++;
	    	arg.output = string( argv[i] );
	  	}
    	else if( param == "-h" or param == "--help")
		{
		    print_help(argv); exit(-1);
		    // fall through
		}
	    else if( param == "<")
	    {
	        break;
	        // skip
	    }
	    else{
	        cerr << "Unknown option. Use -h for help." << endl;
	        exit(-1);
	    }
	    // check mode
	    uint32_t counter = (int)arg.stream + (int)arg.delta + (int)arg.merge;
	    if(counter != 1)
	    {
	    	cerr << "Please select one option out of stream|delta|merge" << endl;
	    	exit(1);
	    }
	}
}

void extend_window(sketch_MT<> * s, vector<uint8_t> * buffer, uint64_t i)
{
	for(uint64_t j=0;j<i;++j)
		s->extend_window((*buffer)[j]);
}

void extend_rlbwt(sketch_MT<> * s, vector<uint8_t> * buffer, uint64_t i)
{
	for(uint64_t j=0;j<i;++j)
	{
		s->extend_rlbwt((*buffer)[j]);
		if(s->is_rlbwt_dropped())
			break;
	}
}

uint64_t loadBuffer(vector<uint8_t>& buffer, const uint64_t K)
{
	uint64_t i = 0;
	
	while(i < K)
	{
		uint8_t c = cin.get();
		if(!cin){ break; }
		buffer[i++] = c;
	}

	return i;
}

void print_stats(sketch<>& s)
{
	cout << "number of sampled lengths : " << s.get_number_of_samples() << endl;
	cout << "stream length = " << s.stream_length() << endl;
	cout << "delta = " << s.estimate_delta() << endl;
	
	ofstream output("delta.txt");
	output << s.estimate_delta();
	output.close();
}

double get_delta(sketch_MT<>& s)
{
	return s.estimate_delta();
}

void load_sketch(sketch<>& s, string sketch)
{
		vector<uint64_t> sampled_lengths;
		sample_kmer_lengths(sampled_lengths,e,u,a);

		ifstream input(sketch); 
		s.load(input,&sampled_lengths);
		input.close();
}

/*
	build sketch on the input stream and save it to outfile, if outfile name is not empty
*/
void stream_delta(uint8_t registers, string outfile = string()){

	//sketch<HyperLogLog> s(324289893284831);	
	uint64_t rand = random(1,255);
	sketch<> s(rand,registers);

	uint64_t i = 0;
	while(cin){
		uint8_t c = cin.get();
		if(cin)
		{
			s.extend_window(c);
			if(!s.is_rlbwt_dropped())
			{
				s.extend_rlbwt(c);
			}
			i++;
			//if(i%10000000==0) cout << "Processed " << i << " characters." << endl;
		}
	}
	// print delta stats
	print_stats(s);
	// store sketches if needed
	if(outfile != string())
	{
		ofstream os(outfile);
		s.store(os);
		os.close();
	}
}

/*
	build sketch (in a parallel way) on the input stream and save it to outfile, if outfile name is not empty
*/
void stream_delta_parallel(uint8_t registers, int32_t threads, string outfile = {}){
	
	uint64_t rand = random(1,255);

	vector<thread> thread_list;
	thread rlbwt_thread;
	vector<sketch_MT<>> sketch_list;

	const uint32_t K = 1000000;
	vector<uint8_t> buffer(K,0);

	vector<uint64_t> sampled_lengths;
	sample_kmer_lengths(sampled_lengths,e,u,a);
	uint64_t window_size = compute_window_size(32);
	uint64_t no_sampled = sampled_lengths.size();

	uint64_t i = 0;
	for(;i<no_sampled;++i)
	{
		if( sampled_lengths[i] > window_size )
			break;
	}

	double no_k_thread = round(i/(double)(threads-1));
	if( no_k_thread*(threads-1) >= i ){ threads--; }

	cout << "Sampled lengths: " << no_sampled << endl;
	cout << "Sampled lengths <= window: " << i << endl;
	cout << "k values per thread: " << no_k_thread << endl;

	// initialize k mer lengths for each thread
	vector< vector<uint64_t> * > length_lists;
	uint64_t y = 0;
	for(uint64_t t=0;t<(threads-2);++t)
	{	
		vector<uint64_t>* l = new vector<uint64_t>(no_k_thread,0);
		for(uint64_t j=0;j<no_k_thread;++j){ (*l)[j] = sampled_lengths[y++]; }
		length_lists.push_back(l);
	}
	vector<uint64_t>* l = new vector<uint64_t>(i-y,0);

	// last list may be larger than previous ones
	for(uint64_t j=y;j<i;++j)
	{
		(*l)[j-y] = sampled_lengths[j];
	}
	length_lists.push_back(l);

	// rlbwt list
	l = new vector<uint64_t>(sampled_lengths.size()-i,0);
	for(uint64_t j=i;j<sampled_lengths.size();++j)
	{
		(*l)[j-i] = sampled_lengths[j];
	}
	length_lists.push_back(l);

	// init threads and sketches
	sketch_list.resize(threads-1);
	thread_list.resize(threads-1);
	for(uint64_t i=0;i<(threads-1);++i)
	{
		// init ith sketch
		uint64_t curr_window = (*length_lists[i])[length_lists[i]->size()-1];
		//uint64_t curr_window = window_size;
		sketch_list[i] = sketch_MT<>(rand,registers,curr_window,length_lists[i]);
	}
	sketch_MT<> rlbwt_s(rand,registers,0,length_lists[threads-1]);

	// main execuction
	i = 0;
	while(true)
	{
		uint64_t loadedChars = loadBuffer(buffer,K);
		if( loadedChars == 0 ){ break; }
		for(uint64_t j=0;j<(threads-1);++j)
		{
			thread_list[j] = thread(extend_window, &sketch_list[j], &buffer, loadedChars);
		}

		if(!rlbwt_s.is_rlbwt_dropped())
		{
			thread t = thread(extend_rlbwt, &rlbwt_s, &buffer, loadedChars);
			t.join();
		}

		// threads sincronization
		for(uint64_t j=0;j<(threads-1);++j)
		{
			thread_list[j].join();
		}

		i += loadedChars;
		//if(i%K==0) cout << "Processed " << i << " characters." << endl;
	}

	double delta = 0;
	for(uint64_t j=0;j<threads-1;++j)
	{
		delta = max(delta,get_delta(sketch_list[j]));
	}
	delta = max(delta,get_delta(rlbwt_s));
	
	cout << "number of sampled lengths : " << no_sampled << endl;
	cout << "stream length = " << i << endl;
	cout << "delta = " << delta << endl;

	ofstream output("delta.txt");
	output << delta;
	output.close();
}

int main(int argc, char* argv[]){

    Args arg;
    parseArgs(argc, argv, arg);

    if(arg.stream)
    {
	    if( arg.threads > 1 )
	    {
	    	stream_delta_parallel(arg.registers,arg.threads);
	    }
	    else
	    {
	    	stream_delta(arg.registers,arg.output);
	    }
	}
	else if(arg.delta)
	{
			sketch<> s;
			load_sketch(s,arg.sketch1);
			print_stats(s);
	}
	else if(arg.merge)
	{
			sketch<> s1, s2;
			load_sketch(s1,arg.sketch1);
			load_sketch(s2,arg.sketch2);
			s1.merge(s2);

			if(arg.output != string())
			{
				ofstream os(arg.output);
				s1.store(os);
				os.close();
			}
	}

    return 0;
}