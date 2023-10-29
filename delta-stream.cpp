// Copyright (c) 2023, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include "internal/delta-stream.hpp"

// function that prints the instructions for using the tool
void print_help(char** argv) { 
	cout << endl <<
		"Usage: " << argv[0] << " [options] [< input stream]" << endl 

		<< "Tool to compute and compare compressibility sketches based on the delta measure." << endl << endl

	 	<< "	-s, --stream" << endl 
		<< "		Computes the sketch of the input stream (pipeline or redirect). Outputs delta(stream). Can be combined with -o to save the sketch." << endl 

	 	<< "	-z, --base" << endl 
		<< "		Set the base for KR-hashing. Range of z: [1,255]. Default: Random number in [1,...,255]" << endl 

	 	<< "	-l, --rlbwt" << endl 
		<< "		Keep a dynamic run-length BWT to compute the fingerpints for k larger than the maximum window size. Default: False" << endl 

		<< "	-d s, --delta s" << endl
		<< "		Given file s containing a sketch, outputs delta(s)." << endl 

		<< "	-o s, --output s" << endl
		<< "		Store the resulting sketch to file s (only for single thread so far)." << endl 
		
		<< "	-m s1 s2, --merge s1 s2" << endl
		<< "		Merge the sketches contained in files s1 and s2. Outputs delta(s1,s2). The two input sketches must have the same parameters. Can be combined with -o to save the merged sketch." << endl 

		<< "	-c s1 s2, --ncd s1 s2" << endl
		<< "		Outputs the normalized compression distance (value in [0,1]) of sketches s1 and s2. The two input sketches must have the same parameters." << endl 

		<< "	-u u, --upper-bound u" << endl
		<< "		Logarithm in base 2 of the upper bound to the stream length (used with -s). The tool uses working space O(2^(u/2) * u) to build the sketch of the stream. Default: u = 32" << endl 

		//<< "	-a a, --sample-rate a" << endl
		//<< "		Sample rate. Samples log_a(stream length) factor _lengths. Must be a double a>1. Default: a = 1.2" << endl << endl

		//<< "	-e E, --exact E" << endl
		//<< "		Factor _lengths {1,2,...,E} are always sampled. Default: E = 20" << endl << endl

		//<< "	-r R, --registers R" << endl
		//<< "		Logarithm in base 2 of the number of registers used by each HLL sketch. Default: R = 10" << ". Range of R: [4,30]." << endl << endl

		<< "	-p P, --precision P" << endl
		<< "		Sketch precision, higher values result in a more precise delta estimation. Range of p [1,7]. Default: p = 3" << endl 

		<< "	-t T, --threads T" << endl
		<< "		Number of threads used to construct the sketches. Default: T = 1" << endl 

		<< endl;
}

// function for parsing the input arguments
void parseArgs(int argc, char** argv, Args& arg) {

	if(argc < 2){ print_help(argv); exit(1); }

	for(size_t i=1;i<argc;++i)
	{
		string param = argv[i];

		if( param == "-s" or param == "--stream" )
		{
			arg.stream = true;
		}
		else if( param == "-l" or param == "--rlbwt" )
		{
			arg.rlbwt = true;
		}
		else if( param == "-z" or param == "--base" )
		{
			i++;
			arg.z = atoi( argv[i] );
			if( arg.z < 1 or arg.z > 255 )
			{
					cerr << "The KR-hasing base range is [1,255]." << endl;
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
		else if( param == "-u" or param == "--U" )
		{
			i++;
			arg.u = atoi( argv[i] );
		}
		else if( param == "-p" or param == "--precision" )
		{
			i++;
			arg.precision = atoi( argv[i] );
			if( arg.precision < 1 and arg.precision > 7 )
			{
				cerr << "Precision has to be in the range [1,7]." << endl;
				exit(1);
			}
		}
		/*
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
		else if( param == "-e" or param == "--E" )
		{
			i++;
			arg.E = atoi( argv[i] );
		}
		else if( param == "-a" or param == "--A" )
		{
			i++;
			arg.a = atof( argv[i] );
		}
		*/
		else if( param == "-d" or param == "--delta" )
		{
			i++;
			arg.sketch1 = string( argv[i] );
			arg.delta = true;
		}
		else if( param == "-c" or param == "--ncd" )
		{
			i++;
			arg.sketch1 = string( argv[i++] );
			arg.sketch2 = string( argv[i] );
			arg.ncd = true;
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
			arg.outfile = string( argv[i] );
		}
		else if( param == "-h" or param == "--help")
		{
				print_help(argv); exit(-1);
		}
		else if( param == "<")
		{	
			break; // skip
		}
		else
		{
			cerr << "Unknown option. Use -h for help." << endl;
			exit(-1);
		}
	}
	// check mode
	uint32_t counter = (int)arg.stream + (int)arg.delta + (int)arg.merge + (int)arg.ncd;
	if(counter != 1)
	{
			cerr << "Please select one option out of stream|delta|merge|ncd" << endl;
			exit(1);
	}
	// check base z
	if( arg.z == 0 )
	{
		arg.z = random(1,255);
	}
	// set parameters from the selected configurations
	arg.e = conf_list[arg.precision-1]._e;
	arg.a = conf_list[arg.precision-1]._a;
	arg.registers = conf_list[arg.precision-1]._r;
	arg.precision = conf_list[arg.precision-1]._p;
}

void extend_window_mt(sketch<> * s, vector<uint8_t> * buffer, uint64_t i)
{
	for(uint64_t j=0;j<i;++j)
		s->extend_window((*buffer)[j]);
}

void extend_rlbwt_mt(sketch<> * s, vector<uint8_t> * buffer, uint64_t i)
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

void print_delta(sketch<>& s)
{
	cout << "Stream length = " << s.stream_length() << endl;
	cout << "delta = " << s.estimate_delta() << endl;
	
	ofstream output("delta.txt");
	output << s.estimate_delta();
	output.close();
}

void print_ncd(sketch<>& s, double mind, double maxd)
{
	double ncd = (s.estimate_delta() - mind)/maxd;
	cout << "delta NCD = " << ncd << endl;
	
	ofstream output("ncd.txt");
	output << ncd;
	output.close();
}

void load_sketch(sketch<>& s, string sketch)
{
	ifstream input(sketch); 
	s.load(input);
	input.close();
}

void compute_delta(string sketch_path)
{
	sketch<> s;
	load_sketch(s,sketch_path);
	print_delta(s);
}

void merge_sketches(string sketch_path1, string sketch_path2, string outfile, bool ncd = false)
{
	sketch<> s1, s2;
	double d1, d2;
	load_sketch(s1,sketch_path1);
	load_sketch(s2,sketch_path2);
	// compute max and min delta if computing ncd
	if( ncd )
	{
		d1 = s1.estimate_delta();
		d2 = s2.estimate_delta();

		if(d1 > d2)
		{
			swap(d1,d2);
		}
	}
	// merge two sketches
	s1.merge(s2);
	// save the merged sketch
	if(outfile != string())
	{
		ofstream os(outfile);
		s1.store(os);
		os.close();
	}
	// compute ncd
	if( ncd )
		print_ncd(s1,d1,d2);
}

/*
	build sketch on the input stream and save it to outfile, if outfile name is not empty
*/
void stream_delta(Args& arg){

	//sketch<HyperLogLog> s(324289893284831);	
	sketch<> s(arg);

	uint64_t i = 0;
	while(cin){
		uint8_t c = cin.get();
		if(cin)
		{
			s.extend_window(c);
			if(arg.rlbwt and !s.is_rlbwt_dropped())
			{
				s.extend_rlbwt(c,true);
			}
			i++;
		}
	}
	// print delta stats
	print_delta(s);
	// store sketches if needed
	if(arg.outfile != string())
	{
		ofstream os(arg.outfile);
		s.store(os);
		os.close();
	}
}

/*
	build sketch (in a parallel way) on the input stream and save it to outfile, if outfile name is not empty
*/
void stream_delta_parallel(Args& arg){

	vector<thread> thread_list;
	thread rlbwt_thread;
	vector<sketch<>> sketch_list;
	sketch<> rlbwt_s;

	const uint32_t K = 1000000;
	vector<uint8_t> buffer(K,0);

	vector<uint64_t> sampled_lengths, rlbwt_lengths;
	//sample_kmer_lengths
	kmer_lengths_sampling(sampled_lengths,arg.e,arg.u,arg.a,arg.precision);
	uint64_t window_size = compute_window_size(arg.u);
	uint64_t no_sampled = sampled_lengths.size();

	uint64_t i = 0;
	for(;i<no_sampled;++i)
	{
			if( sampled_lengths[i] > window_size )
			break;
	}
	if(!arg.rlbwt)
	{
			sampled_lengths.resize(i);
			sampled_lengths.shrink_to_fit();
	}

	uint64_t threads = arg.threads;
	if(arg.rlbwt){ threads--; }
	double no_k_thread = ceil(i/(double)(threads));
	if(no_k_thread < min_k_thread){ no_k_thread = min_k_thread; }
	while(true)
	{
			if( no_k_thread*(threads-1) >= i ){ threads--; }
			else
				break;
	}

	// initialize k mer lengths for each thread
	vector< vector<uint64_t> * > length_lists;
	uint64_t y = 0;
	for(uint64_t t=0;t<(threads-1);++t)
	{	
		vector<uint64_t>* l = new vector<uint64_t>(no_k_thread,0);
		for(uint64_t j=0;j<no_k_thread;++j){ (*l)[j] = sampled_lengths[y++]; }
		length_lists.push_back(l);
	}

	// initialize last list
	vector<uint64_t>* l = new vector<uint64_t>(i-y,0);
	// last list may be larger than previous ones
	for(uint64_t j=y;j<i;++j)
	{
		(*l)[j-y] = sampled_lengths[j];
	}
	length_lists.push_back(l);

	// rlbwt list
	if(arg.rlbwt)
	{
		rlbwt_lengths.resize(sampled_lengths.size()-i);
		//l = new vector<uint64_t>(sampled_lengths.size()-i,0);
		for(uint64_t j=i;j<sampled_lengths.size();++j)
		{
			rlbwt_lengths[j-i] = sampled_lengths[j];
		}
	}

	cout << "Number of sampled lengths: " << i << endl;
	if(arg.rlbwt) cout << "Run-length BWT sampled lengths: " << rlbwt_lengths.size() << endl;
	cout << "Number of threads: " << threads << endl;
	cout << "Lengths per thread: " << no_k_thread << endl;
	//cout << "base z: " << (int)arg.z << endl;

	// init threads and sketches
	sketch_list.resize(threads);
	thread_list.resize(threads);
	for(uint64_t i=0;i<threads;++i)
	{
		// init ith sketch
		uint64_t curr_window = (*length_lists[i])[length_lists[i]->size()-1];
		//uint64_t curr_window = window_size;
		sketch_list[i] = sketch<>(arg,curr_window,length_lists[i]);
	}
	if(arg.rlbwt){ rlbwt_s = sketch<>(arg,0,&rlbwt_lengths); }

	// main execuction
	i = 0;
	while(true)
	{
		uint64_t loadedChars = loadBuffer(buffer,K);
		if( loadedChars == 0 ){ break; }
		for(uint64_t j=0;j<threads;++j)
		{
			thread_list[j] = thread(extend_window_mt, &sketch_list[j], &buffer, loadedChars);
		}

		if(arg.rlbwt and !rlbwt_s.is_rlbwt_dropped())
		{
			thread t = thread(extend_rlbwt_mt, &rlbwt_s, &buffer, loadedChars);
			t.join();
		}

		// threads sincronization
		for(uint64_t j=0;j<threads;++j)
		{
			thread_list[j].join();
		}

		i += loadedChars;
	}

	// merge all sketches
	for(uint64_t j=1;j<threads;++j)
			sketch_list[0].merge_mt(sketch_list[j]);
	// merge rlbwt sketch	
	if(arg.rlbwt)
		sketch_list[0].merge_mt(rlbwt_s);
  	// resize sketches vector
	sketch_list.resize(1); sketch_list.shrink_to_fit();

	// estimate delta
	print_delta(sketch_list[0]);

	// store sketches if needed
	if(arg.outfile != string())
	{
		ofstream os(arg.outfile);
		sketch_list[0].store(os);
		os.close();
	}
}

/* main */
int main(int argc, char* argv[]){

    Args arg;
    parseArgs(argc, argv, arg);

	if(arg.stream)
	{
			if( arg.threads > 1 )
			{
					stream_delta_parallel(arg);
			}
			else
			{
					stream_delta(arg);
			}
	}
	else if(arg.delta)
	{
			compute_delta(arg.sketch1);
	}
	else if(arg.merge)
	{
			merge_sketches(arg.sketch1,arg.sketch2,arg.outfile);
	}
	else if(arg.ncd)
	{
			merge_sketches(arg.sketch1,arg.sketch2,arg.outfile,true);
	}

    return 0;
}