#ifndef DELTA_STRAM_HPP_
#define DELTA_STRAM_HPP_

#include "common.hpp"

//template on the underlying count-distinct sketch and on the char type
template<class hll_t = HyperLogLogHIP, class char_t = uint8_t>
class sketch{

public:

    //total sketch size is <= (_e + (_u/log2(a))) * 2^_r HLL registers 
    //(_u/log2(a) = log_a(U), where U=2^_u is the upper bound to the stream's length)
    static constexpr uint64_t drop_threshold_ = 100000;

    //static constexpr uint64_t q = (uint64_t(1)<<61) - 1; // prime for Karp-Rabin fingerprinting (M61)
    //static constexpr uint64_t q = (uint64_t(1)<<55) + 3; // first prime number after 2^55
    static constexpr uint64_t q = (uint64_t(1)<<54) + (uint64_t(1)<<7) + 31; // first prime number after 2^54

    sketch(){// empty constructor
    }

    // first constructor computing a list of sampled k values
    sketch( Args& arg )
    {
        //TO DO adjust the window size once theory is ready!!
        //for now, sqrt(max stream length)*log(max stream length)
        _window_size = compute_window_size(arg.u);
        // set KR-base
        _z = arg.z;
        // set logarithm base
        _a = arg.a;
        // set sampling parameter
        _e = arg.e;
        // set number of registers
        _r = arg.registers;
        // set logarithm to stream length upper bound
        _u = arg.u;

        // sample a k value for each sketch
        kmer_lengths_sampling(_lengths,_e,_u,_a,arg.precision);
        // check number of lengths smaller than maximum window size
        _k = 0;
        for(uint64_t i=0;i<_lengths.size();++i)
        {
            if(_lengths[i] <= _window_size){ _k++; }
            else
                break;
        }
        // if rlbwt is disabled drop all k values larger than _window_size
        if( !arg.rlbwt )
        {
            _lengths.resize(_k);
            _lengths.shrink_to_fit();
            assert(_lengths.size() == _k);
        }
        cout << "Number of sampled lengths: " << _lengths.size() << endl;
        //cout << "Number of registers: " << (int)_r << endl;

        // sketch initialization
        _hll_sketches = vector<hll_t>(_lengths.size(),{arg.registers});
        // window and window bookmarks initialization
        _window = vector<char_t>(_window_size,0);
        _bookmarks = vector<uint64_t>(_lengths.size(),0);
        // rlbwt and rlbwt bookmarks initialization
        if( arg.rlbwt )
        {
            _rlbwt = new rle_bwt();
            _dropped = false;
        }
        // fingerprints initialization
        _fingerprints = vector<uint64_t>(_lengths.size(),0);
        _exponents = vector<uint64_t>(_lengths.size(),0);

        // exponents initialization
        for(uint64_t i=0;i<_lengths.size();++i){
            assert(i<_exponents.size());
            _exponents[i] = fexp(_z,_lengths[i]-1,q); //fast exponentiation: z^(_lengths[i]-1) mod q
        }
    }

    // 2nd constructor taking a list of sampling lengths
    sketch( Args& arg, 
            uint64_t window_size,
            vector<uint64_t> * kmer_lengths): _window_size(window_size), _lengths(*kmer_lengths)
    {
        // set KR-base
        _z = arg.z;
        // set number of registers
        _r = arg.registers;
        // initialize sketches
        _hll_sketches = vector<hll_t>(_lengths.size(),{_r});
        _k = _lengths.size();

        // window and window bookmarks initialization
        _bookmarks = vector<uint64_t>(_lengths.size(),0);
        if( _window_size > 0)
        {
            _window = vector<char_t>(_window_size,0);
        }
        // rlbwt initialization
        else
        {
            _rlbwt = new rle_bwt();
            _dropped = false;
        }
        // fingerprints and exponents initialization
        _fingerprints = vector<uint64_t>(_lengths.size(),0);
        _exponents = vector<uint64_t>(_lengths.size(),0);

        // exponents initialization
        for(uint64_t i=0;i<_lengths.size();++i){
            assert(i<_exponents.size());
            _exponents[i] = fexp(_z,_lengths[i]-1,q); //fast exponentiation: z^(_lengths[i]-1) mod q
        }
    }
 
    //void update_fingerprint_window(uint64_t i)
    inline void update_fingerprint_window(uint64_t i, char_t c)
    {
        // remove from the active fingerprints the oldest character
        // and move forward their iterators
        
        if(_stream_length >= _lengths[i]){
        
            char_t b = _window[_bookmarks[i]]; //get oldest character

            uint64_t remove = (b * _exponents[i]) % q;
            _fingerprints[i] = ((_fingerprints[i] + q) - remove) % q;

            // shift forward the bookmark
            _bookmarks[i] = (_bookmarks[i]+1) % _window_size;

        }

        //append new character to all fingerprints:
        //since we have already removed the oldest character (done in previous if),
        //we only need to left-shift fingerprint (multiply by _z) and add c
        _fingerprints[i] = (_fingerprints[i] * _z + c) % q;

        //_stream_length increased by 1 in the fingerprints: if the updated stream length
        // is >= the current k-mer length, insert the k-mer fingerprint in the HLL sketch
        if(_stream_length+1 >= _lengths[i])
        {
            _hll_sketches[i].add(reinterpret_cast<char*>(&_fingerprints[i]),sizeof(_fingerprints[i]));
        }
    }

    void extend_window(char_t c)
    {
        for(uint64_t i=0;i<_k;++i)
        {
            update_fingerprint_window(i,c);
        }

        //if empty stream, leave _window_head to 0 (_window_head always points to
        //most recent stream character). Otherwise, increment it.
        _window_head = (_window_head + (_stream_length>0))%_window_size;

        //store new character in updated _window_head
        _window[_window_head] = c;

        //update stream length
        _stream_length++;
    }

    inline void update_fingerprint_rlbwt(uint64_t i, char_t c)
    {
        // remove from the active fingerprints the oldest character
        // and move forward their iterators   
        if(_stream_length >= _lengths[i]){
         
            char_t b = _rlbwt->at(_bookmarks[i]); //get oldest character

            uint64_t remove = (b * _exponents[i]) % q;
            _fingerprints[i] = ((_fingerprints[i] + q) - remove) % q;

            // shift forward the bookmark
            _bookmarks[i] = _rlbwt->LF(_bookmarks[i]);
            // shift the bookmark if we insert a new bwt character
            // in a preeceding position
            if(_bookmarks[i] >= _rlbwt->get_terminator_position())
                _bookmarks[i]++; 
        }

        //append new character to all fingerprints:
        //since we have already removed the oldest character (done in previous if),
        //we only need to left-shift fingerprint (multiply by _z) and add c
        _fingerprints[i] = (_fingerprints[i] * _z + c) % q;

        //_stream_length increased by 1 in the fingerprints: if the updated stream length
        // is >= the current k-mer length, insert the k-mer fingerprint in the HLL sketch
        if(_stream_length+1 >= _lengths[i])
        {
            _hll_sketches[i].add(reinterpret_cast<char*>(&_fingerprints[i]),sizeof(_fingerprints[i]));
        }
    }

    void extend_rlbwt(char_t c, bool d = false)
    {
        // decrease stream length if we have already inserted a character
        // in the window (for single thread execution)
        if(d){ _stream_length --; }

        for(uint64_t i=_k;i<get_number_of_samples();++i)
        {
            update_fingerprint_rlbwt(i,c);
        }

        //add new character to the rlbwt
        _rlbwt->extend(c);
        if( _rlbwt->number_of_runs() >= drop_threshold_ )
        { 
            _dropped = true;
            delete _rlbwt;
        }

        //update stream length
        _stream_length++;
    }

    void decrease_stream_length()
    {
        _stream_length--;
    }

    uint64_t no_runs()
    {
        return _rlbwt->number_of_runs();
    }

    bool is_rlbwt_dropped()
    {
        return _dropped;
    }

    hll_t at_sketch(uint64_t i) const
    {
        assert(i < _lengths.size());
        return _hll_sketches[i];
    }

    uint64_t at_length(uint64_t i) const
    {
        assert(i < _lengths.size());
        return _lengths[i];
    }

    //merge this sketch with s
    void merge(const sketch& s)
    {
        // update stream length
        _stream_length += s.stream_length();
        // merge count-distinct sketches
        for(uint64_t i=0;i<_lengths.size();++i)
        {
            assert(_lengths[i]==s.at_length(i));
            _hll_sketches[i].merge(s.at_sketch(i));
        }
    }

    //get estimate of measure delta
    double estimate_delta() const {
        
        double delta = 0;

        for(uint64_t i=0;i<_lengths.size();++i)
            delta = max(delta,_hll_sketches[i].estimate()/_lengths[i]);

        return delta;
    }

    //store sketch to output stream
    void store(ostream& os) const
    {
        // store sketch parameters
        os << _stream_length << '\n';
        os << _lengths.size() << '\n';
        // store sketches of all k values
        for(uint64_t i=0;i<_hll_sketches.size();++i)
        {
            os << _lengths[i];
            _hll_sketches[i].dump(os);
        }
    }

    //load sketch from input stream
    void load(istream& is)
    {   
        // read stream length
        is >> _stream_length;
        // read input sketch
        int no_sketches = 0;
        is >> no_sketches;
        // initialize sketch and length vectors
        _lengths.resize(no_sketches);
        _hll_sketches.resize(_lengths.size());
        // load sketches of all k values
        for(uint64_t i=0;i<_lengths.size();++i)
        {
            is >> _lengths[i]; 
            _hll_sketches[i].restore(is);
        }
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

    //return log base
    double get_log_base() const {
        return _a;
    }

    //return log base
    uint64_t get_sampling_param() const {
        return _e;
    }

    //number of sampled factor _lengths (for each, we store a HLL sketch)
    uint64_t get_number_of_samples(){
        return _lengths.size();
    }

    void clear_sketch()
    {
        // delete window and bookmarks
        _window.clear(); _bookmarks.clear();
        // delete count-distinct sketches and sampled lengths
        _hll_sketches.clear(); _lengths.clear();
        // delete fingerprints and exponents vectors
        _fingerprints.clear(); _exponents.clear();
    }

    //merge this sketch with s MT
    void merge_mt(sketch& s)
    {
        for(uint64_t i=0;i<s.get_number_of_samples();++i)
        {
            _lengths.push_back(s.at_length(i));
            _hll_sketches.push_back(s.at_sketch(i));
        }
        // clear sketch
        s.clear_sketch();
    }

private:

    /*
     * sketch parameters
    */

    vector<uint64_t> _fingerprints; // Karp-Rabin fingerprints of the windows
    vector<uint64_t> _exponents; // z^_lengths[0], z^_lengths[1], z^_lengths[2], ...

    vector<char_t> _window; // window containing the last _window_size characters of the stream
    uint64_t _window_head = 0; // head of the window
    vector<uint64_t> _bookmarks; // window bookmarks vector
    rle_bwt* _rlbwt; // dynamic run-length bwt
    bool _dropped; // bool variable indicating if the rlbwt has been dropped

    vector<hll_t> _hll_sketches; // the count-distinct sketches
    vector<uint64_t> _lengths;   // sampled k-mer _lengths corresponding to the sketches

    uint64_t _stream_length = 0; //current stream length (or sum of streams _lengths if the sketch is the result of union of streams)
    uint64_t _window_size; // kmer window size

    uint64_t _z; // base of Karp-Rabin hashing (should be a random number in (0,q))
    uint8_t _u;  // log2(upper bound to stream length)
    double _a;   //logarithm base for the sampling of factor _lengths
    uint8_t _r;  // log2(number of registers used by each HLL sketch)
    uint64_t _e; // the sampling patameter
    uint32_t _k; // number of k-mer _lengths smaller than _window_size
};


#endif