#ifndef DELTA_STRAM_HPP_
#define DELTA_STRAM_HPP_

#include "common.hpp"

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
    static constexpr uint64_t drop_threshold_ = 100000;

    //static constexpr uint64_t q = (uint64_t(1)<<61) - 1; // prime for Karp-Rabin fingerprinting (M61)
    static constexpr uint64_t q = (uint64_t(1)<<55) + 3; // first prime number after 2^55

    sketch(){// empty constructor
    }

    sketch( uint64_t z, //random base for KR hashing
            uint8_t r,
            uint8_t u = default_u, 
            double a = default_a,
            uint64_t e = default_e) : _z(z), _u(u), _a(a), _r(r), _e(e)
    {
        //TO DO adjust the window size once theory is ready!!
        //for now, sqrt(max stream length)*log(max stream length)
        _window_size = compute_window_size(_u);

        // sample a k value for each sketch
        sample_kmer_lengths(_lengths,_e,_u,_a);

        //std::cout << "Number of sampled lengths: " << _lengths.size() << "\n";
        //std::cout << "Number of registers: " << (int)_r << "\n";
        // init sketches
        _hll_sketches = vector<hll_t>(_lengths.size(),{_r});

        for(uint64_t i=0;i<_lengths.size();++i)
        {
            if(_lengths[i] <= _window_size){ _k++; }
            else
                break;
        }

        // sketch initialization
        // window and window bookmarks initialization
        _window = vector<char_t>(_window_size,0);
        _bookmarks = vector<uint64_t>(_lengths.size(),0);
        // rlbwt and rlbwt bookmarks initialization
        _rlbwt = new rle_bwt();
        _dropped = false;
        // fingerprints initialization
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
        if(_stream_length_ >= _lengths[i]){
         
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
        if(_stream_length_+1 >= _lengths[i])
        {
            _hll_sketches[i].add(reinterpret_cast<char*>(&_fingerprints[i]),sizeof(_fingerprints[i]));
        }
    }

    void extend_rlbwt(char_t c)
    {
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
        _stream_length_++;
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

    uint64_t at_length(uint64_t i)
    {
        assert(i < _lengths.size());
        return _lengths[i];
    }

    //merge this sketch with s
    void merge(const sketch& s)
    {
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
        // store sketches of all k values
        for(uint64_t i=0;i<_hll_sketches.size();++i)
        {
            _hll_sketches[i].dump(os);
        }
    }

    //load sketch from input stream
    void load(istream& is, vector<uint64_t> * lengths)
    {   
        _lengths = *lengths;
        _hll_sketches.resize(_lengths.size());
        // load sketches of all k values
        for(uint64_t i=0;i<_lengths.size();++i)
        {
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
    vector<uint64_t> _bookmarks;
    // store run-length bwt
    rle_bwt* _rlbwt;
    bool _dropped;

    vector<hll_t> _hll_sketches;    // the count-distinct sketches
    vector<uint64_t> _lengths;  // sampled k-mer _lengths corresponding to the sketches

    //current stream length (or sum of streams _lengths if the sketch is the result of union of streams)
    uint64_t _stream_length = 0; 
    uint64_t _stream_length_ = 0; 
    // k-mer window size
    uint64_t _window_size;

    uint64_t _z; // base of Karp-Rabin hashing (should be a random number in (0,q))
    uint8_t _u;  // log2(upper bound to stream length)
    double _a;   //logarithm base for the sampling of factor _lengths
    uint8_t _r;  // log2(number of registers used by each HLL sketch)
    uint64_t _e; // the first default_e k-mer _lengths are 1..default_e
    uint32_t _k = 0; // number of k-mer _lengths smaller than _window_size
};


#endif