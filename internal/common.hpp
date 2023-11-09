#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <iostream>
#include <vector>
#include <random>
#include <thread>
#include <fstream>

#include "hyperloglog.hpp"
#include "dynamic/dynamic.hpp"

using namespace std;
using namespace hll;
using namespace dyn; 

// struct containing command line parameters and other globals
struct Args
{
  double e = 20;
  uint64_t u = 26;
  double a = 1.2;
  uint8_t registers = 10;
  bool stream = false, rlbwt = false, merge = false;
  bool delta = false, ncd = false;
  bool fast = false;
  int threads = 1;
  string outfile = string();
  string sketch1 = string(), sketch2 = string();
  int precision = 2;
  uint64_t prime = 0;
  uint64_t buffer = 1;
  bool verbose = false; 
};

// configuration object
struct conf{

    conf(){}

    conf(double e, double a, uint8_t r, int p): _e(e), _a(a), _r(r), _p(p)
    {}

    double _e = 10;double _a = 1.4;uint8_t _r = 10;int _p=4;
};
// configurations
//vector<conf> conf_list{ conf(9,1.5,12,4), conf(15,1.3,12,4), conf(10,1.4,14,3), conf(14,1.3,14,3), 
vector<conf> conf_list{ conf(7.5,1.45,12,4), conf(10.5,1.35,14,3), conf(11,1.28,16,2), conf(15,1.2,16,2), conf(15,1.1,18,1) };
// 7.5,1.45
template<typename T, typename M>
T uniform_random(T range_from, T range_to, M seed) {
    //std::random_device                  rand_dev;
    std::mt19937                        generator(seed);
    std::uniform_int_distribution<T>    distr(range_from, range_to);
    return distr(generator);
}

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

// sample kmer lengths
void kmer_lengths_sampling(vector<uint64_t>& lengths, double e, uint64_t u, double a, uint64_t p)
{
    uint64_t sampled_lenghts = e;
    uint64_t U = uint64_t(1)<<u;
    double exp = 1; //exponential
    uint64_t k = 1; //last sampled length
    // up to e*2^1
    uint64_t i=1, j=2;
    if( p%2 != 0 ){ j--; }
    if( p==1 ){ i--; }
    //uint64_t i=0, j=0;
    while( i++ < e*(2) )
    {
        if( (++j)%p == 0 )
            lengths.push_back(i);
    }
    // up to e*2^2
    i--; j = 0;
    while( i++ < (e*(2*2)) )
    {
        if( (++j)%(p+1) == 0 )
            lengths.push_back(i);
    }
    // up to e*2^3
    i--; j = 0;
    while( i++ < (e*(2*2*2)) ) 
    {
        if( (++j)%(p+2) == 0 )
            lengths.push_back(i);
    }
    // exponential selection
    exp = 1;
    k = e*(2*2*2);
    while(exp < U){
        if(ceil(exp) > k and ceil(exp) > e){
            lengths.push_back(exp);
            i++;
            k = ceil(exp);
        }
        exp *= a;
    }
}

static const uint64_t compute_window_size(uint64_t u)
{
    return (uint64_t(1)<<(u/2))*u;
}

static constexpr uint8_t min_k_thread = 1;

#endif