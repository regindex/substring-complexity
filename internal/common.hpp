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

const uint8_t min_k_thread = 4;

// struct containing command line parameters and other globals
struct Args
{
  uint64_t e = 20, u = 32;
  double a = 1.2;
  uint8_t registers = 10, z=0;
  bool stream = false, rlbwt = false, merge = false;
  bool delta = false, ncd = false;
  int threads = 1;
  string outfile = string();
  string sketch1 = string(), sketch2 = string();
  int precision = 3;
  bool verbose = false; 
};

// configuration object
struct conf{

    conf(){}

    conf(uint64_t e, double a, uint8_t r, int p): _e(e), _a(a), _r(r), _p(p)
    {}

    uint64_t _e = 10;double _a = 1.4;uint8_t _r = 10;int _p=4;
};
// configurations
vector<conf> conf_list{ conf(9,1.5,12,4), conf(15,1.3,12,4), conf(10,1.4,14,3), conf(14,1.3,14,3), 
                        conf(11,1.30,16,2), conf(15,1.2,16,2), conf(15,1.1,18,1) };

template<typename T>
T random(T range_from, T range_to) {
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
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
void kmer_lengths_sampling(vector<uint64_t>& lengths, uint64_t e, uint64_t u, double a, uint64_t p)
{
    uint64_t sampled_lenghts = e;
    uint64_t U = uint64_t(1)<<u;
    double exp = 1; //exponential
    uint64_t k = 1; //last sampled length
    // up to e*2^1
    uint64_t i=1, j=2;
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
        if(uint64_t(exp) > k and uint64_t(exp) > e){
            lengths.push_back(exp);
            i++;
            k = exp;
        }
        exp *= a;
    }
}

static const uint64_t compute_window_size(uint64_t u)
{
    return (uint64_t(1)<<(u/2))*u;
}

#endif