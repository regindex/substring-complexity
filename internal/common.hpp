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

const uint64_t e = 20;
const uint64_t u = 32;
const double a = 1.2;

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

void sample_kmer_lengths(vector<uint64_t>& lengths, uint64_t e, uint64_t u, double a)
{
    uint64_t sampled_lenghts = e;
    uint64_t U = uint64_t(1)<<u;
    double exp = 1; //exponential
    uint64_t k = 1; //last sampled length

    while(exp < U){
        if(uint64_t(exp) > k and uint64_t(exp) > e){
            sampled_lenghts++;
            k = exp;
        }
        exp *= a;
    }

    lengths = vector<uint64_t>(sampled_lenghts,0);

    uint64_t i=0;
    while(i++ < e)
        lengths[i-1] = i;
    i--;

    exp = 1;
    k = e;

    while(exp < U){
        if(uint64_t(exp) > k and uint64_t(exp) > e){
            lengths[i++]=exp;
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