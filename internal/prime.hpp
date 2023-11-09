#ifndef PRIME_HPP_
#define PRIME_HPP_

#include "common.hpp"

//{ Simple algorithm running Miller-Rabin primality test
//  source: https://cp-algorithms.com/algebra/primality_tests.html
	uint64_t binpower(uint64_t base, uint64_t e, uint64_t mod) {
	    uint64_t result = 1;
	    base %= mod;
	    while (e) {
	        if (e & 1)
	            result = (__uint128_t)result * base % mod;
	        base = (__uint128_t)base * base % mod;
	        e >>= 1;
	    }
	    return result;
	}

	bool check_composite(uint64_t n, uint64_t a, uint64_t d, int s) {
	    uint64_t x = binpower(a, d, n);
	    if (x == 1 || x == n - 1)
	        return false;
	    for (int r = 1; r < s; r++) {
	        x = (__uint128_t)x * x % n;
	        if (x == n - 1)
	            return false;
	    }
	    return true;
	};

	bool MillerRabin(uint64_t n, uint64_t seed, int iter=5) { // returns true if n is probably prime, else returns false.
		// set seed for random generation
		srand(seed); 

	    if (n < 4)
	        return n == 2 || n == 3;

	    int s = 0;
	    uint64_t d = n - 1;
	    while ((d & 1) == 0) {
	        d >>= 1;
	        s++;
	    }

	    for (int i = 0; i < iter; i++) {
	        int a = 2 + rand() % (n - 3);
	        if (check_composite(n, a, d, s))
	            return false;
	    }
	    return true;
	}

	bool DeterministicMillerRabin(uint64_t n) { // returns true if n is prime, else returns false.
	    if (n < 2)
	        return false;

	    int r = 0;
	    uint64_t d = n - 1;
	    while ((d & 1) == 0) {
	        d >>= 1;
	        r++;
	    }

	    for (int a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
	        if (n == a)
	            return true;
	        if (check_composite(n, a, d, r))
	            return false;
	    }
	    return true;
	}
//}

uint64_t compute_random_prime(uint64_t seed, int bits = 54) 
{
	// compute a random prime number of a certain width
	uint64_t rand_prime = uniform_random((uint64_t(1)<<54),(uint64_t(1)<<55)-1,seed);
	if( rand_prime % 2 == 0 )
		rand_prime++;

	while(true)
	{
		//cout << rand_prime << " " << isprime(rand_prime) << endl;
		if(MillerRabin(rand_prime,seed))
		{
			if(DeterministicMillerRabin(rand_prime))
				break;
		}

		rand_prime += 2;
	}

	return rand_prime;
}

#endif