/*
 *
 * gengraph - generation of random simple connected graphs with prescribed
 *            degree sequence
 *
 * Copyright (C) 2006  Fabien Viger
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef RNG_H
#define RNG_H

#include "igraph_random.h"

namespace KW_RNG {

typedef signed int  sint;
typedef unsigned int uint;
typedef signed long  slong;
typedef unsigned long ulong;

class RNG {
public:
    RNG() { }
    RNG(ulong z_, ulong w_, ulong jsr_, ulong jcong_ ) {
        IGRAPH_UNUSED(z_); IGRAPH_UNUSED(w_); IGRAPH_UNUSED(jsr_);
        IGRAPH_UNUSED(jcong_);
    };
    ~RNG() { }

    void init(ulong z_, ulong w_, ulong jsr_, ulong jcong_ ) {
        IGRAPH_UNUSED(z_); IGRAPH_UNUSED(w_); IGRAPH_UNUSED(jsr_);
        IGRAPH_UNUSED(jcong_);
    }
    long rand_int31() {
        return RNG_INT31();
    }
    double rand_halfopen01() { // (0,1]
        return RNG_UNIF01();
    }
    int binomial(double pp, int n) {
        return RNG_BINOM(n, pp);
    }
};

} // namespace KW_RNG

/* This was the original RNG, but now we use the igraph version */

// __________________________________________________________________________
// random.h   - a Random Number Generator Class
// random.cpp - contains the non-inline class methods

// __________________________________________________________________________
// This C++ code uses the simple, very fast "KISS" (Keep It Simple
// Stupid) random number generator suggested by George Marsaglia in a
// Usenet posting from 1999.  He describes it as "one of my favorite
// generators".  It generates high-quality random numbers that
// apparently pass all commonly used tests for randomness.  In fact, it
// generates random numbers by combining the results of three other good
// random number generators that have different periods and are
// constructed from completely different algorithms.  It does not have
// the ultra-long period of some other generators - a "problem" that can
// be fixed fairly easily - but that seems to be its only potential
// problem.  The period is about 2^123.

// The ziggurat method of Marsaglia is used to generate exponential and
// normal variates.  The method as well as source code can be found in
// the article "The Ziggurat Method for Generating Random Variables" by
// Marsaglia and Tsang, Journal of Statistical Software 5, 2000.

// The method for generating gamma variables appears in "A Simple Method
// for Generating Gamma Variables" by Marsaglia and Tsang, ACM
// Transactions on Mathematical Software, Vol. 26, No 3, Sep 2000, pages
// 363-372.

// The code for Poisson and Binomial random numbers comes from
// Numerical Recipes in C.

// Some of this code is unlikely to work correctly as is on 64 bit
// machines.

// #include <cstdlib>
// #include <ctime>
// #ifdef _WIN32
// #include <process.h>
// #define getpid _getpid
// #else
// #include <unistd.h>
// #endif

// //#ifdef _WIN32
//   static const double PI   =  3.1415926535897932;
//   static const double AD_l =  0.6931471805599453;
//   static const double AD_a =  5.7133631526454228;
//   static const double AD_b =  3.4142135623730950;
//   static const double AD_c = -1.6734053240284925;
//   static const double AD_p =  0.9802581434685472;
//   static const double AD_A =  5.6005707569738080;
//   static const double AD_B =  3.3468106480569850;
//   static const double AD_H =  0.0026106723602095;
//   static const double AD_D =  0.0857864376269050;
// //#endif //_WIN32

// namespace KW_RNG {

// class RNG
// {
// private:
//   ulong z, w, jsr, jcong; // Seeds

//   ulong kn[128], ke[256];
//   double wn[128],fn[128], we[256],fe[256];

// /*
// #ifndef _WIN32
//   static const double PI   =  3.1415926535897932;
//   static const double AD_l =  0.6931471805599453;
//   static const double AD_a =  5.7133631526454228;
//   static const double AD_b =  3.4142135623730950;
//   static const double AD_c = -1.6734053240284925;
//   static const double AD_p =  0.9802581434685472;
//   static const double AD_A =  5.6005707569738080;
//   static const double AD_B =  3.3468106480569850;
//   static const double AD_H =  0.0026106723602095;
//   static const double AD_D =  0.0857864376269050;
// #endif //_WIN32
// */

// public:
//   RNG() { init(); zigset(); }
//   RNG(ulong z_, ulong w_, ulong jsr_, ulong jcong_ ) :
//     z(z_), w(w_), jsr(jsr_), jcong(jcong_) { zigset(); }
//   ~RNG() { }


//   inline ulong znew()
//     { return (z = 36969 * (z & 65535) + (z >> 16)); }
//   inline ulong wnew()
//     { return (w = 18000 * (w & 65535) + (w >> 16)); }
//   inline ulong MWC()
//     { return (((znew() & 65535) << 16) + wnew()); }
//   inline ulong SHR3()
//     { jsr ^= ((jsr & 32767) << 17); jsr ^= (jsr >> 13); return (jsr ^= ((jsr << 5) & 0xFFFFFFFF)); }
//   inline ulong CONG()
//     { return (jcong = (69069 * jcong + 1234567) & 0xFFFFFFFF); }
//   inline double RNOR() {
//     slong h = rand_int32();
//     ulong i = h & 127;
//     return (((ulong) abs((sint) h) < kn[i]) ? h * wn[i] : nfix(h, i));
//   }
//   inline double REXP() {
//     ulong j = rand_int32();
//     ulong i = j & 255;
//     return ((j < ke[i]) ? j * we[i] : efix(j, i));
//   }

//   double nfix(slong h, ulong i);
//   double efix(ulong j, ulong i);
//   void zigset();

//   inline void init()
//     { ulong yo = time(0) + getpid();
//       z = w = jsr = jcong = yo; }
//   inline void init(ulong z_, ulong w_, ulong jsr_, ulong jcong_ )
//     { z = z_; w = w_; jsr = jsr_; jcong = jcong_; }

//   inline ulong rand_int32()         // [0,2^32-1]
//     { return ((MWC() ^ CONG()) + SHR3()) & 0xFFFFFFFF; }
//   inline long rand_int31()          // [0,2^31-1]
//     { return long(rand_int32() >> 1);}
//   inline double rand_closed01()     // [0,1]
//     { return ((double) rand_int32() / 4294967295.0); }
//   inline double rand_open01()       // (0,1)
//     { return (((double) rand_int32() + 0.5) / 4294967296.0); }
//   inline double rand_halfclosed01() // [0,1)
//     { return ((double) rand_int32() / 4294967296.0); }
//   inline double rand_halfopen01()   // (0,1]
//     { return (((double) rand_int32() + 0.5) / 4294967295.5); }

//   // Continuous Distributions
//   inline double uniform(double x = 0.0, double y = 1.0)
//     { return rand_closed01() * (y - x) + x; }
//   inline double normal(double mu = 0.0, double sd = 1.0)
//     { return RNOR() * sd + mu; }
//   inline double exponential(double lambda = 1)
//     { return REXP() / lambda; }
//   double gamma(double shape = 1, double scale = 1);
//   double chi_square(double df)
//     { return gamma(df / 2.0, 0.5); }
//   double beta(double a1, double a2)
//     { double x1 = gamma(a1, 1); return (x1 / (x1 + gamma(a2, 1))); }

//   // Discrete Distributions
//   double poisson(double lambda);
//   int binomial(double pp, int n);

// }; // class RNG

// } // namespace

#endif // RNG_H
