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
#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#ifndef _MSC_VER
    #ifndef register
        #define register
    #endif
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>

namespace gengraph {

// Max line size in files
#define FBUFF_SIZE 1000000

// disable lousy VC++ warnings
#ifdef _ATL_VER_
    #pragma warning(disable : 4127)
#endif //_ATL_VER_

// Verbose
#define VERBOSE_NONE 0
#define VERBOSE_SOME 1
#define VERBOSE_LOTS 2
int VERBOSE();
void SET_VERBOSE(int v);

// Random number generator
void my_srandom(int);
int my_random();
int my_binomial(double pp, int n);
double my_random01(); // (0,1]

#define MY_RAND_MAX 0x7FFFFFFF

// IPv4 address direct translation into 32-bit uint + special IP defs
typedef unsigned int ip_addr;
#define IP_NONE   0x7FFFFFFF
#define IP_STAR   0x00000000
#define IP_MYSELF 0x7F000001

// Compatibility
#ifdef _WIN32
    #define strcasecmp _stricmp
#endif
//inline double round(double x) throw () { return (floor(0.5+x)); }

// No assert
#ifndef _DEBUG
    #ifndef NDEBUG
        #define NDEBUG
    #endif //NDEBUG
#endif //_DEBUG

// Min & Max
#ifndef min
    #define defmin(type) inline type min(type a, type b) { return a<b ? a : b; }
    defmin(int)
    defmin(double)
    defmin(unsigned long)
#endif //min
#ifndef max
    #define defmax(type) inline type max(type a, type b) { return a>b ? a : b; }
    defmax(int)
    defmax(double)
    defmax(unsigned long)
#endif //max

// Traceroute Sampling
#define MODE_USP 0
#define MODE_ASP 1
#define MODE_RSP 2

// Debug definitions
//#define PERFORMANCE_MONITOR
//#define OPT_ISOLATED

// Max Int
#ifndef MAX_INT
    #define MAX_INT 0x7FFFFFFF
#endif //MAX_INT

//Edge type
typedef struct {
    int from;
    int to;
} edge;

// Tag Int
#define TAG_INT 0x40000000

// Oldies ....
#define S_VECTOR_RAW

//*********************
// Routine definitions
//*********************

/* log(1+x)
inline double logp(double x) {
  if(fabs(x)<1e-6) return x+0.5*x*x+0.333333333333333*x*x*x;
  else return log(1.0+x);
}
//*/


//Fast search or replace
inline int* fast_rpl(int *m, const int a, const int b) {
    while (*m != a) {
        m++;
    }
    *m = b;
    return m;
}
inline int* fast_search(int *m, const int size, const int a) {
    int *p = m + size;
    while (m != p--) if (*p == a) {
            return p;
        }
    return NULL;
}

// Lovely percentage print
// inline void print_percent(double yo, FILE *f = stderr) {
//   int arf = int(100.0*yo);
//   if(double(arf)>100.0*yo) arf--;
//   if(arf<100) fprintf(f," ");
//   if(arf<10) fprintf(f," ");
//   fprintf(f,"%d.%d%%",arf,int(1000.0*yo-double(10*arf)));
// }

// Skips non-numerical chars, then numerical chars, then non-numerical chars.
inline char skip_int(char* &c) {
    while (*c < '0' || *c > '9') {
        c++;
    }
    while (*c >= '0' && *c <= '9') {
        c++;
    }
    while (*c != 0 && (*c < '0' || *c > '9')) {
        c++;
    }
    return *c;
}

// distance+1 modulo 255 for breadth-first search
inline unsigned char next_dist(const unsigned char c) {
    return c == 255 ? 1 : c + 1;
}
inline unsigned char prev_dist(const unsigned char c) {
    return c == 1 ? 255 : c - 1;
}

// 1/(RANDMAX+1)
#define inv_RANDMAX (1.0/(1.0+double(MY_RAND_MAX)))

// random number in ]0,1[, _very_ accurate around 0
inline double random_float() {
    int r = my_random();
    double mul = inv_RANDMAX;
    while (r <= 0x7FFFFF) {
        r <<= 8;
        r += (my_random() & 0xFF);
        mul *= (1.0 / 256.0);
    }
    return double(r) * mul;
}

// Return true with probability p. Very accurate when p is small.
#define test_proba(p) (random_float()<(p))

// Random bit generator, sparwise.
static int _random_bits_stored = 0;
static int _random_bits = 0;

inline int random_bit() {
    register int a = _random_bits;
    _random_bits = a >> 1;
    if (_random_bits_stored--) {
        return a & 0x1;
    }
    a = my_random();
    _random_bits = a >> 1;
    _random_bits_stored = 30;
    return a & 0x1;
}

// Hash Profiling (see hash.h)
void _hash_prof();

} // namespace gengraph

#endif //DEFINITIONS_H
