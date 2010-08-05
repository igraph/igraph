#include "gengraph_definitions.h"
#include <cstdlib>
#include <stdio.h>

#include "gengraph_random.h"

namespace gengraph {

static KW_RNG::RNG _my_random;
int my_random() { return _my_random.rand_int31(); }
void my_srandom(int x) { _my_random.init(x,!x*13,x*x+1,(x>>16)+(x<<16)); }
int my_binomial(double pp, int n) { return _my_random.binomial(pp,n); }
double my_random01() { return _my_random.rand_halfopen01(); }

}

namespace gengraph {

static int VERB;
int VERBOSE() { return VERB; }
void SET_VERBOSE(int v) { VERB = v; } 

//Hash profiling
static unsigned long _hash_rm_i   = 0;
static unsigned long _hash_rm_c   = 0;
static unsigned long _hash_add_i  = 0;
static unsigned long _hash_add_c  = 0;
static unsigned long _hash_put_i  = 0;
static unsigned long _hash_put_c  = 0;
static unsigned long _hash_find_i = 0;
static unsigned long _hash_find_c = 0;
static unsigned long _hash_rand_i = 0;
static unsigned long _hash_rand_c = 0;
static unsigned long _hash_expand = 0;
inline void _hash_add_iter()  { _hash_add_i++;  }
inline void _hash_add_call()  { _hash_add_c++;  }
inline void _hash_put_iter()  { _hash_put_i++;  }
inline void _hash_put_call()  { _hash_put_c++;  }
inline void _hash_rm_iter()   { _hash_rm_i++;   }
inline void _hash_rm_call()   { _hash_rm_c++;   }
inline void _hash_find_iter() { _hash_find_i++; }
inline void _hash_find_call() { _hash_find_c++; }
inline void _hash_rand_iter() { _hash_rand_i++; }
inline void _hash_rand_call() { _hash_rand_c++; }
inline void _hash_expand_call() { _hash_expand++; }
void _hash_prof() {
  fprintf(stderr,"HASH_ADD : %lu / %lu\n", _hash_add_c , _hash_add_i);
  fprintf(stderr,"HASH_PUT : %lu / %lu\n", _hash_put_c , _hash_put_i);
  fprintf(stderr,"HASH_FIND: %lu / %lu\n", _hash_find_c, _hash_find_i);
  fprintf(stderr,"HASH_RM  : %lu / %lu\n", _hash_rm_c  , _hash_rm_i);
  fprintf(stderr,"HASH_RAND: %lu / %lu\n", _hash_rand_c, _hash_rand_i);
  fprintf(stderr,"HASH_EXPAND : %lu calls\n", _hash_expand);
}

} // namespace gengraph
