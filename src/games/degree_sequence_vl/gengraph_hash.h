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
#ifndef HASH_H
#define HASH_H

#include <assert.h>
#include "gengraph_definitions.h"

//_________________________________________________________________________
// Hash table profiling... Active only if definition below is uncommented
//_________________________________________________________________________
//#define _HASH_PROFILE

namespace gengraph {

#ifdef _HASH_PROFILE
    void _hash_add_iter();
    void _hash_add_call();
    void _hash_put_iter();
    void _hash_put_call();
    void _hash_rm_iter();
    void _hash_rm_call();
    void _hash_find_iter();
    void _hash_find_call();
    void _hash_rand_iter();
    void _hash_rand_call();
    void _hash_expand_call();
    void _hash_prof();
    #define _HASH_ADD_ITER()  _hash_add_iter()
    #define _HASH_ADD_CALL()  _hash_add_call()
    #define _HASH_PUT_ITER()  _hash_put_iter()
    #define _HASH_PUT_CALL()  _hash_put_call()
    #define _HASH_RM_ITER()   _hash_rm_iter()
    #define _HASH_RM_CALL()   _hash_rm_call()
    #define _HASH_FIND_ITER() _hash_find_iter()
    #define _HASH_FIND_CALL() _hash_find_call()
    #define _HASH_RAND_ITER() _hash_rand_iter()
    #define _HASH_RAND_CALL() _hash_rand_call()
    #define _HASH_EXP_CALL()  _hash_expand_call()
#else
    #define _HASH_ADD_ITER()  {}
    #define _HASH_ADD_CALL()  {}
    #define _HASH_PUT_ITER()  {}
    #define _HASH_PUT_CALL()  {}
    #define _HASH_RM_ITER()   {}
    #define _HASH_RM_CALL()   {}
    #define _HASH_FIND_ITER() {}
    #define _HASH_FIND_CALL() {}
    #define _HASH_RAND_ITER() {}
    #define _HASH_RAND_CALL() {}
    #define _HASH_EXP_CALL()  {}
#endif

//_________________________________________________________________________
// Hash Table properties. Works best when HASH_SIZE_IS_POWER2 is uncommented
// but takes 2.25 times the needed space, in average (from 1.5 to 3)
// If you have memory issues, Try to comment it: tables will take 1.5 times
// the minimal space
//_________________________________________________________________________

#define HASH_SIZE_IS_POWER2
#define MACRO_RATHER_THAN_INLINE

// under HASH_MIN_SIZE, vectors are not hash table (just a simle array)
#define HASH_MIN_SIZE 100
#define IS_HASH(x) ((x)>HASH_MIN_SIZE)
#define HASH_NONE (-1)

#ifdef HASH_SIZE_IS_POWER2
inline igraph_integer_t HASH_EXPAND(igraph_integer_t x) {
    /* Returns pow(2, floor(log2(x)) + 2) if x > 0, 1 otherwise. Works up to
     * x == 2^64, starts to break down afterwards */
    _HASH_EXP_CALL();
    x += x;
    x |= x >> 1;  x |= x >> 2;  x |= x >> 4;  x |= x >> 8;  x |= x >> 16;
#if IGRAPH_INTEGER_SIZE == 64
    x |= x >> 32;
#endif
    return x + 1;
}
#define HASH_KEY(x,size) ((x*2198737)&((size)-1))
#endif //HASH_SIZE_IS_POWER2

#ifdef MACRO_RATHER_THAN_INLINE
#ifndef HASH_SIZE_IS_POWER2
    #define HASH_EXPAND(x) ((x)+((x)>>1))
    #define HASH_UNEXPAND(x) ((((x)<<1)+1)/3)
    #define HASH_KEY(x,size) ((x)%(size))
#endif //HASH_SIZE_IS_POWER2
#define HASH_SIZE(x) (IS_HASH(x) ? HASH_EXPAND(x) : (x) )
#define HASH_REKEY(k,size) ((k)==0 ? (size)-1 : (k)-1)
#else //MACRO_RATHER_THAN_INLINE
#ifndef HASH_SIZE_IS_POWER2
inline igraph_integer_t  HASH_KEY(igraph_integer_t x, igraph_integer_t size) {
    assert(x >= 0);
    return x % size;
};
inline igraph_integer_t  HASH_EXPAND(igraph_integer_t x) {
    _HASH_EXP_CALL();
    return x + (x >> 1);
};
inline int  HASH_UNEXPAND(igraph_integer_t x) {
    return ((x << 1) + 1) / 3;
};
#endif //HASH_SIZE_IS_POWER2
inline int  HASH_REKEY(igraph_integer_t k, igraph_integer_t s) {
    assert(k >= 0);
    if (k == 0) {
        return s - 1;
    } else {
        return k - 1;
    }
};
inline int  HASH_SIZE(igraph_integer_t x) {
    if (IS_HASH(x)) {
        return HASH_EXPAND(x);
    } else {
        return x;
    }
};
#endif //MACRO_RATHER_THAN_INLINE

inline igraph_integer_t HASH_PAIR_KEY(igraph_integer_t x, igraph_integer_t y, igraph_integer_t size) {
    return HASH_KEY(x * 1434879443 + y, size);
}

//_________________________________________________________________________
// Hash-only functions : table must NOT be Raw.
// the argument 'size' is the total size of the hash table
//_________________________________________________________________________

// copy hash table into raw vector
inline void H_copy(igraph_integer_t *mem, igraph_integer_t *h, igraph_integer_t size) {
    for (igraph_integer_t i = HASH_EXPAND(size); i--; h++) {
        if (*h != HASH_NONE) {
            *(mem++) = *h;
        }
    }
}

// Look for the place to add an element. Return NULL if element is already here.
inline igraph_integer_t* H_add(igraph_integer_t* h, igraph_integer_t size, igraph_integer_t a) {
    _HASH_ADD_CALL();
    _HASH_ADD_ITER();
    igraph_integer_t k = HASH_KEY(a, size);
    if (h[k] == HASH_NONE) {
        return h + k;
    }
    while (h[k] != a) {
        _HASH_ADD_ITER();
        k = HASH_REKEY(k, size);
        if (h[k] == HASH_NONE) {
            return h + k;
        }
    }
    return NULL;
}

// would element be well placed in newk ?
inline bool H_better(igraph_integer_t a, igraph_integer_t size, igraph_integer_t currentk, igraph_integer_t newk) {
    igraph_integer_t k = HASH_KEY(a, size);
    if (newk < currentk) {
        return (k < currentk && k >= newk);
    } else {
        return (k < currentk || k >= newk);
    }
}

// removes h[k]
inline void H_rm(igraph_integer_t* h, igraph_integer_t size, igraph_integer_t k) {
    _HASH_RM_CALL();
    igraph_integer_t lasthole = k;
    do {
        _HASH_RM_ITER();
        k = HASH_REKEY(k, size);
        igraph_integer_t next = h[k];
        if (next == HASH_NONE) {
            break;
        }
        if (H_better(next, size, k, lasthole)) {
            h[lasthole] = next;
            lasthole = k;
        }
    } while (true);
    h[lasthole] = HASH_NONE;
}

//put a
inline igraph_integer_t* H_put(igraph_integer_t* h, igraph_integer_t size, igraph_integer_t a) {
    assert(H_add(h, size, a) != NULL);
    _HASH_PUT_CALL();
    _HASH_PUT_ITER();
    igraph_integer_t k = HASH_KEY(a, size);
    while (h[k] != HASH_NONE) {
        k = HASH_REKEY(k, size);
        _HASH_PUT_ITER();
    }
    h[k] = a;
    assert(H_add(h, size, a) == NULL);
    return h + k;
}

// find A
inline igraph_integer_t H_find(igraph_integer_t *h, igraph_integer_t size, igraph_integer_t a) {
    assert(H_add(h, size, a) == NULL);
    _HASH_FIND_CALL();
    _HASH_FIND_ITER();
    igraph_integer_t k = HASH_KEY(a, size);
    while (h[k] != a) {
        k = HASH_REKEY(k, size);
        _HASH_FIND_ITER();
    }
    return k;
}

// Look for the place to add an element. Return NULL if element is already here.
inline bool H_pair_insert(igraph_integer_t* h, igraph_integer_t size, igraph_integer_t a, igraph_integer_t b) {
    _HASH_ADD_CALL();
    _HASH_ADD_ITER();
    igraph_integer_t k = HASH_PAIR_KEY(a, b, size);
    if (h[2 * k] == HASH_NONE) {
        h[2 * k] = a;
        h[2 * k + 1] = b;
        return true;
    }
    while (h[2 * k] != a || h[2 * k + 1] != b) {
        _HASH_ADD_ITER();
        k = HASH_REKEY(k, size);
        if (h[2 * k] == HASH_NONE) {
            h[2 * k] = a;
            h[2 * k + 1] = b;
            return true;
        }
    }
    return false;
}


//_________________________________________________________________________
// Generic functions : table can be either Hash or Raw.
// the argument 'size' is the number of elements
//_________________________________________________________________________

// Look for an element
inline bool H_is(igraph_integer_t *mem, igraph_integer_t size, igraph_integer_t elem) {
    if (IS_HASH(size)) {
        return (H_add(mem, HASH_EXPAND(size), elem) == NULL);
    } else {
        return fast_search(mem, size, elem) != NULL;
    }
}

//pick random location (containing an element)
inline igraph_integer_t* H_random(igraph_integer_t* mem, igraph_integer_t size) {
    if (!IS_HASH(size)) {
        return mem + (my_random() % size);
    }
    _HASH_RAND_CALL();
    size = HASH_EXPAND(size);
    igraph_integer_t* yo;
    do {
        yo = mem + HASH_KEY(my_random(), size);
        _HASH_RAND_ITER();
    } while (*yo == HASH_NONE);
    return yo;
}

// replace *k by b
inline igraph_integer_t* H_rpl(igraph_integer_t *mem, igraph_integer_t size, igraph_integer_t* k, igraph_integer_t b) {
    assert(!H_is(mem, size, b));
    if (!IS_HASH(size)) {
        *k = b;
        return k;
    } else {
        size = HASH_EXPAND(size);
        assert(mem + int(k - mem) == k);
        H_rm(mem, size, int(k - mem));
        return H_put(mem, size, b);
    }
}

// replace a by b
inline igraph_integer_t* H_rpl(igraph_integer_t *mem, igraph_integer_t size, igraph_integer_t a, igraph_integer_t b) {
    assert(H_is(mem, size, a));
    assert(!H_is(mem, size, b));
    if (!IS_HASH(size)) {
        return fast_rpl(mem, a, b);
    } else {
        size = HASH_EXPAND(size);
        H_rm(mem, size, H_find(mem, size, a));
        return H_put(mem, size, b);
    }
}

} // namespace gengraph

#endif //HASH_H
