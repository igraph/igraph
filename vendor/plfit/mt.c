/* mt.c
 *
 * Mersenne Twister random number generator, based on the implementation of
 * Michael Brundage (which has been placed in the public domain).
 *
 * Author: Tamas Nepusz (original by Michael Brundage)
 *
 * See the following URL for the original implementation:
 * http://www.qbrundage.com/michaelb/pubs/essays/random_number_generation.html
 *
 * This file has been placed in the public domain.
 */

#include <stdlib.h>

#include "igraph_random.h"

#include "plfit_mt.h"

static uint16_t get_random_uint16(void) {
    return RNG_INT31() & 0xFFFF;
}

void plfit_mt_init(plfit_mt_rng_t* rng) {
    plfit_mt_init_from_rng(rng, 0);
}

void plfit_mt_init_from_rng(plfit_mt_rng_t* rng, plfit_mt_rng_t* seeder) {
    int i;

    if (seeder == 0) {
        for (i = 0; i < PLFIT_MT_LEN; i++) {
            /* RAND_MAX is guaranteed to be at least 32767, so we can use two
             * calls to rand() to produce a random 32-bit number */
            rng->mt_buffer[i] = (get_random_uint16() << 16) + get_random_uint16();
        }
    } else {
        for (i = 0; i < PLFIT_MT_LEN; i++) {
            rng->mt_buffer[i] = plfit_mt_random(seeder);
        }
    }

    rng->mt_index = 0;
}

#define MT_IA           397
#define MT_IB           (PLFIT_MT_LEN - MT_IA)
#define UPPER_MASK      0x80000000
#define LOWER_MASK      0x7FFFFFFF
#define MATRIX_A        0x9908B0DF
#define TWIST(b,i,j)    ((b)[i] & UPPER_MASK) | ((b)[j] & LOWER_MASK)
#define MAGIC(s)        (((s)&1)*MATRIX_A)

uint32_t plfit_mt_random(plfit_mt_rng_t* rng) {
    uint32_t * b = rng->mt_buffer;
    int idx = rng->mt_index;
    uint32_t s;
    int i;
	
    if (idx == PLFIT_MT_LEN * sizeof(uint32_t)) {
        idx = 0;
        i = 0;
        for (; i < MT_IB; i++) {
            s = TWIST(b, i, i+1);
            b[i] = b[i + MT_IA] ^ (s >> 1) ^ MAGIC(s);
        }
        for (; i < PLFIT_MT_LEN-1; i++) {
            s = TWIST(b, i, i+1);
            b[i] = b[i - MT_IB] ^ (s >> 1) ^ MAGIC(s);
        }
        
        s = TWIST(b, PLFIT_MT_LEN-1, 0);
        b[PLFIT_MT_LEN-1] = b[MT_IA-1] ^ (s >> 1) ^ MAGIC(s);
    }

    rng->mt_index = idx + sizeof(uint32_t);
    return *(uint32_t *)((unsigned char *)b + idx);
    /*
    Matsumoto and Nishimura additionally confound the bits returned to the caller
    but this doesn't increase the randomness, and slows down the generator by
    as much as 25%.  So I omit these operations here.
    
    r ^= (r >> 11);
    r ^= (r << 7) & 0x9D2C5680;
    r ^= (r << 15) & 0xEFC60000;
    r ^= (r >> 18);
    */
}


double plfit_mt_uniform_01(plfit_mt_rng_t* rng) {
    return ((double)plfit_mt_random(rng)) / PLFIT_MT_RAND_MAX;
}
