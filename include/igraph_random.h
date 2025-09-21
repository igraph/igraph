/*
   igraph library.
   Copyright (C) 2003-2025  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef IGRAPH_RANDOM_H
#define IGRAPH_RANDOM_H

#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"

#include <stdint.h>
#include <stdlib.h>

IGRAPH_BEGIN_C_DECLS

/* The new RNG interface is (somewhat) modelled on the GSL */

/* When implementing your own RNG in igraph, the following methods must be
 * supplied in the corresponding igraph_rng_type_t structure:
 *
 * - init()
 * - destroy()
 * - seed()
 * - get()
 *
 * Optionally, you can provide specialized routines for several distributions
 * in the following functions:
 *
 * - get_int()
 * - get_real()
 * - get_norm()
 * - get_geom()
 * - get_binom()
 * - get_exp()
 * - get_gamma()
 * - get_pois()
 *
 * The best is probably to define get() leave the others as NULL; igraph will use
 * default implementations for these.
 *
 * Note that if all that you would do in a get_real() implementation is to
 * generate random bits with get() and divide by the maximum, don't do that;
 * The default implementation takes care of calling get() a sufficient number of
 * times to utilize most of the precision of the igraph_real_t type, and generate
 * accurate variates. Inaccuracies in the output of get_real() can get magnified
 * when using the default generators for non-uniform distributions.
 * When implementing get_real(), the sampling range must be half-open, i.e. [0, 1).
 * If unsure, leave get_real() unimplemented and igraph will provide an implementation
 * in terms of get().
 *
 * When implementing get_int(), you do not need to check whether lo < hi;
 * the caller is responsible for ensuring that this is the case. You can always
 * assume that hi > lo. Note that both endpoints are _inclusive_, and you must
 * make sure that your generation scheme works for both 32-bit and 64-bit
 * versions of igraph_int_t as igraph can be compiled for both cases. If
 * you are unsure, leave get_int() unimplemented and igraph will provide its
 * own implementation based on get().
 */
typedef struct igraph_rng_type_t {
    const char *name;
    uint8_t bits;

    /* Initialization and destruction */
    igraph_error_t (*init)(void **state);
    void (*destroy)(void *state);

    /* Seeding */
    igraph_error_t (*seed)(void *state, igraph_uint_t seed);

    /* Fundamental generator: return as many random bits as the RNG supports in
     * a single round */
    igraph_uint_t (*get)(void *state);

    /* Optional generators; defaults are provided by igraph that rely solely
     * on get() */
    igraph_int_t (*get_int)(void *state, igraph_int_t l, igraph_int_t h);
    igraph_real_t (*get_real)(void *state);
    igraph_real_t (*get_norm)(void *state);
    igraph_real_t (*get_geom)(void *state, igraph_real_t p);
    igraph_real_t (*get_binom)(void *state, igraph_int_t n, igraph_real_t p);
    igraph_real_t (*get_exp)(void *state, igraph_real_t rate);
    igraph_real_t (*get_gamma)(void *state, igraph_real_t shape,
                               igraph_real_t scale);
    igraph_real_t (*get_pois)(void *state, igraph_real_t mu);
} igraph_rng_type_t;

typedef struct igraph_rng_t {
    const igraph_rng_type_t *type;
    void *state;
    igraph_bool_t is_seeded;
} igraph_rng_t;

/* --------------------------------- */

IGRAPH_EXPORT igraph_error_t igraph_rng_init(igraph_rng_t *rng, const igraph_rng_type_t *type);
IGRAPH_EXPORT void igraph_rng_destroy(igraph_rng_t *rng);

IGRAPH_EXPORT igraph_error_t igraph_rng_seed(igraph_rng_t *rng, igraph_uint_t seed);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t igraph_rng_bits(const igraph_rng_t* rng);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_uint_t igraph_rng_max(const igraph_rng_t *rng);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE const char *igraph_rng_name(const igraph_rng_t *rng);

IGRAPH_EXPORT igraph_bool_t igraph_rng_get_bool(igraph_rng_t *rng) ;
IGRAPH_EXPORT igraph_int_t igraph_rng_get_integer(
    igraph_rng_t *rng, igraph_int_t l, igraph_int_t h
);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_normal(
    igraph_rng_t *rng, igraph_real_t m, igraph_real_t s
);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_unif(
    igraph_rng_t *rng, igraph_real_t l, igraph_real_t h
);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_unif01(igraph_rng_t *rng);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_geom(igraph_rng_t *rng, igraph_real_t p);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_binom(
    igraph_rng_t *rng, igraph_int_t n, igraph_real_t p
);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_exp(igraph_rng_t *rng, igraph_real_t rate);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_gamma(
    igraph_rng_t *rng, igraph_real_t shape, igraph_real_t scale
);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_pois(igraph_rng_t *rng, igraph_real_t rate);

/* --------------------------------- */

IGRAPH_EXPORT extern const igraph_rng_type_t igraph_rngtype_glibc2;
IGRAPH_EXPORT extern const igraph_rng_type_t igraph_rngtype_mt19937;
IGRAPH_EXPORT extern const igraph_rng_type_t igraph_rngtype_pcg32;
IGRAPH_EXPORT extern const igraph_rng_type_t igraph_rngtype_pcg64;

IGRAPH_EXPORT igraph_rng_t *igraph_rng_default(void);
IGRAPH_EXPORT igraph_rng_t *igraph_rng_set_default(igraph_rng_t *rng);

/* --------------------------------- */

#define RNG_BOOL()       (igraph_rng_get_bool(igraph_rng_default()))
#define RNG_INTEGER(l,h) (igraph_rng_get_integer(igraph_rng_default(),(l),(h)))
#define RNG_NORMAL(m,s)  (igraph_rng_get_normal(igraph_rng_default(),(m),(s)))
#define RNG_UNIF(l,h)    (igraph_rng_get_unif(igraph_rng_default(),(l),(h)))
#define RNG_UNIF01()     (igraph_rng_get_unif01(igraph_rng_default()))
#define RNG_GEOM(p)      (igraph_rng_get_geom(igraph_rng_default(),(p)))
#define RNG_BINOM(n,p)   (igraph_rng_get_binom(igraph_rng_default(),(n),(p)))
#define RNG_EXP(rate)    (igraph_rng_get_exp(igraph_rng_default(),(rate)))
#define RNG_POIS(rate)   (igraph_rng_get_pois(igraph_rng_default(),(rate)))
#define RNG_GAMMA(shape, scale) \
                         (igraph_rng_get_gamma(igraph_rng_default(), (shape), (scale)))

IGRAPH_END_C_DECLS

#endif
