/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef IGRAPH_RANDOM_H
#define IGRAPH_RANDOM_H

#include "igraph_decls.h"

__BEGIN_DECLS

#include <stdlib.h>
#include <time.h>

#include "igraph_types.h"
#include "igraph_vector.h"

/* The new RNG interface is (somewhat) modelled based on the GSL */

typedef struct igraph_rng_type_t {
    const char *name;
    unsigned long int min; /* at the moment, 'min' must always be set to 0 */
    unsigned long int max;
    int (*init)(void **state);
    void (*destroy)(void *state);
    int (*seed)(void *state, unsigned long int seed);
    unsigned long int (*get)(void *state);
    igraph_real_t (*get_real)(void *state);
    igraph_real_t (*get_norm)(void *state);
    igraph_real_t (*get_geom)(void *state, igraph_real_t p);
    igraph_real_t (*get_binom)(void *state, long int n, igraph_real_t p);
    igraph_real_t (*get_exp)(void *state, igraph_real_t rate);
    igraph_real_t (*get_gamma)(void *state, igraph_real_t shape,
                               igraph_real_t scale);
} igraph_rng_type_t;

typedef struct igraph_rng_t {
    const igraph_rng_type_t *type;
    void *state;
    int def;
} igraph_rng_t;

/* --------------------------------- */

IGRAPH_EXPORT int igraph_rng_init(igraph_rng_t *rng, const igraph_rng_type_t *type);
IGRAPH_EXPORT void igraph_rng_destroy(igraph_rng_t *rng);

IGRAPH_EXPORT int igraph_rng_seed(igraph_rng_t *rng, unsigned long int seed);
IGRAPH_EXPORT unsigned long int igraph_rng_max(igraph_rng_t *rng);
IGRAPH_EXPORT unsigned long int igraph_rng_min(igraph_rng_t *rng);
IGRAPH_EXPORT const char *igraph_rng_name(igraph_rng_t *rng);

IGRAPH_EXPORT long int igraph_rng_get_integer(igraph_rng_t *rng,
                                              long int l, long int h);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_normal(igraph_rng_t *rng,
                                                  igraph_real_t m, igraph_real_t s);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_unif(igraph_rng_t *rng,
                                                igraph_real_t l, igraph_real_t h);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_unif01(igraph_rng_t *rng);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_geom(igraph_rng_t *rng, igraph_real_t p);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_binom(igraph_rng_t *rng, long int n,
                                                 igraph_real_t p);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_exp(igraph_rng_t *rng, igraph_real_t rate);
IGRAPH_EXPORT unsigned long int igraph_rng_get_int31(igraph_rng_t *rng);
IGRAPH_EXPORT igraph_real_t igraph_rng_get_gamma(igraph_rng_t *rng, igraph_real_t shape,
                                                 igraph_real_t scale);
IGRAPH_EXPORT int igraph_rng_get_dirichlet(igraph_rng_t *rng,
                                           const igraph_vector_t *alpha,
                                           igraph_vector_t *result);

/* --------------------------------- */

IGRAPH_EXPORT extern const igraph_rng_type_t igraph_rngtype_glibc2;
IGRAPH_EXPORT extern const igraph_rng_type_t igraph_rngtype_rand;
IGRAPH_EXPORT extern const igraph_rng_type_t igraph_rngtype_mt19937;

IGRAPH_EXPORT igraph_rng_t *igraph_rng_default(void);
IGRAPH_EXPORT void igraph_rng_set_default(igraph_rng_t *rng);

/* --------------------------------- */

#ifdef USING_R

void GetRNGstate(void);
void PutRNGstate(void);
#define RNG_BEGIN()    GetRNGstate()
#define RNG_END()      PutRNGstate()

double Rf_dnorm4(double x, double mu, double sigma, int give_log);
#define igraph_dnorm Rf_dnorm4

#else

#define RNG_BEGIN() \
    if (igraph_rng_default()->def == 1) { \
        igraph_rng_seed(igraph_rng_default(), time(0)); \
        igraph_rng_default()->def=2; \
    }
#define RNG_END()       /* do nothing */

IGRAPH_EXPORT double igraph_dnorm(double x, double mu, double sigma, int give_log);

#endif

#define RNG_INTEGER(l,h) (igraph_rng_get_integer(igraph_rng_default(),(l),(h)))
#define RNG_NORMAL(m,s)  (igraph_rng_get_normal(igraph_rng_default(),(m),(s)))
#define RNG_UNIF(l,h)    (igraph_rng_get_unif(igraph_rng_default(),(l),(h)))
#define RNG_UNIF01()     (igraph_rng_get_unif01(igraph_rng_default()))
#define RNG_GEOM(p)      (igraph_rng_get_geom(igraph_rng_default(),(p)))
#define RNG_BINOM(n,p)   (igraph_rng_get_binom(igraph_rng_default(),(n),(p)))
#define RNG_INT31()      (igraph_rng_get_int31(igraph_rng_default()))

__END_DECLS

#endif
