/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/


#include "igraph_random.h"

#include "igraph_nongraph.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector.h"

#include "core/interruption.h"
#include "core/math.h"
#include "math/safe_intop.h"
#include "random/random_internal.h"

#include "config.h" /* IGRAPH_THREAD_LOCAL, HAVE___UINT128_T, HAVE__UMUL128 */

#if defined(HAVE__UMUL128) || defined(HAVE___UMULH)
#include <intrin.h> /* _umul128() or __umulh() are defined in intrin.h */
#endif

#include <assert.h>
#include <math.h>
#include <float.h> /* DBL_MANT_DIG */

/**
 * \section about_rngs
 *
 * <section id="about-random-numbers-in-igraph">
 * <title>About random numbers in igraph</title>
 *
 * <para>
 * Some algorithms in igraph, such as sampling from random graph models,
 * require random number generators (RNGs). igraph includes a flexible
 * RNG framework that allows hooking up arbitrary random number generators,
 * and comes with several ready-to-use generators. This framework is used
 * in igraph's high-level interfaces to integrate with the host language's
 * own RNG.
 * </para>
 * </section>
 *
 */

/**
 * \section rng_use_cases
 *
 * <section id="random-use-cases"><title>Use cases</title>
 *
 * <section id="random-normal-use"><title>Normal (default) use</title>
 * <para>
 * If the user does not use any of the RNG functions explicitly, but calls
 * some of the randomized igraph functions, then a default RNG is set
 * up the first time an igraph function needs random numbers. The
 * seed of this RNG is the output of the <code>time(0)</code> function
 * call, using the <code>time</code> function from the standard C
 * library. This ensures that igraph creates a different random graph,
 * each time the C program is called.
 * </para>
 *
 * <para>
 * The created default generator is stored internally and can be
 * queried with the \ref igraph_rng_default() function.
 * </para>
 * </section>
 *
 * <section id="random-reproducible-simulations"><title>Reproducible simulations</title>
 * <para>
 * If reproducible results are needed, then the user should set the
 * seed of the default random number generator explicitly, using the
 * \ref igraph_rng_seed() function on the default generator, \ref
 * igraph_rng_default(). When setting the seed to the same number,
 * igraph generates exactly the same random graph (or series of random
 * graphs).
 * </para>
 * </section>
 *
 * <section id="random-changing-default-generator"><title>Changing the default generator</title>
 * <para>
 * By default igraph uses the \ref igraph_rng_default() random number
 * generator. This can be changed any time by calling \ref
 * igraph_rng_set_default(), with an already initialized random number
 * generator. Note that the old (replaced) generator is not
 * destroyed, so no memory is deallocated.
 * </para>
 * </section>
 *
 * <section id="random-using-multiple-generators"><title>Using multiple generators</title>
 * <para>
 * igraph also provides functions to set up multiple random number
 * generators, using the \ref igraph_rng_init() function, and then
 * generating random numbers from them, e.g. with \ref igraph_rng_get_integer()
 * and/or \ref igraph_rng_get_unif() calls.
 * </para>
 *
 * <para>
 * Note that initializing a new random number generator is
 * independent of the generator that the igraph functions themselves
 * use. If you want to replace that, then please use \ref
 * igraph_rng_set_default().
 * </para>
 * </section>
 *
 * <section id="random-example"><title>Example</title>
 * <para>
 * \example examples/simple/random_seed.c
 * </para>
 * </section>
 *
 * </section>
 */

/* ------------------------------------ */

/**
 * \var igraph_i_rng_default
 * The default igraph random number generator
 *
 * This generator is used by all builtin igraph functions that need to
 * generate random numbers; e.g. all random graph generators.
 *
 * You can use \ref igraph_i_rng_default with \ref igraph_rng_seed()
 * to set its seed.
 *
 * You can change the default generator using the \ref
 * igraph_rng_set_default() function.
 */

extern IGRAPH_THREAD_LOCAL igraph_rng_t igraph_i_rng_default; /* defined in rng_pcg32.c */

/**
 * \function igraph_rng_set_default
 * \brief Set the default igraph random number generator.
 *
 * This function \em copies the internal structure of the given \type igraph_rng_t
 * object to igraph's internal default RNG structure. The structure itself
 * contains two pointers only, one to the "methods" of the RNG and one to the
 * memory buffer holding the internal state of the RNG. This means that if you
 * keep on generating random numbers from the RNG after setting it as the
 * default, it will affect the state of the default RNG as well because the two
 * share the same state pointer. However, do \em not expect
 * \ref igraph_rng_default() to return the same pointer as the one you passed
 * in here - the state is shared, but the entire structure is not.
 *
 * \param rng The random number generator to use as default from now
 *    on. Calling \ref igraph_rng_destroy() on it, while it is still
 *    being used as the default will result in crashes and/or
 *    unpredictable results.
 *
 * Time complexity: O(1).
 */

void igraph_rng_set_default(igraph_rng_t *rng) {
    igraph_i_rng_default = (*rng);
}


/* ------------------------------------ */

/**
 * \function igraph_rng_default
 * \brief Query the default random number generator.
 *
 * \return A pointer to the default random number generator.
 *
 * \sa \ref igraph_rng_set_default()
 */

igraph_rng_t *igraph_rng_default(void) {
    return &igraph_i_rng_default;
}

/* ------------------------------------ */

static igraph_uint_t igraph_i_rng_get_random_bits(igraph_rng_t *rng, uint8_t bits);
static uint64_t igraph_i_rng_get_random_bits_uint64(igraph_rng_t *rng, uint8_t bits);

static igraph_uint_t igraph_i_rng_get_uint(igraph_rng_t *rng);
static igraph_uint_t igraph_i_rng_get_uint_bounded(igraph_rng_t *rng, igraph_uint_t range);

static uint32_t igraph_i_rng_get_uint32(igraph_rng_t *rng);
static uint32_t igraph_i_rng_get_uint32_bounded(igraph_rng_t *rng, uint32_t range);

#if IGRAPH_INTEGER_SIZE == 64
static uint64_t igraph_i_rng_get_uint64(igraph_rng_t *rng);
static uint64_t igraph_i_rng_get_uint64_bounded(igraph_rng_t *rng, uint64_t range);
#endif

static double igraph_i_norm_rand(igraph_rng_t *rng);
static double igraph_i_exp_rand(igraph_rng_t *rng);
static double igraph_i_rbinom(igraph_rng_t *rng, igraph_integer_t n, double pp);
static double igraph_i_rexp(igraph_rng_t *rng, double rate);
static double igraph_i_rgamma(igraph_rng_t *rng, double shape, double scale);
static double igraph_i_rpois(igraph_rng_t *rng, double rate);

/**
 * \function igraph_rng_init
 * \brief Initializes a random number generator.
 *
 * This function allocates memory for a random number generator, with
 * the given type, and sets its seed to the default.
 *
 * \param rng Pointer to an uninitialized RNG.
 * \param type The type of the RNG, such as \ref igraph_rngtype_mt19937,
 * \ref igraph_rngtype_glibc2, \ref igraph_rngtype_pcg32 or
 * \ref igraph_rngtype_pcg64.
 * \return Error code.
 */

igraph_error_t igraph_rng_init(igraph_rng_t *rng, const igraph_rng_type_t *type) {
    rng->type = type;
    IGRAPH_CHECK(rng->type->init(&rng->state));
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_rng_destroy
 * \brief Deallocates memory associated with a random number generator.
 *
 * \param rng The RNG to destroy. Do not destroy an RNG that is used
 *    as the default igraph RNG.
 *
 * Time complexity: O(1).
 */

void igraph_rng_destroy(igraph_rng_t *rng) {
    rng->type->destroy(rng->state);
}

/**
 * \function igraph_rng_seed
 * \brief Seeds a random number generator.
 *
 * \param rng The RNG.
 * \param seed The new seed.
 * \return Error code.
 *
 * Time complexity: usually O(1), but may depend on the type of the
 * RNG.
 */
igraph_error_t igraph_rng_seed(igraph_rng_t *rng, igraph_uint_t seed) {
    const igraph_rng_type_t *type = rng->type;
    IGRAPH_CHECK(type->seed(rng->state, seed));
    rng->is_seeded = true;
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_rng_bits
 * \brief The number of random bits that a random number generator can produces in a single round.
 *
 * \param rng The RNG.
 * \return The number of random bits that can be generated in a single round
 *         with the RNG.
 *
 * Time complexity: O(1).
 */
igraph_integer_t igraph_rng_bits(const igraph_rng_t* rng) {
    return rng->type->bits;
}

/**
 * \function igraph_rng_max
 * \brief The maximum possible integer for a random number generator.
 *
 * Note that this number is only for informational purposes; it returns the
 * maximum possible integer that can be generated with the RNG with a single
 * call to its internals. It is derived directly from the number of random
 * \em bits that the RNG can generate in a single round. When this is smaller
 * than what would be needed by other RNG functions like \ref igraph_rng_get_integer(),
 * igraph will call the RNG multiple times to generate more random bits.
 *
 * \param rng The RNG.
 * \return The largest possible integer that can be generated in a single round
 *         with the RNG.
 *
 * Time complexity: O(1).
 */

igraph_uint_t igraph_rng_max(const igraph_rng_t *rng) {
    const igraph_rng_type_t *type = rng->type;
#if IGRAPH_INTEGER_SIZE == 64
    return (type->bits >= 64) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << type->bits) - 1);
#else
    return (type->bits >= 32) ? 0xFFFFFFFFUL : ((1ULL << type->bits) - 1);
#endif
}

/**
 * \function igraph_rng_name
 * \brief The type of a random number generator.
 *
 * \param rng The RNG.
 * \return The name of the type of the generator. Do not deallocate or
 *         change the returned string.
 *
 * Time complexity: O(1).
 */

const char *igraph_rng_name(const igraph_rng_t *rng) {
    const igraph_rng_type_t *type = rng->type;
    return type->name;
}

/**
 * Generates a given number of random bits, possibly invoking the underlying
 * RNG multiple times if needed, and returns the result in an \c igraph_uint_t .
 *
 * \param rng The RNG.
 * \param bits The number of random bits needed. Must be smaller than or equal
 *        to the size of the \c igraph_uint_t data type. Passing a value larger
 *        than the size of \c igraph_uint_t will throw away random bits except
 *        the last few that are needed to fill an \c igraph_uint_t .
 * \return The random bits, packed into the low bits of an \c igraph_uint_t .
 *         The upper, unused bits of \c igraph_uint_t will be set to zero.
 */
static igraph_uint_t igraph_i_rng_get_random_bits(igraph_rng_t *rng, uint8_t bits) {
    const igraph_rng_type_t *type = rng->type;
    igraph_integer_t rng_bitwidth = igraph_rng_bits(rng);
    igraph_uint_t result;

    if (rng_bitwidth >= bits) {
        /* keep the high bits as RNGs sometimes tend to have lower entropy in
         * low bits than in high bits */
        result = type->get(rng->state) >> (rng_bitwidth - bits);
    } else {
        result = 0;
        do {
            result = (result << rng_bitwidth) + type->get(rng->state);
            bits -= rng_bitwidth;
        } while (bits > rng_bitwidth);

        /* and now the last piece */
        result = (result << bits) + (type->get(rng->state) >> (rng_bitwidth - bits));
    }

    return result;
}

/**
 * Generates a given number of random bits, possibly invoking the underlying
 * RNG multiple times if needed, and returns the result in an \c uint64_t .
 *
 * Prefer \c igraph_i_rng_get_random_bits() if you know that you need at most
 * 32 bits due to the type of the return value. This function might perform
 * worse on 32-bit platforms because the result is always 64 bits.
 *
 * \param rng The RNG.
 * \param bits The number of random bits needed. Must be smaller than or equal
 *        to the size of the \c uint64_t data type. Passing a value larger
 *        than the size of \c uint64_t will throw away random bits except
 *        the last few that are needed to fill an \c uint64_t .
 * \return The random bits, packed into the low bits of an \c uint64_t .
 *         The upper, unused bits of \c uint64_t will be set to zero.
 */
static uint64_t igraph_i_rng_get_random_bits_uint64(igraph_rng_t *rng, uint8_t bits) {
    const igraph_rng_type_t *type = rng->type;
    igraph_integer_t rng_bitwidth = igraph_rng_bits(rng);
    uint64_t result;

    if (rng_bitwidth >= bits) {
        /* keep the high bits as RNGs sometimes tend to have lower entropy in
         * low bits than in high bits */
        result = type->get(rng->state) >> (rng_bitwidth - bits);
    } else {
        result = 0;
        do {
            result = (result << rng_bitwidth) + type->get(rng->state);
            bits -= rng_bitwidth;
        } while (bits > rng_bitwidth);

        /* and now the last piece */
        result = (result << bits) + (type->get(rng->state) >> (rng_bitwidth - bits));
    }

    return result;
}

/**
 * Generates a random integer in the full range of the \c igraph_uint_t
 * data type.
 *
 * \param rng The RNG.
 * \return The random integer.
 */
static igraph_uint_t igraph_i_rng_get_uint(igraph_rng_t *rng) {
    return igraph_i_rng_get_random_bits(rng, sizeof(igraph_uint_t) * 8);
}

/**
 * Generates a random integer in the full range of the \c uint32_t
 * data type.
 *
 * \param rng The RNG.
 * \return The random integer.
 */
static uint32_t igraph_i_rng_get_uint32(igraph_rng_t *rng) {
    return igraph_i_rng_get_random_bits(rng, 32);
}

/**
 * Generates a random integer in the range [0; range) (upper bound exclusive),
 * restricted to at most 32 bits.
 *
 * \param rng The RNG.
 * \param range The upper bound (exclusive).
 * \return The random integer.
 */
static uint32_t igraph_i_rng_get_uint32_bounded(igraph_rng_t *rng, uint32_t range) {
    /* Debiased integer multiplication -- Lemire's method
     * from https://www.pcg-random.org/posts/bounded-rands.html */
    uint32_t x, l, t = (-range) % range;
    uint64_t m;
    do {
        x = igraph_i_rng_get_uint32(rng);
        m = (uint64_t)(x) * (uint64_t)(range);
        l = (uint32_t)m;
    } while (l < t);
    return m >> 32;
}

#if IGRAPH_INTEGER_SIZE == 64
/**
 * Generates a random integer in the full range of the \c uint64_t
 * data type.
 *
 * \param rng The RNG.
 * \param range The upper bound (inclusive).
 * \return The random integer.
 */
static uint64_t igraph_i_rng_get_uint64(igraph_rng_t *rng) {
    return igraph_i_rng_get_random_bits(rng, 64);
}

#if !defined(HAVE___UINT128_T)
static uint64_t igraph_i_umul128(uint64_t a, uint64_t b, uint64_t *hi) {
#if defined(HAVE__UMUL128)
    /* MSVC has _umul128() on x64 but not on arm64 */
    return _umul128(a, b, hi);
#elif defined(HAVE___UMULH)
    /* MSVC has __umulh() on arm64 */
    *hi = __umulh(a, b);
    return a*b;
#else
    /* Portable but slow fallback implementation of unsigned
     * 64-bit multiplication obtaining a 128-bit result.
     * Based on https://stackoverflow.com/a/28904636/695132
     */

    uint64_t a_lo = (uint32_t) a;
    uint64_t a_hi = a >> 32;
    uint64_t b_lo = (uint32_t) b;
    uint64_t b_hi = b >> 32;

    uint64_t a_x_b_hi  = a_hi * b_hi;
    uint64_t a_x_b_mid = a_hi * b_lo;
    uint64_t b_x_a_mid = b_hi * a_lo;
    uint64_t a_x_b_lo  = a_lo * b_lo;

    uint64_t carry_bit = ((uint64_t) (uint32_t) a_x_b_mid +
                          (uint64_t) (uint32_t) b_x_a_mid +
                          (a_x_b_lo >> 32) ) >> 32;

    *hi = a_x_b_hi +
          (a_x_b_mid >> 32) + (b_x_a_mid >> 32) +
          carry_bit;

    return a*b;
#endif
}
#endif /* !defined(HAVE___UINT128_T) */

/**
 * Generates a random integer in the range [0; range) (upper bound exclusive),
 * restricted to at most 64 bits.
 *
 * \param rng The RNG.
 * \param range The upper bound (exclusive).
 * \return The random integer.
 */
static uint64_t igraph_i_rng_get_uint64_bounded(igraph_rng_t *rng, uint64_t range) {
    /* Debiased integer multiplication -- Lemire's method
     * from https://www.pcg-random.org/posts/bounded-rands.html */
    uint64_t x, l, t = (-range) % range;
#if defined(HAVE___UINT128_T)
    /* gcc and clang have __uint128_t */
    __uint128_t m;
    do {
        x = igraph_i_rng_get_uint64(rng);
        m = (__uint128_t)(x) * (__uint128_t)(range);
        l = (uint64_t)m;
    } while (l < t);
    return m >> 64;
#else
    uint64_t hi;
    do {
        x = igraph_i_rng_get_uint64(rng);
        l = igraph_i_umul128(x, range, &hi);
    } while (l < t);
    return hi;
#endif
}

#endif /* IGRAPH_INTEGER_SIZE == 64 */

/**
 * Generates a random integer in the range [0; range) (upper bound exclusive).
 *
 * \param rng The RNG.
 * \param range The upper bound (exclusive).
 * \return The random integer.
 */
static igraph_uint_t igraph_i_rng_get_uint_bounded(igraph_rng_t *rng, igraph_uint_t range) {
    /* We must make this function behave the same way for range < 2^32 so igraph
     * behaves the same way on 32-bit and 64-bit platforms as long as we stick
     * to integers less than 2^32. This is to ensure that the unit tests are
     * consistent */

#if IGRAPH_INTEGER_SIZE == 32
    return igraph_i_rng_get_uint32_bounded(rng, range);
#else
    if (range <= UINT32_MAX) {
        return igraph_i_rng_get_uint32_bounded(rng, range);
    } else {
        return igraph_i_rng_get_uint64_bounded(rng, range);
    }
#endif
}

/**
 * \function igraph_rng_get_integer
 * \brief Generate an integer random number from an interval.
 *
 * \param rng Pointer to the RNG to use for the generation. Use \ref
 *        igraph_rng_default() here to use the default igraph RNG.
 * \param l Lower limit, inclusive, it can be negative as well.
 * \param h Upper limit, inclusive, it can be negative as well, but it
 *        should be at least <code>l</code>.
 * \return The generated random integer.
 *
 * Time complexity: O(log2(h-l) / bits) where bits is the value of
 * \ref igraph_rng_bits(rng).
 */

igraph_integer_t igraph_rng_get_integer(
    igraph_rng_t *rng, igraph_integer_t l, igraph_integer_t h
) {
    const igraph_rng_type_t *type = rng->type;
    igraph_uint_t range;

    assert(h >= l);

    if (h == l) {
        return l;
    }

    if (type->get_int) {
        return type->get_int(rng->state, l, h);
    }

    if (IGRAPH_UNLIKELY(l == IGRAPH_INTEGER_MIN && h == IGRAPH_INTEGER_MAX)) {
        /* Full uint range is needed, we can just grab a random number from
         * the uint range and cast it to a signed integer */
        return (igraph_integer_t) igraph_i_rng_get_uint(rng);
    } else if (l >= 0 || h < 0) {
        /* this is okay, (h - l) will not overflow an igraph_integer_t */
        range = (igraph_uint_t)(h - l) + 1;
    } else {
        /* (h - l) could potentially overflow so we need to play it safe. If we
         * are here, l < 0 and h >= 0 so we can cast -l into an igraph_uint_t
         * safely and do the subtraction that way */
        range = ((igraph_uint_t)(h)) + ((igraph_uint_t)(-l)) + 1;
    }

    return l + igraph_i_rng_get_uint_bounded(rng, range);
}

/**
 * \function igraph_rng_get_normal
 * \brief Samples from a normal distribution.
 *
 * Generates random variates from a normal distribution with probability
 * density
 *
 * </para><para>
 * <code>exp( -(x - m)^2 / (2 s^2) )</code>.
 *
 * \param rng Pointer to the RNG to use. Use \ref igraph_rng_default()
 *        here to use the default igraph RNG.
 * \param m The mean.
 * \param s The standard deviation.
 * \return The generated normally distributed random number.
 *
 * Time complexity: depends on the type of the RNG.
 */

igraph_real_t igraph_rng_get_normal(igraph_rng_t *rng,
                                    igraph_real_t m, igraph_real_t s) {
    const igraph_rng_type_t *type = rng->type;
    if (type->get_norm) {
        return type->get_norm(rng->state) * s + m;
    } else {
        return igraph_i_norm_rand(rng) * s + m;
    }
}

/**
 * \function igraph_rng_get_unif
 * \brief Samples real numbers from a given interval.
 *
 * Generates uniformly distributed real numbers from the <code>[l, h)</code>
 * half-open interval.
 *
 * \param rng Pointer to the RNG to use. Use \ref igraph_rng_default()
 *        here to use the default igraph RNG.
 * \param l The lower bound, it can be negative.
 * \param h The upper bound, it can be negative, but it has to be
 *        larger than the lower bound.
 * \return The generated uniformly distributed random number.
 *
 * Time complexity: depends on the type of the RNG.
 */

igraph_real_t igraph_rng_get_unif(igraph_rng_t *rng,
                                  igraph_real_t l, igraph_real_t h) {
    assert(h >= l);

    if (l == h) return h;

    /* Ensure that 'h' is never produced due to numerical roundoff errors, except when l == h. */
    igraph_real_t r;
    do {
        r = igraph_rng_get_unif01(rng) * (h - l) + l;
    } while (IGRAPH_UNLIKELY(r == h));
    return r;
}

/**
 * \function igraph_rng_get_unif01
 * \brief Samples uniformly from the unit interval.
 *
 * Generates uniformly distributed real numbers from the <code>[0, 1)</code>
 * half-open interval.
 *
 * \param rng Pointer to the RNG to use. Use \ref igraph_rng_default()
 *        here to use the default igraph RNG.
 * \return The generated uniformly distributed random number.
 *
 * Time complexity: depends on the type of the RNG.
 */

igraph_real_t igraph_rng_get_unif01(igraph_rng_t *rng) {
    const igraph_rng_type_t *type = rng->type;
    if (type->get_real) {
        return type->get_real(rng->state);
    } else {
        /* We extract 52 random bits from a 64-bit uint and fill that directly
         * into the mantissa of a double, bit-by-bit, clear the sign bit and
         * set the exponent to 2^0. This way we get a 52-bit random double
         * between 1 (inclusive) and 2 (exclusive), uniformly distributed.
         * Then we subtract 1 to arrive at the [0; 1) interval. This is fast
         * but we lose one bit of precision as there are 2^53 possible doubles
         * between 0 and 1. */
        union {
            uint64_t as_uint64_t;
            double as_double;
        } value;
        value.as_uint64_t =
            (igraph_i_rng_get_random_bits_uint64(rng, 52) & 0xFFFFFFFFFFFFFull) | 0x3FF0000000000000ull;
        return value.as_double - 1.0;
    }
}

/**
 * \function igraph_rng_get_geom
 * \brief Samples from a geometric distribution.
 *
 * Generates random variates from a geometric distribution. The number \c k is
 * generated with probability
 *
 * </para><para>
 * <code>(1 - p)^k p</code>, <code>k = 0, 1, 2, ...</code>.
 *
 * \param rng Pointer to the RNG to use. Use \ref igraph_rng_default()
 *        here to use the default igraph RNG.
 * \param p The probability of success in each trial. Must be larger
 *        than zero and smaller or equal to 1.
 * \return The generated geometrically distributed random number.
 *
 * Time complexity: depends on the RNG.
 */

igraph_real_t igraph_rng_get_geom(igraph_rng_t *rng, igraph_real_t p) {
    const igraph_rng_type_t *type = rng->type;
    if (!isfinite(p) || p <= 0 || p > 1) {
        return IGRAPH_NAN;
    }
    if (type->get_geom) {
        return type->get_geom(rng->state, p);
    } else {
        return igraph_rng_get_pois(rng, igraph_i_exp_rand(rng) * ((1 - p) / p));
    }
}

/**
 * \function igraph_rng_get_binom
 * \brief Samples from a binomial distribution.
 *
 * Generates random variates from a binomial distribution. The number \c k is generated
 * with probability
 *
 * </para><para>
 * <code>(n \choose k) p^k (1-p)^(n-k)</code>, <code>k = 0, 1, ..., n</code>.
 *
 * \param rng Pointer to the RNG to use. Use \ref igraph_rng_default()
 *        here to use the default igraph RNG.
 * \param n Number of observations.
 * \param p Probability of an event.
 * \return The generated binomially distributed random number.
 *
 * Time complexity: depends on the RNG.
 */

igraph_real_t igraph_rng_get_binom(igraph_rng_t *rng, igraph_integer_t n, igraph_real_t p) {
    const igraph_rng_type_t *type = rng->type;
    if (type->get_binom) {
        return type->get_binom(rng->state, n, p);
    } else {
        return igraph_i_rbinom(rng, n, p);
    }
}

/**
 * \function igraph_rng_get_gamma
 * \brief Samples from a gamma distribution.
 *
 * Generates random variates from a gamma distribution with probability
 * density proportional to
 *
 * </para><para>
 * <code>x^(shape-1) exp(-x / scale)</code>.
 *
 * \param rng Pointer to the RNG to use. Use \ref igraph_rng_default()
 *        here to use the default igraph RNG.
 * \param shape Shape parameter.
 * \param scale Scale parameter.
 * \return The generated sample.
 *
 * Time complexity: depends on the RNG.
 */

igraph_real_t igraph_rng_get_gamma(igraph_rng_t *rng, igraph_real_t shape,
                                   igraph_real_t scale) {
    const igraph_rng_type_t *type = rng->type;
    if (type->get_gamma) {
        return type->get_gamma(rng->state, shape, scale);
    } else {
        return igraph_i_rgamma(rng, shape, scale);
    }
}

/**
 * \function igraph_rng_get_exp
 * \brief Samples from an exponential distribution.
 *
 * Generates random variates from an exponential distribution with probability
 * density proportional to
 *
 * </para><para>
 * <code>exp(-rate x)</code>.
 *
 * \param rng Pointer to the RNG to use. Use \ref igraph_rng_default()
 *        here to use the default igraph RNG.
 * \param rate Rate parameter.
 * \return The generated sample.
 *
 * Time complexity: depends on the RNG.
 */

igraph_real_t igraph_rng_get_exp(igraph_rng_t *rng, igraph_real_t rate) {
    const igraph_rng_type_t *type = rng->type;
    if (type->get_exp) {
        return type->get_exp(rng->state, rate);
    } else {
        return igraph_i_rexp(rng, rate);
    }
}

/**
 * \function igraph_rng_get_pois
 * \brief Samples from a Poisson distribution.
 *
 * Generates random variates from a Poisson distribution. The number \c k is generated
 * with probability
 *
 * </para><para>
 * <code>rate^k * exp(-rate) / k!</code>, <code>k = 0, 1, 2, ...</code>.
 *
 * \param rng Pointer to the RNG to use. Use \ref igraph_rng_default()
 *        here to use the default igraph RNG.
 * \param rate The rate parameter of the Poisson distribution. Must not be negative.
 * \return The generated geometrically distributed random number.
 *
 * Time complexity: depends on the RNG.
 */

igraph_real_t igraph_rng_get_pois(igraph_rng_t *rng, igraph_real_t rate) {
    const igraph_rng_type_t *type = rng->type;
    if (isnan(rate) || rate < 0) {
        return IGRAPH_NAN;
    } else if (rate == 0) {
        return 0;
    } else if (type->get_pois) {
        return type->get_pois(rng->state, rate);
    } else {
        return igraph_i_rpois(rng, rate);
    }
}


/**
 * \ingroup internal
 *
 * This function appends the rest of the needed random numbers to the
 * result vector. It is Algoirthm A in Vitter's paper.
 */

static void igraph_i_random_sample_alga(igraph_vector_int_t *res,
                                        igraph_integer_t l, igraph_integer_t h,
                                        igraph_integer_t length) {
    /* Vitter: Variables V, quot, Nreal, and top are of type real */

    igraph_integer_t N = h - l + 1;
    igraph_integer_t n = length;

    igraph_real_t top = N - n;
    igraph_real_t Nreal = N;
    igraph_integer_t S = 0;
    igraph_real_t V, quot;

    l = l - 1;

    while (n >= 2) {
        V = RNG_UNIF01();
        S = 1;
        quot = top / Nreal;
        while (quot > V) {
            S += 1;
            top = -1.0 + top;
            Nreal = -1.0 + Nreal;
            quot = (quot * top) / Nreal;
        }
        l += S;
        igraph_vector_int_push_back(res, l); /* allocated */
        Nreal = -1.0 + Nreal; n = -1 + n;
    }

    S = trunc(round(Nreal) * RNG_UNIF01());
    l += S + 1;
    igraph_vector_int_push_back(res, l); /* allocated */
}

/**
 * \ingroup nongraph
 * \function igraph_random_sample
 * \brief Generates an increasing random sequence of integers.
 *
 * This function generates an increasing sequence of random integer
 * numbers from a given interval. The algorithm is taken literally
 * from (Vitter 1987). This method can be used for generating numbers from a
 * \em very large interval. It is primarily created for randomly
 * selecting some edges from the sometimes huge set of possible edges
 * in a large graph.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * J. S. Vitter. An efficient algorithm for sequential random sampling.
 * ACM Transactions on Mathematical Software, 13(1):58--67, 1987.
 * https://doi.org/10.1145/23002.23003
 *
 * \param res Pointer to an initialized vector. This will hold the
 *        result. It will be resized to the proper size.
 * \param l The lower limit of the generation interval (inclusive). This must
 *        be less than or equal to the upper limit, and it must be integral.
 * \param h The upper limit of the generation interval (inclusive). This must
 *        be greater than or equal to the lower limit, and it must be integral.
 * \param length The number of random integers to generate.
 * \return The error code \c IGRAPH_EINVAL is returned in each of the
 *         following cases: (1) The given lower limit is greater than the
 *         given upper limit, i.e. \c l &gt; \c h. (2) Assuming that
 *         \c l &lt; \c h and N is the sample size, the above error code is
 *         returned if N &gt; |\c h - \c l|, i.e. the sample size exceeds the
 *         size of the candidate pool.
 *
 * Time complexity: according to (Vitter 1987), the expected
 * running time is O(length).
 *
 * \example examples/simple/igraph_random_sample.c
 */

igraph_error_t igraph_random_sample(igraph_vector_int_t *res, igraph_integer_t l, igraph_integer_t h,
                         igraph_integer_t length) {
    igraph_integer_t N; /* := h - l + 1 */
    IGRAPH_SAFE_ADD(h, -l, &N);
    IGRAPH_SAFE_ADD(N, 1, &N);

    igraph_integer_t n = length;

    igraph_real_t nreal = length;
    igraph_real_t ninv = (nreal != 0) ? 1.0 / nreal : 0.0;
    igraph_real_t Nreal = N;
    igraph_real_t Vprime;
    igraph_integer_t qu1 = -n + 1 + N;
    igraph_real_t qu1real = -nreal + 1.0 + Nreal;
    igraph_real_t negalphainv = -13;
    igraph_real_t threshold = -negalphainv * n;
    igraph_integer_t S;

    /* getting back some sense of sanity */
    if (l > h) {
        IGRAPH_ERROR("Lower limit is greater than upper limit.", IGRAPH_EINVAL);
    }
    /* now we know that l <= h */
    if (length > N) {
        IGRAPH_ERROR("Sample size exceeds size of candidate pool.", IGRAPH_EINVAL);
    }

    /* treat rare cases quickly */
    if (l == h) {
        IGRAPH_CHECK(igraph_vector_int_resize(res, 1));
        VECTOR(*res)[0] = l;
        return IGRAPH_SUCCESS;
    }
    if (length == 0) {
        igraph_vector_int_clear(res);
        return IGRAPH_SUCCESS;
    }
    if (length == N) {
        IGRAPH_CHECK(igraph_vector_int_resize(res, length));
        for (igraph_integer_t i = 0; i < length; i++) {
            VECTOR(*res)[i] = l++;
        }
        return IGRAPH_SUCCESS;
    }

    igraph_vector_int_clear(res);
    IGRAPH_CHECK(igraph_vector_int_reserve(res, length));

    RNG_BEGIN();

    Vprime = exp(log(RNG_UNIF01()) * ninv);
    l = l - 1;

    while (n > 1 && threshold < N) {
        igraph_real_t X, U;
        igraph_real_t limit, t;
        igraph_real_t negSreal, y1, y2, top, bottom;
        igraph_real_t nmin1inv = 1.0 / (-1.0 + nreal);
        while (1) {
            while (1) {
                X = Nreal * (-Vprime + 1.0);
                S = floor(X);
                /* if (S==0) { S=1; } */
                if (S < qu1) {
                    break;
                }
                Vprime = exp(log(RNG_UNIF01()) * ninv);
            }
            U = RNG_UNIF01();
            negSreal = -S;

            y1 = exp(log(U * Nreal / qu1real) * nmin1inv);
            Vprime = y1 * (-X / Nreal + 1.0) * (qu1real / (negSreal + qu1real));
            if (Vprime <= 1.0) {
                break;
            }

            y2 = 1.0;
            top = -1.0 + Nreal;
            if (-1 + n > S) {
                bottom = -nreal + Nreal;
                limit = -S + N;
            } else {
                bottom = -1.0 + negSreal + Nreal;
                limit = qu1;
            }
            for (t = -1 + N; t >= limit; t--) {
                y2 = (y2 * top) / bottom;
                top = -1.0 + top;
                bottom = -1.0 + bottom;
            }
            if (Nreal / (-X + Nreal) >= y1 * exp(log(y2)*nmin1inv)) {
                Vprime = exp(log(RNG_UNIF01()) * nmin1inv);
                break;
            }
            Vprime = exp(log(RNG_UNIF01()) * ninv);
        }

        l += S + 1;
        igraph_vector_int_push_back(res, l);    /* allocated */
        N = -S + (-1 + N);   Nreal = negSreal + (-1.0 + Nreal);
        n = -1 + n;   nreal = -1.0 + nreal; ninv = nmin1inv;
        qu1 = -S + qu1; qu1real = negSreal + qu1real;
        threshold = threshold + negalphainv;
    }

    if (n > 1) {
        igraph_i_random_sample_alga(res, l + 1, h, n);
    } else {
        S = floor(N * Vprime);
        l += S + 1;
        igraph_vector_int_push_back(res, l);    /* allocated */
    }

    RNG_END();

    return IGRAPH_SUCCESS;
}

static void igraph_i_random_sample_alga_real(igraph_vector_t *res,
                                       igraph_real_t l, igraph_real_t h,
                                       igraph_real_t length) {
    igraph_real_t N = h - l + 1;
    igraph_real_t n = length;

    igraph_real_t top = N - n;
    igraph_real_t Nreal = N;
    igraph_real_t S = 0;
    igraph_real_t V, quot;

    l = l - 1;

    while (n >= 2) {
        V = RNG_UNIF01();
        S = 1;
        quot = top / Nreal;
        while (quot > V) {
            S += 1;
            top = -1.0 + top;
            Nreal = -1.0 + Nreal;
            quot = (quot * top) / Nreal;
        }
        l += S;
        igraph_vector_push_back(res, l); /* allocated */
        Nreal = -1.0 + Nreal; n = -1 + n;
    }

    S = trunc(round(Nreal) * RNG_UNIF01());
    l += S + 1;
    igraph_vector_push_back(res, l); /* allocated */
}

/**
 * \ingroup nongraph
 * \function igraph_random_sample_real
 * \brief Generates an increasing random sequence of integers (igraph_real_t version).
 *
 * This function is the 'real' version of \ref igraph_random_sample(), and was added
 * so \ref igraph_erdos_renyi_game_gnm() and related functions can use a random sample
 * of doubles instead of integers to prevent overflows on systems with 32-bit
 * \type igraph_integer_t.
 *
 * \param res Pointer to an initialized vector. This will hold the
 *        result. It will be resized to the proper size.
 * \param l The lower limit of the generation interval (inclusive). This must
 *        be less than or equal to the upper limit, and it must be integral.
 *        Passing a fractional number here results in undefined behaviour.
 * \param h The upper limit of the generation interval (inclusive). This must
 *        be greater than or equal to the lower limit, and it must be integral.
 *        Passing a fractional number here results in undefined behaviour.
 * \param length The number of random integers to generate.
 * \return The error code \c IGRAPH_EINVAL is returned in each of the
 *         following cases: (1) The given lower limit is greater than the
 *         given upper limit, i.e. \c l &gt; \c h. (2) Assuming that
 *         \c l &lt; \c h and N is the sample size, the above error code is
 *         returned if N &gt; |\c h - \c l|, i.e. the sample size exceeds the
 *         size of the candidate pool.
 */

igraph_error_t igraph_random_sample_real(igraph_vector_t *res, igraph_real_t l,
                    igraph_real_t h, igraph_integer_t length) {
    /* This function is the 'real' version of igraph_random_sample, and was added
     * so erdos_renyi_game_gnm can use a random sample of doubles instead of integers
     * to prevent overflows on systems with 32-bits igraph_integer_t.
     */
    igraph_real_t N = h - l + 1;
    igraph_real_t n = length;

    igraph_real_t nreal = length;
    igraph_real_t ninv = (nreal != 0) ? 1.0 / nreal : 0.0;
    igraph_real_t Nreal = N;
    igraph_real_t Vprime;
    igraph_real_t qu1 = -n + 1 + N;
    igraph_real_t qu1real = -nreal + 1.0 + Nreal;
    igraph_real_t negalphainv = -13;
    igraph_real_t threshold = -negalphainv * n;
    igraph_real_t S;
    int iter = 0;

    /* getting back some sense of sanity */
    if (l > h) {
        IGRAPH_ERROR("Lower limit is greater than upper limit.", IGRAPH_EINVAL);
    }
    /* now we know that l <= h */
    if (length > N) {
        IGRAPH_ERROR("Sample size exceeds size of candidate pool.", IGRAPH_EINVAL);
    }

    /* ensure that we work in the range where igraph_real_t can represent integers exactly */
    if (h > IGRAPH_MAX_EXACT_REAL || l < -IGRAPH_MAX_EXACT_REAL || N > IGRAPH_MAX_EXACT_REAL) {
        IGRAPH_ERROR("Sampling interval too large.", IGRAPH_EOVERFLOW);
    }

    /* treat rare cases quickly */
    if (l == h) {
        IGRAPH_CHECK(igraph_vector_resize(res, 1));
        VECTOR(*res)[0] = l;
        return IGRAPH_SUCCESS;
    }
    if (length == 0) {
        igraph_vector_clear(res);
        return IGRAPH_SUCCESS;
    }
    if (length == N) {
        IGRAPH_CHECK(igraph_vector_resize(res, length));
        for (igraph_integer_t i = 0; i < length; i++) {
            VECTOR(*res)[i] = l++;
        }
        return IGRAPH_SUCCESS;
    }

    igraph_vector_clear(res);
    IGRAPH_CHECK(igraph_vector_reserve(res, length));

    RNG_BEGIN();

    Vprime = exp(log(RNG_UNIF01()) * ninv);
    l = l - 1;

    while (n > 1 && threshold < N) {
        igraph_real_t X, U;
        igraph_real_t limit, t;
        igraph_real_t negSreal, y1, y2, top, bottom;
        igraph_real_t nmin1inv = 1.0 / (-1.0 + nreal);
        while (1) {
            while (1) {
                X = Nreal * (-Vprime + 1.0);
                S = floor(X);
                /* if (S==0) { S=1; } */
                if (S < qu1) {
                    break;
                }
                Vprime = exp(log(RNG_UNIF01()) * ninv);
            }
            U = RNG_UNIF01();
            negSreal = -S;

            y1 = exp(log(U * Nreal / qu1real) * nmin1inv);
            Vprime = y1 * (-X / Nreal + 1.0) * (qu1real / (negSreal + qu1real));
            if (Vprime <= 1.0) {
                break;
            }

            y2 = 1.0;
            top = -1.0 + Nreal;
            if (-1 + n > S) {
                bottom = -nreal + Nreal;
                limit = -S + N;
            } else {
                bottom = -1.0 + negSreal + Nreal;
                limit = qu1;
            }
            for (t = -1 + N; t >= limit; t--) {
                y2 = (y2 * top) / bottom;
                top = -1.0 + top;
                bottom = -1.0 + bottom;
            }
            if (Nreal / (-X + Nreal) >= y1 * exp(log(y2)*nmin1inv)) {
                Vprime = exp(log(RNG_UNIF01()) * nmin1inv);
                break;
            }
            Vprime = exp(log(RNG_UNIF01()) * ninv);
        }

        l += S + 1;
        igraph_vector_push_back(res, l);    /* allocated */
        N = -S + (-1 + N);   Nreal = negSreal + (-1.0 + Nreal);
        n = -1 + n;   nreal = -1.0 + nreal; ninv = nmin1inv;
        qu1 = -S + qu1; qu1real = negSreal + qu1real;
        threshold = threshold + negalphainv;

        if (++iter >= (1 << 14)) {
            iter = 0;
            IGRAPH_ALLOW_INTERRUPTION();
        }
    }

    if (n > 1) {
        igraph_i_random_sample_alga_real(res, l + 1, h, n);
    } else {
        S = floor(N * Vprime);
        l += S + 1;
        igraph_vector_push_back(res, l);    /* allocated */
    }

    RNG_END();

    return IGRAPH_SUCCESS;
}

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000 The R Development Core Team
 *  based on AS 111 (C) 1977 Royal Statistical Society
 *  and   on AS 241 (C) 1988 Royal Statistical Society
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *  SYNOPSIS
 *
 *  double qnorm5(double p, double mu, double sigma,
 *            int lower_tail, int log_p)
 *            {qnorm (..) is synonymous and preferred inside R}
 *
 *  DESCRIPTION
 *
 *  Compute the quantile function for the normal distribution.
 *
 *  For small to moderate probabilities, algorithm referenced
 *  below is used to obtain an initial approximation which is
 *  polished with a final Newton step.
 *
 *  For very large arguments, an algorithm of Wichura is used.
 *
 *  REFERENCE
 *
 *  Beasley, J. D. and S. G. Springer (1977).
 *  Algorithm AS 111: The percentage points of the normal distribution,
 *  Applied Statistics, 26, 118-121.
 *
 *      Wichura, M.J. (1988).
 *      Algorithm AS 241: The Percentage Points of the Normal Distribution.
 *      Applied Statistics, 37, 477-484.
 */

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2004  The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* The ISNAN macro is used in some of the code borrowed from R below. */
#define ISNAN isnan

/* Indicates that we use systems which support NaN values. */
#define IEEE_754 1

/* Private header file for use during compilation of Mathlib */
#ifndef MATHLIB_PRIVATE_H
#define MATHLIB_PRIVATE_H

#define ML_POSINF IGRAPH_INFINITY
#define ML_NEGINF -IGRAPH_INFINITY
#define ML_NAN    IGRAPH_NAN

#define ML_ERROR(x) /* nothing */
#define ML_UNDERFLOW    (DBL_MIN * DBL_MIN)
#define ML_VALID(x) (!ISNAN(x))

#define ME_NONE     0
/*  no error */
#define ME_DOMAIN   1
/*  argument out of domain */
#define ME_RANGE    2
/*  value out of range */
#define ME_NOCONV   4
/*  process did not converge */
#define ME_PRECISION    8
/*  does not have "full" precision */
#define ME_UNDERFLOW    16
/*  and underflow occurred (important for IEEE)*/

#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN); return ML_NAN; }

#endif /* MATHLIB_PRIVATE_H */


/* Utilities for `dpq' handling (density/probability/quantile) */

/* give_log in "d";  log_p in "p" & "q" : */
#define give_log log_p
/* "DEFAULT" */
/* --------- */
#define R_D__0  (log_p ? ML_NEGINF : 0.)        /* 0 */
#define R_D__1  (log_p ? 0. : 1.)           /* 1 */
#define R_DT_0  (lower_tail ? R_D__0 : R_D__1)      /* 0 */
#define R_DT_1  (lower_tail ? R_D__1 : R_D__0)      /* 1 */

#define R_D_Lval(p) (lower_tail ? (p) : (1 - (p)))  /*  p  */
#define R_D_Cval(p) (lower_tail ? (1 - (p)) : (p))  /*  1 - p */

#define R_D_val(x)  (log_p  ? log(x) : (x))     /*  x  in pF(x,..) */
#define R_D_qIv(p)  (log_p  ? exp(p) : (p))     /*  p  in qF(p,..) */
#define R_D_exp(x)  (log_p  ?  (x)   : exp(x))  /* exp(x) */
#define R_D_log(p)  (log_p  ?  (p)   : log(p))  /* log(p) */
#define R_D_Clog(p) (log_p  ? log1p(-(p)) : (1 - (p)))/* [log](1-p) */

/* log(1-exp(x)):  R_D_LExp(x) == (log1p(- R_D_qIv(x))) but even more stable:*/
#define R_D_LExp(x)     (log_p ? R_Log1_Exp(x) : log1p(-x))

/*till 1.8.x:
 * #define R_DT_val(x)  R_D_val(R_D_Lval(x))
 * #define R_DT_Cval(x) R_D_val(R_D_Cval(x)) */
#define R_DT_val(x) (lower_tail ? R_D_val(x)  : R_D_Clog(x))
#define R_DT_Cval(x)    (lower_tail ? R_D_Clog(x) : R_D_val(x))

/*#define R_DT_qIv(p)   R_D_Lval(R_D_qIv(p))         *  p  in qF ! */
#define R_DT_qIv(p) (log_p ? (lower_tail ? exp(p) : - expm1(p)) \
                     : R_D_Lval(p))

/*#define R_DT_CIv(p)   R_D_Cval(R_D_qIv(p))         *  1 - p in qF */
#define R_DT_CIv(p) (log_p ? (lower_tail ? -expm1(p) : exp(p)) \
                     : R_D_Cval(p))

#define R_DT_exp(x) R_D_exp(R_D_Lval(x))        /* exp(x) */
#define R_DT_Cexp(x)    R_D_exp(R_D_Cval(x))        /* exp(1 - x) */

#define R_DT_log(p) (lower_tail? R_D_log(p) : R_D_LExp(p))/* log(p) in qF */
#define R_DT_Clog(p)    (lower_tail? R_D_LExp(p): R_D_log(p))/* log(1-p) in qF*/
#define R_DT_Log(p) (lower_tail? (p) : R_Log1_Exp(p))
/* ==   R_DT_log when we already "know" log_p == TRUE :*/

#define R_Q_P01_check(p)            \
    if ((log_p  && p > 0) ||            \
        (!log_p && (p < 0 || p > 1)) )      \
        ML_ERR_return_NAN

/* additions for density functions (C.Loader) */
#define R_D_fexp(f,x)     (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))
#define R_D_forceint(x)   floor((x) + 0.5)
#define R_D_nonint(x)     (fabs((x) - floor((x)+0.5)) > 1e-7)
/* [neg]ative or [non int]eger : */
#define R_D_negInonint(x) (x < 0. || R_D_nonint(x))

#define R_D_nonint_check(x)                 \
    if (R_D_nonint(x)) {                  \
        MATHLIB_WARNING("non-integer x = %f", x);   \
        return R_D__0;                  \
    }

static double igraph_i_qnorm5(double p, double mu, double sigma, igraph_bool_t lower_tail, igraph_bool_t log_p) {
    double p_, q, r, val;

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma)) {
        return p + mu + sigma;
    }
#endif
    if (p == R_DT_0) {
        return ML_NEGINF;
    }
    if (p == R_DT_1) {
        return ML_POSINF;
    }
    R_Q_P01_check(p);

    if (sigma  < 0) {
        ML_ERR_return_NAN;
    }
    if (sigma == 0) {
        return mu;
    }

    p_ = R_DT_qIv(p);/* real lower_tail prob. p */
    q = p_ - 0.5;

    /*-- use AS 241 --- */
    /* double ppnd16_(double *p, long *ifault)*/
    /*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

            Produces the normal deviate Z corresponding to a given lower
            tail area of P; Z is accurate to about 1 part in 10**16.

            (original fortran code used PARAMETER(..) for the coefficients
             and provided hash codes for checking them...)
    */
    if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
        r = .180625 - q * q;
        val =
            q * (((((((r * 2509.0809287301226727 +
                       33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                   1971.5909503065514427) * r + 133.14166789178437745) * r +
                 3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                   21213.794301586595867) * r + 5394.1960214247511077) * r +
                 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    } else { /* closer than 0.075 from {0,1} boundary */

        /* r = min(p, 1-p) < 0.075 */
        if (q > 0) {
            r = R_DT_CIv(p);    /* 1-p */
        } else {
            r = p_;    /* = R_DT_Iv(p) ^=  p */
        }

        r = sqrt(- ((log_p &&
                     ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ?
                    p : /* else */ log(r)));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

        if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                         .0227238449892691845833) * r + .24178072517745061177) *
                       r + 1.27045825245236838258) * r +
                      3.64784832476320460504) * r + 5.7694972214606914055) *
                    r + 4.6303378461565452959) * r +
                   1.42343711074968357734)
                  / (((((((r *
                           1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                          r + .0151986665636164571966) * r +
                         .14810397642748007459) * r + .68976733498510000455) *
                       r + 1.6763848301838038494) * r +
                      2.05319162663775882187) * r + 1.);
        } else { /* very close to  0 or 1 */
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                         2.71155556874348757815e-5) * r +
                        .0012426609473880784386) * r + .026532189526576123093) *
                      r + .29656057182850489123) * r +
                     1.7848265399172913358) * r + 5.4637849111641143699) *
                   r + 6.6579046435011037772)
                  / (((((((r *
                           2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                          r + 1.8463183175100546818e-5) * r +
                         7.868691311456132591e-4) * r + .0148753612908506148525)
                       * r + .13692988092273580531) * r +
                      .59983220655588793769) * r + 1.);
        }

        if (q < 0.0) {
            val = -val;
        }
        /* return (q >= 0.)? r : -r ;*/
    }
    return mu + sigma * val;
}

static igraph_integer_t imax2(igraph_integer_t x, igraph_integer_t y) {
    return (x < y) ? y : x;
}

static igraph_integer_t imin2(igraph_integer_t x, igraph_integer_t y) {
    return (x < y) ? x : y;
}

static double igraph_i_norm_rand(igraph_rng_t *rng) {
    double r;

    /* Use the inversion method based on uniform variates from (0, 1).
     * We exclude 0.0 as it would lead to generating -infinity.
     * It is assumed that unif01() provides sufficient accuracy.
     * A resolution of 2^-32 may not be sufficient. igraph's default
     * implementaton provides an accuracy of 2^-52.
     */
    do {
        r = igraph_rng_get_unif01(rng);
    } while (r == 0.0);

    return igraph_i_qnorm5(r, 0.0, 1.0, true, false);
}

/*
 * The following function is igraph code (not R / Mathlib).
 *
 * We use simple inverse transform sampling, with the assumption that the
 * quality/resolution of uniform variates is high (52 bits in the default
 * implementation). The quantile function is -log(1 - r) but given that
 * r is sampled uniformly form the unit interval, -log(r) is equivalent.
 * r = 0 is disallowed as it would yield infinity.
 */

static double igraph_i_exp_rand(igraph_rng_t *rng) {
    igraph_real_t r = igraph_rng_get_unif01(rng);
    if (r == 0.0) r = 1.0; /* sample from (0, 1] instead of [0, 1) */
    return -log(r);
}

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2001 The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double rpois(double lambda)
 *
 *  DESCRIPTION
 *
 *    Random variates from the Poisson distribution.
 *
 *  REFERENCE
 *
 *    Ahrens, J.H. and Dieter, U. (1982).
 *    Computer generation of Poisson deviates
 *    from modified normal distributions.
 *    ACM Trans. Math. Software 8, 163-179.
 */

#define a0  -0.5
#define a1   0.3333333
#define a2  -0.2500068
#define a3   0.2000118
#define a4  -0.1661269
#define a5   0.1421878
#define a6  -0.1384794
#define a7   0.1250060

#define one_7   0.1428571428571428571
#define one_12  0.0833333333333333333
#define one_24  0.0416666666666666667

#define repeat for(;;)

#define FALSE 0
#define TRUE  1
#define M_1_SQRT_2PI    0.398942280401432677939946059934     /* 1/sqrt(2pi) */

static double igraph_i_rpois(igraph_rng_t *rng, double mu) {
    /* Factorial Table (0:9)! */
    const double fact[10] = {
        1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.
    };

    /* These are static --- persistent between calls for same mu : */
    static IGRAPH_THREAD_LOCAL int l;
    static IGRAPH_THREAD_LOCAL igraph_integer_t m;

    static IGRAPH_THREAD_LOCAL double b1, b2, c, c0, c1, c2, c3;
    static IGRAPH_THREAD_LOCAL double pp[36], p0, p, q, s, d, omega;
    static IGRAPH_THREAD_LOCAL double big_l;/* integer "w/o overflow" */
    static IGRAPH_THREAD_LOCAL double muprev = 0., muprev2 = 0.;/*, muold    = 0.*/

    /* Local Vars  [initialize some for -Wall]: */
    double del, difmuk = 0., E = 0., fk = 0., fx, fy, g, px, py, t, u = 0., v, x;
    double pois = -1.;
    int k, kflag, big_mu, new_big_mu = FALSE;

    if (!isfinite(mu) || mu < 0) {
        ML_ERR_return_NAN;
    }

    if (mu <= 0.) {
        return 0.;
    }

    big_mu = mu >= 10.;
    if (big_mu) {
        new_big_mu = FALSE;
    }

    if (!(big_mu && mu == muprev)) {/* maybe compute new persistent par.s */

        if (big_mu) {
            new_big_mu = TRUE;
            /* Case A. (recalculation of s,d,l  because mu has changed):
             * The Poisson probabilities pk exceed the discrete normal
             * probabilities fk whenever k >= m(mu).
             */
            muprev = mu;
            s = sqrt(mu);
            d = 6. * mu * mu;
            big_l = floor(mu - 1.1484);
            /* = an upper bound to m(mu) for all mu >= 10.*/
        } else { /* Small mu ( < 10) -- not using normal approx. */

            /* Case B. (start new table and calculate p0 if necessary) */

            /*muprev = 0.;-* such that next time, mu != muprev ..*/
            if (mu != muprev) {
                muprev = mu;
                m = imax2(1, (igraph_integer_t) mu);
                l = 0; /* pp[] is already ok up to pp[l] */
                q = p0 = p = exp(-mu);
            }

            repeat {
                /* Step U. uniform sample for inversion method */
                u = igraph_rng_get_unif01(rng);
                if (u <= p0) {
                    return 0.;
                }

                /* Step T. table comparison until the end pp[l] of the
                   pp-table of cumulative Poisson probabilities
                   (0.458 > ~= pp[9](= 0.45792971447) for mu=10 ) */
                if (l != 0) {
                    for (k = (u <= 0.458) ? 1 : imin2(l, m);  k <= l; k++)
                        if (u <= pp[k]) {
                            return (double)k;
                        }
                    if (l == 35) { /* u > pp[35] */
                        continue;
                    }
                }
                /* Step C. creation of new Poisson
                   probabilities p[l..] and their cumulatives q =: pp[k] */
                l++;
                for (k = l; k <= 35; k++) {
                    p *= mu / k;
                    q += p;
                    pp[k] = q;
                    if (u <= q) {
                        l = k;
                        return (double)k;
                    }
                }
                l = 35;
            } /* end(repeat) */
        }/* mu < 10 */

    } /* end {initialize persistent vars} */

    /* Only if mu >= 10 : ----------------------- */

    /* Step N. normal sample */
    g = mu + s * igraph_i_norm_rand(rng);/* norm_rand() ~ N(0,1), standard normal */

    if (g >= 0.) {
        pois = floor(g);
        /* Step I. immediate acceptance if pois is large enough */
        if (pois >= big_l) {
            return pois;
        }
        /* Step S. squeeze acceptance */
        fk = pois;
        difmuk = mu - fk;
        u = igraph_rng_get_unif01(rng); /* ~ U(0,1) - sample */
        if (d * u >= difmuk * difmuk * difmuk) {
            return pois;
        }
    }

    /* Step P. preparations for steps Q and H.
       (recalculations of parameters if necessary) */

    if (new_big_mu || mu != muprev2) {
        /* Careful! muprev2 is not always == muprev
        because one might have exited in step I or S
        */
        muprev2 = mu;
        omega = M_1_SQRT_2PI / s;
        /* The quantities b1, b2, c3, c2, c1, c0 are for the Hermite
         * approximations to the discrete normal probabilities fk. */

        b1 = one_24 / mu;
        b2 = 0.3 * b1 * b1;
        c3 = one_7 * b1 * b2;
        c2 = b2 - 15. * c3;
        c1 = b1 - 6. * b2 + 45. * c3;
        c0 = 1. - b1 + 3. * b2 - 15. * c3;
        c = 0.1069 / mu; /* guarantees majorization by the 'hat'-function. */
    }

    if (g >= 0.) {
        /* 'Subroutine' F is called (kflag=0 for correct return) */
        kflag = 0;
        goto Step_F;
    }


    repeat {
        /* Step E. Exponential Sample */

        E = igraph_i_exp_rand(rng);/* ~ Exp(1) (standard exponential) */

        /*  sample t from the laplace 'hat'
            (if t <= -0.6744 then pk < fk for all mu >= 10.) */
        u = 2 * igraph_rng_get_unif01(rng) - 1.;
        t = 1.8 + copysign(E, u);
        if (t > -0.6744) {
            pois = floor(mu + s * t);
            fk = pois;
            difmuk = mu - fk;

            /* 'subroutine' F is called (kflag=1 for correct return) */
            kflag = 1;

Step_F: /* 'subroutine' F : calculation of px,py,fx,fy. */

            if (pois < 10) { /* use factorials from table fact[] */
                px = -mu;
                py = pow(mu, pois) / fact[(int)pois];
            } else {
                /* Case pois >= 10 uses polynomial approximation
                   a0-a7 for accuracy when advisable */
                del = one_12 / fk;
                del = del * (1. - 4.8 * del * del);
                v = difmuk / fk;
                if (fabs(v) <= 0.25)
                    px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) *
                                          v + a3) * v + a2) * v + a1) * v + a0)
                    - del;
                else { /* |v| > 1/4 */
                    px = fk * log(1. + v) - difmuk - del;
                }
                py = M_1_SQRT_2PI / sqrt(fk);
            }
            x = (0.5 - difmuk) / s;
            x *= x;/* x^2 */
            fx = -0.5 * x;
            fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
            if (kflag > 0) {
                /* Step H. Hat acceptance (E is repeated on rejection) */
                if (c * fabs(u) <= py * exp(px + E) - fy * exp(fx + E)) {
                    break;
                }
            } else
                /* Step Q. Quotient acceptance (rare case) */
                if (fy - u * fy <= py * exp(px - fx)) {
                    break;
                }
        }/* t > -.67.. */
    }
    return pois;
}

#undef a1
#undef a2
#undef a3
#undef a4
#undef a5
#undef a6
#undef a7

/* This is from nmath/rbinom.c */

#define repeat for(;;)

static double igraph_i_rbinom(igraph_rng_t *rng, igraph_integer_t n, double pp) {

    static IGRAPH_THREAD_LOCAL double c, fm, npq, p1, p2, p3, p4, qn;
    static IGRAPH_THREAD_LOCAL double xl, xll, xlr, xm, xr;

    static IGRAPH_THREAD_LOCAL double psave = -1.0;
    static IGRAPH_THREAD_LOCAL igraph_integer_t nsave = -1;
    static IGRAPH_THREAD_LOCAL igraph_integer_t m;

    double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
    double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
    igraph_integer_t i, ix, k;

    if (!isfinite(pp) ||
        /* n=0, p=0, p=1 are not errors <TSL>*/
        n < 0 || pp < 0. || pp > 1.) {
        ML_ERR_return_NAN;
    }

    if (n == 0 || pp == 0.) {
        return 0;
    }
    if (pp == 1.) {
        return n;
    }

    p = fmin(pp, 1. - pp);
    q = 1. - p;
    np = n * p;
    r = p / q;
    g = r * (n + 1);

    /* Setup, perform only when parameters change [using static (globals): */

    /* FIXING: Want this thread safe
       -- use as little (thread globals) as possible
    */
    if (pp != psave || n != nsave) {
        psave = pp;
        nsave = n;
        if (np < 30.0) {
            /* inverse cdf logic for mean less than 30 */
            qn = pow(q, (double) n);
            goto L_np_small;
        } else {
            ffm = np + p;
            m = ffm;
            fm = m;
            npq = np * q;
            /* Note (igraph): Original code used a cast to (int) for rounding. However,
             * the max npq = n*p*(1-p) value is 0.25*n, thus 2.195 * sqrt(npq) may be
             * as large as 1.0975 * sqrt(n). This is not representable on a 32-bit signed
             * integer when n is a 64-bit signed integer. Thus we use trunc() instead. */
            p1 = trunc(2.195 * sqrt(npq) - 4.6 * q) + 0.5;
            xm = fm + 0.5;
            xl = xm - p1;
            xr = xm + p1;
            c = 0.134 + 20.5 / (15.3 + fm);
            al = (ffm - xl) / (ffm - xl * p);
            xll = al * (1.0 + 0.5 * al);
            al = (xr - ffm) / (xr * q);
            xlr = al * (1.0 + 0.5 * al);
            p2 = p1 * (1.0 + c + c);
            p3 = p2 + c / xll;
            p4 = p3 + c / xlr;
        }
    } else if (n == nsave) {
        if (np < 30.0) {
            goto L_np_small;
        }
    }

    /*-------------------------- np = n*p >= 30 : ------------------- */
    repeat {
        u = igraph_rng_get_unif01(rng) * p4;
        v = igraph_rng_get_unif01(rng);
        /* triangular region */
        if (u <= p1) {
            ix = xm - p1 * v + u;
            goto finis;
        }
        /* parallelogram region */
        if (u <= p2) {
            x = xl + (u - p1) / c;
            v = v * c + 1.0 - fabs(xm - x) / p1;
            if (v > 1.0 || v <= 0.) {
                continue;
            }
            ix = x;
        } else {
            if (u > p3) { /* right tail */
                ix = xr - log(v) / xlr;
                if (ix > n) {
                    continue;
                }
                v = v * (u - p3) * xlr;
            } else {/* left tail */
                ix = xl + log(v) / xll;
                if (ix < 0) {
                    continue;
                }
                v = v * (u - p2) * xll;
            }
        }
        /* determine appropriate way to perform accept/reject test */
        k = imaxabs(ix - m);
        if (k <= 20 || k >= npq / 2 - 1) {
            /* explicit evaluation */
            f = 1.0;
            if (m < ix) {
                for (i = m + 1; i <= ix; i++) {
                    f *= (g / i - r);
                }
            } else if (m != ix) {
                for (i = ix + 1; i <= m; i++) {
                    f /= (g / i - r);
                }
            }
            if (v <= f) {
                goto finis;
            }
        } else {
            /* squeezing using upper and lower bounds on log(f(x)) */
            amaxp = (k / npq) * ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
            ynorm = -k * k / (2.0 * npq);
            alv = log(v);
            if (alv < ynorm - amaxp) {
                goto finis;
            }
            if (alv <= ynorm + amaxp) {
                /* Stirling's formula to machine accuracy */
                /* for the final acceptance/rejection test */
                x1 = ix + 1;
                f1 = fm + 1.0;
                z = n + 1 - fm;
                w = n - ix + 1.0;
                z2 = z * z;
                x2 = x1 * x1;
                f2 = f1 * f1;
                w2 = w * w;
                if (alv <= xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) + (ix - m) * log(w * p / (x1 * q)) + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w / 166320.) {
                    goto finis;
                }
            }
        }
    }

L_np_small:
    /*---------------------- np = n*p < 30 : ------------------------- */

    repeat {
        ix = 0;
        f = qn;
        u = igraph_rng_get_unif01(rng);
        repeat {
            if (u < f) {
                goto finis;
            }
            if (ix > 110) {
                break;
            }
            u -= f;
            ix++;
            f *= (g / ix - r);
        }
    }
finis:
    if (psave > 0.5) {
        ix = n - ix;
    }
    return (double)ix;
}

static igraph_real_t igraph_i_rexp(igraph_rng_t *rng, double rate) {
    igraph_real_t scale = 1.0 / rate;
    if (!isfinite(scale) || scale <= 0.0) {
        if (scale == 0.0) {
            return 0.0;
        }
        return IGRAPH_NAN;
    }
    return scale * igraph_i_exp_rand(rng);
}

/* This is from nmath/rgamma.c */

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000--2008 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double rgamma(double a, double scale);
 *
 *  DESCRIPTION
 *
 *    Random variates from the gamma distribution.
 *
 *  REFERENCES
 *
 *    [1] Shape parameter a >= 1.  Algorithm GD in:
 *
 *    Ahrens, J.H. and Dieter, U. (1982).
 *    Generating gamma variates by a modified
 *    rejection technique.
 *    Comm. ACM, 25, 47-54.
 *
 *
 *    [2] Shape parameter 0 < a < 1. Algorithm GS in:
 *
 *    Ahrens, J.H. and Dieter, U. (1974).
 *    Computer methods for sampling from gamma, beta,
 *    poisson and binomial distributions.
 *    Computing, 12, 223-246.
 *
 *    Input: a = parameter (mean) of the standard gamma distribution.
 *    Output: a variate from the gamma(a)-distribution
 */

static double igraph_i_rgamma(igraph_rng_t *rng, double a, double scale) {
    /* Constants : */
    static const double sqrt32 = 5.656854;
    static const double exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */

    /* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
     * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
     * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
     */
    static const double q1 = 0.04166669;
    static const double q2 = 0.02083148;
    static const double q3 = 0.00801191;
    static const double q4 = 0.00144121;
    static const double q5 = -7.388e-5;
    static const double q6 = 2.4511e-4;
    static const double q7 = 2.424e-4;

    static const double a1 = 0.3333333;
    static const double a2 = -0.250003;
    static const double a3 = 0.2000062;
    static const double a4 = -0.1662921;
    static const double a5 = 0.1423657;
    static const double a6 = -0.1367177;
    static const double a7 = 0.1233795;

    /* State variables: */
    static IGRAPH_THREAD_LOCAL double aa = 0.;
    static IGRAPH_THREAD_LOCAL double aaa = 0.;
    static IGRAPH_THREAD_LOCAL double s, s2, d;    /* no. 1 (step 1) */
    static IGRAPH_THREAD_LOCAL double q0, b, si, c;/* no. 2 (step 4) */

    double e, p, q, r, t, u, v, w, x, ret_val;

    if (!isfinite(a) || !isfinite(scale) || a < 0.0 || scale <= 0.0) {
        if (scale == 0.) {
            return 0.;
        }
        ML_ERR_return_NAN;
    }

    if (a < 1.) { /* GS algorithm for parameters a < 1 */
        if (a == 0) {
            return 0.;
        }
        e = 1.0 + exp_m1 * a;
        repeat {
            p = e * igraph_rng_get_unif01(rng);
            if (p >= 1.0) {
                x = -log((e - p) / a);
                if (igraph_i_exp_rand(rng) >= (1.0 - a) * log(x)) {
                    break;
                }
            } else {
                x = exp(log(p) / a);
                if (igraph_i_exp_rand(rng) >= x) {
                    break;
                }
            }
        }
        return scale * x;
    }

    /* --- a >= 1 : GD algorithm --- */

    /* Step 1: Recalculations of s2, s, d if a has changed */
    if (a != aa) {
        aa = a;
        s2 = a - 0.5;
        s = sqrt(s2);
        d = sqrt32 - s * 12.0;
    }
    /* Step 2: t = standard normal deviate,
               x = (s,1/2) -normal deviate. */

    /* immediate acceptance (i) */
    t = igraph_i_norm_rand(rng);
    x = s + 0.5 * t;
    ret_val = x * x;
    if (t >= 0.0) {
        return scale * ret_val;
    }

    /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
    u = igraph_rng_get_unif01(rng);
    if (d * u <= t * t * t) {
        return scale * ret_val;
    }

    /* Step 4: recalculations of q0, b, si, c if necessary */

    if (a != aaa) {
        aaa = a;
        r = 1.0 / a;
        q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
               + q2) * r + q1) * r;

        /* Approximation depending on size of parameter a */
        /* The constants in the expressions for b, si and c */
        /* were established by numerical experiments */

        if (a <= 3.686) {
            b = 0.463 + s + 0.178 * s2;
            si = 1.235;
            c = 0.195 / s - 0.079 + 0.16 * s;
        } else if (a <= 13.022) {
            b = 1.654 + 0.0076 * s2;
            si = 1.68 / s + 0.275;
            c = 0.062 / s + 0.024;
        } else {
            b = 1.77;
            si = 0.75;
            c = 0.1515 / s;
        }
    }
    /* Step 5: no quotient test if x not positive */

    if (x > 0.0) {
        /* Step 6: calculation of v and quotient q */
        v = t / (s + s);
        if (fabs(v) <= 0.25)
            q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
                                      + a3) * v + a2) * v + a1) * v;
        else {
            q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
        }


        /* Step 7: quotient acceptance (q) */
        if (log(1.0 - u) <= q) {
            return scale * ret_val;
        }
    }

    repeat {
        /* Step 8: e = standard exponential deviate
         *  u =  0,1 -uniform deviate
         *  t = (b,si)-double exponential (laplace) sample */
        e = igraph_i_exp_rand(rng);
        u = igraph_rng_get_unif01(rng);
        u = u + u - 1.0;
        if (u < 0.0) {
            t = b - si * e;
        } else {
            t = b + si * e;
        }
        /* Step  9:  rejection if t < tau(1) = -0.71874483771719 */
        if (t >= -0.71874483771719) {
            /* Step 10:  calculation of v and quotient q */
            v = t / (s + s);
            if (fabs(v) <= 0.25)
                q = q0 + 0.5 * t * t *
                ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
                  + a2) * v + a1) * v;
            else {
                q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
            }
            /* Step 11:  hat acceptance (h) */
            /* (if q not positive go to step 8) */
            if (q > 0.0) {
                w = expm1(q);
                /*  ^^^^^ original code had approximation with rel.err < 2e-7 */
                /* if t is rejected sample again at step 8 */
                if (c * fabs(u) <= w * exp(e - 0.5 * t * t)) {
                    break;
                }
            }
        }
    } /* repeat .. until  `t' is accepted */
    x = s + 0.5 * t;
    return scale * x * x;
}

igraph_error_t igraph_rng_get_dirichlet(igraph_rng_t *rng,
                             const igraph_vector_t *alpha,
                             igraph_vector_t *result) {

    igraph_integer_t len = igraph_vector_size(alpha);
    igraph_real_t sum = 0.0;

    if (len < 2) {
        IGRAPH_ERROR("Dirichlet parameter vector too short, must have at least two entries.",
                     IGRAPH_EINVAL);
    }
    if (igraph_vector_min(alpha) <= 0) {
        IGRAPH_ERROR("Dirichlet concentration parameters must be positive.",
                     IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_resize(result, len));

    for (igraph_integer_t i = 0; i < len; i++) {
        VECTOR(*result)[i] = igraph_rng_get_gamma(rng, VECTOR(*alpha)[i], 1.0);
        sum += VECTOR(*result)[i];
    }
    for (igraph_integer_t i = 0; i < len; i++) {
        VECTOR(*result)[i] /= sum;
    }

    return IGRAPH_SUCCESS;
}
