/* mt.h
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

#ifndef __MT_H__
#define __MT_H__

#ifdef _MSC_VER
#  define uint32_t __int32
#else
#  include <stdint.h>
#endif

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#define MT_LEN       624

/**
 * \def MT_RAND_MAX
 *
 * The maximum random number that \c mt_random() can generate.
 */
#define MT_RAND_MAX 0xFFFFFFFF

/**
 * Struct that stores the internal state of a Mersenne Twister random number
 * generator.
 */
typedef struct {
    int mt_index;
    uint32_t mt_buffer[MT_LEN];
} mt_rng_t;

/**
 * \brief Initializes a Mersenne Twister random number generator.
 *
 * The random number generator is seeded with random 32-bit numbers obtained
 * from the \em built-in random number generator using consecutive calls to
 * \c rand().
 *
 * \param  rng  the random number generator to initialize
 */
void mt_init(mt_rng_t* rng);

/**
 * \brief Initializes a Mersenne Twister random number generator, seeding it
 *        from another one.
 *
 * The random number generator is seeded with random 32-bit numbers obtained
 * from another, initialized Mersenne Twister random number generator.
 *
 * \param  rng     the random number generator to initialize
 * \param  seeder  the random number generator that will seed the one being
 *                 initialized. When null, the random number generator will
 *                 be initialized from the built-in RNG as if \ref mt_init()
 *                 was called.
 */
void mt_init_from_rng(mt_rng_t* rng, mt_rng_t* seeder);

/**
 * \brief Returns the next 32-bit random number from the given Mersenne Twister
 * random number generator.
 *
 * \param  rng  the random number generator to use
 * \return the next 32-bit random number from the generator
 */
uint32_t mt_random(mt_rng_t* rng);

/**
 * \brief Returns a uniformly distributed double from the interval [0;1)
 * based on the next value of the given Mersenne Twister random number
 * generator.
 *
 * \param  rng  the random number generator to use
 * \return a uniformly distributed random number from the interval [0;1)
 */
double mt_uniform_01(mt_rng_t* rng);

__END_DECLS

#endif


