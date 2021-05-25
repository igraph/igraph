
#ifndef CLIQUERCONF_H
#define CLIQUERCONF_H

/*
 * setelement is the basic memory type used in sets.  It is often fastest
 * to be as large as can fit into the CPU registers.
 *
 * ELEMENTSIZE is the size of one setelement, measured in bits.  It must
 * be either 16, 32 or 64  (otherwise additional changes must be made to
 * the source).
 *
 * The default is to use "unsigned long int" and attempt to guess the
 * size using <limits.h>, which should work pretty well.  Check functioning
 * with "make test".
 */

/* typedef unsigned long int setelement; */
/* #define ELEMENTSIZE 64 */


/*
 * INLINE is a command prepended to function declarations to instruct the
 * compiler to inline the function.  If inlining is not desired, define blank.
 *
 * The default is to use "inline", which is recognized by most compilers.
 */

/* #define INLINE */
/* #define INLINE __inline__ */
#if __STDC_VERSION__ >= 199901L
 #define INLINE inline
#else
 #if defined(_MSC_VER)
  #define INLINE __inline
 #elif defined(__GNUC__)
  #define INLINE __inline__
 #else
  #define INLINE
 #endif
#endif


/*
 * Set handling functions are defined as static functions in set.h for
 * performance reasons.  This may cause unnecessary warnings from the
 * compiler.  Some compilers (such as GCC) have the possibility to turn
 * off the warnings on a per-function basis using a flag prepended to
 * the function declaration.
 *
 * The default is to use the correct attribute when compiling with GCC,
 * or no flag otherwise.
 */

/* #define UNUSED_FUNCTION __attribute__((unused)) */
/* #define UNUSED_FUNCTION */


/*
 * Uncommenting the following will disable all assertions  (checks that
 * function arguments and other variables are correct).  This is highly
 * discouraged, as it allows bugs to go unnoticed easier.  The assertions
 * are set so that they do not slow down programs notably.
 */

/* #define ASSERT(x) */

#endif /* !CLIQUERCONF_H */
