/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.in by autoheader.  */

/* Define to 1 if you have the <dlfcn.h> header file. */
/* #undef HAVE_DLFCN_H */

/* Define to 1 if you have the `expm1' function. */
/* #undef HAVE_EXPM1 */

/* Define to 1 if you have the `fabsl' function. */
#define HAVE_FABSL 1

/* Define to 1 if you have the `finite' function. */
#define HAVE_FINITE 1

/* Define to 1 if you have the `fmin' function. */
/* #undef HAVE_FMIN */

/* Define to 1 if you have the `ftruncate' function. */
/* #undef HAVE_FTRUNCATE */

/* Define to 1 if you have the GLPK library */
/* #undef HAVE_GLPK */

/* Define to 1 if you have the GMP library */
/* #undef HAVE_GMP */

/* Define to 1 if you have the <inttypes.h> header file. */
/* #undef HAVE_INTTYPES_H */

/* Define to 1 if you have the `isnan' function. */
#define HAVE_ISNAN 1

/* Define to 1 if you have the `arpack' library (-larpack). */
/* #undef HAVE_LIBARPACK */

/* Define to 1 if you have the `blas' library (-lblas). */
/* #undef HAVE_LIBBLAS */

/* Define to 1 if you have the `f2c' library (-lf2c). */
/* #undef HAVE_LIBF2C */

/* Define to 1 if you have the `lapack' library (-llapack). */
/* #undef HAVE_LIBLAPACK */

/* Define to 1 if you have the libxml2 libraries installed */
#define HAVE_LIBXML 0

/* Define to 1 if you have the `log1p' function. */
/* #undef HAVE_LOG1P */

/* Define to 1 if you have the `log2' function. */
/* #undef HAVE_LOG2 */

/* Define to 1 if you have the `logbl' function. */
/* #undef HAVE_LOGBL */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `rint' function. */
/* #undef HAVE_RINT */

/* Define to 1 if you have the `rintf' function. */
/* #undef HAVE_RINTF */

/* Define to 1 if you have the `round' function. */
/* #undef HAVE_ROUND */

/* Define to 1 if you have the `snprintf' function. */
#define HAVE_SNPRINTF 1

/* Define to 1 if you have the <stdarg.h> header file. */
#define HAVE_STDARG_H 1

/* Define to 1 if you have the <stdint.h> header file. */
/* #undef HAVE_STDINT_H */

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strcasecmp' function. */
#define HAVE_STRCASECMP 1

/* Define to 1 if you have the `strdup' function. */
/* #undef HAVE_STRDUP */

/* Define to 1 if you have the <strings.h> header file. */
/* #undef HAVE_STRINGS_H */

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/int_types.h> header file. */
/* #undef HAVE_SYS_INT_TYPES_H */

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the sys/times.h header */
/* #undef HAVE_TIMES_H */

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* Define to 1 if you want to use thread-local storage for global igraph
   structures */
/* #undef HAVE_TLS */

/* Define to 1 if you have the <unistd.h> header file. */
/* #undef HAVE_UNISTD_H */

/* Define to 1 if you have the `_strdup' function. */
#define HAVE__STRDUP 1

/* Keyword for thread local storage, or just static if not available */
#define IGRAPH_F77_SAVE static

/* Keyword for thread local storage, or empty if not available */
#define IGRAPH_THREAD_LOCAL

/* Define to 1 if you use the internal ARPACK library */
#define INTERNAL_ARPACK 1

/* Define to 1 if you use the internal BLAS library */
#define INTERNAL_BLAS 1

/* Define to 1 if you use the internal F2C library */
#define INTERNAL_F2C 1

/* Define to 1 if you use the internal GLPK library */
#define INTERNAL_GLPK 1
/* Define to 1 if you use the internal LAPACK library */
#define INTERNAL_LAPACK 1

/* Name of package */
#define PACKAGE "igraph"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "csardi.gabor@gmail.com"

/* Define to the full name of this package. */
#define PACKAGE_NAME "igraph"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "igraph 0.6"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "igraph"

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.6"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "0.6"

/* Define to 1 if `lex' declares `yytext' as a `char *' by default, not a
   `char[]'. */
#define YYTEXT_POINTER 1

#pragma warning (disable:4244)

#define strcasecmp _stricmp

#define isnan _isnan
#define finite _finite
#define hypot _hypot
#include <float.h>

#define snprintf igraph_i_snprintf

/* To turn off some warnings about fscanf, strcpy etc */
#define _CRT_SECURE_NO_WARNINGS 1

