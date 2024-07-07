#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
    #define __BEGIN_DECLS extern "C" {
    #define __END_DECLS }
#else
    #define __BEGIN_DECLS /* empty */
    #define __END_DECLS /* empty */
#endif

/* This is to eliminate gcc warnings about unused parameters */
#define IGRAPH_UNUSED(x) (void)(x)

/* The pure function attribute of GCC-compatible compilers indicates
 * that the function does not have side-effects, i.e. it does not
 * modify global memory. This enables additional compiler optimizations
 * such as common subexpression elimination.
 *
 * The const attribute is similar but with much more stringent requirements.
 * The function must also not read global memory. Generally, const functions
 * should not take pointers, and must compute the return value solely based
 * on their input.
 */
#ifdef __GNUC__
#define IGRAPH_FUNCATTR_PURE __attribute__((__pure__))
#define IGRAPH_FUNCATTR_CONST __attribute__((__const__))
#else
#define IGRAPH_FUNCATTR_PURE
#define IGRAPH_FUNCATTR_CONST
#endif

/* IGRAPH_ASSUME() provides hints to the compiler about conditions
 * that are true yet the compiler cannot deduce. Use with great care.
 * Assuming a condition that is not actually true leads to undefined behaviour. */
#if defined(__clang__)
   /* For Clang, see https://clang.llvm.org/docs/LanguageExtensions.html */
#  if __has_builtin(__builtin_assume)
#    define IGRAPH_ASSUME(expr) __builtin_assume(expr)
#  else
#    define IGRAPH_ASSUME(expr) /* empty */
#  endif
#elif defined(__GNUC__) && !defined(__ICC)
   /* Introduced in GCC 4.5, https://gcc.gnu.org/gcc-4.5/changes.html */
#  define IGRAPH_ASSUME(expr) do { if (expr) {} else { __builtin_unreachable(); } } while (0)
#elif defined(_MSC_VER) || defined(__ICC)
#  define IGRAPH_ASSUME(expr) __assume(expr)
#else
#  define IGRAPH_ASSUME(expr) /* empty */
#endif

/* IGRAPH_I_STRINGIFY(X) evaluates X and converts the result to a string. */
#define IGRAPH_I_STRINGIFY_I(X) #X
#define IGRAPH_I_STRINGIFY(X) IGRAPH_I_STRINGIFY_I(X)

/* Include the definition of macros controlling symbol visibility */
#include "igraph_export.h"

/* Used instead of IGRAPH_EXPORT with functions that need to be tested,
 * but are not part of the public API. */
#define IGRAPH_PRIVATE_EXPORT IGRAPH_EXPORT
