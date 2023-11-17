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

/* Include the definition of macros controlling symbol visibility */
#include "igraph_export.h"

/* Used instead of IGRAPH_EXPORT with functions that need to be tested,
 * but are not part of the public API. */
#define IGRAPH_PRIVATE_EXPORT IGRAPH_EXPORT

/* The pure function attribute of GCC-compatible compilers indicates
 * that the function does not have side-effects, i.e. it does not
 * modify global memory. This enables additional compiler optimizations
 * such as common subexpression elimination. */
#ifdef __GNUC__
#define IGRAPH_FUNCATTR_PURE __attribute__((__pure__))
#else
#define IGRAPH_FUNCATTR_PURE
#endif
