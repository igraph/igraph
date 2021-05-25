/* amd_internal.h */

/* Written by Andrew Makhorin <mao@gnu.org>. */

#ifndef AMD_INTERNAL_H
#define AMD_INTERNAL_H

/* AMD will be exceedingly slow when running in debug mode. */
#if 1
#ifndef NDEBUG
#define NDEBUG
#endif
#endif

#include "amd.h"
#define _GLPSTD_STDIO
#include "glpenv.h"

#define Int int
#define ID "%d"
#define Int_MAX INT_MAX

#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t)(-1))
#endif

#define EMPTY (-1)
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) ((i < EMPTY) ? FLIP (i) : (i))

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define IMPLIES(p, q) (!(p) || (q))

#define GLOBAL

#define AMD_order amd_order
#define AMD_defaults amd_defaults
#define AMD_control amd_control
#define AMD_info amd_info
#define AMD_1 amd_1
#define AMD_2 amd_2
#define AMD_valid amd_valid
#define AMD_aat amd_aat
#define AMD_postorder amd_postorder
#define AMD_post_tree amd_post_tree
#define AMD_dump amd_dump
#define AMD_debug amd_debug
#define AMD_debug_init amd_debug_init
#define AMD_preprocess amd_preprocess

#define amd_malloc xmalloc
#if 0 /* 24/V-2009 */
#define amd_free xfree
#else
#define amd_free(ptr) { if ((ptr) != NULL) xfree(ptr); } 
#endif
#define amd_printf xprintf

#define PRINTF(params) { amd_printf params; }

#ifndef NDEBUG
#define ASSERT(expr) xassert(expr)
#define AMD_DEBUG0(params) { PRINTF(params); }
#define AMD_DEBUG1(params) { if (AMD_debug >= 1) PRINTF(params); }
#define AMD_DEBUG2(params) { if (AMD_debug >= 2) PRINTF(params); }
#define AMD_DEBUG3(params) { if (AMD_debug >= 3) PRINTF(params); }
#define AMD_DEBUG4(params) { if (AMD_debug >= 4) PRINTF(params); }
#else
#define ASSERT(expression)
#define AMD_DEBUG0(params)
#define AMD_DEBUG1(params)
#define AMD_DEBUG2(params)
#define AMD_DEBUG3(params)
#define AMD_DEBUG4(params)
#endif

#define amd_aat _glp_amd_aat
size_t AMD_aat(Int n, const Int Ap[], const Int Ai[], Int Len[],
      Int Tp[], double Info[]);

#define amd_1 _glp_amd_1
void AMD_1(Int n, const Int Ap[], const Int Ai[], Int P[], Int Pinv[],
      Int Len[], Int slen, Int S[], double Control[], double Info[]);

#define amd_postorder _glp_amd_postorder
void AMD_postorder(Int nn, Int Parent[], Int Npiv[], Int Fsize[],
      Int Order[], Int Child[], Int Sibling[], Int Stack[]);

#define amd_post_tree _glp_amd_post_tree
#ifndef NDEBUG
Int AMD_post_tree(Int root, Int k, Int Child[], const Int Sibling[],
      Int Order[], Int Stack[], Int nn);
#else
Int AMD_post_tree(Int root, Int k, Int Child[], const Int Sibling[],
      Int Order[], Int Stack[]);
#endif

#define amd_preprocess _glp_amd_preprocess
void AMD_preprocess(Int n, const Int Ap[], const Int Ai[], Int Rp[],
      Int Ri[], Int W[], Int Flag[]);

#define amd_debug _glp_amd_debug
extern Int AMD_debug;

#define amd_debug_init _glp_amd_debug_init
void AMD_debug_init(char *s);

#define amd_dump _glp_amd_dump
void AMD_dump(Int n, Int Pe[], Int Iw[], Int Len[], Int iwlen,
      Int pfree, Int Nv[], Int Next[], Int Last[], Int Head[],
      Int Elen[], Int Degree[], Int W[], Int nel);

#endif

/* eof */
