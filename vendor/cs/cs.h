/* This is a MODIFIED version of the original CXSparse/Include/cs.h file from
 * SuiteSparse 5.12.0 (CXSparse version 3.2.0). The modifications are outlined
 * here:
 * 
 * - Dependency on SuiteSparse_long was removed
 * - CXSparse is configured to use igraph_integer_t as cs_long_t
 * - CXSparse function prefix is set to cs_igraph instead of cs_igraph
 * - Unneeded CXSparse function variants are removed
 *
 * The remaining comments below are from the original cs.h header */

/* ========================================================================== */
/* CXSparse/Include/cs.h file */
/* ========================================================================== */

/* This is the CXSparse/Include/cs.h file.  It has the same name (cs.h) as
   the CSparse/Include/cs.h file.  The 'make install' for SuiteSparse installs
   CXSparse, and this file, instead of CSparse.  The two packages have the same
   cs.h include filename, because CXSparse is a superset of CSparse.  Any user
   program that uses CSparse can rely on CXSparse instead, with no change to the
   user code.  The #include "cs.h" line will work for both versions, in user
   code, and the function names and user-visible typedefs from CSparse all
   appear in CXSparse.  For experimenting and changing the package itself, I
   recommend using CSparse since it's simpler and easier to modify.  For
   using the package in production codes, I recommend CXSparse since it has
   more features (support for complex matrices, and both int and long
   versions).
 */

/* ========================================================================== */

#ifndef _CXS_H
#define _CXS_H
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

#include "igraph_types.h"

#ifdef __cplusplus
#ifndef NCOMPLEX
#include <complex>
typedef std::complex<double> cs_complex_t ;
#endif
extern "C" {
#else
#ifndef NCOMPLEX
#include <complex.h>
#define cs_complex_t double _Complex
#endif
#endif

#define CS_VER 3                    /* CXSparse Version */
#define CS_SUBVER 2
#define CS_SUBSUB 0
#define CS_DATE "Sept 12, 2017"       /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006-2016"
#define CXSPARSE

#define cs_long_t       igraph_integer_t
#define cs_long_t_id    "%" IGRAPH_PRId
#define cs_long_t_max   IGRAPH_INTEGER_MAX

/* -------------------------------------------------------------------------- */
/* double/cs_long_t version of CXSparse */
/* -------------------------------------------------------------------------- */

/* --- primary CSparse routines and data structures ------------------------- */

typedef struct cs_igraph_sparse  /* matrix in compressed-column or triplet form */
{
    cs_long_t nzmax ; /* maximum number of entries */
    cs_long_t m ;     /* number of rows */
    cs_long_t n ;     /* number of columns */
    cs_long_t *p ;    /* column pointers (size n+1) or col indlces (size nzmax) */
    cs_long_t *i ;    /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    cs_long_t nz ;    /* # of entries in triplet matrix, -1 for compressed-col */
} cs_igraph ;

cs_igraph *cs_igraph_add (const cs_igraph *A, const cs_igraph *B, double alpha, double beta) ;
cs_long_t cs_igraph_cholsol (cs_long_t order, const cs_igraph *A, double *b) ;
cs_long_t cs_igraph_dupl (cs_igraph *A) ;
cs_long_t cs_igraph_entry (cs_igraph *T, cs_long_t i, cs_long_t j, double x) ;
cs_long_t cs_igraph_lusol (cs_long_t order, const cs_igraph *A, double *b, double tol) ;
cs_long_t cs_igraph_gaxpy (const cs_igraph *A, const double *x, double *y) ;
cs_igraph *cs_igraph_multiply (const cs_igraph *A, const cs_igraph *B) ;
cs_long_t cs_igraph_qrsol (cs_long_t order, const cs_igraph *A, double *b) ;
cs_igraph *cs_igraph_transpose (const cs_igraph *A, cs_long_t values) ;
cs_igraph *cs_igraph_compress (const cs_igraph *T) ;
double cs_igraph_norm (const cs_igraph *A) ;
/*cs_long_t cs_igraph_print (const cs_igraph *A, cs_long_t brief) ;*/
cs_igraph *cs_igraph_load (FILE *f) ;

/* utilities */
void *cs_igraph_calloc (cs_long_t n, size_t size) ;
void *cs_igraph_free (void *p) ;
void *cs_igraph_realloc (void *p, cs_long_t n, size_t size, cs_long_t *ok) ;
cs_igraph *cs_igraph_spalloc (cs_long_t m, cs_long_t n, cs_long_t nzmax, cs_long_t values,
    cs_long_t t) ;
cs_igraph *cs_igraph_spfree (cs_igraph *A) ;
cs_long_t cs_igraph_sprealloc (cs_igraph *A, cs_long_t nzmax) ;
void *cs_igraph_malloc (cs_long_t n, size_t size) ;

/* --- secondary CSparse routines and data structures ----------------------- */

typedef struct cs_igraph_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    cs_long_t *pinv ;     /* inverse row perm. for QR, fill red. perm for Chol */
    cs_long_t *q ;        /* fill-reducing column permutation for LU and QR */
    cs_long_t *parent ;   /* elimination tree for Cholesky and QR */
    cs_long_t *cp ;       /* column pointers for Cholesky, row counts for QR */
    cs_long_t *leftmost ; /* leftmost[i] = min(find(A(i,:))), for QR */
    cs_long_t m2 ;        /* # of rows for QR, after adding fictitious rows */
    double lnz ;        /* # entries in L for LU or Cholesky; in V for QR */
    double unz ;        /* # entries in U for LU; in R for QR */
} cs_igraphs ;

typedef struct cs_igraph_numeric   /* numeric Cholesky, LU, or QR factorization */
{
    cs_igraph *L ;      /* L for LU and Cholesky, V for QR */
    cs_igraph *U ;      /* U for LU, r for QR, not used for Cholesky */
    cs_long_t *pinv ; /* partial pivoting for LU */
    double *B ;     /* beta [0..n-1] for QR */
} cs_igraphn ;

typedef struct cs_igraph_dmperm_results    /* cs_igraph_dmperm or cs_igraph_scc output */
{
    cs_long_t *p ;    /* size m, row permutation */
    cs_long_t *q ;    /* size n, column permutation */
    cs_long_t *r ;    /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
    cs_long_t *s ;    /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
    cs_long_t nb ;    /* # of blocks in fine dmperm decomposition */
    cs_long_t rr [5] ;    /* coarse row decomposition */
    cs_long_t cc [5] ;    /* coarse column decomposition */
} cs_igraphd ;

cs_long_t *cs_igraph_amd (cs_long_t order, const cs_igraph *A) ;
cs_igraphn *cs_igraph_chol (const cs_igraph *A, const cs_igraphs *S) ;
cs_igraphd *cs_igraph_dmperm (const cs_igraph *A, cs_long_t seed) ;
cs_long_t cs_igraph_droptol (cs_igraph *A, double tol) ;
cs_long_t cs_igraph_dropzeros (cs_igraph *A) ;
cs_long_t cs_igraph_happly (const cs_igraph *V, cs_long_t i, double beta, double *x) ;
cs_long_t cs_igraph_ipvec (const cs_long_t *p, const double *b, double *x, cs_long_t n) ;
cs_long_t cs_igraph_lsolve (const cs_igraph *L, double *x) ;
cs_long_t cs_igraph_ltsolve (const cs_igraph *L, double *x) ;
cs_igraphn *cs_igraph_lu (const cs_igraph *A, const cs_igraphs *S, double tol) ;
cs_igraph *cs_igraph_permute (const cs_igraph *A, const cs_long_t *pinv, const cs_long_t *q,
    cs_long_t values) ;
cs_long_t *cs_igraph_pinv (const cs_long_t *p, cs_long_t n) ;
cs_long_t cs_igraph_pvec (const cs_long_t *p, const double *b, double *x, cs_long_t n) ;
cs_igraphn *cs_igraph_qr (const cs_igraph *A, const cs_igraphs *S) ;
cs_igraphs *cs_igraph_schol (cs_long_t order, const cs_igraph *A) ;
cs_igraphs *cs_igraph_sqr (cs_long_t order, const cs_igraph *A, cs_long_t qr) ;
cs_igraph *cs_igraph_symperm (const cs_igraph *A, const cs_long_t *pinv, cs_long_t values) ;
cs_long_t cs_igraph_usolve (const cs_igraph *U, double *x) ;
cs_long_t cs_igraph_utsolve (const cs_igraph *U, double *x) ;
cs_long_t cs_igraph_updown (cs_igraph *L, cs_long_t sigma, const cs_igraph *C,
    const cs_long_t *parent) ;

/* utilities */
cs_igraphs *cs_igraph_sfree (cs_igraphs *S) ;
cs_igraphn *cs_igraph_nfree (cs_igraphn *N) ;
cs_igraphd *cs_igraph_dfree (cs_igraphd *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */

cs_long_t *cs_igraph_counts (const cs_igraph *A, const cs_long_t *parent,
    const cs_long_t *post, cs_long_t ata) ;
double cs_igraph_cumsum (cs_long_t *p, cs_long_t *c, cs_long_t n) ;
cs_long_t cs_igraph_dfs (cs_long_t j, cs_igraph *G, cs_long_t top, cs_long_t *xi,
    cs_long_t *pstack, const cs_long_t *pinv) ;
cs_long_t *cs_igraph_etree (const cs_igraph *A, cs_long_t ata) ;
cs_long_t cs_igraph_fkeep (cs_igraph *A,
    cs_long_t (*fkeep) (cs_long_t, cs_long_t, double, void *), void *other) ;
double cs_igraph_house (double *x, double *beta, cs_long_t n) ;
cs_long_t *cs_igraph_maxtrans (const cs_igraph *A, cs_long_t seed) ;
cs_long_t *cs_igraph_post (const cs_long_t *parent, cs_long_t n) ;
cs_igraphd *cs_igraph_scc (cs_igraph *A) ;
cs_long_t cs_igraph_scatter (const cs_igraph *A, cs_long_t j, double beta, cs_long_t *w,
    double *x, cs_long_t mark,cs_igraph *C, cs_long_t nz) ;
cs_long_t cs_igraph_tdfs (cs_long_t j, cs_long_t k, cs_long_t *head, const cs_long_t *next,
    cs_long_t *post, cs_long_t *stack) ;
cs_long_t cs_igraph_leaf (cs_long_t i, cs_long_t j, const cs_long_t *first,
    cs_long_t *maxfirst, cs_long_t *prevleaf, cs_long_t *ancestor, cs_long_t *jleaf) ;
cs_long_t cs_igraph_reach (cs_igraph *G, const cs_igraph *B, cs_long_t k, cs_long_t *xi,
    const cs_long_t *pinv) ;
cs_long_t cs_igraph_spsolve (cs_igraph *L, const cs_igraph *B, cs_long_t k, cs_long_t *xi,
    double *x, const cs_long_t *pinv, cs_long_t lo) ;
cs_long_t cs_igraph_ereach (const cs_igraph *A, cs_long_t k, const cs_long_t *parent,
    cs_long_t *s, cs_long_t *w) ;
cs_long_t *cs_igraph_randperm (cs_long_t n, cs_long_t seed) ;

/* utilities */
cs_igraphd *cs_igraph_dalloc (cs_long_t m, cs_long_t n) ;
cs_igraph *cs_igraph_done (cs_igraph *C, void *w, void *x, cs_long_t ok) ;
cs_long_t *cs_igraph_idone (cs_long_t *p, cs_igraph *C, void *w, cs_long_t ok) ;
cs_igraphn *cs_igraph_ndone (cs_igraphn *N, cs_igraph *C, void *w, void *x, cs_long_t ok) ;
cs_igraphd *cs_igraph_ddone (cs_igraphd *D, cs_igraph *C, void *w, cs_long_t ok) ;

/* -------------------------------------------------------------------------- */
/* Macros for constructing each version of CSparse */
/* -------------------------------------------------------------------------- */

#define CS_INT cs_long_t
#define CS_INT_MAX cs_long_t_max
#define CS_ID cs_long_t_id
#define CS_ENTRY double
#define CS_NAME(nm) cs_igraph ## nm
#define cs cs_igraph

#define CS_REAL(x) (x)
#define CS_IMAG(x) (0.)
#define CS_CONJ(x) (x)
#define CS_ABS(x) fabs(x)

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define CS_CSC(A) (A && (A->nz == -1))
#define CS_TRIPLET(A) (A && (A->nz >= 0))

/* --- primary CSparse routines and data structures ------------------------- */

#define cs_add CS_NAME (_add)
#define cs_cholsol CS_NAME (_cholsol)
#define cs_dupl CS_NAME (_dupl)
#define cs_entry CS_NAME (_entry)
#define cs_lusol CS_NAME (_lusol)
#define cs_gaxpy CS_NAME (_gaxpy)
#define cs_multiply CS_NAME (_multiply)
#define cs_qrsol CS_NAME (_qrsol)
#define cs_transpose CS_NAME (_transpose)
#define cs_compress CS_NAME (_compress)
#define cs_norm CS_NAME (_norm)
/*#define cs_print CS_NAME (_print)*/
#define cs_load CS_NAME (_load)

/* utilities */
#define cs_calloc CS_NAME (_calloc)
#define cs_free CS_NAME (_free)
#define cs_realloc CS_NAME (_realloc)
#define cs_spalloc CS_NAME (_spalloc)
#define cs_spfree CS_NAME (_spfree)
#define cs_sprealloc CS_NAME (_sprealloc)
#define cs_malloc CS_NAME (_malloc)

/* --- secondary CSparse routines and data structures ----------------------- */
#define css CS_NAME (s)
#define csn CS_NAME (n)
#define csd CS_NAME (d)

#define cs_amd CS_NAME (_amd)
#define cs_chol CS_NAME (_chol)
#define cs_dmperm CS_NAME (_dmperm)
#define cs_droptol CS_NAME (_droptol)
#define cs_dropzeros CS_NAME (_dropzeros)
#define cs_happly CS_NAME (_happly)
#define cs_ipvec CS_NAME (_ipvec)
#define cs_lsolve CS_NAME (_lsolve)
#define cs_ltsolve CS_NAME (_ltsolve)
#define cs_lu CS_NAME (_lu)
#define cs_permute CS_NAME (_permute)
#define cs_pinv CS_NAME (_pinv)
#define cs_pvec CS_NAME (_pvec)
#define cs_qr CS_NAME (_qr)
#define cs_schol CS_NAME (_schol)
#define cs_sqr CS_NAME (_sqr)
#define cs_symperm CS_NAME (_symperm)
#define cs_usolve CS_NAME (_usolve)
#define cs_utsolve CS_NAME (_utsolve)
#define cs_updown CS_NAME (_updown)

/* utilities */
#define cs_sfree CS_NAME (_sfree)
#define cs_nfree CS_NAME (_nfree)
#define cs_dfree CS_NAME (_dfree)

/* --- tertiary CSparse routines -------------------------------------------- */
#define cs_counts CS_NAME (_counts)
#define cs_cumsum CS_NAME (_cumsum)
#define cs_dfs CS_NAME (_dfs)
#define cs_etree CS_NAME (_etree)
#define cs_fkeep CS_NAME (_fkeep)
#define cs_house CS_NAME (_house)
#define cs_invmatch CS_NAME (_invmatch)
#define cs_maxtrans CS_NAME (_maxtrans)
#define cs_post CS_NAME (_post)
#define cs_scc CS_NAME (_scc)
#define cs_scatter CS_NAME (_scatter)
#define cs_tdfs CS_NAME (_tdfs)
#define cs_reach CS_NAME (_reach)
#define cs_spsolve CS_NAME (_spsolve)
#define cs_ereach CS_NAME (_ereach)
#define cs_randperm CS_NAME (_randperm)
#define cs_leaf CS_NAME (_leaf)

/* utilities */
#define cs_dalloc CS_NAME (_dalloc)
#define cs_done CS_NAME (_done)
#define cs_idone CS_NAME (_idone)
#define cs_ndone CS_NAME (_ndone)
#define cs_ddone CS_NAME (_ddone)

#ifdef __cplusplus
}
#endif
#endif
