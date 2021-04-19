/*
 * CXSPARSE: a Concise Sparse Matrix package - Extended.
 * Copyright (c) 2006-2009, Timothy A. Davis.
 * http://www.cise.ufl.edu/research/sparse/CXSparse
 * 
 * CXSparse is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * CXSparse is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "cs.h"
/* allocate a sparse matrix (triplet form or compressed-column form) */
cs *cs_spalloc (CS_INT m, CS_INT n, CS_INT nzmax, CS_INT values, CS_INT triplet)
{
    cs *A = cs_calloc (1, sizeof (cs)) ;    /* allocate the cs struct */
    if (!A) return (NULL) ;                 /* out of memory */
    A->m = m ;                              /* define dimensions and nzmax */
    A->n = n ;
    A->nzmax = nzmax = CS_MAX (nzmax, 1) ;
    A->nz = triplet ? 0 : -1 ;              /* allocate triplet or comp.col */
    A->p = cs_malloc (triplet ? nzmax : n+1, sizeof (CS_INT)) ;
    A->i = cs_malloc (nzmax, sizeof (CS_INT)) ;
    A->x = values ? cs_malloc (nzmax, sizeof (CS_ENTRY)) : NULL ;
    return ((!A->p || !A->i || (values && !A->x)) ? cs_spfree (A) : A) ;
}

/* change the max # of entries sparse matrix */
CS_INT cs_sprealloc (cs *A, CS_INT nzmax)
{
    CS_INT ok, oki, okj = 1, okx = 1 ;
    if (!A) return (0) ;
    if (nzmax <= 0) nzmax = (CS_CSC (A)) ? (A->p [A->n]) : A->nz ;
    A->i = cs_realloc (A->i, nzmax, sizeof (CS_INT), &oki) ;
    if (CS_TRIPLET (A)) A->p = cs_realloc (A->p, nzmax, sizeof (CS_INT), &okj) ;
    if (A->x) A->x = cs_realloc (A->x, nzmax, sizeof (CS_ENTRY), &okx) ;
    ok = (oki && okj && okx) ;
    if (ok) A->nzmax = nzmax ;
    return (ok) ;
}

/* free a sparse matrix */
cs *cs_spfree (cs *A)
{
    if (!A) return (NULL) ;     /* do nothing if A already NULL */
    cs_free (A->p) ;
    cs_free (A->i) ;
    cs_free (A->x) ;
    return (cs_free (A)) ;      /* free the cs struct and return NULL */
}

/* free a numeric factorization */
csn *cs_nfree (csn *N)
{
    if (!N) return (NULL) ;     /* do nothing if N already NULL */
    cs_spfree (N->L) ;
    cs_spfree (N->U) ;
    cs_free (N->pinv) ;
    cs_free (N->B) ;
    return (cs_free (N)) ;      /* free the csn struct and return NULL */
}

/* free a symbolic factorization */
css *cs_sfree (css *S)
{
    if (!S) return (NULL) ;     /* do nothing if S already NULL */
    cs_free (S->pinv) ;
    cs_free (S->q) ;
    cs_free (S->parent) ;
    cs_free (S->cp) ;
    cs_free (S->leftmost) ;
    return (cs_free (S)) ;      /* free the css struct and return NULL */
}

/* allocate a cs_dmperm or cs_scc result */
csd *cs_dalloc (CS_INT m, CS_INT n)
{
    csd *D ;
    D = cs_calloc (1, sizeof (csd)) ;
    if (!D) return (NULL) ;
    D->p = cs_malloc (m, sizeof (CS_INT)) ;
    D->r = cs_malloc (m+6, sizeof (CS_INT)) ;
    D->q = cs_malloc (n, sizeof (CS_INT)) ;
    D->s = cs_malloc (n+6, sizeof (CS_INT)) ;
    return ((!D->p || !D->r || !D->q || !D->s) ? cs_dfree (D) : D) ;
}

/* free a cs_dmperm or cs_scc result */
csd *cs_dfree (csd *D)
{
    if (!D) return (NULL) ;     /* do nothing if D already NULL */
    cs_free (D->p) ;
    cs_free (D->q) ;
    cs_free (D->r) ;
    cs_free (D->s) ;
    return (cs_free (D)) ;
}

/* free workspace and return a sparse matrix result */
cs *cs_done (cs *C, void *w, void *x, CS_INT ok)
{
    cs_free (w) ;                       /* free workspace */
    cs_free (x) ;
    return (ok ? C : cs_spfree (C)) ;   /* return result if OK, else free it */
}

/* free workspace and return CS_INT array result */
CS_INT *cs_idone (CS_INT *p, cs *C, void *w, CS_INT ok)
{
    cs_spfree (C) ;                     /* free temporary matrix */
    cs_free (w) ;                       /* free workspace */
    return (ok ? p : cs_free (p)) ;     /* return result if OK, else free it */
}

/* free workspace and return a numeric factorization (Cholesky, LU, or QR) */
csn *cs_ndone (csn *N, cs *C, void *w, void *x, CS_INT ok)
{
    cs_spfree (C) ;                     /* free temporary matrix */
    cs_free (w) ;                       /* free workspace */
    cs_free (x) ;
    return (ok ? N : cs_nfree (N)) ;    /* return result if OK, else free it */
}

/* free workspace and return a csd result */
csd *cs_ddone (csd *D, cs *C, void *w, CS_INT ok)
{
    cs_spfree (C) ;                     /* free temporary matrix */
    cs_free (w) ;                       /* free workspace */
    return (ok ? D : cs_dfree (D)) ;    /* return result if OK, else free it */
}
