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
/* ordering and symbolic analysis for a Cholesky factorization */
css *cs_schol (CS_INT order, const cs *A)
{
    CS_INT n, *c, *post, *P ;
    cs *C ;
    css *S ;
    if (!CS_CSC (A)) return (NULL) ;        /* check inputs */
    n = A->n ;
    S = cs_calloc (1, sizeof (css)) ;       /* allocate result S */
    if (!S) return (NULL) ;                 /* out of memory */
    P = cs_amd (order, A) ;                 /* P = amd(A+A'), or natural */
    S->pinv = cs_pinv (P, n) ;              /* find inverse permutation */
    cs_free (P) ;
    if (order && !S->pinv) return (cs_sfree (S)) ;
    C = cs_symperm (A, S->pinv, 0) ;        /* C = spones(triu(A(P,P))) */
    S->parent = cs_etree (C, 0) ;           /* find etree of C */
    post = cs_post (S->parent, n) ;         /* postorder the etree */
    c = cs_counts (C, S->parent, post, 0) ; /* find column counts of chol(C) */
    cs_free (post) ;
    cs_spfree (C) ;
    S->cp = cs_malloc (n+1, sizeof (CS_INT)) ; /* allocate result S->cp */
    S->unz = S->lnz = cs_cumsum (S->cp, c, n) ; /* find column pointers for L */
    cs_free (c) ;
    return ((S->lnz >= 0) ? S : cs_sfree (S)) ;
}
