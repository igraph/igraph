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
/* [L,U,pinv]=lu(A, [q lnz unz]). lnz and unz can be guess */
csn *cs_lu (const cs *A, const css *S, double tol)
{
    cs *L, *U ;
    csn *N ;
    CS_ENTRY pivot, *Lx, *Ux, *x ;
    double a, t ;
    CS_INT *Lp, *Li, *Up, *Ui, *pinv, *xi, *q, n, ipiv, k, top, p, i, col, lnz,unz;
    if (!CS_CSC (A) || !S) return (NULL) ;          /* check inputs */
    n = A->n ;
    q = S->q ; lnz = S->lnz ; unz = S->unz ;
    x = cs_malloc (n, sizeof (CS_ENTRY)) ;            /* get CS_ENTRY workspace */
    xi = cs_malloc (2*n, sizeof (CS_INT)) ;            /* get CS_INT workspace */
    N = cs_calloc (1, sizeof (csn)) ;               /* allocate result */
    if (!x || !xi || !N) return (cs_ndone (N, NULL, xi, x, 0)) ;
    N->L = L = cs_spalloc (n, n, lnz, 1, 0) ;       /* allocate result L */
    N->U = U = cs_spalloc (n, n, unz, 1, 0) ;       /* allocate result U */
    N->pinv = pinv = cs_malloc (n, sizeof (CS_INT)) ;  /* allocate result pinv */
    if (!L || !U || !pinv) return (cs_ndone (N, NULL, xi, x, 0)) ;
    Lp = L->p ; Up = U->p ;
    for (i = 0 ; i < n ; i++) x [i] = 0 ;           /* clear workspace */
    for (i = 0 ; i < n ; i++) pinv [i] = -1 ;       /* no rows pivotal yet */
    for (k = 0 ; k <= n ; k++) Lp [k] = 0 ;         /* no cols of L yet */
    lnz = unz = 0 ;
    for (k = 0 ; k < n ; k++)       /* compute L(:,k) and U(:,k) */
    {
        /* --- Triangular solve --------------------------------------------- */
        Lp [k] = lnz ;              /* L(:,k) starts here */
        Up [k] = unz ;              /* U(:,k) starts here */
        if ((lnz + n > L->nzmax && !cs_sprealloc (L, 2*L->nzmax + n)) ||
            (unz + n > U->nzmax && !cs_sprealloc (U, 2*U->nzmax + n)))
        {
            return (cs_ndone (N, NULL, xi, x, 0)) ;
        }
        Li = L->i ; Lx = L->x ; Ui = U->i ; Ux = U->x ;
        col = q ? (q [k]) : k ;
        top = cs_spsolve (L, A, col, xi, x, pinv, 1) ;  /* x = L\A(:,col) */
        /* --- Find pivot --------------------------------------------------- */
        ipiv = -1 ;
        a = -1 ;
        for (p = top ; p < n ; p++)
        {
            i = xi [p] ;            /* x(i) is nonzero */
            if (pinv [i] < 0)       /* row i is not yet pivotal */
            {
                if ((t = CS_ABS (x [i])) > a)
                {
                    a = t ;         /* largest pivot candidate so far */
                    ipiv = i ;
                }
            }
            else                    /* x(i) is the entry U(pinv[i],k) */
            {
                Ui [unz] = pinv [i] ;
                Ux [unz++] = x [i] ;
            }
        }
        if (ipiv == -1 || a <= 0) return (cs_ndone (N, NULL, xi, x, 0)) ;
        if (pinv [col] < 0 && CS_ABS (x [col]) >= a*tol) ipiv = col ;
        /* --- Divide by pivot ---------------------------------------------- */
        pivot = x [ipiv] ;          /* the chosen pivot */
        Ui [unz] = k ;              /* last entry in U(:,k) is U(k,k) */
        Ux [unz++] = pivot ;
        pinv [ipiv] = k ;           /* ipiv is the kth pivot row */
        Li [lnz] = ipiv ;           /* first entry in L(:,k) is L(k,k) = 1 */
        Lx [lnz++] = 1 ;
        for (p = top ; p < n ; p++) /* L(k+1:n,k) = x / pivot */
        {
            i = xi [p] ;
            if (pinv [i] < 0)       /* x(i) is an entry in L(:,k) */
            {
                Li [lnz] = i ;      /* save unpermuted row in L */
                Lx [lnz++] = x [i] / pivot ;    /* scale pivot column */
            }
            x [i] = 0 ;             /* x [0..n-1] = 0 for next k */
        }
    }
    /* --- Finalize L and U ------------------------------------------------- */
    Lp [n] = lnz ;
    Up [n] = unz ;
    Li = L->i ;                     /* fix row indices of L for final pinv */
    for (p = 0 ; p < lnz ; p++) Li [p] = pinv [Li [p]] ;
    cs_sprealloc (L, 0) ;           /* remove extra space from L and U */
    cs_sprealloc (U, 0) ;
    return (cs_ndone (N, NULL, xi, x, 1)) ;     /* success */
}
