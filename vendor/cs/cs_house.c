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
/* create a Householder reflection [v,beta,s]=house(x), overwrite x with v,
 * where (I-beta*v*v')*x = s*e1 and e1 = [1 0 ... 0]'.
 * Note that this CXSparse version is different than CSparse.  See Higham,
 * Accuracy & Stability of Num Algorithms, 2nd ed, 2002, page 357. */
CS_ENTRY cs_house (CS_ENTRY *x, double *beta, CS_INT n)
{
    CS_ENTRY s = 0 ;
    CS_INT i ;
    if (!x || !beta) return (-1) ;          /* check inputs */
    /* s = norm(x) */
    for (i = 0 ; i < n ; i++) s += x [i] * CS_CONJ (x [i]) ;
    s = sqrt (s) ;
    if (s == 0)
    {
        (*beta) = 0 ;
        x [0] = 1 ;
    }
    else
    {
        /* s = sign(x[0]) * norm (x) ; */
        if (x [0] != 0)
        {
            s *= x [0] / CS_ABS (x [0]) ;
        }
        x [0] += s ;
        (*beta) = 1. / CS_REAL (CS_CONJ (s) * x [0]) ;
    }
    return (-s) ;
}
