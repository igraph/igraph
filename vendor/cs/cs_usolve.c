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
/* solve Ux=b where x and b are dense.  x=b on input, solution on output. */
CS_INT cs_usolve (const cs *U, CS_ENTRY *x)
{
    CS_INT p, j, n, *Up, *Ui ;
    CS_ENTRY *Ux ;
    if (!CS_CSC (U) || !x) return (0) ;                     /* check inputs */
    n = U->n ; Up = U->p ; Ui = U->i ; Ux = U->x ;
    for (j = n-1 ; j >= 0 ; j--)
    {
        x [j] /= Ux [Up [j+1]-1] ;
        for (p = Up [j] ; p < Up [j+1]-1 ; p++)
        {
            x [Ui [p]] -= Ux [p] * x [j] ;
        }
    }
    return (1) ;
}
