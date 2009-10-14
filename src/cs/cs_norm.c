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
/* 1-norm of a sparse matrix = max (sum (abs (A))), largest column sum */
double cs_norm (const cs *A)
{
    CS_INT p, j, n, *Ap ;
    CS_ENTRY *Ax ;
    double norm = 0, s ;
    if (!CS_CSC (A) || !A->x) return (-1) ;             /* check inputs */
    n = A->n ; Ap = A->p ; Ax = A->x ;
    for (j = 0 ; j < n ; j++)
    {
        for (s = 0, p = Ap [j] ; p < Ap [j+1] ; p++) s += CS_ABS (Ax [p]) ;
        norm = CS_MAX (norm, s) ;
    }
    return (norm) ;
}
