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

#include "igraph_random.h"

#include "cs.h"
/* return a random permutation vector, the identity perm, or p = n-1:-1:0.
 * seed = -1 means p = n-1:-1:0.  seed = 0 means p = identity.  otherwise
 * p = random permutation.  */
CS_INT *cs_randperm (CS_INT n, CS_INT seed)
{
    CS_INT *p, k, j, t ;
    if (seed == 0) return (NULL) ;      /* return p = NULL (identity) */
    p = cs_malloc (n, sizeof (CS_INT)) ;   /* allocate result */
    if (!p) return (NULL) ;             /* out of memory */
    for (k = 0 ; k < n ; k++) p [k] = n-k-1 ;
    if (seed == -1) return (p) ;        /* return reverse permutation */
    /* srand (seed) ;                      /\* get new random number seed *\/ */
    RNG_BEGIN();
    for (k = 0 ; k < n ; k++)
    {
        /* j = k + (rand ( ) % (n-k)) ;    /\* j = rand CS_INT in range k to n-1 *\/ */
      j = k + RNG_INTEGER(k, n-1) ;
        t = p [j] ;                     /* swap p[k] and p[j] */
        p [j] = p [k] ;
        p [k] = t ;
    }
    RNG_END();
    return (p) ;
}
