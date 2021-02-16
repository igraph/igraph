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
/* load a triplet matrix from a file */
cs *cs_load (FILE *f)
{
    CS_INT i, j ;
    double x ;
#ifdef CS_COMPLEX
    double xi ;
#endif
    cs *T ;
    if (!f) return (NULL) ;                             /* check inputs */
    T = cs_spalloc (0, 0, 1, 1, 1) ;                    /* allocate result */
#ifdef CS_COMPLEX
    while (fscanf (f, ""CS_ID" "CS_ID" %lg %lg\n", &i, &j, &x, &xi) == 4)
#else
    while (fscanf (f, ""CS_ID" "CS_ID" %lg\n", &i, &j, &x) == 3)
#endif
    {
#ifdef CS_COMPLEX
        if (!cs_entry (T, i, j, x + xi*I)) return (cs_spfree (T)) ;
#else
        if (!cs_entry (T, i, j, x)) return (cs_spfree (T)) ;
#endif
    }
    return (T) ;
}
