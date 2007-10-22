/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "config.h"
#include "error.h"

/* This is from GNU R's print.c */

/* Fortran-callable error routine for lapack */

void F77_FUNC(igraphxerbla,IGRAPHXERBLA)(const char *srname, int *info)
{
   /* srname is not null-terminated.  It should be 6 characters. */
    char buf[7], buf2[200];
    strncpy(buf, srname, 6);
    buf[6] = '\0';
    snprintf(buf2, sizeof(buf2)/sizeof(char)-1, 
	     "BLAS/LAPACK routine '%6s' gave error code %d", buf, -*info);
    igraph_error(buf2, __FILE__, __LINE__, IGRAPH_FAILURE);
}
