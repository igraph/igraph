/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#ifndef REST_MEMORY_H
#define REST_MEMORY_H

#include <stdlib.h>

#define Calloc(n,t)    (t*) calloc( (size_t)(n), sizeof(t) )
#define Realloc(p,n,t) (t*) realloc((void*)(p), (size_t)((n)*sizeof(t)))
#define Free(p)        (free( (void *)(p) ), (p) = NULL)

#endif
