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

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

#if defined(BASE_IGRAPH_REAL)
#define BASE igraph_real_t
#define SHORT
#define ATOMIC
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1

#elif defined(BASE_BOOL)
#define BASE igraph_bool_t
#define SHORT bool
#define ATOMIC igraph_bool_t
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1

#else
#error unknown BASE_ directive
#endif

#if defined(BASE_IGRAPH_REAL)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) CONCAT2(dir,t)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT3(dir,SHORT,t)
#endif

#if defined(HEAP_TYPE_MIN)
#define HEAPMORE <
#define HEAPMOREEQ <=
#define HEAPLESS >
#define HEAPLESSEQ >=
#undef FUNCTION
#undef TYPE
#if defined(BASE_IGRAPH_REAL)
#define FUNCTION(dir,name) CONCAT3(dir,min,name)
#define TYPE(dir) CONCAT3(dir,min,t)
#else
#define FUNCTION(a,c) CONCAT4(a,min,SHORT,c)
#define TYPE(dir) CONCAT4(dir,min,SHORT,t)
#endif
#endif

#if defined(HEAP_TYPE_MAX)
#define HEAPMORE >
#define HEAPMOREEQ >=
#define HEAPLESS <
#define HEAPLESSEQ <=
#endif

