/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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
    #define OUT_FORMAT "%G"
    #define PRINTFUNC(val) igraph_real_printf(val)
    #define FPRINTFUNC(file, val) igraph_real_fprintf(file, val)
    #define ZERO 0.0
    #define ONE 1.0
    #define MULTIPLICITY 1

#elif defined(BASE_FLOAT)
    #define BASE float
    #define SHORT float
    #define OUT_FORMAT "%f"
    #define ZERO 0.0F
    #define ONE 1.0F
    #define MULTIPLICITY 1

#elif defined(BASE_LONG)
    #define BASE long
    #define SHORT long
    #define OUT_FORMAT "%ld"
    #define ZERO 0L
    #define ONE 1L
    #define MULTIPLICITY 1

#elif defined(BASE_CHAR)
    #define BASE char
    #define SHORT char
    #define OUT_FORMAT "%d"
    #define ZERO 0
    #define ONE 1
    #define MULTIPLICITY 1

#elif defined(BASE_BOOL)
    #define BASE igraph_bool_t
    #define SHORT bool
    #define OUT_FORMAT "%d"
    #define ZERO 0
    #define ONE 1
    #define MULTIPLICITY 1

#elif defined(BASE_INT)
    #define BASE int
    #define SHORT int
    #define OUT_FORMAT "%d"
    #define ZERO 0
    #define ONE 1
    #define MULTIPLICITY 1

#elif defined(BASE_LIMB)
    #define BASE limb_t
    #define SHORT limb
    #define ZERO 0
    #define ONE 1
    #define MULTIPLICITY 1
    #define UNSIGNED 1

#elif defined(BASE_PTR)
    #define BASE void*
    #define SHORT ptr
    #define ZERO 0
    #define MULTIPLICITY 1

#elif defined(BASE_COMPLEX)
    #undef complex
    #define BASE igraph_complex_t
    #define SHORT complex
    #define ZERO igraph_complex(0,0)
    #define ONE {{1.0,0.0}}
    #define MULTIPLICITY 2
    #define NOTORDERED 1
    #define NOABS 1
    #define SUM(a,b,c) ((a) = igraph_complex_add((b),(c)))
    #define DIFF(a,b,c) ((a) = igraph_complex_sub((b),(c)))
    #define PROD(a,b,c) ((a) = igraph_complex_mul((b),(c)))
    #define DIV(a,b,c) ((a) = igraph_complex_div((b),(c)))
    #define EQ(a,b) IGRAPH_COMPLEX_EQ((a),(b))
    #define SQ(a) IGRAPH_REAL(igraph_complex_mul((a),(a)))

#else
    #error unknown BASE_ directive
#endif

#if defined(BASE_IGRAPH_REAL)
    #define FUNCTION(dir,name) CONCAT2(dir,name)
    #define TYPE(dir) CONCAT2(dir,t)
#elif defined(BASE_BOOL)
    /* Special case because stdbool.h defines bool as a macro to _Bool which would
    * screw things up */
    #define FUNCTION(a,c) CONCAT3x(a,bool,c)
    #define TYPE(dir) CONCAT3x(dir,bool,t)
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

