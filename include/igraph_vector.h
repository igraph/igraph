/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_VECTOR_H
#define IGRAPH_VECTOR_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_complex.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Flexible vector                                    */
/* -------------------------------------------------- */

#define BASE_IGRAPH_REAL
#include "igraph_pmt.h"
#include "igraph_vector_type.h"
#include "igraph_pmt_off.h"
#undef BASE_IGRAPH_REAL

#define BASE_FLOAT
#include "igraph_pmt.h"
#include "igraph_vector_type.h"
#include "igraph_pmt_off.h"
#undef BASE_FLOAT

#define BASE_LONG
#include "igraph_pmt.h"
#include "igraph_vector_type.h"
#include "igraph_pmt_off.h"
#undef BASE_LONG

#define BASE_CHAR
#include "igraph_pmt.h"
#include "igraph_vector_type.h"
#include "igraph_pmt_off.h"
#undef BASE_CHAR

#define BASE_BOOL
#include "igraph_pmt.h"
#include "igraph_vector_type.h"
#include "igraph_pmt_off.h"
#undef BASE_BOOL

#define BASE_INT
#include "igraph_pmt.h"
#include "igraph_vector_type.h"
#include "igraph_pmt_off.h"
#undef BASE_INT

#define BASE_COMPLEX
#include "igraph_pmt.h"
#include "igraph_vector_type.h"
#include "igraph_pmt_off.h"
#undef BASE_COMPLEX

#define BASE_IGRAPH_REAL
#include "igraph_pmt.h"
#include "igraph_vector_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_IGRAPH_REAL

#define BASE_FLOAT
#include "igraph_pmt.h"
#include "igraph_vector_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_FLOAT

#define BASE_LONG
#include "igraph_pmt.h"
#include "igraph_vector_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_LONG

#define BASE_CHAR
#include "igraph_pmt.h"
#include "igraph_vector_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_CHAR

#define BASE_BOOL
#include "igraph_pmt.h"
#include "igraph_vector_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_BOOL

#define BASE_INT
#include "igraph_pmt.h"
#include "igraph_vector_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_INT

#define BASE_COMPLEX
#include "igraph_pmt.h"
#include "igraph_vector_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_COMPLEX

/* -------------------------------------------------- */
/* Helper macros                                      */
/* -------------------------------------------------- */

#ifndef IGRAPH_VECTOR_NULL
    #define IGRAPH_VECTOR_NULL { 0,0,0 }
#endif

#ifndef IGRAPH_VECTOR_INIT_FINALLY
#define IGRAPH_VECTOR_INIT_FINALLY(v, size) \
    do { IGRAPH_CHECK(igraph_vector_init(v, size)); \
        IGRAPH_FINALLY(igraph_vector_destroy, v); } while (0)
#endif
#ifndef IGRAPH_VECTOR_BOOL_INIT_FINALLY
#define IGRAPH_VECTOR_BOOL_INIT_FINALLY(v, size) \
    do { IGRAPH_CHECK(igraph_vector_bool_init(v, size)); \
        IGRAPH_FINALLY(igraph_vector_bool_destroy, v); } while (0)
#endif
#ifndef IGRAPH_VECTOR_CHAR_INIT_FINALLY
#define IGRAPH_VECTOR_CHAR_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(igraph_vector_char_init(v, size)); \
  IGRAPH_FINALLY(igraph_vector_char_destroy, v); } while (0)
#endif
#ifndef IGRAPH_VECTOR_INT_INIT_FINALLY
#define IGRAPH_VECTOR_INT_INIT_FINALLY(v, size) \
    do { IGRAPH_CHECK(igraph_vector_int_init(v, size)); \
        IGRAPH_FINALLY(igraph_vector_int_destroy, v); } while (0)
#endif
#ifndef IGRAPH_VECTOR_LONG_INIT_FINALLY
#define IGRAPH_VECTOR_LONG_INIT_FINALLY(v, size) \
    do { IGRAPH_CHECK(igraph_vector_long_init(v, size)); \
        IGRAPH_FINALLY(igraph_vector_long_destroy, v); } while (0)
#endif

/* -------------------------------------------------- */
/* Type-specific vector functions                     */
/* -------------------------------------------------- */

IGRAPH_EXPORT int igraph_vector_floor(const igraph_vector_t *from, igraph_vector_long_t *to);
IGRAPH_EXPORT int igraph_vector_round(const igraph_vector_t *from, igraph_vector_long_t *to);

IGRAPH_EXPORT igraph_bool_t igraph_vector_e_tol(const igraph_vector_t *lhs,
                                                const igraph_vector_t *rhs,
                                                igraph_real_t tol);

IGRAPH_EXPORT int igraph_vector_zapsmall(igraph_vector_t *v, igraph_real_t tol);

IGRAPH_EXPORT int igraph_vector_order(const igraph_vector_t* v, const igraph_vector_t *v2,
                                      igraph_vector_t* res, igraph_real_t maxval);
IGRAPH_EXPORT int igraph_vector_order1(const igraph_vector_t* v,
                                       igraph_vector_t* res, igraph_real_t maxval);
IGRAPH_EXPORT int igraph_vector_order1_int(const igraph_vector_t* v,
                                           igraph_vector_int_t* res, igraph_real_t maxval);
IGRAPH_EXPORT int igraph_vector_order2(igraph_vector_t *v);
IGRAPH_EXPORT int igraph_vector_rank(const igraph_vector_t *v, igraph_vector_t *res,
                                     long int nodes);
IGRAPH_EXPORT int igraph_vector_is_nan(const igraph_vector_t *v,
                                       igraph_vector_bool_t *is_nan);
IGRAPH_EXPORT igraph_bool_t igraph_vector_is_any_nan(const igraph_vector_t *v);

__END_DECLS

#endif
