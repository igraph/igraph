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

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

#include "igraph_types.h"
#include "igraph_complex.h"

#ifdef HAVE_STDINT_H
#  include <stdint.h>
#else
#  if HAVE_SYS_INT_TYPES_H
#    include <sys/int_types.h>    /* for Solaris */
#  endif
#endif

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

#define BASE_TIME
#include "igraph_pmt.h"
#include "igraph_vector_type.h"
#include "igraph_pmt_off.h"
#undef BASE_TIME

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

#define BASE_TIME
#include "igraph_pmt.h"
#include "igraph_vector_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_TIME

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
#ifndef IGRAPH_VECTOR_LONG_INIT_FINALLY
#define IGRAPH_VECTOR_LONG_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(igraph_vector_long_init(v, size)); \
  IGRAPH_FINALLY(igraph_vector_long_destroy, v); } while (0)
#endif
#ifndef IGRAPH_VECTOR_INT_INIT_FINALLY
#define IGRAPH_VECTOR_INT_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(igraph_vector_int_init(v, size));		\
  IGRAPH_FINALLY(igraph_vector_int_destroy, v); } while (0)
#endif
#ifndef IGRAPH_VECTOR_TIME_INIT_FINALLY
#define IGRAPH_VECTOR_TIME_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(igraph_vector_time_init(v, size));		\
  IGRAPH_FINALLY(igraph_vector_time_destroy, v); } while (0)
#endif

#define IGRAPH_VECTOR_CONSTANT(name, ...)				\
  igraph_real_t IGRAPH_UNIQUE(name,1)[] = { __VA_ARGS__ };		\
  igraph_vector_t (name);						\
  const igraph_vector_t *IGRAPH_UNIQUE(name,2) =			\
    igraph_vector_view(&(name), IGRAPH_UNIQUE(name,1),			\
		       sizeof(IGRAPH_UNIQUE(name,1))/			\
		       sizeof(igraph_real_t))

#define IGRAPH_VECTOR_FLOAT_CONSTANT(name, ...)				\
  float IGRAPH_UNIQUE(name,1)[] = { __VA_ARGS__ };			\
  igraph_vector_float_t (name);						\
  const igraph_vector_float_t *IGRAPH_UNIQUE(name,2) =			\
    igraph_vector_float_view(&(name), IGRAPH_UNIQUE(name,1),		\
			     sizeof(IGRAPH_UNIQUE(name,1))/sizeof(float))

#define IGRAPH_VECTOR_LONG_CONSTANT(name, ...)				\
  long IGRAPH_UNIQUE(name,1)[] = { __VA_ARGS__ };			\
  igraph_vector_long_t (name);						\
  const igraph_vector_long_t *IGRAPH_UNIQUE(name,2) =			\
    igraph_vector_long_view(&(name), IGRAPH_UNIQUE(name,1),		\
			    sizeof(IGRAPH_UNIQUE(name,1))/sizeof(long))

#define IGRAPH_VECTOR_CHAR_CONSTANT(name, ...)				\
  char IGRAPH_UNIQUE(name,1)[] = { __VA_ARGS__ };			\
  igraph_vector_char_t (name);						\
  const igraph_vector_char_t *IGRAPH_UNIQUE(name,2) =			\
    igraph_vector_char_view(&(name), IGRAPH_UNIQUE(name,1),		\
			    sizeof(IGRAPH_UNIQUE(name,1))/sizeof(char))

#define IGRAPH_VECTOR_BOOL_CONSTANT(name, ...)				\
  igraph_bool_t IGRAPH_UNIQUE(name,1)[] = { __VA_ARGS__ };		\
  igraph_vector_bool_t (name);						\
  const igraph_vector_bool_t *IGRAPH_UNIQUE(name,2) =			\
    igraph_vector_bool_view(&(name), IGRAPH_UNIQUE(name,1),		\
			    sizeof(IGRAPH_UNIQUE(name,1))/		\
			    sizeof(igraph_bool_t))

#define IGRAPH_VECTOR_INT_CONSTANT(name, ...)				\
  int IGRAPH_UNIQUE(name,1)[] = { __VA_ARGS__ };			\
  igraph_vector_int_t (name);						\
  const igraph_vector_int_t *IGRAPH_UNIQUE(name,2) =			\
    igraph_vector_int_view(&(name), IGRAPH_UNIQUE(name,1),		\
			   sizeof(IGRAPH_UNIQUE(name,1))/sizeof(int))

#define IGRAPH_VECTOR_TIME_CONSTANT(name, ...)				\
  igraph_time_t IGRAPH_UNIQUE(name,1)[] = { __VA_ARGS__ };		\
  igraph_vector_time_t (name);						\
  const igraph_vector_time_t *IGRAPH_UNIQUE(name,2) =			\
    igraph_vector_time_view(&(name), IGRAPH_UNIQUE(name,1),		\
			    sizeof(IGRAPH_UNIQUE(name,1))/		\
			    sizeof(igraph_time_t))

/* -------------------------------------------------- */
/* Type-specific vector functions                     */
/* -------------------------------------------------- */

int igraph_vector_floor(const igraph_vector_t *from, igraph_vector_long_t *to);
int igraph_vector_round(const igraph_vector_t *from, igraph_vector_long_t *to);

igraph_bool_t igraph_vector_e_tol(const igraph_vector_t *lhs,
				  const igraph_vector_t *rhs,
				  igraph_real_t tol);

/* These are for internal use only */
int igraph_vector_order(const igraph_vector_t* v, const igraph_vector_t *v2,
			igraph_vector_t* res, igraph_real_t maxval);
int igraph_vector_order3(const igraph_vector_t* v,
			 const igraph_vector_time_t *vt,
			 const igraph_vector_t *v2,
			 igraph_vector_t* res, igraph_real_t nodes,
			 igraph_time_t time_steps);
int igraph_vector_order1(const igraph_vector_t* v, 
			 igraph_vector_t* res, igraph_real_t maxval);
int igraph_vector_order1_int(const igraph_vector_t* v,
			 igraph_vector_int_t* res, igraph_real_t maxval);
int igraph_vector_order2(igraph_vector_t *v);
int igraph_vector_rank(const igraph_vector_t *v, igraph_vector_t *res, 
		       long int nodes);

__END_DECLS

#endif
