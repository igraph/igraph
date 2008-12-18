/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2008  Gabor Csardi <csardi@rmki.kfki.hu>
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

#ifndef IGRAPH_COMPLEX_H
#define IGRAPH_COMPLEX_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#include "types.h"

typedef struct igraph_complex_t {
  igraph_real_t real, imag;
} igraph_complex_t;

#define REALPART(x) ((x).real)
#define IMAGPART(x) ((x).imag)

igraph_real_t igraph_complex_real(igraph_complex_t cmplx);
igraph_real_t igraph_complex_imag(igraph_complex_t cmplx);
igraph_real_t igraph_complex_mod(igraph_complex_t cmplx);
igraph_real_t igraph_complex_arg(igraph_complex_t cmplx);
igraph_complex_t igraph_complex_conj(igraph_complex_t cmplx);
igraph_complex_t igraph_complex_mul(igraph_complex_t lhs, igraph_complex_t rhs);
igraph_complex_t igraph_complex_div(igraph_complex_t lhs, igraph_complex_t rhs);

typedef struct igraph_vector_complex_t {
  igraph_complex_t *stor_begin;
  igraph_complex_t *stor_end;
  igraph_complex_t *end;
} igraph_vector_complex_t;

/* Note the VECTOR() still works for this!!! */

int igraph_vector_complex_real(const igraph_vector_complex_t *cv, igraph_vector_t *res);
int igraph_vector_complex_imag(const igraph_vector_complex_t *cv, igraph_vector_t *res);
int igraph_vector_complex_mod(const igraph_vector_complex_t *cv, igraph_vector_t *res);
int igraph_vector_complex_arg(const igraph_vector_complex_t *cv, igraph_vector_t *res);
int igraph_vector_complex_conj(const igraph_vector_complex_t *cv, igraph_vector_complex_t *res);

int igraph_vector_complex_init(igraph_vector_complex_t* v, long int size);
int igraph_vector_complex_copy(igraph_vector_complex_t *to, 
			       const igraph_vector_complex_t *from);
void igraph_vector_complex_destroy(igraph_vector_complex_t* v);

int igraph_vector_complex_init_realimag(igraph_vector_complex_t *v,
					const igraph_vector_t *real,
					const igraph_vector_t *imag);

igraph_complex_t igraph_vector_complex_e(const igraph_vector_complex_t* v, long int pos);
igraph_complex_t* igraph_vector_complex_e_ptr(const igraph_vector_complex_t* v, long int pos);
void igraph_vector_complex_set(igraph_vector_complex_t* v, long int pos, igraph_complex_t value);
igraph_complex_t igraph_vector_complex_tail(const igraph_vector_complex_t *v);

void igraph_vector_complex_null(igraph_vector_complex_t* v);
void igraph_vector_complex_fill(igraph_vector_complex_t* v, igraph_complex_t e);

const igraph_vector_complex_t *igraph_vector_complex_view(const igraph_vector_complex_t *v,
							  const igraph_complex_t *data, 
							  long int length);
void igraph_vector_complex_copy_to(const igraph_vector_complex_t *v, igraph_complex_t* to);
int igraph_vector_complex_update(igraph_vector_complex_t *to, 
				 const igraph_vector_complex_t *from);
int igraph_vector_complex_append(igraph_vector_complex_t *to, 
				 const igraph_vector_complex_t *from);
int igraph_vector_complex_swap(igraph_vector_complex_t *v1, igraph_vector_complex_t *v2);

int igraph_vector_complex_swap_elements(igraph_vector_complex_t *v,
					long int i, long int j);
int igraph_vector_complex_reverse(igraph_vector_complex_t *v);

void igraph_vector_complex_add_constant(igraph_vector_complex_t *v, igraph_complex_t plus);
void igraph_vector_complex_scale(igraph_vector_complex_t *v, igraph_complex_t by);
int igraph_vector_complex_add(igraph_vector_complex_t *v1, 
				const igraph_vector_complex_t *v2);
int igraph_vector_complex_sub(igraph_vector_complex_t *v1, 
				const igraph_vector_complex_t *v2);
int igraph_vector_complex_mul(igraph_vector_complex_t *v1, 
				const igraph_vector_complex_t *v2);
int igraph_vector_complex_div(igraph_vector_complex_t *v1, 
				const igraph_vector_complex_t *v2);

igraph_bool_t igraph_vector_complex_empty (const igraph_vector_complex_t* v);
long int igraph_vector_complex_size(const igraph_vector_complex_t* v);
igraph_bool_t igraph_vector_complex_isnull(const igraph_vector_complex_t *v);
igraph_complex_t igraph_vector_complex_sum(const igraph_vector_complex_t *v);
igraph_complex_t igraph_vector_complex_prod(const igraph_vector_complex_t *v);
igraph_bool_t igraph_vector_complex_is_equal(const igraph_vector_complex_t *lhs, 
					     const igraph_vector_complex_t *rhs);

igraph_bool_t igraph_vector_complex_contains(const igraph_vector_complex_t *v, igraph_complex_t e);
igraph_bool_t igraph_vector_complex_search(const igraph_vector_complex_t *v,
					     long int from, igraph_complex_t what, 
					     long int *pos);

void igraph_vector_complex_clear(igraph_vector_complex_t* v);
int igraph_vector_complex_resize(igraph_vector_complex_t* v, long int newsize);
int igraph_vector_complex_reserve(igraph_vector_complex_t* v, long int size);
int igraph_vector_complex_push_back(igraph_vector_complex_t* v, igraph_complex_t e);
igraph_complex_t igraph_vector_complex_pop_back(igraph_vector_complex_t* v);
int igraph_vector_complex_insert(igraph_vector_complex_t *v, long int pos, igraph_complex_t value);
void igraph_vector_complex_remove(igraph_vector_complex_t *v, long int elem);
void igraph_vector_complex_remove_section(igraph_vector_complex_t *v, 
					    long int from, long int to);


__END_DECLS

#endif
