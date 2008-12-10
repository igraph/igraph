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

#include "igraph_complex.h"
#include "memory.h"
#include <string.h>
#include <math.h>

igraph_real_t igraph_complex_real(igraph_complex_t cmplx) {
  return REALPART(cmplx);
}

igraph_real_t igraph_complex_imag(igraph_complex_t cmplx) {
  return IMAGPART(cmplx);
}

igraph_real_t igraph_complex_mod(igraph_complex_t cmplx) {
  return sqrt(REALPART(cmplx)*REALPART(cmplx) + 
	      IMAGPART(cmplx)*IMAGPART(cmplx));
}

igraph_real_t igraph_complex_arg(igraph_complex_t cmplx) {
  return atan2(IMAGPART(cmplx), REALPART(cmplx));
}

igraph_complex_t igraph_complex_conj(igraph_complex_t cmplx) {
  igraph_complex_t ret;
  REALPART(ret)=REALPART(cmplx);
  IMAGPART(ret)=-IMAGPART(cmplx);
  return ret;
}

igraph_complex_t igraph_complex_mul(igraph_complex_t lhs, igraph_complex_t rhs) {
  igraph_complex_t res;
  REALPART(res) = REALPART(lhs) * REALPART(rhs) - IMAGPART(lhs) * IMAGPART(rhs);
  IMAGPART(res) = REALPART(lhs) * IMAGPART(rhs) + IMAGPART(lhs) * REALPART(rhs);
  return res;
}

igraph_complex_t igraph_complex_div(igraph_complex_t lhs, igraph_complex_t rhs) {
  igraph_real_t ratio, den;
  igraph_real_t abr, abi;
  igraph_complex_t ret;

  if( (abr = REALPART(rhs)) < 0)
    abr = - abr;
  if( (abi = IMAGPART(rhs)) < 0)
    abi = - abi;
  if( abr <= abi ) {
    ratio = REALPART(rhs) / IMAGPART(rhs) ;
    den = IMAGPART(rhs) * (1 + ratio*ratio);
    REALPART(ret) = (REALPART(lhs)*ratio + IMAGPART(lhs)) / den;
    IMAGPART(ret) = (IMAGPART(lhs)*ratio - REALPART(lhs)) / den;
  } else {
    ratio = IMAGPART(rhs) / REALPART(rhs) ;
    den = REALPART(rhs) * (1 + ratio*ratio);
    REALPART(ret) = (REALPART(lhs) + IMAGPART(lhs)*ratio) / den;
    IMAGPART(ret) = (IMAGPART(lhs) - REALPART(lhs)*ratio) / den;
  }
  return ret;
}

int igraph_vector_complex_init(igraph_vector_complex_t* v, long int size) {
  long int alloc_size= size > 0 ? size : 1;
  if (size < 0) { size=0; }
  v->stor_begin=igraph_Calloc(alloc_size, igraph_complex_t);
  if (v->stor_begin==0) {
    IGRAPH_ERROR("cannot init vector", IGRAPH_ENOMEM);
  }
  v->stor_end=v->stor_begin + alloc_size;
  v->end=v->stor_begin+size;
  
  return 0;
}

int igraph_vector_complex_copy(igraph_vector_complex_t *to, 
			       const igraph_vector_complex_t *from) {
  to->stor_begin=igraph_Calloc(igraph_vector_complex_size(from), igraph_complex_t);
  if (to->stor_begin==0) {
    IGRAPH_ERROR("canot copy vector", IGRAPH_ENOMEM);
  }
  to->stor_end=to->stor_begin+igraph_vector_complex_size(from);
  to->end=to->stor_end;
  memcpy(to->stor_begin, from->stor_begin, 
	 igraph_vector_complex_size(from)*sizeof(igraph_complex_t));
  return 0;
}

void igraph_vector_complex_destroy(igraph_vector_complex_t* v) {
  if (v->stor_begin != 0) {
    igraph_Free(v->stor_begin);
    v->stor_begin = NULL;
  }
}


igraph_complex_t igraph_vector_complex_e(const igraph_vector_complex_t* v, long int pos) {
  return * (v->stor_begin + pos);
}


igraph_complex_t* igraph_vector_complex_e_ptr(const igraph_vector_complex_t* v, long int pos) {
  return v->stor_begin+pos;
}


void igraph_vector_complex_set(igraph_vector_complex_t* v, long int pos, igraph_complex_t value) {
  *(v->stor_begin + pos) = value;
}


igraph_complex_t igraph_vector_complex_tail(const igraph_vector_complex_t *v) {
  return *((v->end)-1);  
}



void igraph_vector_complex_null(igraph_vector_complex_t* v) {
  if (igraph_vector_complex_size(v)>0) {
    memset(v->stor_begin, 0, 
	   sizeof(igraph_complex_t)*igraph_vector_complex_size(v));
  }
}


void igraph_vector_complex_fill(igraph_vector_complex_t* v, igraph_complex_t e) {
  igraph_complex_t *ptr;
  for (ptr = v->stor_begin; ptr < v->end; ptr++) {
    *ptr = e;
  }  
}



const igraph_vector_complex_t *igraph_vector_complex_view(const igraph_vector_complex_t *v,
							  const igraph_complex_t *data, 
							  long int length) {
  igraph_vector_complex_t *v2=(igraph_vector_complex_t*)v;
  v2->stor_begin=(igraph_complex_t*)data;
  v2->stor_end=(igraph_complex_t*)data+length;
  v2->end=v->stor_end;
  return v;
}


void igraph_vector_complex_copy_to(const igraph_vector_complex_t *v, igraph_complex_t* to) {
  if (v->end != v->stor_begin) {
    memcpy(to, v->stor_begin, sizeof(igraph_complex_t) * (v->end - v->stor_begin));
  }
}


int igraph_vector_complex_update(igraph_vector_complex_t *to, 
				 const igraph_vector_complex_t *from) {
  long int n=igraph_vector_complex_size(from);
  igraph_vector_complex_resize(to, n);
  memcpy(to->stor_begin, from->stor_begin, sizeof(igraph_complex_t)*n);
  return 0;
}


int igraph_vector_complex_append(igraph_vector_complex_t *to, 
				 const igraph_vector_complex_t *from) {
  long int tosize, fromsize;
  
  tosize=igraph_vector_complex_size(to);
  fromsize=igraph_vector_complex_size(from);
  IGRAPH_CHECK(igraph_vector_complex_resize(to, tosize+fromsize));
  memcpy(to->stor_begin+tosize, from->stor_begin, 
	 sizeof(igraph_real_t)*fromsize);
  to->end=to->stor_begin+tosize+fromsize;
  
  return 0;
}


int igraph_vector_complex_swap(igraph_vector_complex_t *v1, igraph_vector_complex_t *v2) {
  long int i, n1=igraph_vector_complex_size(v1);
  long int n2=igraph_vector_complex_size(v2);
  if (n1 != n2) {
    IGRAPH_ERROR("Vectors must have the same number of elements for swapping",
		 IGRAPH_EINVAL);
  }
  
  for (i=0; i<n1; i++) {
    igraph_complex_t tmp;
    tmp=VECTOR(*v1)[i];
    VECTOR(*v1)[i]=VECTOR(*v2)[i];
    VECTOR(*v2)[i]=tmp;
  }
  return 0;
}


int igraph_vector_complex_swap_elements(igraph_vector_complex_t *v,
					long int i, long int j) {
  igraph_complex_t tmp=VECTOR(*v)[i];
  VECTOR(*v)[i]=VECTOR(*v)[j];
  VECTOR(*v)[j]=tmp;
  return 0;
}


int igraph_vector_complex_reverse(igraph_vector_complex_t *v) {
  long int n=igraph_vector_complex_size(v), n2=n/2;
  long int i,j;
  for (i=0, j=n-1; i<n2; i++, j--) {
    igraph_complex_t tmp;
    tmp=VECTOR(*v)[i];
    VECTOR(*v)[i]=VECTOR(*v)[j];
    VECTOR(*v)[j]=tmp;
  }
  return 0;
}



void igraph_vector_complex_add_constant(igraph_vector_complex_t *v, igraph_complex_t plus) {
  long int i, n=igraph_vector_complex_size(v);
  for (i=0; i<n; i++) {
    REALPART(VECTOR(*v)[i]) += REALPART(plus);
    IMAGPART(VECTOR(*v)[i]) += IMAGPART(plus);
  }
}


void igraph_vector_complex_scale(igraph_vector_complex_t *v, igraph_complex_t by) {
  long int i;
  for (i=0; i<igraph_vector_complex_size(v); i++) {
    VECTOR(*v)[i] = igraph_complex_mul(VECTOR(*v)[i], by);
  }
}


int igraph_vector_complex_add(igraph_vector_complex_t *v1, 
				const igraph_vector_complex_t *v2) {
  long int n1=igraph_vector_complex_size(v1);
  long int n2=igraph_vector_complex_size(v2);
  long int i;
  if (n1 != n2) {
    IGRAPH_ERROR("Vectors must have the same number of elements for swapping",
		 IGRAPH_EINVAL);
  }
  
  for (i=0; i<n1; i++) {
    REALPART(VECTOR(*v1)[i]) += REALPART(VECTOR(*v2)[i]);
    IMAGPART(VECTOR(*v1)[i]) += IMAGPART(VECTOR(*v2)[i]);
  }
  
  return 0;
}


int igraph_vector_complex_sub(igraph_vector_complex_t *v1, 
				const igraph_vector_complex_t *v2) {
  long int n1=igraph_vector_complex_size(v1);
  long int n2=igraph_vector_complex_size(v2);
  long int i;
  if (n1 != n2) {
    IGRAPH_ERROR("Vectors must have the same number of elements for swapping",
		 IGRAPH_EINVAL);
  }
  
  for (i=0; i<n1; i++) {
    REALPART(VECTOR(*v1)[i]) -= REALPART(VECTOR(*v2)[i]);
    IMAGPART(VECTOR(*v1)[i]) -= IMAGPART(VECTOR(*v2)[i]);
  }
  
  return 0;
}


int igraph_vector_complex_mul(igraph_vector_complex_t *v1, 
				const igraph_vector_complex_t *v2) {
  long int n1=igraph_vector_complex_size(v1);
  long int n2=igraph_vector_complex_size(v2);
  long int i;
  if (n1 != n2) {
    IGRAPH_ERROR("Vectors must have the same number of elements for swapping",
		 IGRAPH_EINVAL);
  }
  
  for (i=0; i<n1; i++) {
    VECTOR(*v1)[i] = igraph_complex_mul(VECTOR(*v1)[i], VECTOR(*v2)[i]);
  }
  
  return 0;
}


int igraph_vector_complex_div(igraph_vector_complex_t *v1, 
				const igraph_vector_complex_t *v2) {
  long int n1=igraph_vector_complex_size(v1);
  long int n2=igraph_vector_complex_size(v2);
  long int i;
  if (n1 != n2) {
    IGRAPH_ERROR("Vectors must have the same number of elements for swapping",
		 IGRAPH_EINVAL);
  }
  
  for (i=0; i<n1; i++) {
    VECTOR(*v1)[i] = igraph_complex_div(VECTOR(*v1)[i], VECTOR(*v2)[i]);
  }
  
  return 0;
}



igraph_bool_t igraph_vector_complex_empty (const igraph_vector_complex_t* v) {
  return v->stor_begin == v->end;
}


long int igraph_vector_complex_size(const igraph_vector_complex_t* v) {
  return v->end - v->stor_begin;
}


igraph_bool_t igraph_vector_complex_isnull(const igraph_vector_complex_t *v) {
  long int n=igraph_vector_complex_size(v);
  long int i=0;
  
  while ( i<n && REALPART(VECTOR(*v)[i])==0 && IMAGPART(VECTOR(*v)[i])==0) {
    i++;
  }
  
  return i==n;
}


igraph_complex_t igraph_vector_complex_sum(const igraph_vector_complex_t *v) {
  igraph_complex_t res={0,0};
  igraph_complex_t *p;
  for (p=v->stor_begin; p<v->end; p++) {
    REALPART(res) += REALPART(*p);
    IMAGPART(res) += IMAGPART(*p);
  }
  return res;
}


igraph_complex_t igraph_vector_complex_prod(const igraph_vector_complex_t *v) {
  igraph_complex_t res={1,0};
  igraph_complex_t *p;

  for (p=v->stor_begin; p<v->end; p++) {
    igraph_complex_t c1=res;
    REALPART(res) = REALPART(c1) * REALPART(*p) - IMAGPART(c1) * IMAGPART(*p);
    IMAGPART(res) = REALPART(c1) * IMAGPART(*p) + IMAGPART(c1) * REALPART(*p);
  }
  return res;
}


igraph_bool_t igraph_vector_complex_is_equal(const igraph_vector_complex_t *lhs, 
					     const igraph_vector_complex_t *rhs) {
  long int i, s;
  
  s=igraph_vector_complex_size(lhs);
  if (s != igraph_vector_complex_size(rhs)) {
    return 0;
  } else {
    for (i=0; i<s; i++) {
      if (REALPART(VECTOR(*lhs)[i]) != REALPART(VECTOR(*rhs)[i]) ||
	  IMAGPART(VECTOR(*lhs)[i]) != IMAGPART(VECTOR(*rhs)[i])) {
	return 0;
      }
    }
    return 1;
  }
}


igraph_bool_t igraph_vector_complex_contains(const igraph_vector_complex_t *v, igraph_complex_t e) {
  igraph_complex_t *p=v->stor_begin;
  while (p<v->end) {
    if (REALPART(*p)==REALPART(e) &&
	IMAGPART(*p)==IMAGPART(e)) { 
      return 1;
    }
    p++;
  }
  return 0;
}


igraph_bool_t igraph_vector_complex_search(const igraph_vector_complex_t *v,
					     long int from, igraph_complex_t what, 
					     long int *pos) {
  long int i, n=igraph_vector_complex_size(v);  
  for (i=from; i<n; i++) {
    if (REALPART(VECTOR(*v)[i])==REALPART(what) && 
	IMAGPART(VECTOR(*v)[i])==IMAGPART(what)) break;
  }
  
  if (i<n) {
    if (pos != 0) {
      *pos=i;
    }
    return 1;
  } else {
    return 0;
  }
}



void igraph_vector_complex_clear(igraph_vector_complex_t* v) {
  v->end = v->stor_begin;
}


int igraph_vector_complex_resize(igraph_vector_complex_t* v, long int newsize) {
  IGRAPH_CHECK(igraph_vector_complex_reserve(v, newsize));
  v->end = v->stor_begin+newsize;
  return 0;
}


int igraph_vector_complex_reserve(igraph_vector_complex_t* v, long int size) {
  long int actual_size=igraph_vector_complex_size(v);
  igraph_complex_t *tmp;

  if (size <= igraph_vector_complex_size(v)) { return 0; }
  
  tmp=igraph_Realloc(v->stor_begin, size, igraph_complex_t);
  if (tmp==0) {
    IGRAPH_ERROR("cannot reserve space for vector", IGRAPH_ENOMEM);
  }
  v->stor_begin=tmp;
  v->stor_end=v->stor_begin + size;
  v->end=v->stor_begin+actual_size;
  
  return 0;
}


int igraph_vector_complex_push_back(igraph_vector_complex_t* v, igraph_complex_t e) {
	
  /* full, allocate more storage */
  if (v->stor_end == v->end) {
    long int new_size = igraph_vector_complex_size(v) * 2;
    if (new_size == 0) { new_size = 1; }
    IGRAPH_CHECK(igraph_vector_complex_reserve(v, new_size));
  }
  
  *(v->end) = e;
  v->end += 1;
  
  return 0;
}


igraph_complex_t igraph_vector_complex_pop_back(igraph_vector_complex_t* v) {
  igraph_complex_t tmp;
  tmp=igraph_vector_complex_e(v, igraph_vector_complex_size(v)-1);
  v->end -= 1;
  return tmp;
}


int igraph_vector_complex_insert(igraph_vector_complex_t *v, long int pos, igraph_complex_t value) {
  long int size = igraph_vector_complex_size(v);
  IGRAPH_CHECK(igraph_vector_complex_resize(v, size+1));
  if (pos<size) {
    memmove(v->stor_begin+pos+1, v->stor_begin+pos, 
	    sizeof(igraph_complex_t)*(size-pos));
  }
  v->stor_begin[pos] = value;
  return 0;
}


void igraph_vector_complex_remove(igraph_vector_complex_t *v, long int elem) {
  igraph_vector_complex_remove_section(v, elem, elem+1);
}


void igraph_vector_complex_remove_section(igraph_vector_complex_t *v, 
					    long int from, long int to) {
  memmove(v->stor_begin+from, v->stor_begin+to,
	  sizeof(igraph_complex_t)*(v->end-v->stor_begin-to));
  v->end -= (to-from);
}



