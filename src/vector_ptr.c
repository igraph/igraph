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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "types.h"
#include "memory.h"
#include "random.h"
#include "error.h"

#include <assert.h>
#include <string.h> 		/* memcpy & co. */
#include <stdlib.h>

/**
 * \ingroup vectorptr
 * \brief Initialize a poiter vector (constructor).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int igraph_vector_ptr_init      (igraph_vector_ptr_t* v, int long size) {	
        long int alloc_size= size > 0 ? size : 1;
	assert(v != NULL);
	if (size < 0) { size=0; }
	v->stor_begin=Calloc(alloc_size, void*);
	if (v->stor_begin==0) {
	  IGRAPH_ERROR("vector ptr init failed", IGRAPH_ENOMEM);
	}
	v->stor_end=v->stor_begin + alloc_size;
	v->end=v->stor_begin+size;

	return 0;
}

/**
 */ 

const igraph_vector_ptr_t *igraph_vector_ptr_view (const igraph_vector_ptr_t *v, void *const *data, 
				     long int length) {
  igraph_vector_ptr_t *v2=(igraph_vector_ptr_t*) v;
  v2->stor_begin=(void **)data;
  v2->stor_end=(void**)data+length;
  v2->end=v2->stor_end;
  return v;
}

/**
 * \ingroup vectorptr
 * \brief Destroys a pointer vector.
 */

void igraph_vector_ptr_destroy   (igraph_vector_ptr_t* v) {
  assert(v != 0);
  if (v->stor_begin != 0) {
    Free(v->stor_begin);
    v->stor_begin = NULL;
  }
}

/**
 * \ingroup vectorptr
 * \brief Calls free() on all elements of a pointer vector.
 */

void igraph_vector_ptr_free_all   (igraph_vector_ptr_t* v) {
  void **ptr;
  assert(v != 0);
  assert(v->stor_begin != 0);
  for (ptr=v->stor_begin; ptr<v->end; ptr++) {
    Free(*ptr);
  }
}

/**
 * \ingroup vectorptr
 * \brief Calls free() on all elements and destroys the pointer vector.
 */

void igraph_vector_ptr_destroy_all   (igraph_vector_ptr_t* v) { 
  assert(v != 0);
  assert(v->stor_begin != 0);
  igraph_vector_ptr_free_all(v);
  igraph_vector_ptr_destroy(v);
}

/**
 * \ingroup vectorptr
 * \brief Reserves memory for a pointer vector for later use.
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int igraph_vector_ptr_reserve   (igraph_vector_ptr_t* v, long int size) {
	long int actual_size=igraph_vector_ptr_size(v);
	void **tmp;
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	
	if (size <= igraph_vector_ptr_size(v)) { return 0; }

	tmp=Realloc(v->stor_begin, size, void*);
	if (tmp==0) {
	  IGRAPH_ERROR("vector ptr reserve failed", IGRAPH_ENOMEM);
	}
	v->stor_begin=tmp;
	v->stor_end=v->stor_begin + size;
	v->end=v->stor_begin+actual_size;
	
	return 0;
}

/**
 * \ingroup vectorptr
 * \brief Decides whether the pointer vector is empty.
 */

igraph_bool_t igraph_vector_ptr_empty     (const igraph_vector_ptr_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);	
	return v->stor_begin == v->end;
}

/**
 * \ingroup vectorptr
 * \brief Gives the number of elements in the pointer vector.
 * \todo why does the assert fail???
 */

long int igraph_vector_ptr_size      (const igraph_vector_ptr_t* v) {
	assert(v != NULL);
/* 	assert(v->stor_begin != NULL);		 */
	return v->end - v->stor_begin;
}

/**
 * \ingroup vectorptr
 * \brief Removes all elements in a pointer vector.
 */

void igraph_vector_ptr_clear     (igraph_vector_ptr_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);		
	v->end = v->stor_begin;	
}

/**
 * \ingroup vectorptr
 * \brief Appends an elements to the back of a pointer vector.
 */

int igraph_vector_ptr_push_back (igraph_vector_ptr_t* v, void* e) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);	

	/* full, allocate more storage */
	if (v->stor_end == v->end) {
		long int new_size = igraph_vector_ptr_size(v) * 2;
		if (new_size == 0) { new_size = 1; }
		IGRAPH_CHECK(igraph_vector_ptr_reserve(v, new_size));
	}
	
	*(v->end) = e;
	v->end += 1;
	
	return 0;
}

/**
 * \ingroup vectorptr
 * \brief Access an element of a pointer vector.
 */

void* igraph_vector_ptr_e         (const igraph_vector_ptr_t* v, long int pos) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);	
	return * (v->stor_begin + pos);
}

/**
 * \ingroup vectorptr
 * \brief Assign to an element of a pointer vector.
 */

void igraph_vector_ptr_set       (igraph_vector_ptr_t* v, long int pos, void* value) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);		
	*(v->stor_begin + pos) = value;
}

/**
 * \ingroup vectorptr
 * \brief Set all elements of a pointer vector to the NULL pointer.
 */

void igraph_vector_ptr_null      (igraph_vector_ptr_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);		
	if (igraph_vector_ptr_size(v)>0) {
		memset(v->stor_begin, 0, sizeof(void*)*igraph_vector_ptr_size(v));
	}
}

/**
 * \ingroup vectorptr
 * \brief Resizes a pointer vector.
 */

int igraph_vector_ptr_resize(igraph_vector_ptr_t* v, long int newsize) {
  IGRAPH_CHECK(igraph_vector_ptr_reserve(v, newsize));
  v->end = v->stor_begin+newsize;
  return 0;
}

/**
 * \ingroup vectorptr
 * \brief Initializes a pointer vector from an array (constructor).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int igraph_vector_ptr_init_copy(igraph_vector_ptr_t *v, void* *data, long int length) {
  v->stor_begin=Calloc(length, void*);
  if (v->stor_begin==0) {
    IGRAPH_ERROR("cannot init ptr vector from array", IGRAPH_ENOMEM);
  }
  v->stor_end=v->stor_begin+length;
  v->end=v->stor_end;
  memcpy(v->stor_begin, data, length*sizeof(void*));
  
  return 0;
}

/**
 * \ingroup vectorptr
 * \brief Copy the contents of a pointer vector to a regular C array.
 */

void igraph_vector_ptr_copy_to(const igraph_vector_ptr_t *v, void** to) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);		
  if (v->end != v->stor_begin) {
    memcpy(to, v->stor_begin, sizeof(void*) * (v->end - v->stor_begin));
  }
}

/**
 * \ingroup vectorptr
 * \brief Copy a pointer vector (constructor).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 * \todo why does the assert fail???
 */

int igraph_vector_ptr_copy(igraph_vector_ptr_t *to, const igraph_vector_ptr_t *from) {
  assert(from != NULL);
/*   assert(from->stor_begin != NULL); */
  to->stor_begin=Calloc(igraph_vector_ptr_size(from), void*);
  if (to->stor_begin==0) {
    IGRAPH_ERROR("cannot copy ptr vector", IGRAPH_ENOMEM);
  }
  to->stor_end=to->stor_begin+igraph_vector_ptr_size(from);
  to->end=to->stor_end;
  memcpy(to->stor_begin, from->stor_begin, igraph_vector_ptr_size(from)*sizeof(void*));
  
  return 0;
}

/**
 * \ingroup vectorptr
 * \brief Remove an element from a pointer vector.
 */

void igraph_vector_ptr_remove(igraph_vector_ptr_t *v, long int pos) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  if (pos+1<igraph_vector_ptr_size(v)) { /* TOOD: why is this needed */
    memcpy(v->stor_begin+pos, v->stor_begin+pos+1,
	   sizeof(void*)*(igraph_vector_ptr_size(v)-pos));
  }
  v->end--;
}
