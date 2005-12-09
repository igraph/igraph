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

int vector_ptr_init      (vector_ptr_t* v, int long size) {	
        long int alloc_size= size > 0 ? size : 1;
	assert(v != NULL);
	if (size < 0) { size=0; }
	v->stor_begin=Calloc(alloc_size, void*);
	if (v->stor_begin==0) {
	  IGRAPH_FERROR("vector ptr init failed", IGRAPH_ENOMEM);
	}
	v->stor_end=v->stor_begin + alloc_size;
	v->end=v->stor_begin+size;

	return 0;
}

/**
 * \ingroup vectorptr
 * \brief Destroys a pointer vector.
 */

void vector_ptr_destroy   (vector_ptr_t* v) {
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

void vector_ptr_free_all   (vector_ptr_t* v) {
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

void vector_ptr_destroy_all   (vector_ptr_t* v) { 
  assert(v != 0);
  assert(v->stor_begin != 0);
  vector_ptr_free_all(v);
  vector_ptr_destroy(v);
}

/**
 * \ingroup vectorptr
 * \brief Reserves memory for a pointer vector for later use.
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int vector_ptr_reserve   (vector_ptr_t* v, long int size) {
	long int actual_size=vector_ptr_size(v);
	void **tmp;
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	
	if (size <= vector_ptr_size(v)) { return 0; }

	tmp=Realloc(v->stor_begin, size, void*);
	if (tmp==0) {
	  IGRAPH_FERROR("vector ptr reserve failed", IGRAPH_ENOMEM);
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

bool_t vector_ptr_empty     (vector_ptr_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);	
	return v->stor_begin == v->end;
}

/**
 * \ingroup vectorptr
 * \brief Gives the number of elements in the pointer vector.
 * \todo why does the assert fail???
 */

long int vector_ptr_size      (vector_ptr_t* v) {
	assert(v != NULL);
/* 	assert(v->stor_begin != NULL);		 */
	return v->end - v->stor_begin;
}

/**
 * \ingroup vectorptr
 * \brief Removes all elements in a pointer vector.
 */

void vector_ptr_clear     (vector_ptr_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);		
	v->end = v->stor_begin;	
}

/**
 * \ingroup vectorptr
 * \brief Appends an elements to the back of a pointer vector.
 */

int vector_ptr_push_back (vector_ptr_t* v, void* e) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);	

	/* full, allocate more storage */
	if (v->stor_end == v->end) {
		long int new_size = vector_ptr_size(v) * 2;
		if (new_size == 0) { new_size = 1; }
		IGRAPH_CHECK(vector_ptr_reserve(v, new_size));
	}
	
	*(v->end) = e;
	v->end += 1;
	
	return 0;
}

/**
 * \ingroup vectorptr
 * \brief Access an element of a pointer vector.
 */

void* vector_ptr_e         (vector_ptr_t* v, long int pos) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);	
	return * (v->stor_begin + pos);
}

/**
 * \ingroup vectorptr
 * \brief Assign to an element of a pointer vector.
 */

void vector_ptr_set       (vector_ptr_t* v, long int pos, void* value) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);		
	*(v->stor_begin + pos) = value;
}

/**
 * \ingroup vectorptr
 * \brief Set all elements of a pointer vector to the NULL pointer.
 */

void vector_ptr_null      (vector_ptr_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);		
	if (vector_ptr_size(v)>0) {
		memset(v->stor_begin, 0, sizeof(void*)*vector_ptr_size(v));
	}
}

/**
 * \ingroup vectorptr
 * \brief Resizes a pointer vector.
 */

int vector_ptr_resize(vector_ptr_t* v, long int newsize) {
  int ret;
  IGRAPH_CHECK(vector_ptr_reserve(v, newsize));
  v->end = v->stor_begin+newsize;
  return 0;
}

/**
 * \ingroup vectorptr
 * \brief Handle a C array as a pointer vector temporarily. Don't ever 
 * call vector_ptr_destroy() on these vectors.
 */

vector_ptr_t vector_ptr_as_vector(void** data, long int length) {
  vector_ptr_t v;
  v.stor_begin=data;
  v.stor_end=v.end=data+length;
  
  return v;
}

/**
 * \ingroup vectorptr
 * \brief Initializes a pointer vector from an array (constructor).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int vector_ptr_init_copy(vector_ptr_t *v, void* *data, long int length) {
  v->stor_begin=Calloc(length, void*);
  if (v->stor_begin==0) {
    IGRAPH_FERROR("cannot init ptr vector from array", IGRAPH_ENOMEM);
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

void vector_ptr_copy_to(vector_ptr_t *v, void** to) {
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

int vector_ptr_copy(vector_ptr_t *to, vector_ptr_t *from) {
  assert(from != NULL);
/*   assert(from->stor_begin != NULL); */
  to->stor_begin=Calloc(vector_ptr_size(from), void*);
  if (to->stor_begin==0) {
    IGRAPH_FERROR("cannot copy ptr vector", IGRAPH_ENOMEM);
  }
  to->stor_end=to->stor_begin+vector_ptr_size(from);
  to->end=to->stor_end;
  memcpy(to->stor_begin, from->stor_begin, vector_ptr_size(from)*sizeof(void*));
  
  return 0;
}

/**
 * \ingroup vectorptr
 * \brief Remove an element from a pointer vector.
 */

void vector_ptr_remove(vector_ptr_t *v, long int pos) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  if (pos+1<vector_ptr_size(v)) { /* TOOD: why is this needed */
    memcpy(v->stor_begin+pos, v->stor_begin+pos+1,
	   sizeof(void*)*(vector_ptr_size(v)-pos));
  }
  v->end--;
}
