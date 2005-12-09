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
 * \ingroup vector
 * \brief Initializes a vector object (constructor).
 * 
 * Every vector needs to be initialized before it can be used, and
 * there are a number initialization functions, the others are
 * vector_init_copy(), vector_init_seq() and vector_copy().
 * 
 * Every vector object initialized by this function should be
 * destroyed (ie. the memory allocated for it should be freed) when it
 * is not needed any more, the vector_destroy() function is
 * responsible for this.
 * @param v Pointer to a not yet initialized vector object.
 * @param size The size of the vector.
 * @return error code:
 *         - <b>IGRAPH_ENOMEM</b> if there is not enough memory.
 * 
 * Time complexity: operating system dependent.
 */

int vector_init      (vector_t* v, int long size) {	
        long int alloc_size= size > 0 ? size : 1;
	if (size < 0) { size=0; }
	v->stor_begin=Calloc(alloc_size, real_t);
	if (v->stor_begin==0) {
	  IGRAPH_FERROR("cannot init vector", IGRAPH_ENOMEM);
	}
	v->stor_end=v->stor_begin + alloc_size;
	v->end=v->stor_begin+size;

	return 0;
}

/**
 * \ingroup vector
 * \brief Destroys a vector object
 *
 * All vectors initialized by vector_init() should be properly
 * destroyed by this function. A destroyed vector needs to be
 * reinitialized by vector_init(), vector_init_copy() or
 * vector_as_vector() before using it again.
 * @param v Pointer to the (previously initialized) vector object to
 *        destroy. 
 * 
 * Time complexity: operating system dependent.
 */

void vector_destroy   (vector_t* v) {
  assert(v != 0);
  if (v->stor_begin != 0) {
    Free(v->stor_begin);
    v->stor_begin = NULL;
  }
}

/**
 * \ingroup vector
 * \brief Reserves memory for a vector
 * 
 * \a igraph vectors are flexible, they can grow and shrink. Growing
 * however occasionally needs the data in the vector to be copyed.
 * In order to avoid you can call this function to reserve space for
 * future growth of the vector. 
 * 
 * Note that this function does <em>not</em> change the size of the
 * vector, let us see a small example to clarify things: if you
 * reserve space for 100 elements and the size of your
 * vector was (and still is) 60, then you can surely add additional 40
 * elements to your vector before it will be copied.
 * @param v The vector object.
 * @param size The new <em>allocated</em> size of the vector.
 * @return Error code:
 *         - <b>IGRPAH_ENOMEM</b> if there is not enough memory.
 *
 * Time complexity: operating system dependent, should be around
 * <code>O(n)</code>, <code>n</code> is the new allocated size of the
 * vector.
 */

int vector_reserve   (vector_t* v, long int size) {
	long int actual_size=vector_size(v);
	real_t *tmp;
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	if (size <= vector_size(v)) { return 0; }

	tmp=Realloc(v->stor_begin, size, real_t);
	if (tmp==0) {
	  IGRAPH_FERROR("cannot reserve space for vector", IGRAPH_ENOMEM);
	}
	v->stor_begin=tmp;
	v->stor_end=v->stor_begin + size;
	v->end=v->stor_begin+actual_size;
	
	return 0;
}

/**
 * \ingroup vector
 * \brief Decides whether the size of the vector is zero.
 *
 * @param v The vector object.
 * @return Non-zero number if the size of the vector is not zero and
 *         zero otherwise.
 * 
 * Time complexity: <code>O(1)</code>.
 */

bool_t vector_empty     (vector_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	return v->stor_begin == v->end;
}

/**
 * \ingroup vector
 * \brief Gives the size (=length) of the vector.
 * 
 * @param v The vector object
 * @return The size of the vector.
 *
 * Time complexity: <code>O(1)</code>. 
 */

long int vector_size      (vector_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	return v->end - v->stor_begin;
}

/**
 * \ingroup vector
 * \brief Removes all elements from a vector.
 * 
 * This function simply sets the size of the vector to zero, it does
 * not free any allocated memory. For that you have to call
 * vector_destroy().
 * @param v The vector object.
 * 
 * Time complexity: <code>O(1)</code>.
 */

void vector_clear     (vector_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	v->end = v->stor_begin;
}

/**
 * \ingroup vector
 * \brief Appends one element to a vector.
 * 
 * This function resizes the vector to be one element longer and
 * sets the very last element in the vector to <code>e</code>.
 * @param v The vector object.
 * @param e The element to append to the vector.
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: not enough memory.
 * 
 * Time complexity: operating system dependent. What is important that
 * it can be assumed that <code>n</code> subsequent calls to this
 * function has time complexity <code>O(n)</code>, even if there
 * hadn't been any space reserved for the new elements by
 * vector_reserve(). This is 
 * implemented by a trick similar to the C++ <code>vector</code>
 * class: each time more memory is allocated for a vector, the size of
 * the additionally allocated memory is the same as the vector's
 * current length.
 */

int vector_push_back (vector_t* v, real_t e) {
  	assert(v != NULL);
	assert(v->stor_begin != NULL);
	
	/* full, allocate more storage */
	if (v->stor_end == v->end) {
		long int new_size = vector_size(v) * 2;
		if (new_size == 0) { new_size = 1; }
		IGRAPH_CHECK(vector_reserve(v, new_size));
	}
	
	*(v->end) = e;
	v->end += 1;
	
	return 0;
}

/**
 * \ingroup vector
 * \brief Access an element of a vector.
 */

real_t vector_e         (vector_t* v, long int pos) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	return * (v->stor_begin + pos);
}

/**
 * \ingroup vector
 * \brief Get the address of an element of a vector
 */

real_t*vector_e_ptr  (vector_t* v, long int pos) {
  assert(v!=NULL);
  assert(v->stor_begin != NULL);
  return v->stor_begin+pos;
}

/**
 * \ingroup vector
 * \brief Assignment to an element of a vector.
 */

void vector_set       (vector_t* v, long int pos, real_t value) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);	
	*(v->stor_begin + pos) = value;
}

/**
 * \ingroup vector
 * \brief Sets each element in the vector to zero.
 * 
 * @param v The vector object.
 *
 * Time complexity: <code>O(n)</code>, the size of the vector.
 */

void vector_null      (vector_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	if (vector_size(v)>0) {
		memset(v->stor_begin, 0, sizeof(real_t)*vector_size(v));
	}
}

/**
 * \ingroup vector
 * \brief Returns the last element in a vector.
 *
 * It is an error to call this function on an empty vector.
 * @param v The vector object.
 * @return The last element.
 * 
 * Time complexity: <code>O(1)</code>.
 */

real_t vector_tail(vector_t *v) {
  assert(v!=NULL);
  assert(v->stor_begin != NULL);
  return *((v->end)-1);
}

/**
 * \ingroup vector
 * \brief Removes and returns the last element of a vector.
 */

real_t vector_pop_back(vector_t* v) {
  real_t tmp;
  assert(v!=NULL);
  assert(v->stor_begin != NULL);
  assert(v->end != v->stor_begin);
  tmp=vector_e(v, vector_size(v)-1);
  v->end -= 1;
  return tmp;
}

/**
 * \ingroup vector
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int vector_order(vector_t* v, vector_t* res, integer_t nodes) {
  long int edges=vector_size(v);
  vector_t ptr;
  vector_t rad;
  long int i, j;
  int ret1;

  assert(v!=NULL);
  assert(v->stor_begin != NULL);

  VECTOR_INIT_FINALLY(&ptr, nodes+1);
  VECTOR_INIT_FINALLY(&rad, edges);
  IGRAPH_CHECK(vector_resize(res, edges));
  
  for (i=0; i<edges; i++) {
    long int radix=v->stor_begin[i];
    if (VECTOR(ptr)[radix]!=0) {
      VECTOR(rad)[i]=VECTOR(ptr)[radix];
    }
    VECTOR(ptr)[radix]=i+1;
  }
  
  j=0;
  for (i=0; i<nodes+1; i++) {
    if (VECTOR(ptr)[i] != 0) {
      long int next=VECTOR(ptr)[i]-1;
      res->stor_begin[j++]=next;
      while (VECTOR(rad)[next] != 0) {
	next=VECTOR(rad)[next]-1;
	res->stor_begin[j++]=next;
      }
    }
  }
  
  vector_destroy(&ptr);
  vector_destroy(&rad);
  IGRAPH_FINALLY_CLEAN(2);
  
  return 0;
}

/**
 * \ingroup vector
 * \brief Internal comparision function of vector elements, used by 
 * vector_sort().
 */

int vector_sort_cmp(const void *a, const void *b) {
  const real_t *da = (const real_t *) a;
  const real_t *db = (const real_t *) b;

  return (*da > *db) - (*da < *db);
}

/**
 * \ingroup vector
 * \brief Sorts the elements of the vector into ascending order.
 * 
 * This function uses the built-in sort function of the C library.
 * @param v Pointer to an initialized vector object.
 *
 * Time complexity: should be <code>O(nlogn)</code> for <code>n</code>
 * elements.
 */

void vector_sort(vector_t *v) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  qsort(v->stor_begin, vector_size(v), sizeof(real_t), vector_sort_cmp);
}

/**
 * \ingroup vector
 * \brief Resize the vector.
 *
 * Note that this function does not free any memory, just sets the
 * size of the vector to the given one. It can on the other hand 
 * allocate more memory if the new size is larger than the previous
 * one. In this case the newly appeared elements in the vectors are
 * <em>not</em> set to zero, they are uninitialized.
 * @param v The vector object
 * @param newsize The new size of the vector.
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: not enough memory.
 * 
 * Time complexity: <code>O(1)</code> if the new size is smaller,
 * operating system dependent if it is larger. In the latter case it
 * is usually around <code>O(n)</code>, <code>n</code> is the new size
 * of the vector.
 */

int vector_resize(vector_t* v, long int newsize) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  IGRAPH_CHECK(vector_reserve(v, newsize));
  v->end = v->stor_begin+newsize;
  return 0;
}

/**
 * \ingroup vector
 * \brief Gives the maximum element of the vector.
 *
 * If the size of the vector is zero, an arbitrary number is
 * returned.
 * @param v The vector object.
 * @return The maximum element.
 *
 * Time complexity: <code>O(n)</code>, <code>n</code> is the size of
 * the vector. 
 */

real_t vector_max(vector_t* v) {
  real_t max;
  real_t *ptr;
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  max=*(v->stor_begin);
  ptr=v->stor_begin+1;
  while (ptr < v->end) {
    if ((*ptr) > max) {
      max=*ptr;
    }
    ptr++;
  }
  return max;
}

/**
 * \ingroup vector
 * \brief Handle a regular C array as a vector, don't ever call 
 * vector_destroy() on this!
 */

vector_t vector_as_vector(real_t* data, long int length) {
  vector_t v;
  v.stor_begin=data;
  v.stor_end=v.end=data+length;
  
  return v;
}

/**
 * \ingroup vector
 * \brief Initializes a vector from an ordinary C array (constructor).
 * 
 * @param v Pointer to an uninitialized vector object.
 * @param data A regular C array.
 * @param length The length of the C array.
 * @return Error code: 
 *         - <b>IGRAPH_ENOMEM</b> if there is not enough memory.
 * 
 * Time complexity: operating system specific, usually
 * <code>O(length)</code>.
 */

int vector_init_copy(vector_t *v, real_t *data, long int length) {
  v->stor_begin=Calloc(length, real_t);
  if (v->stor_begin==0) {
    IGRAPH_FERROR("cannot init vector from array", IGRAPH_ENOMEM);
  }
  v->stor_end=v->stor_begin+length;
  v->end=v->stor_end;
  memcpy(v->stor_begin, data, length*sizeof(real_t));
  
  return 0;
}

/**
 * \ingroup vector
 * \brief Copies the contents of a vector to a C array.
 * 
 * The C array should have sufficient length.
 * @param v The vector object.
 * @param to The C array.
 * 
 * Time complexity: <code>O(n)</code>, <code>n</code> is the size of
 * the vector.
 */

void vector_copy_to(vector_t *v, real_t* to) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  if (v->end != v->stor_begin) {
    memcpy(to, v->stor_begin, sizeof(real_t) * (v->end - v->stor_begin));
  }
}

/**
 * \ingroup vector
 * \brief Initializes a vector from another vector object (constructor).
 * 
 * The contents of the existing vector object will be copied to
 * the new one.
 * @param to Pointer to a not yet initialized vector object.
 * @param from The original vector object to copy.
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b> if there is not enough memory.
 * 
 * Time complexity: operating system dependent, usually
 * <code>O(n)</code>, <code>n</code> is the size of the vector.
 */

int vector_copy(vector_t *to, vector_t *from) {
  assert(from != NULL);
  assert(from->stor_begin != NULL);
  to->stor_begin=Calloc(vector_size(from), real_t);
  if (to->stor_begin==0) {
    IGRAPH_FERROR("canot copy vector", IGRAPH_ENOMEM);
  }
  to->stor_end=to->stor_begin+vector_size(from);
  to->end=to->stor_end;
  memcpy(to->stor_begin, from->stor_begin, vector_size(from)*sizeof(real_t));
  
  return 0;
}

/**
 * \ingroup vector
 * \brief Calculates the sum of the elements in the vector.
 *
 * For the empty vector 0.0 is returned.
 * @param v The vector object.
 * @return The sum of the elements.
 * 
 * 
 * Time complexity: <code>O(n)</code>, the size of the vector.
 */

real_t vector_sum(vector_t *v) {
  real_t res=0;
  real_t *p;
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  for (p=v->stor_begin; p<v->end; p++) {
    res += *p;
  }
  return res;
}

/**
 * \ingroup vector
 * \brief Calculates the product of the elements in the vector.
 * 
 * For the empty vector one (1) is returned.
 * @param v The vector object.
 * @return The product of the elements.
 * 
 * Time complexity: <code>O(n)</code>, the size of the vector.
 */

real_t vector_prod(vector_t *v) {
  real_t res=1;
  real_t *p;
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  for (p=v->stor_begin; p<v->end; p++) {
    res *= *p;
  }
  return res;
}

/**
 * \ingroup vector
 * \brief Initializes a vector with a sequence.
 * 
 * The vector will contain the numbers <code>from</code>,
 * <code>from+1</code>, ..., <code>to</code>.
 * @param v Pointer to an uninitialized vector object.
 * @param from The lower limit in the sequence (inclusive).
 * @param to The upper limit in the sequence (inclusive).
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory.
 *
 * Time complexity: <code>O(n)</code>, the number of elements in the
 * vector. 
 */

int vector_init_seq(vector_t *v, real_t from, real_t to) {
  real_t *p;
  IGRAPH_CHECK(vector_init(v, to-from+1));

  for (p=v->stor_begin; p<v->end; p++) {
    *p = from++;
  }
  
  return 0;
}

/**
 * \ingroup vector
 * \brief Deletes a section from a vector.
 */

void vector_remove_section(vector_t *v, long int from, long int to) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  memmove(v->stor_begin+from, v->stor_begin+to,
	  sizeof(real_t)*(v->end-v->stor_begin-to));
  v->end -= (to-from);
}

/**
 * \ingroup vector
 * \brief Removes a single element from a vector.
 */

void vector_remove(vector_t *v, long int elem) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  vector_remove_section(v, elem, elem+1);
}

/**
 * \ingroup vector
 * \brief Copies a section of a vector.
 */

int vector_move_interval(vector_t *v, long int begin, long int end, 
			 long int to) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  memcpy(v->stor_begin+to, v->stor_begin+begin, 
	 sizeof(real_t)*(end-begin));

  return 0;
}

/**
 * \ingroup vector
 * \brief Remove elements of a vector (for internal use).
 */

void vector_permdelete(vector_t *v, long int *index, long int nremove) {
  long int i;
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  for (i=0; i<vector_size(v); i++) {
    if (index[i] != 0) {
      VECTOR(*v)[ index[i]-1 ] = VECTOR(*v)[i];
    }
  }
  v->end -= nremove;
}

/**
 * \ingroup vector
 * \brief Remove elements of a vector (for internal use).
 */

void vector_remove_negidx(vector_t *v, vector_t *neg, long int nremove) {
  long int i, idx=0;
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  for (i=0; i<vector_size(v); i++) {
    VECTOR(*v)[idx++] = VECTOR(*v)[i];
  }
  v->end -= nremove;
}

/**
 * \ingroup vector
 * \brief Checks if all elements of a vector are in the given interval.
 */

bool_t vector_isininterval(vector_t *v, real_t low, real_t high) {
  real_t *ptr;
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  for (ptr=v->stor_begin; ptr<v->end; ptr++) {
    if (*ptr < low || *ptr >high) {
      return 0;
    }
  }
  return 1;
}

/**
 * \ingroup vector
 * \brief Checks if all elements of a vector are smaller than a limit.
 */

bool_t vector_any_smaller(vector_t *v, real_t limit) {
  real_t *ptr;
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  for (ptr=v->stor_begin; ptr<v->end; ptr++) {
    if (*ptr < limit) {
      return 1;
    }
  }
  return 0;
}

/**
 * \ingroup vecgtor
 * \brief Decides whether two vectors contain exactly the same elements.
 */

bool_t vector_is_equal(vector_t *lhs, vector_t *rhs) {
  long int i, s;
  assert(lhs != 0);
  assert(rhs != 0);
  assert(lhs->stor_begin != 0);
  assert(rhs->stor_begin != 0);
  
  s=vector_size(lhs);
  if (s != vector_size(rhs)) {
    return 0;
  } else {
    for (i=0; i<s; i++) {
      if (VECTOR(*lhs)[i] != VECTOR(*rhs)[i]) {
	return 0;
      }
    }
    return 1;
  }
}
