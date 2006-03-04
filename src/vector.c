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
#include <stdarg.h>		/* va_start & co */

/**
 * \ingroup vector
 * \section about_igraph_vector_t_objects About \type igraph_vector_t objects
 * 
 * <para>The \type igraph_vector_t data type is a simple and efficient
 * interface to arrays containing real numbers. It is something
 * similar as (but much simpler than) the \type vector template
 * in the C++ standard library.</para>
 *
 * <para>Vectors are used extensively in \a igraph, all
 * functions which expect or return a list of numbers use
 * igraph_vector_t to achieve this.</para>
 *
 * <para>The \type igraph_vector_t type usually uses
 * O(n) space
 * to store n elements. Sometimes it
 * uses more, this is because vectors can shrink, but even if they
 * shrink, the current implementation does not free a single bit of
 * memory.</para>
 * 
 * <para>The position of the elements in a \type igraph_vector_t
 * object is numbered from zero, as this is the usual C
 * standard.</para> 
 */

/**
 * \ingroup vector
 * \section igraph_vector_constructors_and_destructors Constructors and
 * Destructors
 * 
 * <para>\type igraph_vector_t objects have to be initialized before using
 * them, this is analogous to calling a constructor on them. There are a
 * number of \type igraph_vector_t constructors, for your
 * convenience. \ref igraph_vector_init() is the basic constructor, it
 * creates a vector of the given length, filled with zeros.
 * \ref igraph_vector_init_real(), \ref igraph_vector_init_real_end(), \ref
 * igraph_vector_init_int() and \ref igraph_vector_init_int_end() are convenience
 * constructors, these create a vector with the elements given as
 * their parameters. \ref igraph_vector_copy() creates a new identical copy
 * of an already existing and initialized vector. \ref
 * igraph_vector_init_copy() creates a vector by copying a regular C array. 
 * \ref igraph_vector_init_seq() creates a vector containing a regular
 * sequence with increment one.</para>
 * 
 * <para>\ref igraph_vector_view() is a special constructor, it allows you to
 * handle a regular C array as a \type vector without copying
 * its elements.
 * </para> 
 *
 * <para>If a \type igraph_vector_t object is not needed any more, it
 * should be destroyed to free its allocated memory by calling the
 * \type igraph_vector_t destructor, \ref igraph_vector_destroy().</para>
 * 
 * <para> Note that vectors created by \ref igraph_vector_view() are special,
 * you mustn't call \ref igraph_vector_destroy() on these.</para>
 */

/**
 * \ingroup vector
 * \function igraph_vector_init
 * \brief Initializes a vector object (constructor).
 * 
 * Every vector needs to be initialized before it can be used, and
 * there are a number of initialization functions or otherwise called
 * constructors. 
 * 
 * Every vector object initialized by this function should be
 * destroyed (ie. the memory allocated for it should be freed) when it
 * is not needed anymore, the \ref igraph_vector_destroy() function is
 * responsible for this.
 * \param v Pointer to a not yet initialized vector object.
 * \param size The size of the vector.
 * \return error code:
 *       \c IGRAPH_ENOMEM if there is not enough memory.
 * 
 * Time complexity: operating system dependent, the amount of 
 * \quote time \endquote required to allocate
 * O(n) elements,
 * n is the number of elements. 
 */

int igraph_vector_init      (igraph_vector_t* v, int long size) {	
        long int alloc_size= size > 0 ? size : 1;
	if (size < 0) { size=0; }
	v->stor_begin=Calloc(alloc_size, real_t);
	if (v->stor_begin==0) {
	  IGRAPH_ERROR("cannot init vector", IGRAPH_ENOMEM);
	}
	v->stor_end=v->stor_begin + alloc_size;
	v->end=v->stor_begin+size;

	return 0;
}

/**
 * \ingroup vector
 * \function igraph_vector_view
 * \brief Handle a regular C array as a \type igraph_vector_t.
 * 
 * This is a special \type igraph_vector_t constructor. It allows to
 * handle a regular C array as a \type igraph_vector_t temporarily.
 * Be sure that you \em don't ever call the destructor (\ref
 * igraph_vector_destroy()) on objects created by this constructor.
 * \param v Pointer to an uninitialized \type igraph_vector_t object.
 * \param data Pointer, the C array.
 * \param length The length of the C array.
 * \return Pointer to the vector object, the same as the 
 *     \p v parameter, for convenience.
 * 
 * Time complexity: O(1)
 */ 

const igraph_vector_t *igraph_vector_view (const igraph_vector_t *v, const real_t *data, 
			     long int length) {
  igraph_vector_t *v2=(igraph_vector_t*) v;
  v2->stor_begin=(real_t*)data;
  v2->stor_end=(real_t*)data+length;
  v2->end=v2->stor_end;
  return v;
}

/**
 * \ingroup vector
 * \function igraph_vector_init_real
 * \brief Create an \type igraph_vector_t from the parameters.
 * 
 * Because of how C and the C library handles variable length argument
 * lists, it is required that you supply real constants to this
 * function. This means that
 * \verbatim igraph_vector_t v;
 * igraph_vector_init_real(&amp;v, 5, 1,2,3,4,5); \endverbatim
 * is an error at runtime and the results are undefined. This is
 * the proper way:
 * \verbatim igraph_vector_t v;
 * igraph_vector_init_real(&amp;v, 5, 1.0,2.0,3.0,4.0,5.0); \endverbatim
 * \param v Pointer to an uninitialized \type igraph_vector_t object.
 * \param no Positive integer, the number of \type real_t
 *    parameters to follow.
 * \param ... The elements of the vector.
 * \return Error code, this can be \c IGRAPH_ENOMEM
 *     if there isn't enough memory to allocate the vector.
 *
 * \sa \ref igraph_vector_init_real_end(), \ref igraph_vector_init_int() for similar
 * functions.
 *
 * Time complexity: depends on the time required to allocate memory,
 * but at least O(n), the number of
 * elements in the vector.
 */

int igraph_vector_init_real(igraph_vector_t *v, int no, ...) {
  int i=0;
  va_list ap;
  IGRAPH_CHECK(igraph_vector_init(v, no));

  va_start(ap, no);
  for (i=0; i<no; i++) {
    VECTOR(*v)[i]=(real_t) va_arg(ap, double);
  }
  va_end(ap);
  
  return 0;
}

/**
 * \ingroup vector
 * \function igraph_vector_init_real_end
 * \brief Create an \type igraph_vector_t from the parameters.
 * 
 * This constructor is similar to \ref igraph_vector_init_real(), the only
 * difference is that instead of giving the number of elements in the
 * vector, a special marker element follows the last real vector
 * element.
 * \param v Pointer to an uninitialized \type igraph_vector_t object.
 * \param endmark This element will signal the end of the vector. It
 *    will \em not be part of the vector.
 * \param ... The elements of the vector.
 * \return Error code, \c IGRAPH_ENOMEM if there
 *    isn't enough memory.
 * 
 * \sa \ref igraph_vector_init_real() and \ref igraph_vector_init_int_end() for
 * similar functions.
 * 
 * Time complexity: at least O(n) for 
 * n elements plus the time
 * complexity of the memory allocation.
 */

int igraph_vector_init_real_end(igraph_vector_t *v, real_t endmark, ...) {
  int i=0, n=0;
  va_list ap;

  va_start(ap, endmark);
  while (1) {
    real_t num = va_arg(ap, double);
    if (num == endmark) {
      break;
    }
    n++;
  }
  va_end(ap);

  IGRAPH_VECTOR_INIT_FINALLY(v, n);
  
  va_start(ap, endmark);
  for (i=0; i<n; i++) {
    VECTOR(*v)[i]=(real_t) va_arg(ap, double);
  }
  va_end(ap);
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \ingroup vector
 * \function igraph_vector_init_int
 * \brief Create an \type igraph_vector_t containing the parameters.
 * 
 * This function is similar to \ref igraph_vector_init_real(), but it expects 
 * \type int parameters. It is important that all parameters
 * should be of this type, otherwise the result of the function call
 * is undefined.
 * \param v Pointer to an uninitialized \type igraph_vector_t object.
 * \param no The number of \type int parameters to follow.
 * \param ... The elements of the vector.
 * \return Error code, \c IGRAPH_ENOMEM if there is
 *    not enough memory.
 * \sa \ref igraph_vector_init_real() and igraph_vector_init_int_end(), these are
 *    similar functions.
 *
 * Time complexity: at least O(n) for 
 * n elements plus the time
 * complexity of the memory allocation.
 */

int igraph_vector_init_int(igraph_vector_t *v, int no, ...) {
  int i=0;
  va_list ap;
  IGRAPH_CHECK(igraph_vector_init(v, no));

  va_start(ap, no);
  for (i=0; i<no; i++) {
    VECTOR(*v)[i]=(real_t) va_arg(ap, int);
  }
  va_end(ap);
  
  return 0;
}

/**
 * \ingroup vector
 * \function igraph_vector_init_int_end
 * \brief Create an \type igraph_vector_t from the parameters.
 * 
 * This constructor is similar to \ref igraph_vector_init_int(), the only
 * difference is that instead of giving the number of elements in the
 * vector, a special marker element follows the last real vector
 * element.
 * \param v Pointer to an uninitialized \type igraph_vector_t object.
 * \param endmark This element will signal the end of the vector. It
 *    will \em not be part of the vector.
 * \param ... The elements of the vector.
 * \return Error code, \c IGRAPH_ENOMEM if there
 *    isn't enough memory.
 * 
 * \sa \ref igraph_vector_init_int() and \ref igraph_vector_init_real_end() for
 * similar functions.
 *
 * Time complexity: at least O(n) for 
 * n elements plus the time
 * complexity of the memory allocation.
 */

int igraph_vector_init_int_end(igraph_vector_t *v, int endmark, ...) {
  int i=0, n=0;
  va_list ap;

  va_start(ap, endmark);
  while (1) {
    int num = va_arg(ap, int);
    if (num == endmark) {
      break;
    }
    n++;
  }
  va_end(ap);

  IGRAPH_VECTOR_INIT_FINALLY(v, n);
  
  va_start(ap, endmark);
  for (i=0; i<n; i++) {
    VECTOR(*v)[i]=(real_t) va_arg(ap, int);
  }
  va_end(ap);
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \ingroup vector
 * \function igraph_vector_destroy
 * \brief Destroys a vector object.
 *
 * All vectors initialized by \ref igraph_vector_init() should be properly
 * destroyed by this function. A destroyed vector needs to be
 * reinitialized by \ref igraph_vector_init(), \ref igraph_vector_init_copy() or
 * another constructor.
 * \param v Pointer to the (previously initialized) vector object to
 *        destroy. 
 *
 * Time complexity: operating system dependent.
 */

void igraph_vector_destroy   (igraph_vector_t* v) {
  assert(v != 0);
  if (v->stor_begin != 0) {
    Free(v->stor_begin);
    v->stor_begin = NULL;
  }
}

/**
 * \ingroup vector
 * \function igraph_vector_reserve
 * \brief Reserves memory for a vector.
 * 
 * \a igraph vectors are flexible, they can grow and
 * shrink. Growing 
 * however occasionally needs the data in the vector to be copyed.
 * In order to avoid you can call this function to reserve space for
 * future growth of the vector. 
 * 
 * Note that this function does \em not change the size of the
 * vector. Let us see a small example to clarify things: if you
 * reserve space for 100 elements and the size of your
 * vector was (and still is) 60, then you can surely add additional 40
 * elements to your vector before it will be copied.
 * \param v The vector object.
 * \param size The new \em allocated size of the vector.
 * \return Error code:
 *         \c IGRPAH_ENOMEM if there is not enough memory.
 *
 * Time complexity: operating system dependent, should be around
 * O(n), n 
 * is the new allocated size of the vector.
 */

int igraph_vector_reserve   (igraph_vector_t* v, long int size) {
	long int actual_size=igraph_vector_size(v);
	real_t *tmp;
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	if (size <= igraph_vector_size(v)) { return 0; }

	tmp=Realloc(v->stor_begin, size, real_t);
	if (tmp==0) {
	  IGRAPH_ERROR("cannot reserve space for vector", IGRAPH_ENOMEM);
	}
	v->stor_begin=tmp;
	v->stor_end=v->stor_begin + size;
	v->end=v->stor_begin+actual_size;
	
	return 0;
}

/**
 * \ingroup vector
 * \function igraph_vector_empty
 * \brief Decides whether the size of the vector is zero.
 *
 * \param v The vector object.
 * \return Non-zero number if the size of the vector is not zero and
 *         zero otherwise.
 * 
 * Time complexity: O(1).
 */

bool_t igraph_vector_empty     (const igraph_vector_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	return v->stor_begin == v->end;
}

/**
 * \ingroup vector
 * \function igraph_vector_size
 * \brief Gives the size (=length) of the vector.
 * 
 * \param v The vector object
 * \return The size of the vector.
 *
 * Time complexity: O(1). 
 */

long int igraph_vector_size      (const igraph_vector_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	return v->end - v->stor_begin;
}

/**
 * \ingroup vector
 * \function igraph_vector_clear
 * \brief Removes all elements from a vector.
 * 
 * This function simply sets the size of the vector to zero, it does
 * not free any allocated memory. For that you have to call
 * \ref igraph_vector_destroy().
 * \param v The vector object.
 * 
 * Time complexity: O(1).
 */

void igraph_vector_clear     (igraph_vector_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	v->end = v->stor_begin;
}

/**
 * \ingroup vector
 * \function igraph_vector_push_back
 * \brief Appends one element to a vector.
 * 
 * This function resizes the vector to be one element longer and
 * sets the very last element in the vector to \p e.
 * \param v The vector object.
 * \param e The element to append to the vector.
 * \return Error code:
 *         \c IGRAPH_ENOMEM: not enough memory.
 * 
 * Time complexity: operating system dependent. What is important that
 * a sequence of n
 * subsequent calls to this function has time complexity
 * O(n), even if there 
 * hadn't been any space reserved for the new elements by
 * \ref igraph_vector_reserve(). This is implemented by a trick similar to the C++
 * \type vector class: each time more memory is allocated for a
 * vector, the size of the additionally allocated memory is the same
 * as the vector's current length. (We assume here that the time
 * complexity of memory allocation is at most linear.)
 */

int igraph_vector_push_back (igraph_vector_t* v, real_t e) {
  	assert(v != NULL);
	assert(v->stor_begin != NULL);
	
	/* full, allocate more storage */
	if (v->stor_end == v->end) {
		long int new_size = igraph_vector_size(v) * 2;
		if (new_size == 0) { new_size = 1; }
		IGRAPH_CHECK(igraph_vector_reserve(v, new_size));
	}
	
	*(v->end) = e;
	v->end += 1;
	
	return 0;
}

/**
 * \ingroup vector
 * \section igraph_vector_accessing_elements Accessing elements of a
 * \type igraph_vector_t.
 * 
 * <para>The simplest way to access an element of a vector is to use the
 * \ref VECTOR macro. This macro can be used both for querying and setting
 * \type igraph_vector_t elements. If you need a function, \ref
 * igraph_vector_e() queries and \ref igraph_vector_set() sets an element of a
 * vector. \ref igraph_vector_e_ptr() returns the address of an element.</para>
 * 
 * <para>\ref igraph_vector_tail() returns the last element of a non-empty
 * vector. There is no <function>igraph_vector_head()</function> function
 * however, as it is easy to write <code>VECTOR(v)[0]</code>
 * instead.</para>
 */

/**
 * \ingroup vector
 * \function igraph_vector_e
 * \brief Access an element of a vector.
 * \param v The \type igraph_vector_t object.
 * \param pos The position of the element, the index of the first
 *    element is zero.
 * \return The desired element.
 * \sa \ref igraph_vector_e_ptr() and the \ref VECTOR macro.
 * 
 * Time complexity: O(1).
 */

real_t igraph_vector_e         (const igraph_vector_t* v, long int pos) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	return * (v->stor_begin + pos);
}

/**
 * \ingroup vector
 * \function igraph_vector_e_ptr
 * \brief Get the address of an element of a vector
 * \param v The \type igraph_vector_t object.
 * \param pos The position of the element, the position of the first
 *   element is zero.
 * \return Pointer to the desired element.
 * \sa \ref igraph_vector_e() and the \ref VECTOR macro.
 * 
 * Time complexity: O(1).
 */

real_t*igraph_vector_e_ptr  (const igraph_vector_t* v, long int pos) {
  assert(v!=NULL);
  assert(v->stor_begin != NULL);
  return v->stor_begin+pos;
}

/**
 * \ingroup vector
 * \function igraph_vector_set
 * \brief Assignment to an element of a vector.
 * \param v The \type igraph_vector_t element.
 * \param pos Position of the element to set.
 * \param value New value of the element.
 * \sa \ref igraph_vector_e().
 */

void igraph_vector_set       (igraph_vector_t* v, long int pos, real_t value) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);	
	*(v->stor_begin + pos) = value;
}

/**
 * \ingroup vector
 * \function igraph_vector_null
 * \brief Sets each element in the vector to zero.
 * 
 * Note that \ref igraph_vector_init() sets the elements to zero as well, so
 * it makes no sense to call this function on a just initialized
 * vector. 
 * \param v The vector object.
 *
 * Time complexity: O(n), the size of
 * the vector. 
 */

void igraph_vector_null      (igraph_vector_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	if (igraph_vector_size(v)>0) {
		memset(v->stor_begin, 0, sizeof(real_t)*igraph_vector_size(v));
	}
}

/**
 * \ingroup vector
 * \function igraph_vector_tail
 * \brief Returns the last element in a vector.
 *
 * It is an error to call this function on an empty vector, the result
 * is undefined.
 * \param v The vector object.
 * \return The last element.
 * 
 * Time complexity: O(1).
 */

real_t igraph_vector_tail(const igraph_vector_t *v) {
  assert(v!=NULL);
  assert(v->stor_begin != NULL);
  return *((v->end)-1);
}

/**
 * \ingroup vector
 * \function igraph_vector_pop_back
 * \brief Removes and returns the last element of a vector.
 *
 * It is an error to call this function with an empty vector.
 * \param v The vector object.
 * \return The removed last element.
 * 
 * Time complexity: O(1).
 */

real_t igraph_vector_pop_back(igraph_vector_t* v) {
  real_t tmp;
  assert(v!=NULL);
  assert(v->stor_begin != NULL);
  assert(v->end != v->stor_begin);
  tmp=igraph_vector_e(v, igraph_vector_size(v)-1);
  v->end -= 1;
  return tmp;
}

/**
 * \ingroup vector
 * \function igraph_vector_order
 * \brief Calculate the order of the elements in a vector.
 *
 * The smallest element will have order zero, the second smallest
 * order one, etc. 
 * \param v The original \type igraph_vector_t object.
 * \param res An initialized \type igraph_vector_t object, it will be
 *    resized to match the size of \p v. The
 *    result of the computation will be stored here.
 * \param nodes Hint, the largest element in \p v.
 * \return Error code:
 *         \c IGRAPH_ENOMEM: out of memory
 * 
 * Time complexity: O()
 */

int igraph_vector_order(const igraph_vector_t* v, igraph_vector_t* res, integer_t nodes) {
  long int edges=igraph_vector_size(v);
  igraph_vector_t ptr;
  igraph_vector_t rad;
  long int i, j;

  assert(v!=NULL);
  assert(v->stor_begin != NULL);

  IGRAPH_VECTOR_INIT_FINALLY(&ptr, nodes+1);
  IGRAPH_VECTOR_INIT_FINALLY(&rad, edges);
  IGRAPH_CHECK(igraph_vector_resize(res, edges));
  
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
  
  igraph_vector_destroy(&ptr);
  igraph_vector_destroy(&rad);
  IGRAPH_FINALLY_CLEAN(2);
  
  return 0;
}

int igraph_vector_order2(igraph_vector_t *v) {

  igraph_indheap_t heap;

  igraph_indheap_init_array(&heap, VECTOR(*v), igraph_vector_size(v));
  IGRAPH_FINALLY(igraph_indheap_destroy, &heap);

  igraph_vector_clear(v);
  while (!igraph_indheap_empty(&heap)) {
    IGRAPH_CHECK(igraph_vector_push_back(v, 
					 igraph_indheap_max_index(&heap)-1));
    igraph_indheap_delete_max(&heap);
  }
  
  igraph_indheap_destroy(&heap);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \ingroup vector
 * \function igraph_vector_sort_cmp
 * \brief Internal comparision function of vector elements, used by 
 * \ref igraph_vector_sort().
 */

int igraph_vector_sort_cmp(const void *a, const void *b) {
  const real_t *da = (const real_t *) a;
  const real_t *db = (const real_t *) b;

  return (*da > *db) - (*da < *db);
}

/**
 * \ingroup vector
 * \function igraph_vector_sort
 * \brief Sorts the elements of the vector into ascending order.
 * 
 * This function uses the built-in sort function of the C library.
 * \param v Pointer to an initialized vector object.
 *
 * Time complexity: should be
 * O(nlogn) for
 * n 
 * elements.
 */

void igraph_vector_sort(igraph_vector_t *v) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  qsort(v->stor_begin, igraph_vector_size(v), sizeof(real_t), igraph_vector_sort_cmp);
}

/**
 * \ingroup vector
 * \function igraph_vector_resize
 * \brief Resize the vector.
 *
 * Note that this function does not free any memory, just sets the
 * size of the vector to the given one. It can on the other hand 
 * allocate more memory if the new size is larger than the previous
 * one. In this case the newly appeared elements in the vector are
 * \em not set to zero, they are uninitialized.
 * \param v The vector object
 * \param newsize The new size of the vector.
 * \return Error code, 
 *         \c IGRAPH_ENOMEM if there is not enough
 *         memory. Note that this function \em never returns an error
 *         if the vector is made smaller.
 * \sa \ref igraph_vector_reserve() for allocating memory for future
 * extensions of a vector.
 * 
 * Time complexity: O(1) if the new
 * size is smaller, operating system dependent if it is larger. In the
 * latter case it is usually around
 * O(n),
 * n is the new size of the vector. 
 */

int igraph_vector_resize(igraph_vector_t* v, long int newsize) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  IGRAPH_CHECK(igraph_vector_reserve(v, newsize));
  v->end = v->stor_begin+newsize;
  return 0;
}

/**
 * \ingroup vector
 * \function igraph_vector_max
 * \brief Gives the maximum element of the vector.
 *
 * If the size of the vector is zero, an arbitrary number is
 * returned.
 * \param v The vector object.
 * \return The maximum element.
 *
 * Time complexity: O(n),
 * n is the size of the vector. 
 */

real_t igraph_vector_max(const igraph_vector_t* v) {
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
 * \function igraph_vector_which_max
 * \brief Gives the position of the maximum element of the vector.
 *
 * If the size of the vector is zero, -1 is 
 * returned.
 * \param v The vector object.
 * \return The position of the first maximum element.
 *
 * Time complexity: O(n),
 * n is the size of the vector. 
 */

long int igraph_vector_which_max(const igraph_vector_t* v) {
  long int which=-1;
  if (!igraph_vector_empty(v)) {
    real_t max;
    real_t *ptr;
    long int pos;
    assert(v != NULL);
    assert(v->stor_begin != NULL);
    max=*(v->stor_begin); which=0;
    ptr=v->stor_begin+1; pos=1;
    while (ptr < v->end) {
      if ((*ptr) > max) {
	max=*ptr;
	which=pos;
      }
      ptr++; pos++;
    }
  }
  return which;
}

/**
 * \ingroup vector
 * \function igraph_vector_init_copy
 * \brief Initializes a vector from an ordinary C array (constructor).
 * 
 * \param v Pointer to an uninitialized vector object.
 * \param data A regular C array.
 * \param length The length of the C array.
 * \return Error code: 
 *         \c IGRAPH_ENOMEM if there is not enough memory.
 * 
 * Time complexity: operating system specific, usually
 * O(\p length).
 */

int igraph_vector_init_copy(igraph_vector_t *v, real_t *data, long int length) {
  v->stor_begin=Calloc(length, real_t);
  if (v->stor_begin==0) {
    IGRAPH_ERROR("cannot init vector from array", IGRAPH_ENOMEM);
  }
  v->stor_end=v->stor_begin+length;
  v->end=v->stor_end;
  memcpy(v->stor_begin, data, length*sizeof(real_t));
  
  return 0;
}

/**
 * \ingroup vector
 * \function igraph_vector_copy_to
 * \brief Copies the contents of a vector to a C array.
 * 
 * The C array should have sufficient length.
 * \param v The vector object.
 * \param to The C array.
 * 
 * Time complexity: O(n),
 * n is the size of the vector.
 */

void igraph_vector_copy_to(const igraph_vector_t *v, real_t* to) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  if (v->end != v->stor_begin) {
    memcpy(to, v->stor_begin, sizeof(real_t) * (v->end - v->stor_begin));
  }
}

/**
 * \ingroup vector
 * \function igraph_vector_copy
 * \brief Initializes a vector from another vector object (constructor).
 * 
 * The contents of the existing vector object will be copied to
 * the new one.
 * \param to Pointer to a not yet initialized vector object.
 * \param from The original vector object to copy.
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory.
 * 
 * Time complexity: operating system dependent, usually
 * O(n),
 * n is the size of the vector. 
 */

int igraph_vector_copy(igraph_vector_t *to, const igraph_vector_t *from) {
  assert(from != NULL);
  assert(from->stor_begin != NULL);
  to->stor_begin=Calloc(igraph_vector_size(from), real_t);
  if (to->stor_begin==0) {
    IGRAPH_ERROR("canot copy vector", IGRAPH_ENOMEM);
  }
  to->stor_end=to->stor_begin+igraph_vector_size(from);
  to->end=to->stor_end;
  memcpy(to->stor_begin, from->stor_begin, igraph_vector_size(from)*sizeof(real_t));
  
  return 0;
}

/**
 * \ingroup vector
 * \function igraph_vector_sum
 * \brief Calculates the sum of the elements in the vector.
 *
 * For the empty vector 0.0 is returned.
 * \param v The vector object.
 * \return The sum of the elements.
 * 
 * 
 * Time complexity: O(n), the size of
 * the vector. 
 */

real_t igraph_vector_sum(const igraph_vector_t *v) {
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
 * \function igraph_vector_prod
 * \brief Calculates the product of the elements in the vector.
 * 
 * For the empty vector one (1) is returned.
 * \param v The vector object.
 * \return The product of the elements.
 * 
 * Time complexity: O(n), the size of
 * the vector. 
 */

real_t igraph_vector_prod(const igraph_vector_t *v) {
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
 * \function igraph_vector_init_seq
 * \brief Initializes a vector with a sequence.
 * 
 * The vector will contain the numbers \p from,
 * \p from+1, ..., \p to.
 * \param v Pointer to an uninitialized vector object.
 * \param from The lower limit in the sequence (inclusive).
 * \param to The upper limit in the sequence (inclusive).
 * \return Error code:
 *         \c IGRAPH_ENOMEM: out of memory.
 *
 * Time complexity: O(n), the number
 * of elements in the vector. 
 */

int igraph_vector_init_seq(igraph_vector_t *v, real_t from, real_t to) {
  real_t *p;
  IGRAPH_CHECK(igraph_vector_init(v, to-from+1));

  for (p=v->stor_begin; p<v->end; p++) {
    *p = from++;
  }
  
  return 0;
}

/**
 * \ingroup vector
 * \function igraph_vector_remove_section
 * \brief Deletes a section from a vector.
 * 
 * Note that this function does not do range checking. The result is
 * undefined if you supply invalid limits.
 * \param v The vector object.
 * \param from The position of the first element to remove.
 * \param to The position of the first element \em not to remove.
 *
 * Time complexity: O(n-from),
 * n is the number of elements in the
 * vector. 
 */

void igraph_vector_remove_section(igraph_vector_t *v, long int from, long int to) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  memmove(v->stor_begin+from, v->stor_begin+to,
	  sizeof(real_t)*(v->end-v->stor_begin-to));
  v->end -= (to-from);
}

/**
 * \ingroup vector
 * \function igraph_vector_remove
 * \brief Removes a single element from a vector.
 *
 * Note that this function does not do range checking.
 * \param v The vector object.
 * \param elem The position of the element to remove.
 * 
 * Time complexity: O(n-elem),
 * n is the number of elements in the
 * vector. 
 */

void igraph_vector_remove(igraph_vector_t *v, long int elem) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  igraph_vector_remove_section(v, elem, elem+1);
}

/**
 * \ingroup vector
 * \function igraph_vector_move_interval
 * \brief Copies a section of a vector.
 *
 * The result of this function is undefined if the source and target
 * intervals overlap.
 * \param v The vector object.
 * \param begin The position of the first element to move.
 * \param end The position of the first element \em not to move.
 * \param to The target position.
 * \return Error code, the current implementation always returns with
 *    success. 
 *
 * Time complexity: O(end-begin).
 */

int igraph_vector_move_interval(igraph_vector_t *v, long int begin, long int end, 
			 long int to) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  memcpy(v->stor_begin+to, v->stor_begin+begin, 
	 sizeof(real_t)*(end-begin));

  return 0;
}

/**
 * \ingroup vector
 * \function igraph_vector_permdelete
 * \brief Remove elements of a vector (for internal use).
 */

void igraph_vector_permdelete(igraph_vector_t *v, long int *index, long int nremove) {
  long int i;
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  for (i=0; i<igraph_vector_size(v); i++) {
    if (index[i] != 0) {
      VECTOR(*v)[ index[i]-1 ] = VECTOR(*v)[i];
    }
  }
  v->end -= nremove;
}

/**
 * \ingroup vector
 * \function igraph_vector_remove_negidx
 * \brief Remove elements of a vector (for internal use).
 */

void igraph_vector_remove_negidx(igraph_vector_t *v, const igraph_vector_t *neg, long int nremove) {
  long int i, idx=0;
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  for (i=0; i<igraph_vector_size(v); i++) {
    VECTOR(*v)[idx++] = VECTOR(*v)[i];
  }
  v->end -= nremove;
}

/**
 * \ingroup vector
 * \function igraph_vector_isininterval
 * \brief Checks if all elements of a vector are in the given
 * interval.
 * 
 * \param v The vector object.
 * \param low The lower limit of the interval (inclusive).
 * \param high The higher limit of the interval (inclusive).
 * \return True (positive integer) if all vector elements are in the
 *   interval, false (zero) otherwise.
 *
 * Time complexity: O(n), the number
 * of elements in the vector.
 */

bool_t igraph_vector_isininterval(const igraph_vector_t *v, real_t low, real_t high) {
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
 * \function igraph_vector_any_smaller
 * \brief Checks if any element of a vector is smaller than a limit.
 * 
 * \param v The \type igraph_vector_t object.
 * \param limit The limit.
 * \return True (positive integer) if the vector contains at least one
 *   smaller element than \p limit, false (zero)
 *   otherwise. 
 * 
 * Time complexity: O(n), the number
 * of elements in the vector.
 */

bool_t igraph_vector_any_smaller(const igraph_vector_t *v, real_t limit) {
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
 * \ingroup vector
 * \function igraph_vector_is_equal
 * \brief Decides whether two vectors contain exactly the same elements
 * (in the same order).
 * 
 * \param lhs The first vector.
 * \param rhs The second vector.
 * \return Positive integer if the two vectors are equal element by
 * element or zero if they are not.
 * 
 * Time complexity: O(n), the length
 * of the vectors.
 */

bool_t igraph_vector_is_equal(const igraph_vector_t *lhs, const igraph_vector_t *rhs) {
  long int i, s;
  assert(lhs != 0);
  assert(rhs != 0);
  assert(lhs->stor_begin != 0);
  assert(rhs->stor_begin != 0);
  
  s=igraph_vector_size(lhs);
  if (s != igraph_vector_size(rhs)) {
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

/**
 * \ingroup vector
 * \function igraph_vector_binsearch
 * \brief Finds an element by binary searching a sorted vector.
 * 
 * It is assumed that the vector is sorted. If the specified element
 * (\p what) is not in the vector, then the
 * position of where it should be inserted (to keep the vector sorted)
 * is returned.
 * \param v The \type igraph_vector_t object.
 * \param what The element to search for.
 * \param pos Pointer to a \type long int. This is set to the
 *   position of an instance of \p what in the
 *   vector if it is present. If \p v does not
 *   contain \p what then
 *   \p pos is set to the position to which it
 *   should be inserted (to keep the the vector sorted of course).
 * \return Positive integer (true) if \p what is
 *   found in the vector, zero (false) otherwise.
 * 
 * Time complexity: O(log(n)),
 * n is the number of elements in
 * \p v.
 */

bool_t igraph_vector_binsearch(const igraph_vector_t *v, real_t what, long int *pos) {
  long int left=0;
  long int right=igraph_vector_size(v)-1;

  while (left < right-1) {
    long int middle=(left+right)/2;
    if (VECTOR(*v)[middle] > what) {
      right=middle;
    } else if (VECTOR(*v)[middle] < what) {
      left=middle;
    } else {
      left=middle;
      break;
    }
  }

  if (VECTOR(*v)[left] != what && VECTOR(*v)[right]==what) {
    left=right;
  }
  
  if (pos != 0) {
    *pos=left;
  }
  return VECTOR(*v)[left]==what;
}

/**
 * \function igraph_vector_multiply
 * \brief Multiply all elements of a vector by a constant
 * 
 * \param v The vector.
 * \param by The constant.
 * \return Error code. The current implementation always returns with success.
 * 
 * Added in version 0.2.</para><para>
 * 
 * Time complexity: O(n), the number of elements in a vector.
 */

void igraph_vector_multiply(igraph_vector_t *v, real_t by) {
  long int i;
  for (i=0; i<igraph_vector_size(v); i++) {
    VECTOR(*v)[i] *= by;
  }
}

bool_t igraph_vector_search(igraph_vector_t *v, long int from, real_t what, 
			    long int *pos) {
  long int i, n=igraph_vector_size(v);  
  for (i=from; i<n; i++) {
    if (VECTOR(*v)[i]==what) break;
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

int igraph_vector_filter_smaller(igraph_vector_t *v, real_t elem) {
  long int i=0, n=igraph_vector_size(v);
  long int s;
  while (i<n && VECTOR(*v)[i]<elem) {
    i++;
  }
  s=i;
  
  while (s<n && VECTOR(*v)[s]==elem) {
    s++;
  }
  
  igraph_vector_remove_section(v, 0, i+(s-i)/2);
  return 0;
}

int igraph_vector_append(igraph_vector_t *to, const igraph_vector_t *from) {
  long int tosize, fromsize;
  
  tosize=igraph_vector_size(to);
  fromsize=igraph_vector_size(from);
  IGRAPH_CHECK(igraph_vector_resize(to, tosize+fromsize));
  memcpy(to->stor_begin+tosize, from->stor_begin, 
	 sizeof(real_t)*fromsize);
  to->end=to->stor_begin+tosize+fromsize;
  
  return 0;
}
