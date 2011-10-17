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
#include "config.h"

#include <assert.h>
#include <string.h> 		/* memcpy & co. */
#include <stdlib.h>

/**
 * \section about_igraph_vector_ptr_objects Pointer vectors
 * (<type>igraph_vector_ptr_t</type>)
 * 
 * <para>The \type igraph_vector_ptr_t data type is very similar to
 * the \type igraph_vector_t type, but it stores generic pointers instead of
 * real numbers.</para>
 * 
 * <para>This type has the same space complexity as \type
 * igraph_vector_t, and most implemented operations work the same way
 * as for \type igraph_vector_t. </para>
 * 
 * <para>This type is mostly used to pass to or receive from a set of
 * graphs to some \a igraph functions, such as \ref
 * igraph_decompose(), which decomposes a graph to connected
 * components.</para>
 * 
 * <para>The same \ref VECTOR macro used for ordinary vectors can be
 * used for pointer vectors as well, please note that a typeless
 * generic pointer will be provided by this macro and you may need to
 * cast it to a specific pointer before starting to work with it.</para>
 */


/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_init
 * \brief Initialize a pointer vector (constructor).
 *
 * </para><para>
 * This is the constructor of the pointer vector data type. All
 * pointer vectors constructed this way should be destroyed via
 * calling \ref igraph_vector_ptr_destroy().
 * \param v Pointer to an uninitialized
 *        <type>igraph_vector_ptr_t</type> object, to be created.
 * \param size Integer, the size of the pointer vector.
 * \return Error code:
 *         \c IGRAPH_ENOMEM if out of memory
 * 
 * Time complexity: operating system dependent, the amount of \quote
 * time \endquote required to allocate \p size elements.
 */

int igraph_vector_ptr_init      (igraph_vector_ptr_t* v, int long size) {	
        long int alloc_size= size > 0 ? size : 1;
	assert(v != NULL);
	if (size < 0) { size=0; }
	v->stor_begin=igraph_Calloc(alloc_size, void*);
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
 * \function igraph_vector_ptr_destroy
 * \brief Destroys a pointer vector.
 * 
 * </para><para>
 * The destructor for pointer vectors.
 * \param v Pointer to the pointer vector to destroy.
 * 
 * Time complexity: operating system dependend, the \quote time
 * \endquote required to deallocate O(n) bytes, n is the number of
 * elements allocated for the pointer vector (not neccessarily the
 * number of elements in the vector).
 */

void igraph_vector_ptr_destroy   (igraph_vector_ptr_t* v) {
  assert(v != 0);
  if (v->stor_begin != 0) {
    igraph_Free(v->stor_begin);
    v->stor_begin = NULL;
  }
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_free_all
 * \brief Calls free() on all elements of a pointer vector.
 *
 * \param v Pointer to the pointer vector whose elements will be freed.
 * 
 * Time complexity: operating system dependent, the \quote time
 * \endquote required to deallocate O(n) pointers, each pointing to
 * a memory area of arbitrary size. n is the number of
 * elements in the pointer vector.
 */

void igraph_vector_ptr_free_all   (igraph_vector_ptr_t* v) {
  void **ptr;
  assert(v != 0);
  assert(v->stor_begin != 0);
  for (ptr=v->stor_begin; ptr<v->end; ptr++) {
    igraph_Free(*ptr);
  }
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_destroy_all
 * \brief Calls free() on all elements and destroys the pointer vector.
 *
 * \param v Pointer to the pointer vector to destroy.
 *
 * Time complexity: operating system dependent, the \quote time
 * \endquote required to deallocate O(n) pointers, each pointing to
 * a memory area of arbitrary size, plus the \quote time \endquote
 * required to deallocate O(n) bytes, n being the number of elements
 * allocated for the pointer vector (not necessarily the number of
 * elements in the vector).
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

	tmp=igraph_Realloc(v->stor_begin, size, void*);
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
 * \function igraph_vector_ptr_size
 * \brief Gives the number of elements in the pointer vector.
 * 
 * \param v The pointer vector object.
 * \return The size of the object, ie. the number of pointers stored.
 * 
 * Time complexity: O(1).
 */

long int igraph_vector_ptr_size      (const igraph_vector_ptr_t* v) {
	assert(v != NULL);
/* 	assert(v->stor_begin != NULL);		 */ /* TODO */
	return v->end - v->stor_begin;
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_clear
 * \brief Removes all elements from a pointer vector.
 * 
 * </para><para>
 * This function resizes a pointer to vector to zero length. Note that
 * the pointed objects are \em not deallocated, you should call
 * free() on them, or make sure that their allocated memory is freed
 * in some other way, you'll get memory leaks otherwise.
 * 
 * </para><para>
 * Note that the current implementation of this function does
 * \em not deallocate the memory required for storing the
 * pointers, so making a pointer vector smaller this way does not give
 * back any memory. This behavior might change in the future.
 * \param v The pointer vector to clear.
 * 
 * Time complexity: O(1).
 */

void igraph_vector_ptr_clear     (igraph_vector_ptr_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);		
	v->end = v->stor_begin;	
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_push_back
 * \brief Appends an elements to the back of a pointer vector.
 * 
 * \param v The pointer vector.
 * \param e The new element to include in the pointer vector.
 * \return Error code.
 * \sa igraph_vector_push_back() for the corresponding operation of
 * the ordinary vector type.
 * 
 * Time complexity: O(1) or O(n), n is the number of elements in the
 * vector. The pointer vector implementation ensures that n subsequent
 * push_back operations need O(n) time to complete.
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

void *igraph_vector_ptr_pop_back (igraph_vector_ptr_t *v) {
	void *tmp;
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	assert(v->stor_begin != v->end);
	tmp=*(v->end);
	v->end -= 1;
	  
	return tmp;
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_insert
 * \brief Inserts a single element into a pointer vector.
 *
 * Note that this function does not do range checking. Insertion will shift the
 * elements from the position given to the end of the vector one position to the
 * right, and the new element will be inserted in the empty space created at
 * the given position. The size of the vector will increase by one.
 *
 * \param v The pointer vector object.
 * \param pos The position where the new element is inserted.
 * \param e The inserted element
 */
int igraph_vector_ptr_insert(igraph_vector_ptr_t* v, long int pos, void* e) {
  long int size = igraph_vector_ptr_size(v);
  IGRAPH_CHECK(igraph_vector_ptr_resize(v, size+1));
  if (pos<size) {
    memmove(v->stor_begin+pos+1, v->stor_begin+pos, 
	    sizeof(void*)*(size-pos));
  }
  v->stor_begin[pos] = e;
  return 0;
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_e
 * \brief Access an element of a pointer vector.
 * 
 * \param v Pointer to a pointer vector.
 * \param pos The index of the pointer to return.
 * \return The pointer at \p pos position.
 * 
 * Time complexity: O(1).
 */

void* igraph_vector_ptr_e         (const igraph_vector_ptr_t* v, long int pos) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);	
	return * (v->stor_begin + pos);
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_set
 * \brief Assign to an element of a pointer vector.
 * 
 * \param v Pointer to a pointer vector.
 * \param pos The index of the pointer to update.
 * \param value The new pointer to set in the vector.
 *
 * Time complexity: O(1).
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
 * \function igraph_vector_ptr_resize
 * \brief Resizes a pointer vector.
 * 
 * </para><para>
 * Note that if a vector is made smaller the pointed object are not
 * deallocated by this function.
 * \param v A pointer vector.
 * \param newsize The new size of the pointer vector.
 * \return Error code.
 * 
 * Time complexity: O(1) if the vector if made smaller. Operating
 * system dependent otherwise, the amount of \quote time \endquote
 * needed to allocate the memory for the vector elements.
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
 * \return Error code:
 *         \c IGRAPH_ENOMEM if out of memory
 */

int igraph_vector_ptr_init_copy(igraph_vector_ptr_t *v, void* *data, long int length) {
  v->stor_begin=igraph_Calloc(length, void*);
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
 * \function igraph_vector_ptr_copy
 * \brief Copy a pointer vector (constructor).
 *
 * </para><para>
 * This function creates a pointer vector by copying another one. This
 * is shallow copy, only the pointers in the vector will be copyed.
 * \param to Pointer to an uninitialized pointer vector object.
 * \param from A pointer vector object.
 * \return Error code:
 *         \c IGRAPH_ENOMEM if out of memory
 * 
 * Time complexity: O(n) if allocating memory for n elements can be
 * done in O(n) time.
 */

int igraph_vector_ptr_copy(igraph_vector_ptr_t *to, const igraph_vector_ptr_t *from) {
  assert(from != NULL);
/*   assert(from->stor_begin != NULL); */ /* TODO */
  to->stor_begin=igraph_Calloc(igraph_vector_ptr_size(from), void*);
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
    memmove(v->stor_begin+pos, v->stor_begin+pos+1,
	    sizeof(void*)*(igraph_vector_ptr_size(v)-pos-1));
  }
  v->end--;
}

/**
 * \ingroup vectorptr
 * \brief Sort the pointer vector based on an external comparison function
 *
 * Sometimes it is necessary to sort the pointers in the vector based on
 * the property of the element being referenced by the pointer. This
 * function allows us to sort the vector based on an arbitrary external
 * comparison function which accepts two \c void* pointers \c p1 and \c p2
 * and returns an integer less than, equal to or greater than zero if the
 * first argument is considered to be respectively less than, equal to, or
 * greater than the second. \c p1 and \c p2 will point to the pointer in the
 * vector, so they have to be double-dereferenced if one wants to get access
 * to the underlying object the address of which is stored in \c v .
 */
void igraph_vector_ptr_sort(igraph_vector_ptr_t *v, int (*compar)(const void*, const void*)) {
  qsort(v->stor_begin, igraph_vector_ptr_size(v), sizeof(void*), compar);
}

