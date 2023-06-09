/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_vector_ptr.h"

#include "igraph_types.h"
#include "igraph_memory.h"
#include "igraph_qsort.h"

#include <string.h>         /* memcpy & co. */
#include <stdint.h>         /* uintptr_t */
#include <stdlib.h>

/**
 * \section about_igraph_vector_ptr_objects Pointer vectors
 * (<type>igraph_vector_ptr_t</type>)
 *
 * <para>The \type igraph_vector_ptr_t data type is very similar to
 * the \ref igraph_vector_t type, but it stores generic pointers instead of
 * real numbers.</para>
 *
 * <para>This type has the same space complexity as \ref
 * igraph_vector_t, and most implemented operations work the same way
 * as for \ref igraph_vector_t.</para>
 *
 * <para>The same \ref VECTOR macro used for ordinary vectors can be
 * used for pointer vectors as well, please note that a typeless
 * generic pointer will be provided by this macro and you may need to
 * cast it to a specific pointer before starting to work with it.</para>
 *
 * <para>Pointer vectors may have an associated item destructor function
 * which takes a pointer and returns nothing. The item destructor will
 * be called on each item in the pointer vector when it is destroyed by
 * \ref igraph_vector_ptr_destroy() or \ref igraph_vector_ptr_destroy_all(),
 * or when its elements are freed by \ref igraph_vector_ptr_free_all().
 * Note that the semantics of an item destructor does not coincide with
 * C++ destructors; for instance, when a pointer vector is resized to a
 * smaller size, the extra items will \em not be destroyed automatically!
 * Nevertheless, item destructors may become handy in many cases; for
 * instance, a vector of graphs generated by some function can
 * be destroyed with a single call to \ref igraph_vector_ptr_destroy_all()
 * if the item destructor is set to \ref igraph_destroy().</para>
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

igraph_error_t igraph_vector_ptr_init(igraph_vector_ptr_t* v, igraph_integer_t size) {
    igraph_integer_t alloc_size = size > 0 ? size : 1;
    IGRAPH_ASSERT(v != NULL);
    if (size < 0) {
        size = 0;
    }
    v->stor_begin = IGRAPH_CALLOC(alloc_size, void*);
    if (v->stor_begin == 0) {
        IGRAPH_ERROR("vector ptr init failed", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    v->stor_end = v->stor_begin + alloc_size;
    v->end = v->stor_begin + size;
    v->item_destructor = 0;

    return IGRAPH_SUCCESS;
}

/**
 */

const igraph_vector_ptr_t *igraph_vector_ptr_view(
    const igraph_vector_ptr_t *v, void *const *data, igraph_integer_t length
) {
    igraph_vector_ptr_t *v2 = (igraph_vector_ptr_t*) v;
    v2->stor_begin = (void **)data;
    v2->stor_end = (void**)data + length;
    v2->end = v2->stor_end;
    v2->item_destructor = 0;
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
 * Time complexity: operating system dependent, the \quote time
 * \endquote required to deallocate O(n) bytes, n is the number of
 * elements allocated for the pointer vector (not necessarily the
 * number of elements in the vector).
 */

void igraph_vector_ptr_destroy(igraph_vector_ptr_t* v) {
    IGRAPH_ASSERT(v != 0);
    if (v->stor_begin != 0) {
        IGRAPH_FREE(v->stor_begin);
        v->stor_begin = NULL;
    }
}

static void igraph_i_vector_ptr_call_item_destructor_all(igraph_vector_ptr_t* v) {
    void **ptr;

    if (v->item_destructor != 0) {
        for (ptr = v->stor_begin; ptr < v->end; ptr++) {
            if (*ptr != 0) {
                v->item_destructor(*ptr);
            }
        }
    }
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_free_all
 * \brief Frees all the elements of a pointer vector.
 *
 * If an item destructor is set for this pointer vector, this function will
 * first call the destructor on all elements of the vector and then
 * free all the elements using \ref igraph_free(). If an item destructor is not set,
 * the elements will simply be freed.
 *
 * \param v Pointer to the pointer vector whose elements will be freed.
 *
 * Time complexity: operating system dependent, the \quote time
 * \endquote required to call the destructor n times and then
 * deallocate O(n) pointers, each pointing to a memory area of
 * arbitrary size. n is the number of elements in the pointer vector.
 */

void igraph_vector_ptr_free_all(igraph_vector_ptr_t* v) {
    void **ptr;
    IGRAPH_ASSERT(v != 0);
    IGRAPH_ASSERT(v->stor_begin != 0);

    igraph_i_vector_ptr_call_item_destructor_all(v);
    for (ptr = v->stor_begin; ptr < v->end; ptr++) {
        IGRAPH_FREE(*ptr);
    }
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_destroy_all
 * \brief Frees all the elements and destroys the pointer vector.
 *
 * This function is equivalent to \ref igraph_vector_ptr_free_all()
 * followed by \ref igraph_vector_ptr_destroy().
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

void igraph_vector_ptr_destroy_all(igraph_vector_ptr_t* v) {
    IGRAPH_ASSERT(v != 0);
    IGRAPH_ASSERT(v->stor_begin != 0);
    igraph_vector_ptr_free_all(v);
    igraph_vector_ptr_set_item_destructor(v, 0);
    igraph_vector_ptr_destroy(v);
}

/**
 * \ingroup vectorptr
 * \brief Reserves memory for a pointer vector for later use.
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

igraph_error_t igraph_vector_ptr_reserve(igraph_vector_ptr_t* v, igraph_integer_t capacity) {
    igraph_integer_t actual_size = igraph_vector_ptr_size(v);
    void **tmp;
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);
    IGRAPH_ASSERT(capacity >= 0);

    if (capacity <= igraph_vector_ptr_size(v)) {
        return IGRAPH_SUCCESS;
    }

    tmp = IGRAPH_REALLOC(v->stor_begin, (size_t) capacity, void*);
    IGRAPH_CHECK_OOM(tmp, "Cannot reserve space for pointer vector.");

    v->stor_begin = tmp;
    v->stor_end = v->stor_begin + capacity;
    v->end = v->stor_begin + actual_size;

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vectorptr
 * \brief Decides whether the pointer vector is empty.
 */

igraph_bool_t igraph_vector_ptr_empty(const igraph_vector_ptr_t* v) {
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);
    return v->stor_begin == v->end;
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_size
 * \brief Gives the number of elements in the pointer vector.
 *
 * \param v The pointer vector object.
 * \return The size of the object, i.e. the number of pointers stored.
 *
 * Time complexity: O(1).
 */

igraph_integer_t igraph_vector_ptr_size(const igraph_vector_ptr_t* v) {
    IGRAPH_ASSERT(v != NULL);
    /*  IGRAPH_ASSERT(v->stor_begin != NULL);       */ /* TODO */
    return v->end - v->stor_begin;
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_capacity
 * \brief Returns the allocated capacity of the pointer vector.
 *
 * Note that this might be different from the size of the vector (as
 * queried by \ref igraph_vector_ptr_size()), and specifies how many elements
 * the vector can hold, without reallocation.
 *
 * \param v Pointer to the (previously initialized) pointer vector object to query.
 * \return The allocated capacity.
 *
 * \sa \ref igraph_vector_ptr_size().
 *
 * Time complexity: O(1).
 */

igraph_integer_t igraph_vector_ptr_capacity(const igraph_vector_ptr_t* v) {
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);
    return v->stor_end - v->stor_begin;
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_clear
 * \brief Removes all elements from a pointer vector.
 *
 * </para><para>
 * This function resizes a pointer to vector to zero length. Note that
 * the pointed objects are \em not deallocated, you should call
 * \ref igraph_free() on them, or make sure that their allocated memory is freed
 * in some other way, you'll get memory leaks otherwise. If you have
 * set up an item destructor earlier, the destructor will be called
 * on every element.
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

void igraph_vector_ptr_clear(igraph_vector_ptr_t* v) {
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);
    igraph_i_vector_ptr_call_item_destructor_all(v);
    v->end = v->stor_begin;
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_push_back
 * \brief Appends an element to the back of a pointer vector.
 *
 * \param v The pointer vector.
 * \param e The new element to include in the pointer vector.
 * \return Error code.
 * \sa \ref igraph_vector_push_back() for the corresponding operation of
 * the ordinary vector type.
 *
 * Time complexity: O(1) or O(n), n is the number of elements in the
 * vector. The pointer vector implementation ensures that n subsequent
 * push_back operations need O(n) time to complete.
 */

igraph_error_t igraph_vector_ptr_push_back(igraph_vector_ptr_t* v, void* e) {
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);

    /* full, allocate more storage */
    if (v->stor_end == v->end) {
        igraph_integer_t new_size = igraph_vector_ptr_size(v) * 2;
        if (new_size == 0) {
            new_size = 1;
        }
        IGRAPH_CHECK(igraph_vector_ptr_reserve(v, new_size));
    }

    *(v->end) = e;
    v->end += 1;

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_pop_back
 * \brief Removes and returns the last element of a pointer vector.
 *
 * </para><para>
 * It is an error to call this function with an empty vector.
 *
 * \param v The pointer vector.
 * \return The removed last element.
 *
 * Time complexity: O(1).
 */

void *igraph_vector_ptr_pop_back(igraph_vector_ptr_t *v) {
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);
    IGRAPH_ASSERT(v->stor_begin != v->end);
    v->end -= 1;
    return *(v->end);
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
igraph_error_t igraph_vector_ptr_insert(igraph_vector_ptr_t* v, igraph_integer_t pos, void* e) {
    igraph_integer_t size = igraph_vector_ptr_size(v);
    IGRAPH_CHECK(igraph_vector_ptr_resize(v, size + 1));
    if (pos < size) {
        memmove(v->stor_begin + pos + 1, v->stor_begin + pos,
                sizeof(void*) * (size_t) (size - pos));
    }
    v->stor_begin[pos] = e;
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_get
 * \brief Access an element of a pointer vector.
 *
 * \param v Pointer to a pointer vector.
 * \param pos The index of the pointer to return.
 * \return The pointer at \p pos position.
 *
 * Time complexity: O(1).
 */

void *igraph_vector_ptr_get(const igraph_vector_ptr_t* v, igraph_integer_t pos) {
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);
    return *(v->stor_begin + pos);
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_e
 * \brief Access an element of a pointer vector (deprecated alias).
 *
 * \deprecated-by igraph_vector_ptr_get 0.10.0
 */

void *igraph_vector_ptr_e(const igraph_vector_ptr_t* v, igraph_integer_t pos) {
    return igraph_vector_ptr_get(v, pos);
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

void igraph_vector_ptr_set(igraph_vector_ptr_t* v, igraph_integer_t pos, void* value) {
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);
    *(v->stor_begin + pos) = value;
}

/**
 * \ingroup vectorptr
 * \brief Set all elements of a pointer vector to the NULL pointer.
 */

void igraph_vector_ptr_null(igraph_vector_ptr_t* v) {
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);
    if (igraph_vector_ptr_size(v) > 0) {
        memset(v->stor_begin, 0, sizeof(void*) *
               (size_t) igraph_vector_ptr_size(v));
    }
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_resize
 * \brief Resizes a pointer vector.
 *
 * </para><para>
 * Note that if a vector is made smaller the pointed object are not
 * deallocated by this function and the item destructor is not called
 * on the extra elements.
 *
 * \param v A pointer vector.
 * \param newsize The new size of the pointer vector.
 * \return Error code.
 *
 * Time complexity: O(1) if the vector if made smaller. Operating
 * system dependent otherwise, the amount of \quote time \endquote
 * needed to allocate the memory for the vector elements.
 */

igraph_error_t igraph_vector_ptr_resize(igraph_vector_ptr_t* v, igraph_integer_t newsize) {
    IGRAPH_CHECK(igraph_vector_ptr_reserve(v, newsize));
    v->end = v->stor_begin + newsize;
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vectorptr
 * \brief Initializes a pointer vector from an array (constructor).
 *
 * \param v Pointer to an uninitialized
 *        <type>igraph_vector_ptr_t</type> object to be initialized.
 * \param data The array of pointers that serves as the initial contents of the
 *        pointer vector.
 * \param length Integer, the length of the array.
 * \return Error code:
 *         \c IGRAPH_ENOMEM if out of memory
 */

igraph_error_t igraph_vector_ptr_init_array(igraph_vector_ptr_t *v, void *const *data, igraph_integer_t length) {
    v->stor_begin = IGRAPH_CALLOC(length, void*);
    if (v->stor_begin == 0) {
        IGRAPH_ERROR("Cannot initialize pointer vector from array", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    v->stor_end = v->stor_begin + length;
    v->end = v->stor_end;
    v->item_destructor = 0;
    memcpy(v->stor_begin, data, (size_t) length * sizeof(void*));

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vectorptr
 * \brief Copy the contents of a pointer vector to a regular C array.
 */

void igraph_vector_ptr_copy_to(const igraph_vector_ptr_t *v, void** to) {
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);
    if (v->end != v->stor_begin) {
        memcpy(to, v->stor_begin, sizeof(void*) *
               (size_t) (v->end - v->stor_begin));
    }
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_init_copy
 * \brief Initializes a pointer vector from another one (constructor).
 *
 * </para><para>
 * This function creates a pointer vector by copying another one. This
 * is shallow copy, only the pointers in the vector will be copied.
 *
 * </para><para>
 * It is potentially dangerous to copy a pointer vector with an associated
 * item destructor. The copied vector will inherit the item destructor,
 * which may cause problems when both vectors are destroyed as the items
 * might get destroyed twice. Make sure you know what you are doing when
 * copying a pointer vector with an item destructor, or unset the item
 * destructor on one of the vectors later.
 *
 * \param to Pointer to an uninitialized pointer vector object.
 * \param from A pointer vector object.
 * \return Error code:
 *         \c IGRAPH_ENOMEM if out of memory
 *
 * Time complexity: O(n) if allocating memory for n elements can be
 * done in O(n) time.
 */

igraph_error_t igraph_vector_ptr_init_copy(igraph_vector_ptr_t *to, const igraph_vector_ptr_t *from) {
    igraph_integer_t from_size;

    IGRAPH_ASSERT(from != NULL);
    /*   IGRAPH_ASSERT(from->stor_begin != NULL); */ /* TODO */

    from_size = igraph_vector_ptr_size(from);

    to->stor_begin = IGRAPH_CALLOC(from_size, void*);
    if (to->stor_begin == 0) {
        IGRAPH_ERROR("Cannot copy pointer vector", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    to->stor_end = to->stor_begin + igraph_vector_ptr_size(from);
    to->end = to->stor_end;
    to->item_destructor = from->item_destructor;
    memcpy(to->stor_begin, from->stor_begin,
           (size_t) igraph_vector_ptr_size(from)*sizeof(void*));

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_copy
 * \brief Initializes a pointer vector from another one (deprecated alias).
 *
 * \deprecated-by igraph_vector_ptr_init_copy 0.10
 */

igraph_error_t igraph_vector_ptr_copy(igraph_vector_ptr_t *to, const igraph_vector_ptr_t *from) {
    return igraph_vector_ptr_init_copy(to, from);
}

/**
 * \ingroup vectorptr
 * \brief Remove an element from a pointer vector.
 */

void igraph_vector_ptr_remove(igraph_vector_ptr_t *v, igraph_integer_t pos) {
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);
    if (pos + 1 < igraph_vector_ptr_size(v)) { /* No need to move data when removing the last element. */
        memmove(v->stor_begin + pos, v->stor_begin + pos + 1,
                sizeof(void*) * (size_t) (igraph_vector_ptr_size(v) - pos - 1));
    }
    v->end--;
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_sort
 * \brief Sorts the pointer vector based on an external comparison function.
 *
 * Sometimes it is necessary to sort the pointers in the vector based on
 * the property of the element being referenced by the pointer. This
 * function allows us to sort the vector based on an arbitrary external
 * comparison function which accepts two <type>void *</type> pointers \c p1 and \c p2
 * and returns an integer less than, equal to or greater than zero if the
 * first argument is considered to be respectively less than, equal to, or
 * greater than the second. \c p1 and \c p2 will point to the pointer in the
 * vector, so they have to be double-dereferenced if one wants to get access
 * to the underlying object the address of which is stored in \c v.
 *
 * \param v The pointer vector to be sorted.
 * \param compar A qsort-compatible comparison function. It must take pointers to the
 *    elements of the pointer vector. For example, if the pointer vector contains
 *    <code>igraph_vector_t *</code> pointers, then the comparison function must
 *    interpret its arguments as <code>igraph_vector_t **</code>.
 */
void igraph_vector_ptr_sort(igraph_vector_ptr_t *v, int (*compar)(const void*, const void*)) {
    igraph_qsort(v->stor_begin, (size_t) igraph_vector_ptr_size(v), sizeof(void*),
          compar);
}

igraph_error_t igraph_vector_ptr_append(igraph_vector_ptr_t *to, const igraph_vector_ptr_t *from) {
    igraph_integer_t origsize = igraph_vector_ptr_size(to);
    igraph_integer_t othersize = igraph_vector_ptr_size(from);
    igraph_integer_t i;

    IGRAPH_CHECK(igraph_vector_ptr_resize(to, origsize + othersize));
    for (i = 0; i < othersize; i++, origsize++) {
        to->stor_begin[origsize] = from->stor_begin[i];
    }

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_set_item_destructor
 * \brief Sets the item destructor for this pointer vector.
 *
 * The item destructor is a function which will be called on every non-null
 * pointer stored in this vector when \ref igraph_vector_ptr_destroy(),
 * igraph_vector_ptr_destroy_all() or \ref igraph_vector_ptr_free_all()
 * is called.
 *
 * \return The old item destructor.
 *
 * Time complexity: O(1).
 */
igraph_finally_func_t* igraph_vector_ptr_set_item_destructor(
    igraph_vector_ptr_t *v, igraph_finally_func_t *func) {
    igraph_finally_func_t* result = v->item_destructor;

    v->item_destructor = func;

    return result;
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_get_item_destructor
 * \brief Gets the current item destructor for this pointer vector.
 *
 * The item destructor is a function which will be called on every non-null
 * pointer stored in this vector when \ref igraph_vector_ptr_destroy(),
 * igraph_vector_ptr_destroy_all() or \ref igraph_vector_ptr_free_all()
 * is called.
 *
 * \return The current item destructor.
 *
 * Time complexity: O(1).
 */
igraph_finally_func_t* igraph_vector_ptr_get_item_destructor(const igraph_vector_ptr_t *v) {
    IGRAPH_ASSERT(v != 0);
    return v->item_destructor;
}

typedef int cmp_t (const void *, const void *);

/**
 * Comparison function passed to qsort_r from  igraph_vector_ptr_sort_ind
 */
static int igraph_vector_ptr_i_sort_ind_cmp(void *thunk, const void *p1, const void *p2) {
    cmp_t *cmp = (cmp_t *) thunk;
    uintptr_t *pa = (uintptr_t*) p1;
    uintptr_t *pb = (uintptr_t*) p2;
    void **item_a_ptr = (void**) *pa;
    void **item_b_ptr = (void**) *pb;
    return cmp(*item_a_ptr, *item_b_ptr);
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_sort_ind
 * \brief Returns a permutation of indices that sorts a vector of pointers.
 *
 * Takes an unsorted array \c v as input and computes an array of
 * indices inds such that v[ inds[i] ], with i increasing from 0, is
 * an ordered array (either ascending or descending, depending on
 * \v order). The order of indices for identical elements is not
 * defined.
 *
 * \param v the array to be sorted
 * \param inds the output array of indices. This must be initialized,
 *         but will be resized
 * \param cmp a comparator function that takes two elements of the pointer
 *        vector being sorted (these are constant pointers on their own)
 *        and returns a negative value if the item \em "pointed to" by the
 *        first pointer is smaller than the item \em "pointed to" by the
 *        second pointer, a positive value if it is larger, or zero if the
 *        two items are equal
 * \return Error code.
 *
 * This routine uses the C library qsort routine.
 * Algorithm: 1) create an array of pointers to the elements of v. 2)
 * Pass this array to qsort. 3) after sorting the difference between
 * the pointer value and the first pointer value gives its original
 * position in the array. Use this to set the values of inds.
 */

igraph_error_t igraph_vector_ptr_sort_ind(igraph_vector_ptr_t *v,
        igraph_vector_int_t *inds, cmp_t *cmp) {
    igraph_integer_t i;
    uintptr_t *vind, first;
    igraph_integer_t n = igraph_vector_ptr_size(v);

    IGRAPH_CHECK(igraph_vector_int_resize(inds, n));
    if (n == 0) {
        return IGRAPH_SUCCESS;
    }

    vind = IGRAPH_CALLOC(n, uintptr_t);
    if (vind == 0) {
        IGRAPH_ERROR("igraph_vector_ptr_sort_ind failed", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }

    for (i = 0; i < n; i++) {
        vind[i] = (uintptr_t) &VECTOR(*v)[i];
    }

    first = vind[0];

    igraph_qsort_r(vind, n, sizeof(vind[0]), (void*)cmp, igraph_vector_ptr_i_sort_ind_cmp);

    for (i = 0; i < n; i++) {
        VECTOR(*inds)[i] = (vind[i] - first) / sizeof(uintptr_t);
    }

    IGRAPH_FREE(vind);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup vectorptr
 * \function igraph_vector_ptr_permute
 * \brief Permutes the elements of a pointer vector in place according to an index vector.
 *
 * </para><para>
 * This function takes a vector \c v and a corresponding index vector \c ind,
 * and permutes the elements of \c v such that \c v[ind[i]] is moved to become
 * \c v[i] after the function is executed.
 *
 * </para><para>
 * It is an error to call this function with an index vector that does not
 * represent a valid permutation. Each element in the index vector must be
 * between 0 and the length of the vector minus one (inclusive), and each such
 * element must appear only once. The function does not attempt to validate the
 * index vector.
 *
 * </para><para>
 * The index vector that this function takes is compatible with the index vector
 * returned from \ref igraph_vector_ptr_sort_ind(); passing in the index vector
 * from \ref igraph_vector_ptr_sort_ind() will sort the original vector.
 *
 * </para><para>
 * As a special case, this function allows the index vector to be \em shorter
 * than the vector being permuted, in which case the elements whose indices do
 * not occur in the index vector will be removed from the vector.
 *
 * \param v    the vector to permute
 * \param ind  the index vector
 *
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory.
 *
 * Time complexity: O(n), the size of the vector.
 */
igraph_error_t igraph_vector_ptr_permute(igraph_vector_ptr_t* v, const igraph_vector_int_t* index) {
    IGRAPH_ASSERT(v != NULL);
    IGRAPH_ASSERT(v->stor_begin != NULL);
    IGRAPH_ASSERT(index != NULL);
    IGRAPH_ASSERT(index->stor_begin != NULL);
    IGRAPH_ASSERT(igraph_vector_ptr_size(v) >= igraph_vector_int_size(index));

    igraph_vector_ptr_t v_copy;
    void** v_ptr;
    igraph_integer_t *ind_ptr;

    /* There is a more space-efficient algorithm that needs O(1) space only,
     * but it messes up the index vector, which we don't want */

    IGRAPH_CHECK(igraph_vector_ptr_init(&v_copy, igraph_vector_int_size(index)));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &v_copy);

    for (
        v_ptr = v_copy.stor_begin, ind_ptr = index->stor_begin;
        ind_ptr < index->end;
        v_ptr++, ind_ptr++
    ) {
        *v_ptr = VECTOR(*v)[*ind_ptr];
    }

    IGRAPH_CHECK(igraph_vector_ptr_resize(v, igraph_vector_int_size(index)));
    igraph_vector_ptr_copy_to(&v_copy, VECTOR(*v));

    igraph_vector_ptr_destroy(&v_copy);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_vector_ptr_index(
    const igraph_vector_ptr_t *v, igraph_vector_ptr_t *newv,
    const igraph_vector_int_t *idx
) {
    igraph_integer_t i, j, newlen = igraph_vector_int_size(idx);
    IGRAPH_CHECK(igraph_vector_ptr_resize(newv, newlen));

    for (i = 0; i < newlen; i++) {
        j = VECTOR(*idx)[i];
        VECTOR(*newv)[i] = VECTOR(*v)[j];
    }

    return IGRAPH_SUCCESS;
}
