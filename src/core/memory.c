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

#include "igraph_memory.h"

/**
 * \section about-alloc-funcs About allocation functions
 *
 * Some igraph functions return a pointer vector (igraph_vector_ptr_t)
 * containing pointers to other igraph or other data types. These data
 * types are dynamically allocated and have to be deallocated
 * manually when the user does not need them any more. \c igraph_vector_ptr_t
 * has functions to deallocate the contained pointers on its own, but in this
 * case it has to be ensured that these pointers are allocated by a function
 * that corresponding to the deallocator function that igraph uses.
 *
 * </para><para>
 * To this end, igraph exports the memory allocation functions that are used
 * internally so the user of the library can ensure that the proper functions
 * are used when pointers are moved between the code written by the user and
 * the code of the igraph library.
 *
 * </para><para>
 * Additionally, the memory allocator functions used by igraph work around the
 * quirk of classical \c malloc(), \c realloc() and \c calloc() implementations
 * where the behaviour of allocating zero bytes is undefined. igraph allocator
 * functions will always allocate at least one byte.
 */

/**
 * \function igraph_free
 * \brief Deallocate memory that was allocated by igraph functions.
 *
 * This function exposes the \c free() function used internally by igraph.
 *
 * \param ptr Pointer to the piece of memory to be deallocated.
 *
 * Time complexity: platform dependent, ideally it should be O(1).
 *
 * \sa \ref igraph_calloc(), \ref igraph_malloc(), \ref igraph_realloc()
 */

void igraph_free(void *ptr) {
    IGRAPH_FREE(ptr);
}


/**
 * \function igraph_calloc
 * \brief Allocate memory that can be safely deallocated by igraph functions.
 *
 * This function behaves like \c calloc(), but it ensures that at least one
 * byte is allocated even when the caller asks for zero bytes.
 *
 * \param count Number of items to be allocated.
 * \param size Size of a single item to be allocated.
 * \return Pointer to the piece of allocated memory; \c NULL if the allocation
 * failed.
 *
 * \sa \ref igraph_malloc(), \ref igraph_realloc(), \ref igraph_free()
 */

void *igraph_calloc(size_t count, size_t size) {
    return (void *) IGRAPH_CALLOC(count * size, char);
}


/**
 * \function igraph_malloc
 * \brief Allocate memory that can be safely deallocated by igraph functions.
 *
 * This function behaves like \c malloc(), but it ensures that at least one
 * byte is allocated even when the caller asks for zero bytes.
 *
 * \param size Number of bytes to be allocated. Zero is treated as one byte.
 * \return Pointer to the piece of allocated memory; \c NULL if the allocation
 * failed.
 *
 * \sa \ref igraph_calloc(), \ref igraph_realloc(), \ref igraph_free()
 */

void *igraph_malloc(size_t size) {
    return IGRAPH_MALLOC(size);
}


/**
 * \function igraph_realloc
 * \brief Reallocate memory that can be safely deallocated by igraph functions.
 *
 * This function behaves like \c realloc(), but it ensures that at least one
 * byte is allocated even when the caller asks for zero bytes.
 *
 * \param ptr The pointer to reallocate.
 * \param size Number of bytes to be allocated.
 * \return Pointer to the piece of allocated memory; \c NULL if the allocation
 * failed.
 *
 * \sa \ref igraph_free(), \ref igraph_malloc()
 */

void *igraph_realloc(void* ptr, size_t size) {
    return (void*) IGRAPH_REALLOC(ptr, size, char);
}
