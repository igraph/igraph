/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2022  The igraph development team

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

#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector_list.h"

#define VECTOR_LIST

#define BASE_IGRAPH_REAL
#include "igraph_pmt.h"
#include "typed_list.pmt"
#include "igraph_pmt_off.h"
#undef BASE_IGRAPH_REAL

#define BASE_INT
#include "igraph_pmt.h"
#include "typed_list.pmt"
#include "igraph_pmt_off.h"
#undef BASE_INT

#undef VECTOR_LIST

/**
 * \ingroup vector
 * \section about_igraph_vector_list_t_objects About \type igraph_vector_list_t objects
 *
 * <para>The \type igraph_vector_list_t data type is essentially a list of
 * \type igraph_vector_t objects with automatic memory management. It is something
 * similar to (but much simpler than) the \type vector template in the C++
 * standard library where the elements are vectors themselves.</para>
 *
 * <para>There are multiple variants of \type igraph_vector_list_t; the basic variant
 * stores vectors of doubles (i.e. each item is an \ref igraph_vector_t), but
 * there is also \type igraph_vector_int_list_t for integers (where each item is
 * an \type igraph_vector_int_t). Lists of vectors are used in \a igraph in many
 * cases, e.g., when returning lists of paths, cliques or vertex sets.
 * Functions that expect or return a list of numeric vectors typically use
 * \type igraph_vector_list_t or \type igraph_vector_int_list_t to achieve this.
 * Lists of integer vectors are used when the vectors in the list are supposed
 * to hold vertex or edge identifiers, while lists of floating-point vectors
 * are used when the vectors are expected to hold fractional numbers or
 * infinities.</para>
 *
 * <para>The elements in an \type igraph_vector_list_t object and its variants are
 * indexed from zero, we follow the usual C convention here.</para>
 *
 * <para>Almost all of the functions described below for \type igraph_vector_list_t
 * also exist for all the other vector list variants. These variants are not
 * documented separately; you can simply replace \c vector_list with, say,
 * \c vector_int_list if you need a function for another variant. For instance,
 * to initialize a list of integer vectors, you need to use
 * \c igraph_vector_int_list_init() and not \ref igraph_vector_list_init().</para>
 *
 * <para>Before diving into a detailed description of the functions related to
 * lists of vectors, we must also talk about the \em "ownership rules" of these
 * objects. The most important rule is that the vectors in the list are
 * \em owned by the list itself, meaning that the user is \em not responsible
 * for allocating memory for the vectors or for freeing the memory associated
 * to the vectors. It is the responsibility of the list to allocate and initialize
 * the vectors when new items are created in the list, and it is also the
 * responsibility of the list to destroy the items when they are removed from
 * the list without passing on their ownership to the user. As a consequence,
 * the list may not contain "uninitialized" or "null" items; each item is
 * initialized when it comes to existence. If you create a list containing
 * one million vectors, you are not only allocating memory for one million
 * \ref igraph_vector_t object but you are also initializing one million
 * vectors. Also, if you have a list containing one million vectors and you
 * clear the list by calling \c igraph_vector_list_clear(), the list will
 * implicitly destroy these lists, and any pointers that you may hold to the
 * items become invalid at once.</para>
 *
 * <para>Speaking about pointers, the typical way of working with vectors in
 * a list is to obtain a pointer to one of the items via the
 * \c igraph_vector_list_get_ptr() method and then passing this pointer
 * onwards to functions that manipulate \ref igraph_vector_t objects. However,
 * note that the pointers are \em ephemeral in the sense that they may be
 * invalidated any time when the list is modified because a modification may
 * involve the re-allocation of the internal storage of the list if more space
 * is needed, and the pointers that you obtained will not follow the
 * reallocation. This limitation does not appear often in real-world usage of
 * \c igraph_vector_list_t and in general, the advantages of the automatic
 * memory management outweigh this limitation.</para>
 */

/**
 * \ingroup vector
 * \section igraph_vector_list_constructors_and_destructors Constructors and
 * Destructors
 *
 * <para>\type igraph_vector_list_t objects have to be initialized before using
 * them, this is analogous to calling a constructor on them.
 * \ref igraph_vector_list_init() is the basic constructor; it creates a list
 * of the given length and also initializes each vector in the newly created
 * list to zero length.</para>
 *
 * <para>If an \type igraph_vector_list_t object is not needed any more, it
 * should be destroyed to free its allocated memory by calling the
 * \type igraph_vector_list_t destructor, \ref igraph_vector_list_destroy().
 * Calling the destructor also destroys all the vectors inside the vector
 * list due to the ownership rules. If you want to keep a few of the vectors
 * in the vector list, you need to copy them with \ref igraph_vector_copy() or
 * \ref igraph_vector_update(), or you need to remove them from the list and
 * take ownership by calling \c igraph_vector_list_remove() or
 * \c igraph_vector_list_remove_fast() .</para>
 */
