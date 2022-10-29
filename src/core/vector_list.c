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
 * \ingroup vector_list
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
 * an \type igraph_vector_int_t), \type igraph_matrix_list_t for matrices of
 * doubles and so on. The following list summarizes the variants that are
 * currently available in the library:</para>
 *
 * \ilist
 * \ili \type igraph_vector_list_t for lists of vectors of floating-point numbers
 *      (\type igraph_vector_t)
 * \ili \type igraph_vector_int_list_t for lists of integer vectors
 *      (\type igraph_vector_int_t)
 * \ili \type igraph_matrix_list_t for lists of matrices of floating-point numbers
 *      (\type igraph_matrix_t)
 * \ili \type igraph_graph_list_t for lists of graphs (\type igraph_t)
 * \endilist
 *
 * <para>Lists of vectors are used in \a igraph in many
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
 * lists of vectors, we must also talk about the \em ownership rules of these
 * objects. The most important rule is that the vectors in the list are
 * owned by the list itself, meaning that the user is \em not responsible
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
 * clear the list by calling \ref igraph_vector_list_clear(), the list will
 * implicitly destroy these lists, and any pointers that you may hold to the
 * items become invalid at once.</para>
 *
 * <para>Speaking about pointers, the typical way of working with vectors in
 * a list is to obtain a pointer to one of the items via the
 * \ref igraph_vector_list_get_ptr() method and then passing this pointer
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
 * \ingroup vector_list
 * \section igraph_vector_list_constructors_and_destructors Constructors and
 * destructors
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
 * in the vector list, you need to copy them with \ref igraph_vector_init_copy() or
 * \ref igraph_vector_update(), or you need to remove them from the list and
 * take ownership by calling \ref igraph_vector_list_pop_back(),
 * \ref igraph_vector_list_remove() or \ref igraph_vector_list_remove_fast() .</para>
 */


/**
 * \ingroup vector_list
 * \section igraph_vector_list_accessing_elements Accessing elements
 *
 * <para>Elements of a vector list may be accessed with the
 * \ref igraph_vector_list_get_ptr() function. The function returns a \em pointer
 * to the vector with a given index inside the list, and you may then pass
 * this pointer onwards to other functions that can query or manipulate
 * vectors. The pointer itself is guaranteed to stay valid as long as the
 * list itself is not modified; however, \em any modification to the list
 * will invalidate the pointer, even modifications that are seemingly unrelated
 * to the vector that the pointer points to (such as adding a new vector at
 * the end of the list). This is because the list data structure may be forced
 * to re-allocate its internal storage if a new element does not fit into the
 * already allocated space, and there are no guarantees that the re-allocated
 * block remains at the same memory location (typically it gets moved elsewhere).
 * </para>
 *
 * <para>Note that the standard \ref VECTOR macro that works for ordinary vectors
 * does not work for lists of vectors to access the i-th element (but of course
 * you can use it to index into an existing vector that you retrieved from the
 * vector list with \ref igraph_vector_list_get_ptr() ). This is because the
 * macro notation would allow one to overwrite the vector in the list with
 * another one without the list knowing about it, so the list would not be able
 * to destroy the vector that was overwritten by a new one.
 * </para>
 *
 * <para> \ref igraph_vector_list_tail_ptr() returns a pointer to the last
 * vector in the list, or \c NULL if the list is empty. There is no
 * <function>igraph_vector_list_head_ptr()</function>, however, as it is easy to
 * write <code>igraph_vector_list_get_ptr(v, 0)</code> instead.</para>
 */
