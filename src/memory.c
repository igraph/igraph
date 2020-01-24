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
#include "config.h"

/**
 * \function igraph_free
 * Deallocate memory that was allocated by igraph functions
 *
 * Some igraph functions return a pointer vector (igraph_vector_ptr_t)
 * containing pointers to other igraph or other data types. These data
 * types are dynamically allocated and have to be deallocated
 * manually, if the user does not need them any more. This can be done
 * by calling igraph_free on them.
 *
 * </para><para>
 * Here is a complete example on how to use \c igraph_free properly.
 * <programlisting>
 * <![CDATA[#include <igraph.h>
 *
 * int main(void)
 * {
 *    igraph_t graph;
 *    igraph_vector_ptr_t seps;
 *    long int i;
 *
 *    igraph_famous(&graph, "tutte");
 *    igraph_vector_ptr_init(&seps, 0);
 *    igraph_minimum_size_separators(&graph, &seps);
 *
 *    for (i=0; i<igraph_vector_ptr_size(&seps); i++) {
 *      igraph_vector_t *v=VECTOR(seps)[i];
 *      igraph_vector_print(v);
 *      igraph_vector_destroy(v);
 *      igraph_free(v);
 *    }
 *
 *    igraph_vector_ptr_destroy(&seps);
 *    igraph_destroy(&graph);
 *    return 0;
 * }]]>
 * </programlisting>
 *
 *
 * \param p Pointer to the piece of memory to be deallocated.
 * \return Error code, currently always zero, meaning success.
 *
 * Time complexity: platform dependent, ideally it should be O(1).
 *
 * \sa \ref igraph_malloc()
 */

int igraph_free(void *p) {
    igraph_Free(p);
    return 0;
}


/**
 * \function igraph_malloc
 * Allocate memory that can be safely deallocated by igraph functions
 *
 * Some igraph functions, such as \ref igraph_vector_ptr_free_all() and
 * \ref igraph_vector_ptr_destroy_all() can free memory that may have been
 * allocated by the user.  \c igraph_malloc() works exactly like \c malloc()
 * from the C standard library, but it is guaranteed that it can be safely
 * paired with the \c free() function used by igraph internally (which is
 * also user-accessible through \ref igraph_free()).
 *
 * \param n Number of bytes to be allocated.
 * \return Pointer to the piece of allocated memory.
 *
 * \sa \ref igraph_free()
 */

void *igraph_malloc(size_t n) {
    return malloc(n);
}
