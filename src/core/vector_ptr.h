/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2021  The igraph development team

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

#ifndef IGRAPH_CORE_VECTOR_PTR_H
#define IGRAPH_CORE_VECTOR_PTR_H

#include "igraph_decls.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

IGRAPH_PRIVATE_EXPORT igraph_finally_func_t* igraph_i_vector_ptr_get_item_destructor(const igraph_vector_ptr_t *v);
IGRAPH_PRIVATE_EXPORT igraph_finally_func_t* igraph_i_vector_ptr_set_item_destructor(
    igraph_vector_ptr_t *v, igraph_finally_func_t *func);

/**
 * \define IGRAPH_I_VECTOR_PTR_SET_ITEM_DESTRUCTOR
 * \brief Sets the item destructor for this pointer vector (macro version).
 *
 * This macro is expanded to \ref igraph_i_vector_ptr_set_item_destructor(), the
 * only difference is that the second argument is automatically cast to an
 * \c igraph_finally_func_t*. The cast is necessary in most cases as the
 * destructor functions we use (such as \ref igraph_vector_destroy()) take a
 * pointer to some concrete igraph data type, while \c igraph_finally_func_t
 * expects \c void*
 */
#define IGRAPH_I_VECTOR_PTR_SET_ITEM_DESTRUCTOR(v, func) \
    igraph_i_vector_ptr_set_item_destructor((v), (igraph_finally_func_t*)(func))

__END_DECLS

#endif
