/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_ESTACK_H
#define IGRAPH_ESTACK_H

#include "igraph_decls.h"
#include "igraph_stack.h"
#include "igraph_vector.h"

__BEGIN_DECLS

typedef struct igraph_estack_t {
    igraph_stack_int_t stack;
    igraph_vector_bool_t isin;
} igraph_estack_t;

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_estack_init(
    igraph_estack_t *s, igraph_integer_t setsize, igraph_integer_t stacksize);
IGRAPH_PRIVATE_EXPORT void igraph_estack_destroy(igraph_estack_t *s);

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_estack_push(igraph_estack_t *s, igraph_integer_t elem);
IGRAPH_PRIVATE_EXPORT igraph_integer_t igraph_estack_pop(igraph_estack_t *s);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_estack_iselement(const igraph_estack_t *s,
                                                            igraph_integer_t elem);
IGRAPH_PRIVATE_EXPORT igraph_integer_t igraph_estack_size(const igraph_estack_t *s);

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_estack_print(const igraph_estack_t *s);

__END_DECLS

#endif
