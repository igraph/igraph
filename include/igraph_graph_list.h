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

#ifndef IGRAPH_GRAPH_LIST_H
#define IGRAPH_GRAPH_LIST_H

#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* List of graphs                                     */
/* -------------------------------------------------- */

#define GRAPH_LIST
#define BASE_GRAPH
#define EXTRA_TYPE_FIELDS igraph_bool_t directed;
#include "igraph_pmt.h"
#include "igraph_typed_list_pmt.h"
#include "igraph_pmt_off.h"
#undef EXTRA_TYPE_FIELDS
#undef BASE_GRAPH
#undef GRAPH_LIST

void igraph_graph_list_set_directed(igraph_graph_list_t* list, igraph_bool_t directed);

/* -------------------------------------------------- */
/* Helper macros                                      */
/* -------------------------------------------------- */

#define IGRAPH_GRAPH_LIST_INIT_FINALLY(v, size) \
    do { IGRAPH_CHECK(igraph_graph_list_init(v, size)); \
        IGRAPH_FINALLY(igraph_graph_list_destroy, v); } while (0)

__END_DECLS

#endif
