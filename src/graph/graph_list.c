/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2022  The igraph development team
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_graph_list.h"

#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_types.h"

#define GRAPH_LIST
#define BASE_GRAPH
#define CUSTOM_INIT_DESTROY
#include "igraph_pmt.h"
#include "../core/typed_list.pmt"
#include "igraph_pmt_off.h"
#undef CUSTOM_INIT_DESTROY
#undef BASE_GRAPH
#undef GRAPH_LIST

void igraph_graph_list_set_directed(
    igraph_graph_list_t* list, igraph_bool_t directed
) {
    IGRAPH_ASSERT(list != 0);
    list->directed = directed;
}

static igraph_error_t igraph_i_graph_list_init_item(
    const igraph_graph_list_t* list, igraph_t* item
) {
    return igraph_empty(item, 0, list->directed);
}

static igraph_error_t igraph_i_graph_list_copy_item(
    igraph_t* dest, const igraph_t* source
) {
    return igraph_copy(dest, source);
}

static void igraph_i_graph_list_destroy_item(igraph_t* item) {
    igraph_destroy(item);
}
