/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_INTERFACE_H
#define IGRAPH_INTERFACE_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_datatype.h"
#include "igraph_iterators.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Interface                                          */
/* -------------------------------------------------- */

DECLDIR int igraph_empty(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed);
DECLDIR int igraph_empty_attrs(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, void *attr);
DECLDIR void igraph_destroy(igraph_t *graph);
DECLDIR int igraph_copy(igraph_t *to, const igraph_t *from);
DECLDIR int igraph_add_edges(igraph_t *graph, const igraph_vector_t *edges,
                             void *attr);
DECLDIR int igraph_add_vertices(igraph_t *graph, igraph_integer_t nv,
                                void *attr);
DECLDIR int igraph_delete_edges(igraph_t *graph, igraph_es_t edges);
DECLDIR int igraph_delete_vertices(igraph_t *graph, const igraph_vs_t vertices);
DECLDIR int igraph_delete_vertices_idx(igraph_t *graph, const igraph_vs_t vertices,
                                       igraph_vector_t *idx,
                                       igraph_vector_t *invidx);
DECLDIR igraph_integer_t igraph_vcount(const igraph_t *graph);
DECLDIR igraph_integer_t igraph_ecount(const igraph_t *graph);
DECLDIR int igraph_neighbors(const igraph_t *graph, igraph_vector_t *neis, igraph_integer_t vid,
                             igraph_neimode_t mode);
DECLDIR igraph_bool_t igraph_is_directed(const igraph_t *graph);
DECLDIR int igraph_degree(const igraph_t *graph, igraph_vector_t *res,
                          const igraph_vs_t vids, igraph_neimode_t mode,
                          igraph_bool_t loops);
DECLDIR int igraph_edge(const igraph_t *graph, igraph_integer_t eid,
                        igraph_integer_t *from, igraph_integer_t *to);
DECLDIR int igraph_edges(const igraph_t *graph, igraph_es_t eids,
                         igraph_vector_t *edges);
DECLDIR int igraph_get_eid(const igraph_t *graph, igraph_integer_t *eid,
                           igraph_integer_t from, igraph_integer_t to,
                           igraph_bool_t directed, igraph_bool_t error);
DECLDIR int igraph_get_eids(const igraph_t *graph, igraph_vector_t *eids,
                            const igraph_vector_t *pairs,
                            const igraph_vector_t *path,
                            igraph_bool_t directed, igraph_bool_t error);
DECLDIR int igraph_get_eids_multi(const igraph_t *graph, igraph_vector_t *eids,
                                  const igraph_vector_t *pairs,
                                  const igraph_vector_t *path,
                                  igraph_bool_t directed, igraph_bool_t error);
DECLDIR int igraph_adjacent(const igraph_t *graph, igraph_vector_t *eids, igraph_integer_t vid,
                            igraph_neimode_t mode);          /* deprecated */
DECLDIR int igraph_incident(const igraph_t *graph, igraph_vector_t *eids, igraph_integer_t vid,
                            igraph_neimode_t mode);

#define IGRAPH_FROM(g,e) ((igraph_integer_t)(VECTOR((g)->from)[(long int)(e)]))
#define IGRAPH_TO(g,e)   ((igraph_integer_t)(VECTOR((g)->to)  [(long int)(e)]))
#define IGRAPH_OTHER(g,e,v) \
    ((igraph_integer_t)(IGRAPH_TO(g,(e))==(v) ? IGRAPH_FROM((g),(e)) : IGRAPH_TO((g),(e))))

__END_DECLS

#endif
