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

#ifndef IGRAPH_CONSTRUCTORS_H
#define IGRAPH_CONSTRUCTORS_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_types.h"
#include "igraph_matrix.h"
#include "igraph_datatype.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Constructors, deterministic                        */
/* -------------------------------------------------- */

DECLDIR int igraph_create(igraph_t *graph, const igraph_vector_t *edges, igraph_integer_t n,
                          igraph_bool_t directed);
DECLDIR int igraph_small(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed,
                         ...);
DECLDIR int igraph_adjacency(igraph_t *graph, igraph_matrix_t *adjmatrix,
                             igraph_adjacency_t mode);
DECLDIR int igraph_weighted_adjacency(igraph_t *graph, igraph_matrix_t *adjmatrix,
                                      igraph_adjacency_t mode, const char* attr,
                                      igraph_bool_t loops);
DECLDIR int igraph_star(igraph_t *graph, igraph_integer_t n, igraph_star_mode_t mode,
                        igraph_integer_t center);
DECLDIR int igraph_lattice(igraph_t *graph, const igraph_vector_t *dimvector, igraph_integer_t nei,
                           igraph_bool_t directed, igraph_bool_t mutual, igraph_bool_t circular);
DECLDIR int igraph_ring(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed,
                        igraph_bool_t mutual, igraph_bool_t circular);
DECLDIR int igraph_tree(igraph_t *graph, igraph_integer_t n, igraph_integer_t children,
                        igraph_tree_mode_t type);
DECLDIR int igraph_from_prufer(igraph_t *graph, const igraph_vector_int_t *prufer);
DECLDIR int igraph_full(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, igraph_bool_t loops);
DECLDIR int igraph_full_citation(igraph_t *graph, igraph_integer_t n,
                                 igraph_bool_t directed);
DECLDIR int igraph_atlas(igraph_t *graph, int number);
DECLDIR int igraph_extended_chordal_ring(igraph_t *graph, igraph_integer_t nodes,
        const igraph_matrix_t *W, igraph_bool_t directed);
DECLDIR int igraph_connect_neighborhood(igraph_t *graph, igraph_integer_t order,
                                        igraph_neimode_t mode);
DECLDIR int igraph_linegraph(const igraph_t *graph, igraph_t *linegraph);

DECLDIR int igraph_de_bruijn(igraph_t *graph, igraph_integer_t m, igraph_integer_t n);
DECLDIR int igraph_kautz(igraph_t *graph, igraph_integer_t m, igraph_integer_t n);
DECLDIR int igraph_famous(igraph_t *graph, const char *name);
DECLDIR int igraph_lcf_vector(igraph_t *graph, igraph_integer_t n,
                              const igraph_vector_t *shifts,
                              igraph_integer_t repeats);
DECLDIR int igraph_lcf(igraph_t *graph, igraph_integer_t n, ...);
DECLDIR int igraph_realize_degree_sequence(igraph_t *graph,
        const igraph_vector_t *outdeg, const igraph_vector_t *indeg,
        igraph_realize_degseq_t method);

__END_DECLS

#endif
