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

#ifndef IGRAPH_HRG_H
#define IGRAPH_HRG_H

#include "igraph_decls.h"

#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_graph_list.h"
#include "igraph_vector.h"

__BEGIN_DECLS

/**
 * \struct igraph_hrg_t
 * \brief Data structure to store a hierarchical random graph.
 *
 * A hierarchical random graph (HRG) can be given as a binary tree,
 * where the internal vertices are labeled with real numbers.
 *
 * </para><para>Note that you don't necessarily have to know this
 * internal representation for using the HRG functions, just pass the
 * HRG objects created by one igraph function, to another igraph
 * function.
 *
 * </para><para>
 * It has the following members:
 *
 * \member left Vector that contains the left children of the internal
 *    tree vertices. The first vertex is always the root vertex, so
 *    the first element of the vector is the left child of the root
 *    vertex. Internal vertices are denoted with negative numbers,
 *    starting from -1 and going down, i.e. the root vertex is
 *    -1. Leaf vertices are denoted by non-negative number, starting
 *    from zero and up.
 * \member right Vector that contains the right children of the
 *    vertices, with the same encoding as the \c left vector.
 * \member prob The connection probabilities attached to the internal
 *    vertices, the first number belongs to the root vertex
 *    (i.e. internal vertex -1), the second to internal vertex -2,
 *    etc.
 * \member edges The number of edges in the subtree below the given
 *    internal vertex.
 * \member vertices The number of vertices in the subtree below the
 *    given internal vertex, including itself.
 */

typedef struct igraph_hrg_t {
    igraph_vector_int_t left;
    igraph_vector_int_t right;
    igraph_vector_t prob;
    igraph_vector_int_t vertices;
    igraph_vector_int_t edges;
} igraph_hrg_t;

IGRAPH_EXPORT igraph_error_t igraph_hrg_init(igraph_hrg_t *hrg, igraph_integer_t n);
IGRAPH_EXPORT void igraph_hrg_destroy(igraph_hrg_t *hrg);
IGRAPH_EXPORT igraph_integer_t igraph_hrg_size(const igraph_hrg_t *hrg);
IGRAPH_EXPORT igraph_error_t igraph_hrg_resize(igraph_hrg_t *hrg, igraph_integer_t newsize);

IGRAPH_EXPORT igraph_error_t igraph_hrg_fit(
    const igraph_t *graph, igraph_hrg_t *hrg, igraph_bool_t start,
    igraph_integer_t steps
);

IGRAPH_EXPORT igraph_error_t igraph_hrg_sample(
    const igraph_hrg_t *hrg, igraph_t *sample
);

IGRAPH_EXPORT igraph_error_t igraph_hrg_sample_many(
    const igraph_hrg_t *hrg, igraph_graph_list_t *samples,
    igraph_integer_t num_samples
);

IGRAPH_EXPORT igraph_error_t igraph_hrg_game(
    igraph_t *graph, const igraph_hrg_t *hrg
);

IGRAPH_EXPORT igraph_error_t igraph_from_hrg_dendrogram(
    igraph_t *graph, const igraph_hrg_t *hrg, igraph_vector_t *prob
);

IGRAPH_EXPORT igraph_error_t igraph_hrg_consensus(const igraph_t *graph,
                                       igraph_vector_int_t *parents,
                                       igraph_vector_t *weights,
                                       igraph_hrg_t *hrg,
                                       igraph_bool_t start,
                                       igraph_integer_t num_samples);

IGRAPH_EXPORT igraph_error_t igraph_hrg_predict(const igraph_t *graph,
                                     igraph_vector_int_t *edges,
                                     igraph_vector_t *prob,
                                     igraph_hrg_t *hrg,
                                     igraph_bool_t start,
                                     igraph_integer_t num_samples,
                                     igraph_integer_t num_bins);

IGRAPH_EXPORT igraph_error_t igraph_hrg_create(igraph_hrg_t *hrg,
                                    const igraph_t *graph,
                                    const igraph_vector_t *prob);

/* Deprecated functions: */

IGRAPH_DEPRECATED IGRAPH_EXPORT igraph_error_t igraph_hrg_dendrogram(
    igraph_t *graph, const igraph_hrg_t *hrg
);

__END_DECLS

#endif  /* IGRAPH_HRG_H */
