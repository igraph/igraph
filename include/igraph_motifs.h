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

#ifndef IGRAPH_MOTIFS_H
#define IGRAPH_MOTIFS_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_iterators.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Graph motifs                                       */
/* -------------------------------------------------- */

/**
 * \typedef igraph_motifs_handler_t
 * \brief Callback type for \c igraph_motifs_randesu_callback.
 *
 * \ref igraph_motifs_randesu_callback() calls a specified callback
 * function whenever a new motif is found during a motif search. This
 * callback function must be of type \c igraph_motifs_handler_t. It has
 * the following arguments:
 *
 * \param graph The graph that that algorithm is working on. Of course
 *   this must not be modified.
 * \param vids The IDs of the vertices in the motif that has just been
 *   found. This vector is owned by the motif search algorithm, so do not
 *   modify or destroy it; make a copy of it if you need it later.
 * \param isoclass The isomorphism class of the motif that has just been
 *   found. Use \ref igraph_graph_count() to find the maximum possible
 *   isoclass for graphs of a given size. See \ref igraph_isoclass and
 *   \ref igraph_isoclass_subgraph for more information.
 * \param extra The extra argument that was passed to \ref
 *   igraph_motifs_randesu_callback().
 * \return \c IGRAPH_SUCCESS to continue the motif search,
 *    \c IGRAPH_STOP to stop the motif search and return to the caller
 *    normally. Any other return value is interpreted as an igraph error code,
 *    which will terminate the search and return the same error code to the
 *    caller.
 *
 * \sa \ref igraph_motifs_randesu_callback()
 */

typedef igraph_error_t igraph_motifs_handler_t(const igraph_t *graph,
        igraph_vector_int_t *vids,
        igraph_integer_t isoclass,
        void* extra);

IGRAPH_EXPORT igraph_error_t igraph_motifs_randesu(const igraph_t *graph, igraph_vector_t *hist,
                                        igraph_integer_t size, const igraph_vector_t *cut_prob);

IGRAPH_EXPORT igraph_error_t igraph_motifs_randesu_callback(const igraph_t *graph, igraph_integer_t size,
                                                 const igraph_vector_t *cut_prob,
                                                 igraph_motifs_handler_t *callback,
                                                 void* extra);

IGRAPH_EXPORT igraph_error_t igraph_motifs_randesu_estimate(const igraph_t *graph, igraph_integer_t *est,
                                                 igraph_integer_t size, const igraph_vector_t *cut_prob,
                                                 igraph_integer_t sample_size,
                                                 const igraph_vector_int_t *sample);
IGRAPH_EXPORT igraph_error_t igraph_motifs_randesu_no(const igraph_t *graph, igraph_integer_t *no,
                                           igraph_integer_t size, const igraph_vector_t *cut_prob);

IGRAPH_EXPORT igraph_error_t igraph_dyad_census(const igraph_t *graph, igraph_real_t *mut,
                                     igraph_real_t *asym, igraph_real_t *null);
IGRAPH_EXPORT igraph_error_t igraph_triad_census(const igraph_t *igraph, igraph_vector_t *res);

IGRAPH_EXPORT igraph_error_t igraph_adjacent_triangles(const igraph_t *graph,
                                            igraph_vector_t *res,
                                            const igraph_vs_t vids);

IGRAPH_EXPORT igraph_error_t igraph_list_triangles(const igraph_t *graph,
                                        igraph_vector_int_t *res);

__END_DECLS

#endif
