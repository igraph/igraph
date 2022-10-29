/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2021 The igraph development team

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

#include "igraph_constructors.h"

#include "igraph_operators.h"

#include "math/safe_intop.h"

/**
 * \function igraph_lcf_vector
 * \brief Creates a graph from LCF notation.
 *
 * This function is essentially the same as \ref igraph_lcf(), only
 * the way for giving the arguments is different. See \ref
 * igraph_lcf() for details.
 * \param graph Pointer to an uninitialized graph object.
 * \param n Integer constant giving the number of vertices.
 * \param shifts A vector giving the shifts.
 * \param repeats An integer constant giving the number of repeats
 *        for the shifts.
 * \return Error code.
 *
 * \sa \ref igraph_lcf(), \ref igraph_extended_chordal_ring()
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices plus
 * the number of edges.
 */
igraph_error_t igraph_lcf_vector(igraph_t *graph, igraph_integer_t n,
                      const igraph_vector_int_t *shifts,
                      igraph_integer_t repeats) {

    igraph_vector_int_t edges;
    igraph_integer_t no_of_shifts = igraph_vector_int_size(shifts);
    igraph_integer_t ptr = 0, i, sptr = 0;
    igraph_integer_t no_of_nodes = n;
    igraph_integer_t no_of_edges = n + no_of_shifts * repeats;
    igraph_integer_t no_of_edges2;

    if (repeats < 0) {
        IGRAPH_ERROR("Number of repeats must not be negative.", IGRAPH_EINVAL);
    }

    /* no_of_edges = n + no_of_shifts * repeats */
    IGRAPH_SAFE_MULT(no_of_shifts, repeats, &no_of_edges);
    IGRAPH_SAFE_ADD(no_of_edges, n, &no_of_edges);
    IGRAPH_SAFE_MULT(no_of_edges, 2, &no_of_edges2);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges2);

    if (no_of_nodes > 0) {
        /* Create a ring first */
        for (i = 0; i < no_of_nodes; i++) {
            VECTOR(edges)[ptr++] = i;
            VECTOR(edges)[ptr++] = i + 1;
        }
        VECTOR(edges)[ptr - 1] = 0;
    }

    /* Then add the rest */
    while (ptr < 2 * no_of_edges) {
        igraph_integer_t sh = VECTOR(*shifts)[sptr % no_of_shifts];
        igraph_integer_t from = sptr % no_of_nodes;
        igraph_integer_t to = (no_of_nodes + sptr + sh) % no_of_nodes;
        VECTOR(edges)[ptr++] = from;
        VECTOR(edges)[ptr++] = to;
        sptr++;
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, IGRAPH_UNDIRECTED));
    IGRAPH_CHECK(igraph_simplify(graph, true, true, NULL));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_lcf
 * \brief Creates a graph from LCF notation.
 *
 * </para><para>
 * LCF is short for Lederberg-Coxeter-Frucht, it is a concise notation for
 * 3-regular Hamiltonian graphs. It consists of three parameters: the
 * number of vertices in the graph, a list of shifts giving additional
 * edges to a cycle backbone, and another integer giving how many times
 * the shifts should be performed. See
 * http://mathworld.wolfram.com/LCFNotation.html for details.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n Integer, the number of vertices in the graph.
 * \param ... The shifts and the number of repeats for the shifts,
 *        plus an additional 0 to mark the end of the arguments.
 * \return Error code.
 *
 * \sa See \ref igraph_lcf_vector() for a similar function using a
 * vector_t instead of the variable length argument list.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges.
 *
 * \example examples/simple/igraph_lcf.c
 */
igraph_error_t igraph_lcf(igraph_t *graph, igraph_integer_t n, ...) {
    igraph_vector_int_t shifts;
    igraph_integer_t repeats;
    va_list ap;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&shifts, 0);

    va_start(ap, n);
    while (1) {
        igraph_error_t err;
        int num = va_arg(ap, int);
        if (num == 0) {
            break;
        }
        err = igraph_vector_int_push_back(&shifts, num);
        if (err != IGRAPH_SUCCESS) {
            va_end(ap);
            IGRAPH_ERROR("", err);
        }
    }
    va_end(ap);
    if (igraph_vector_int_size(&shifts) == 0) {
        repeats = 0;
    } else {
        repeats = igraph_vector_int_pop_back(&shifts);
    }

    IGRAPH_CHECK(igraph_lcf_vector(graph, n, &shifts, repeats));
    igraph_vector_int_destroy(&shifts);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
