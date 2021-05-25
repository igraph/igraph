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

#include "igraph_interface.h"

/**
 * \ingroup generators
 * \function igraph_from_prufer
 * \brief Generates a tree from a Pr&uuml;fer sequence.
 *
 * A Pr&uuml;fer sequence is a unique sequence of integers associated
 * with a labelled tree. A tree on n vertices can be represented by a
 * sequence of n-2 integers, each between 0 and n-1 (inclusive).
 *
 * The algorithm used by this function is based on
 * Paulius Micikevi&ccaron;ius, Saverio Caminiti, Narsingh Deo:
 * Linear-time Algorithms for Encoding Trees as Sequences of Node Labels
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param prufer The Pr&uuml;fer sequence
 * \return Error code:
 *          \clist
 *          \cli IGRAPH_ENOMEM
 *             there is not enough memory to perform the operation.
 *          \cli IGRAPH_EINVAL
 *             invalid Pr&uuml;fer sequence given
 *          \endclist
 *
 * \sa \ref igraph_to_prufer(), \ref igraph_tree(), \ref igraph_tree_game()
 *
 */
int igraph_from_prufer(igraph_t *graph, const igraph_vector_int_t *prufer) {
    igraph_vector_int_t degree;
    igraph_vector_t edges;
    long n;
    long i, k;
    long u, v; /* vertices */
    long ec;

    n = igraph_vector_int_size(prufer) + 2;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&degree, n); /* initializes vector to zeros */
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 2 * (n - 1));

    /* build out-degree vector (i.e. number of child vertices) and verify Prufer sequence */
    for (i = 0; i < n - 2; ++i) {
        long u = VECTOR(*prufer)[i];
        if (u >= n || u < 0) {
            IGRAPH_ERROR("Invalid Prufer sequence", IGRAPH_EINVAL);
        }
        VECTOR(degree)[u] += 1;
    }

    v = 0;  /* initialize v now, in case Prufer sequence is empty */
    k = 0;  /* index into the Prufer vector */
    ec = 0; /* index into the edges vector */
    for (i = 0; i < n; ++i) {
        u = i;

        while (k < n - 2 && u <= i && (VECTOR(degree)[u] == 0)) {
            /* u is a leaf here */

            v = VECTOR(*prufer)[k]; /* parent of u */

            /* add edge */
            VECTOR(edges)[ec++] = v;
            VECTOR(edges)[ec++] = u;

            k += 1;

            VECTOR(degree)[v] -= 1;

            u = v;
        }

        if (k == n - 2) {
            break;
        }
    }

    /* find u for last edge, v is already set */
    for (u = i + 1; u < n; ++u)
        if ((VECTOR(degree)[u] == 0) && u != v) {
            break;
        }

    /* add last edge */
    VECTOR(edges)[ec++] = v;
    VECTOR(edges)[ec++] = u;

    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) n, /* directed = */ 0));

    igraph_vector_destroy(&edges);
    igraph_vector_int_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
