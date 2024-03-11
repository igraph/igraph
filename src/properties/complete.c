/*
  IGraph library.
  Copyright (C) 2024 The igraph development team <igraph@igraph.org>

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  this program. If not, see <https://www.gnu.org/licenses/>.
*/


#include "igraph_interface.h"
#include "igraph_structural.h"

#include "graph/internal.h"

/**
 * \ingroup structural
 * \function igraph_is_complete
 * \brief Decides whether the graph is complete.
 *
 * \experimental
 *
 * A graph is considered complete if all pairs of different vertices are
 * adjacent.
 *
 * </para><para>
 * The null graph and the singleton graph are considered complete.
 *
 * \param graph The graph object to analyze.
 * \param res Pointer to a logical variable, the result will be stored here.
 *
 * \return Error code.
 *
 * Time complexity: O(|V| + |E|) at worst.
 */

igraph_error_t igraph_is_complete(const igraph_t *graph, igraph_bool_t *res)
{
    const igraph_integer_t vcount = igraph_vcount(graph);
    const igraph_integer_t ecount = igraph_ecount(graph);
    igraph_integer_t complete_ecount;
    igraph_bool_t simple, directed = igraph_is_directed(graph);
    igraph_vector_int_t neighbours;

    /* If the graph is the null graph or the singleton graph, return early */
    if (vcount == 0 || vcount == 1) {
        *res = true;
        return IGRAPH_SUCCESS;
    }

    /* Compute the amount of edges a complete graph of vcount vertices would
       have */

    /* Depends on whether the graph is directed */

    /* We have to take care of integer overflowing */

#if IGRAPH_INTEGER_SIZE == 32
    if (directed) {
        /* Highest x s.t. x² - x < 2^31 - 1 */
        if (vcount > 46341) {
            *res = false;
            return IGRAPH_SUCCESS;
        } else {
            complete_ecount = vcount * (vcount - 1);
        }
    } else {
        /* Highest x s.t. (x² - x) / 2 < 2^31 - 1 */
        if (vcount > 65536) {
            *res = false;
            return IGRAPH_SUCCESS;
        } else {
            complete_ecount = vcount % 2 == 0 ?
                              (vcount / 2) * (vcount - 1) :
                              vcount * ((vcount - 1) / 2);
        }
    }
#elif IGRAPH_INTEGER_SIZE == 64
    if (directed) {
        /* Highest x s.t. x² - x < 2^63 - 1 */
        if (vcount > 3037000500) {
            *res = false;
            return IGRAPH_SUCCESS;
        } else {
            complete_ecount = vcount * (vcount - 1);
        }
    } else {
        /* Highest x s.t. (x² - x) / 2 < 2^63 - 1 */
        if (vcount > 4294967296) {
            *res = false;
            return IGRAPH_SUCCESS;
        } else {
            complete_ecount = vcount % 2 == 0 ?
                              (vcount / 2) * (vcount - 1) :
                              vcount * ((vcount - 1) / 2);
        }
    }
#else
    /* If values other than 32 or 64 become allowed,
     * this code will need to be updated. */
#  error "Unexpected IGRAPH_INTEGER_SIZE value."
#endif

    /* If the amount of edges is strictly lower than what it should be for a
       complete graph, return early */

    if (ecount < complete_ecount) {
        *res = false;
        return IGRAPH_SUCCESS;
    }

    /* If the graph is simple, compare and conclude */
    IGRAPH_CHECK(igraph_is_simple(graph, &simple));

    if (simple) {
        *res = (ecount == complete_ecount);
        return IGRAPH_SUCCESS;
    }

    /* Allocate memory for vector of size v */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neighbours, vcount);

    for (igraph_integer_t i = 0; i < vcount; ++i) {

        igraph_vector_int_clear(&neighbours);

        IGRAPH_CHECK(igraph_i_neighbors(graph, &neighbours, i, IGRAPH_OUT,
                                        IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));

        if ((igraph_vector_int_size(&neighbours) < vcount - 1)) {
            *res = false;
            goto cleanup;
        }
    }

    /* If we arrive here, we have found no neighbour vector of size strictly
       less than vcount - 1. The graph is therefore complete */

    *res = true;

cleanup:

    igraph_vector_int_destroy(&neighbours);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
