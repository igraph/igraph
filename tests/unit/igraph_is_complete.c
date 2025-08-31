/*
  igraph library.
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

#include <igraph.h>

#include "test_utilities.h"

/* Check that the null graph is complete */
static igraph_error_t check_null(void)
{
    igraph_t null;
    igraph_bool_t complete;

    igraph_full(&null, 0, false, false);
    IGRAPH_ASSERT(!igraph_vcount(&null));
    IGRAPH_ASSERT(!igraph_ecount(&null));
    igraph_is_complete(&null, &complete);
    IGRAPH_ASSERT(complete);

    igraph_destroy(&null);

    return IGRAPH_SUCCESS;
}

/* Check that the singleton graph is complete */
static igraph_error_t check_singleton(void)
{
    igraph_t singleton;
    igraph_bool_t complete;

    igraph_full(&singleton, 1, false, false);
    IGRAPH_ASSERT(igraph_vcount(&singleton) == 1);
    IGRAPH_ASSERT(!igraph_ecount(&singleton));
    igraph_is_complete(&singleton, &complete);
    IGRAPH_ASSERT(complete);

    igraph_destroy(&singleton);

    return IGRAPH_SUCCESS;
}

/* Run checks on other non-trivial graphs */
/* - Generating a full graph without any loops, adding the mutuals if needed */
/* - Checking that it is complete */
/* - Removing an edge */
/* - Checking that it is not complete anymore */
/* - Then doing the same for multiple graphs */
static igraph_error_t check(const igraph_bool_t directed, const igraph_int_t vertices)
{
    igraph_t simple, multiple;
    igraph_bool_t complete = false;
    igraph_int_t eid = -1;
    igraph_int_t size = -1;
    igraph_es_t es;
    igraph_vector_bool_t mutuals;
    igraph_vector_int_t toadd;
    igraph_int_t i;

    /* generic, simple graph */
    igraph_full(&simple, vertices, directed, false);

    /* igraph_full does not produce complete directed graphs, merely
     * **tournaments**. So the mutuals must be added manually in order to get a
     * complete graph. */

    /* So we must find all edges for which the mutual is missing, and add it
       ourselves */

    if (directed) {

        igraph_vector_bool_init(&mutuals, vertices * vertices);
        igraph_vector_bool_clear(&mutuals);

        igraph_is_mutual(&simple, &mutuals, igraph_ess_all(IGRAPH_EDGEORDER_ID),
                         false);

        igraph_vector_int_init(&toadd, vertices * vertices);
        igraph_vector_int_clear(&toadd);

        for (i = 0 ; i < igraph_vector_bool_size(&mutuals); ++i) {
            /* mutual is missing, adding it */
            if (VECTOR(mutuals)[i] == false) {
                igraph_vector_int_push_back(&toadd, IGRAPH_TO(&simple, i));
                igraph_vector_int_push_back(&toadd, IGRAPH_FROM(&simple, i));
            }
        }

        igraph_add_edges(&simple, &toadd, NULL);
    }
    igraph_is_complete(&simple, &complete);
    IGRAPH_ASSERT(complete);

    /* Remove a specific edge */
    igraph_delete_edges(&simple, igraph_ess_1(42));
    igraph_is_complete(&simple, &complete);
    IGRAPH_ASSERT(!complete);

    /* generic, multiple graph */
    igraph_full(&multiple, vertices, directed, false);

    /* igraph_full does not produce complete directed graphs, merely
     * **tournaments**. So the mutuals must be added manually in order to get a
     * complete graph. */

    /* So we must find all edges for which the mutual is missing, and add it
       ourselves */
    if (directed) {

        igraph_vector_bool_clear(&mutuals);

        igraph_is_mutual(&multiple, &mutuals,
                         igraph_ess_all(IGRAPH_EDGEORDER_ID), false);

        igraph_vector_int_clear(&toadd);

        for (i = 0 ; i < igraph_vector_bool_size(&mutuals); ++i) {
            /* mutual is missing, adding it */
            if (VECTOR(mutuals)[i] == false) {
                igraph_vector_int_push_back(&toadd, IGRAPH_TO(&multiple, i));
                igraph_vector_int_push_back(&toadd, IGRAPH_FROM(&multiple, i));
            }
        }

        igraph_add_edges(&multiple, &toadd, NULL);
    }

    /* Add two loops */
    igraph_add_edge(&multiple, 41, 41);
    igraph_add_edge(&multiple, 42, 42);

    igraph_is_complete(&multiple, &complete);
    IGRAPH_ASSERT(complete);

    igraph_invalidate_cache(&multiple);

    /* Add two multiple edges */
    igraph_add_edge(&multiple, IGRAPH_FROM(&multiple, 42),
                    IGRAPH_TO(&multiple, 42));
    igraph_add_edge(&multiple, IGRAPH_FROM(&multiple, 42),
                    IGRAPH_TO(&multiple, 42));

    igraph_is_complete(&multiple, &complete);
    IGRAPH_ASSERT(complete);

    igraph_rng_seed(igraph_rng_default(), 42);

    /* Here, we cannot just simply remove any edge anymore */
    /* We must make sure that it is not a loop or an edge that has parallel
       edges */
    do {
        eid = igraph_rng_get_integer(igraph_rng_default(), 0,
                                     igraph_ecount(&multiple) - 1);
        igraph_es_pairs_small(&es, false, IGRAPH_FROM(&multiple, eid),
                              IGRAPH_TO(&multiple, eid), -1);

        igraph_es_size(&multiple, &es, &size);

    } while (size != 1 ||
             IGRAPH_FROM(&multiple, eid) == IGRAPH_TO(&multiple, eid));

    /* Remove the edge we found. There is one necessarily, provided vertices is
       high enough */
    /* We only added two loops, and two parallels */
    igraph_delete_edges(&multiple, es);

    igraph_is_complete(&multiple, &complete);
    IGRAPH_ASSERT(!complete);

    igraph_destroy(&simple);
    igraph_destroy(&multiple);
    igraph_es_destroy(&es);

    if (directed) {
        igraph_vector_bool_destroy(&mutuals);
        igraph_vector_int_destroy(&toadd);
    }

    return IGRAPH_SUCCESS;
}

int main(void)
{
    check_null();
    check_singleton();
    check(false, 50);
    check(false, 51);
    check(true, 50);
    check(true, 51);

    VERIFY_FINALLY_STACK();

    return 0;
}
