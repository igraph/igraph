#include "prpack_igraph_graph.h"
#include <climits>
#include <cstring>
#include <new>

#include "igraph_interface.h"

using namespace prpack;
using namespace std;

#ifdef PRPACK_IGRAPH_SUPPORT

igraph_error_t prpack_igraph_graph::convert_from_igraph(
        const igraph_t *g, const igraph_vector_t *weights, bool directed) {

    const bool treat_as_directed = igraph_is_directed(g) && directed;
    const igraph_integer_t vcount = igraph_vcount(g);
    const igraph_integer_t ecount = igraph_ecount(g);
    double *p_weight;
    int *p_head;

    if (vcount > INT_MAX) {
        IGRAPH_ERROR("Too many vertices for PRPACK.", IGRAPH_EINVAL);
    }
    if (ecount > (treat_as_directed ? INT_MAX : INT_MAX/2)) {
        IGRAPH_ERROR("Too many edges for PRPACK.", IGRAPH_EINVAL);
    }

    if (weights && igraph_vector_size(weights) != ecount) {
        IGRAPH_ERROR("Weight vector length must agree with number of edges.", IGRAPH_EINVAL);
    }

    // Get the number of vertices and edges. For undirected graphs, we add
    // an edge in both directions.
    num_vs = (int) vcount;
    num_es = (int) ecount;
    num_self_es = 0;
    if (!treat_as_directed) {
        num_es *= 2;
    }

    // Allocate memory for heads and tails
    p_head = heads = new int[num_es];
    tails = new int[num_vs];
    memset(tails, 0, num_vs * sizeof(tails[0]));

    // Allocate memory for weights if needed
    if (weights) {
        p_weight = vals = new double[num_es];
    }

    // Count the number of ignored edges (those with negative or zero weight)
    int num_ignored_es = 0;

    if (treat_as_directed) {
        // Use of igraph "finally" stack is safe in this block
        // since no exceptions can be thrown from here.

        // Select all the edges and iterate over them by the source vertices
        // Add the edges
        igraph_eit_t eit;
        IGRAPH_CHECK(igraph_eit_create(g, igraph_ess_all(IGRAPH_EDGEORDER_TO), &eit));
        IGRAPH_FINALLY(igraph_eit_destroy, &eit);
        while (!IGRAPH_EIT_END(eit)) {
            igraph_integer_t eid = IGRAPH_EIT_GET(eit);
            IGRAPH_EIT_NEXT(eit);

            // Handle the weight
            if (weights != NULL) {
                // Does this edge have zero or negative weight?
                if (VECTOR(*weights)[eid] < 0) {
                    // Negative weights are disallowed.
                    IGRAPH_ERROR("Edge weights must not be negative.", IGRAPH_EINVAL);
                } else if (isnan(VECTOR(*weights)[eid])) {
                    IGRAPH_ERROR("Edge weights must not be NaN.", IGRAPH_EINVAL);
                } else if (VECTOR(*weights)[eid] == 0) {
                    // Edges with zero weight are ignored.
                    num_ignored_es++;
                    continue;
                }

                *p_weight = VECTOR(*weights)[eid];
                ++p_weight;
            }

            *p_head = IGRAPH_FROM(g, eid);
            ++p_head;
            ++tails[IGRAPH_TO(g, eid)];

            if (IGRAPH_FROM(g, eid) == IGRAPH_TO(g, eid)) {
                ++num_self_es;
            }
        }
        igraph_eit_destroy(&eit);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        // Use of igraph "finally" stack is safe in this block
        // since no exceptions can be thrown from here.

        // Select all the edges and iterate over them by the target vertices
        igraph_vector_int_t neis;
        IGRAPH_CHECK(igraph_vector_int_init(&neis, 0));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &neis);

        for (int i = 0; i < num_vs; i++) {
            IGRAPH_CHECK(igraph_incident(g, &neis, i, IGRAPH_ALL));

            int temp = igraph_vector_int_size(&neis);

            // TODO: should loop edges be added in both directions?
            int *p_head_copy = p_head;
            for (int j = 0; j < temp; j++) {
                if (weights != NULL) {
                    if (VECTOR(*weights)[VECTOR(neis)[j]] <= 0) {
                        // Ignore
                        num_ignored_es++;
                        continue;
                    }

                    *p_weight = VECTOR(*weights)[VECTOR(neis)[j]];
                    ++p_weight;
                }

                *p_head = IGRAPH_OTHER(g, VECTOR(neis)[j], i);
                if (i == *p_head) {
                    num_self_es++;
                }
                ++p_head;
            }
            tails[i] = p_head - p_head_copy;
        }

        igraph_vector_int_destroy(&neis);
        IGRAPH_FINALLY_CLEAN(1);
    }

    // Decrease num_es by the number of ignored edges
    num_es -= num_ignored_es;

    // Finalize the tails vector
    for (int i = 0, sum = 0; i < num_vs; ++i) {
        int temp = sum;
        sum += tails[i];
        tails[i] = temp;
    }

    // Normalize the weights
    normalize_weights();

    // Debug
    /*
    printf("Heads:");
    for (i = 0; i < num_es; ++i) {
        printf(" %d", heads[i]);
    }
    printf("\n");
    printf("Tails:");
    for (i = 0; i < num_vs; ++i) {
        printf(" %d", tails[i]);
    }
    printf("\n");
    if (vals) {
        printf("Vals:");
        for (i = 0; i < num_es; ++i) {
            printf(" %.4f", vals[i]);
        }
        printf("\n");
    }
    printf("===========================\n");
    */

    return IGRAPH_SUCCESS;
}

// PRPACK_IGRAPH_SUPPORT
#endif
