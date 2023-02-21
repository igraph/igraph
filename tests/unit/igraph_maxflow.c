/* IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>
#include "test_utilities.h"

int check_flow(int errorinc,
               const igraph_t *graph, igraph_real_t flow_value,
               const igraph_vector_t *flow, const igraph_vector_int_t *cut,
               const igraph_vector_int_t *partition,
               const igraph_vector_int_t *partition2,
               igraph_integer_t source, igraph_integer_t target,
               const igraph_vector_t *capacity,
               igraph_bool_t print) {

    igraph_integer_t i, n = igraph_vcount(graph), m = igraph_ecount(graph);
    igraph_integer_t nc = igraph_vector_int_size(cut);
    igraph_vector_int_t inedges, outedges;
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_real_t cutsize;
    igraph_t graph_copy;
    igraph_matrix_t sp;
    igraph_vector_int_t cut_int;

    igraph_vector_int_init(&cut_int, 0);

    if (print) {
        printf("flow value: %g\n", (double) flow_value);
        printf("flow: ");
        igraph_vector_print(flow);
        printf("first partition:  ");
        igraph_vector_int_print(partition);
        printf("second partition: ");
        igraph_vector_int_print(partition2);
        printf("edges in the cut: ");
        for (i = 0; i < nc; i++) {
            igraph_integer_t edge = VECTOR(*cut)[i];
            igraph_integer_t from = IGRAPH_FROM(graph, edge);
            igraph_integer_t to  = IGRAPH_TO  (graph, edge);
            if (!directed && from > to) {
                igraph_integer_t tmp = from;
                from = to;
                to = tmp;
            }
            printf("%" IGRAPH_PRId "-%" IGRAPH_PRId " (%g), ", from, to, VECTOR(*capacity)[edge]);
        }
        printf("\n");
    }
    fflush(stdout);

    /* Always less than the capacity */
    for (i = 0; i < m; i++) {
        if (VECTOR(*flow)[i] > VECTOR(*capacity)[i]) {
            return errorinc + 3;
        }
    }

    /* What comes in goes out, but only in directed graphs,
       there is no in and out in undirected ones...
     */
    if (igraph_is_directed(graph)) {
        igraph_vector_int_init(&inedges, 0);
        igraph_vector_int_init(&outedges, 0);

        for (i = 0; i < n; i++) {
            igraph_integer_t n1, n2, j;
            igraph_real_t in_flow = 0.0, out_flow = 0.0;
            igraph_incident(graph, &inedges,  i, IGRAPH_IN);
            igraph_incident(graph, &outedges, i, IGRAPH_OUT);
            n1 = igraph_vector_int_size(&inedges);
            n2 = igraph_vector_int_size(&outedges);
            for (j = 0; j < n1; j++) {
                igraph_integer_t e = VECTOR(inedges)[j];
                in_flow += VECTOR(*flow)[e];
            }
            for (j = 0; j < n2; j++) {
                igraph_integer_t e = VECTOR(outedges)[j];
                out_flow += VECTOR(*flow)[e];
            }
            if (i == source) {
                if (in_flow > 0) {
                    return errorinc + 4;
                }
                if (out_flow != flow_value) {
                    return errorinc + 5;
                }
            } else if (i == target) {
                if (out_flow > 0) {
                    return errorinc + 6;
                }
                if (in_flow != flow_value) {
                    return errorinc + 7;
                }

            } else {
                if (in_flow != out_flow) {
                    return errorinc + 8;
                }
            }
        }

        igraph_vector_int_destroy(&inedges);
        igraph_vector_int_destroy(&outedges);
    }

    /* Check the minimum cut size*/
    for (i = 0, cutsize = 0.0; i < nc; i++) {
        igraph_integer_t edge = VECTOR(*cut)[i];
        cutsize += VECTOR(*capacity)[edge];
    }
    if (fabs(cutsize - flow_value) > 1e-14) {
        return errorinc + 9;
    }

    /* Check that the cut indeed cuts */
    igraph_copy(&graph_copy, graph);

    n = igraph_vector_int_size(cut);
    igraph_vector_int_resize(&cut_int, n);
    for (i = 0; i < n; i++) {
        VECTOR(cut_int)[i] = VECTOR(*cut)[i];
    }
    igraph_delete_edges(&graph_copy, igraph_ess_vector(&cut_int));
    igraph_matrix_init(&sp, 1, 1);
    igraph_distances(&graph_copy, &sp, /*from=*/ igraph_vss_1(source),
                     /*to=*/ igraph_vss_1(target), IGRAPH_OUT);
    if (MATRIX(sp, 0, 0) != IGRAPH_INFINITY) {
        return errorinc + 10;
    }
    igraph_matrix_destroy(&sp);
    igraph_destroy(&graph_copy);

    igraph_vector_int_destroy(&cut_int);

    return 0;
}

int main(void) {

    igraph_t g;
    igraph_real_t flow_value;
    igraph_vector_int_t cut;
    igraph_vector_t capacity;
    igraph_vector_int_t partition, partition2;
    igraph_vector_t flow;
    igraph_integer_t i, n;
    igraph_integer_t source, target;
    FILE *infile;
    igraph_real_t flow_value2 = 0.0;
    int check;
    igraph_maxflow_stats_t stats;

    igraph_vector_init(&capacity, 0);

    /***************/
    infile = fopen("ak-4102.max", "r");
    igraph_read_graph_dimacs_flow(
        &g, infile, 0, 0, &source, &target, &capacity, IGRAPH_DIRECTED
    );
    fclose(infile);

    igraph_vector_int_init(&cut, 0);
    igraph_vector_int_init(&partition, 0);
    igraph_vector_int_init(&partition2, 0);
    igraph_vector_init(&flow, 0);

    igraph_maxflow(&g, &flow_value, &flow, &cut, &partition,
                   &partition2, source, target, &capacity, &stats);

    if (flow_value != 8207) {
        return 1;
    }

    n = igraph_vector_int_size(&cut);
    for (i = 0; i < n; i++) {
        igraph_integer_t e = VECTOR(cut)[i];
        flow_value2 += VECTOR(capacity)[e];
    }
    if (flow_value != flow_value2) {
        return 2;
    }

    /* Check the flow */
    if ( (check = check_flow(0, &g, flow_value, &flow, &cut, &partition,
                             &partition2, source, target, &capacity,
                             /*print=*/ 0))) {
        return check;
    }

    igraph_destroy(&g);
    igraph_vector_destroy(&capacity);
    igraph_vector_int_destroy(&cut);
    igraph_vector_int_destroy(&partition);
    igraph_vector_int_destroy(&partition2);
    igraph_vector_destroy(&flow);

    /* ------------------------------------- */

    igraph_small(&g, 4, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 1, 2, 1, 3, 2, 3, -1);
    igraph_vector_init_int_end(&capacity, -1, 4, 2, 10, 2, 2, -1);
    igraph_vector_int_init(&cut, 0);
    igraph_vector_int_init(&partition, 0);
    igraph_vector_int_init(&partition2, 0);
    igraph_vector_init(&flow, 0);

    igraph_maxflow(&g, &flow_value, &flow, &cut, &partition, &partition2,
                   /*source=*/ 0, /*target=*/ 3, &capacity, &stats);

    if ( (check = check_flow(20, &g, flow_value, &flow, &cut, &partition,
                             &partition2, 0, 3, &capacity,
                             /*print=*/ 1))) {
        return check;
    }

    igraph_vector_int_destroy(&cut);
    igraph_vector_int_destroy(&partition2);
    igraph_vector_int_destroy(&partition);
    igraph_vector_destroy(&capacity);
    igraph_vector_destroy(&flow);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
