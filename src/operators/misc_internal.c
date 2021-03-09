/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2020 The igraph development team

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

#include "operators/misc_internal.h"

#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_qsort.h"

void igraph_i_union_intersection_destroy_vectors(igraph_vector_ptr_t *v) {
    long int i, n = igraph_vector_ptr_size(v);
    for (i = 0; i < n; i++) {
        if (VECTOR(*v)[i] != 0) {
            igraph_vector_destroy(VECTOR(*v)[i]);
            IGRAPH_FREE(VECTOR(*v)[i]);
        }
    }
    igraph_vector_ptr_destroy(v);
}

void igraph_i_union_intersection_destroy_vector_longs(igraph_vector_ptr_t *v) {
    long int i, n = igraph_vector_ptr_size(v);
    for (i = 0; i < n; i++) {
        if (VECTOR(*v)[i] != 0) {
            igraph_vector_long_destroy(VECTOR(*v)[i]);
            IGRAPH_FREE(VECTOR(*v)[i]);
        }
    }
    igraph_vector_ptr_destroy(v);
}

int igraph_i_order_edgelist_cmp(void *edges, const void *e1, const void *e2) {
    igraph_vector_t *edgelist = edges;
    long int edge1 = (*(const long int*) e1) * 2;
    long int edge2 = (*(const long int*) e2) * 2;
    long int from1 = VECTOR(*edgelist)[edge1];
    long int from2 = VECTOR(*edgelist)[edge2];
    if (from1 < from2) {
        return -1;
    } else if (from1 > from2) {
        return 1;
    } else {
        long int to1 = VECTOR(*edgelist)[edge1 + 1];
        long int to2 = VECTOR(*edgelist)[edge2 + 1];
        if (to1 < to2) {
            return -1;
        } else if (to1 > to2) {
            return 1;
        } else {
            return 0;
        }
    }
}

int igraph_i_merge(igraph_t *res, int mode,
                   const igraph_t *left, const igraph_t *right,
                   igraph_vector_t *edge_map1, igraph_vector_t *edge_map2) {

    long int no_of_nodes_left = igraph_vcount(left);
    long int no_of_nodes_right = igraph_vcount(right);
    long int no_of_nodes;
    long int no_edges_left = igraph_ecount(left);
    long int no_edges_right = igraph_ecount(right);
    igraph_bool_t directed = igraph_is_directed(left);
    igraph_vector_t edges;
    igraph_vector_t edges1, edges2;
    igraph_vector_long_t order1, order2;
    long int i, j, eptr = 0;
    long int idx1, idx2, edge1 = -1, edge2 = -1, from1 = -1, from2 = -1, to1 = -1, to2 = -1;
    igraph_bool_t l;

    if (directed != igraph_is_directed(right)) {
        IGRAPH_ERROR("Cannot make union or intersection of directed "
                     "and undirected graph", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&edges1, no_edges_left * 2);
    IGRAPH_VECTOR_INIT_FINALLY(&edges2, no_edges_right * 2);
    IGRAPH_CHECK(igraph_vector_long_init(&order1, no_edges_left));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &order1);
    IGRAPH_CHECK(igraph_vector_long_init(&order2, no_edges_right));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &order2);

    if (edge_map1) {
        switch (mode) {
        case IGRAPH_MERGE_MODE_UNION:
            IGRAPH_CHECK(igraph_vector_resize(edge_map1, no_edges_left));
            break;
        case IGRAPH_MERGE_MODE_INTERSECTION:
            igraph_vector_clear(edge_map1);
            break;
        }
    }
    if (edge_map2) {
        switch (mode) {
        case IGRAPH_MERGE_MODE_UNION:
            IGRAPH_CHECK(igraph_vector_resize(edge_map2, no_edges_right));
            break;
        case IGRAPH_MERGE_MODE_INTERSECTION:
            igraph_vector_clear(edge_map2);
            break;
        }
    }

    no_of_nodes = no_of_nodes_left > no_of_nodes_right ?
                  no_of_nodes_left : no_of_nodes_right;

    /* We merge the two edge lists. We need to sort them first.
       For undirected graphs, we also need to make sure that
       for every edge, that larger (non-smaller) vertex id is in the
       second column. */

    IGRAPH_CHECK(igraph_get_edgelist(left, &edges1, /*bycol=*/ 0));
    IGRAPH_CHECK(igraph_get_edgelist(right, &edges2, /*bycol=*/ 0));
    if (!directed) {
        for (i = 0, j = 0; i < no_edges_left; i++, j += 2) {
            if (VECTOR(edges1)[j] > VECTOR(edges1)[j + 1]) {
                long int tmp = VECTOR(edges1)[j];
                VECTOR(edges1)[j] = VECTOR(edges1)[j + 1];
                VECTOR(edges1)[j + 1] = tmp;
            }
        }
        for (i = 0, j = 0; i < no_edges_right; i++, j += 2) {
            if (VECTOR(edges2)[j] > VECTOR(edges2)[j + 1]) {
                long int tmp = VECTOR(edges2)[j];
                VECTOR(edges2)[j] = VECTOR(edges2)[j + 1];
                VECTOR(edges2)[j + 1] = tmp;
            }
        }
    }

    for (i = 0; i < no_edges_left; i++) {
        VECTOR(order1)[i] = i;
    }
    for (i = 0; i < no_edges_right; i++) {
        VECTOR(order2)[i] = i;
    }

    igraph_qsort_r(VECTOR(order1), no_edges_left, sizeof(VECTOR(order1)[0]),
                   &edges1, igraph_i_order_edgelist_cmp);
    igraph_qsort_r(VECTOR(order2), no_edges_right, sizeof(VECTOR(order2)[0]),
                   &edges2, igraph_i_order_edgelist_cmp);

#define INC1() if ( (++idx1) < no_edges_left) {          \
        edge1 = VECTOR(order1)[idx1];                    \
        from1 = VECTOR(edges1)[2*edge1];                 \
        to1 = VECTOR(edges1)[2*edge1+1];                 \
    }
#define INC2() if ( (++idx2) < no_edges_right) {         \
        edge2 = VECTOR(order2)[idx2];                    \
        from2 = VECTOR(edges2)[2*edge2];                 \
        to2 = VECTOR(edges2)[2*edge2+1];                 \
    }

    idx1 = idx2 = -1;
    INC1();
    INC2();

#define CONT() switch (mode) {                               \
    case IGRAPH_MERGE_MODE_UNION:                            \
        l = idx1 < no_edges_left || idx2 < no_edges_right;   \
        break;                                               \
    case IGRAPH_MERGE_MODE_INTERSECTION:                     \
        l = idx1 < no_edges_left && idx2 < no_edges_right;   \
        break;                                               \
    }

    CONT();
    while (l) {
        if (idx2 >= no_edges_right ||
            (idx1 < no_edges_left && from1 < from2) ||
            (idx1 < no_edges_left && from1 == from2 && to1 < to2)) {
            /* Edge from first graph */
            if (mode == IGRAPH_MERGE_MODE_UNION) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, from1));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, to1));
                if (edge_map1) {
                    VECTOR(*edge_map1)[edge1] = eptr;
                }
                eptr++;
            }
            INC1();
        } else if (idx1 >= no_edges_left ||
                   (idx2 < no_edges_right && from2 < from1) ||
                   (idx2 < no_edges_right && from1 == from2 && to2 < to1)) {
            /* Edge from second graph */
            if (mode == IGRAPH_MERGE_MODE_UNION) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, from2));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, to2));
                if (edge_map2) {
                    VECTOR(*edge_map2)[edge2] = eptr;
                }
                eptr++;
            }
            INC2();
        } else {
            /* Edge from both */
            IGRAPH_CHECK(igraph_vector_push_back(&edges, from1));
            IGRAPH_CHECK(igraph_vector_push_back(&edges, to1));
            if (mode == IGRAPH_MERGE_MODE_UNION) {
                if (edge_map1) {
                    VECTOR(*edge_map1)[edge1] = eptr;
                }
                if (edge_map2) {
                    VECTOR(*edge_map2)[edge2] = eptr;
                }
            } else if (mode == IGRAPH_MERGE_MODE_INTERSECTION) {
                if (edge_map1) {
                    IGRAPH_CHECK(igraph_vector_push_back(edge_map1, edge1));
                }
                if (edge_map2) {
                    IGRAPH_CHECK(igraph_vector_push_back(edge_map2, edge2));
                }
            }
            eptr++;
            INC1();
            INC2();
        }
        CONT();
    }

#undef INC1
#undef INC2

    igraph_vector_long_destroy(&order2);
    igraph_vector_long_destroy(&order1);
    igraph_vector_destroy(&edges2);
    igraph_vector_destroy(&edges1);
    IGRAPH_FINALLY_CLEAN(4);

    IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}
