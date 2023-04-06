/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2003-2021 The igraph development team

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

#include "igraph_games.h"

#include "igraph_conversion.h"
#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_qsort.h"
#include "igraph_random.h"
#include "igraph_structural.h"

#include "core/interruption.h"

/* The "code" of an edge is a single index representing its location in the adjacency matrix,
 * More specifically, the relevant parts of the adjacency matrix (i.e. non-diagonal in directed,
 * upper triangular in undirected) are column-wise concatenated into an array. The "code" is
 * the index in this array. We use floating point numbers for the code, as it can easily
 * exceed integers representable on 32 bits.
 */
#define D_CODE(f,t) (((t)==no_of_nodes-1 ? (f) : (t)) * no_of_nodes + (f))
#define U_CODE(f,t) ((t) * ((t)-1) / 2 + (f))
#define CODE(f,t) (directed ? D_CODE((double)(f),(double)(t)) : U_CODE((double)(f),(double)(t)))

/* TODO: Slight speedup may be possible if repeated vertex count queries are avoided. */
static int code_cmp(void *graph, const void *va, const void *vb) {
    const igraph_integer_t *a = (const igraph_integer_t *) va;
    const igraph_integer_t *b = (const igraph_integer_t *) vb;
    const igraph_integer_t no_of_nodes = igraph_vcount((igraph_t *) graph);
    const igraph_bool_t directed = igraph_is_directed((igraph_t *) graph);
    const igraph_real_t codea = CODE(a[0], a[1]);
    const igraph_real_t codeb = CODE(b[0], b[1]);
    if (codea < codeb) {
        return -1;
    } else if (codea > codeb) {
        return 1;
    } else {
        return 0;
    }
}

/* Sort an edge vector by edge codes. */
static void sort_edges(igraph_vector_int_t *edges, const igraph_t *graph) {
    igraph_qsort_r(VECTOR(*edges), igraph_vector_int_size(edges) / 2, 2*sizeof(igraph_integer_t), (void *) graph, code_cmp);
}

/**
 * \function igraph_correlated_game
 * \brief Generates a random graph correlated to an existing graph.
 *
 * Sample a new graph by perturbing the adjacency matrix of a
 * given simple graph and shuffling its vertices.
 *
 * \param old_graph The original graph, it must be simple.
 * \param new_graph The new graph will be stored here.
 * \param corr A scalar in the unit interval [0,1], the target Pearson
 *        correlation between the adjacency matrices of the original and the
 *        generated graph (the adjacency matrix being used as a vector).
 * \param p A numeric scalar, the probability of an edge between two
 *        vertices, it must in the open (0,1) interval. Typically,
 *        the density of \p old_graph.
 * \param permutation A permutation to apply to the vertices of the
 *        generated graph. It can also be a null pointer, in which case
 *        the vertices will not be permuted.
 * \return Error code
 *
 * \sa \ref igraph_correlated_pair_game() for generating a pair
 * of correlated random graphs in one go.
 */
igraph_error_t igraph_correlated_game(const igraph_t *old_graph, igraph_t *new_graph,
                           igraph_real_t corr, igraph_real_t p,
                           const igraph_vector_int_t *permutation) {

    igraph_integer_t no_of_nodes = igraph_vcount(old_graph);
    igraph_integer_t no_of_edges = igraph_ecount(old_graph);
    igraph_bool_t directed = igraph_is_directed(old_graph);
    igraph_real_t no_of_all = directed ? ((igraph_real_t) no_of_nodes) * (no_of_nodes - 1) :
                              ((igraph_real_t) no_of_nodes) * (no_of_nodes - 1) / 2;
    igraph_real_t no_of_missing = no_of_all - no_of_edges;
    igraph_real_t q = p + corr * (1 - p);
    igraph_real_t p_del = 1 - q;
    igraph_real_t p_add = ((1 - q) * (p / (1 - p)));
    igraph_vector_t add, delete;
    igraph_vector_int_t edges, newedges;
    igraph_real_t last;
    igraph_integer_t p_e = 0, p_a = 0, p_d = 0;
    igraph_integer_t no_add, no_del;
    igraph_real_t next_e, next_a, next_d;
    igraph_integer_t i, newec;
    igraph_bool_t simple;

    if (corr < 0 || corr > 1) {
        IGRAPH_ERRORF("Correlation must be in [0,1] in correlated Erdos-Renyi game, got %g.",
                      IGRAPH_EINVAL, corr);
    }
    if (p <= 0 || p >= 1) {
        IGRAPH_ERRORF("Edge probability must be in (0,1) in correlated Erdos-Renyi game, got %g.",
                      IGRAPH_EINVAL, p);
    }
    if (permutation) {
        if (igraph_vector_int_size(permutation) != no_of_nodes) {
            IGRAPH_ERROR("Invalid permutation length in correlated Erdos-Renyi game.",
                         IGRAPH_EINVAL);
        }
    }
    IGRAPH_CHECK(igraph_is_simple(old_graph, &simple));
    if (! simple) {
        IGRAPH_ERROR("The original graph must be simple for correlated Erdos-Renyi game.",
                     IGRAPH_EINVAL);
    }

    /* Special cases */

    if (corr == 0) {
        return igraph_erdos_renyi_game_gnp(new_graph, no_of_nodes, p, directed, IGRAPH_NO_LOOPS);
    }
    if (corr == 1) {
        /* We don't copy, because we don't need the attributes.... */
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges * 2);
        IGRAPH_CHECK(igraph_get_edgelist(old_graph, &edges, /* bycol= */ 0));
        if (permutation) {
            newec = igraph_vector_int_size(&edges);
            for (i = 0; i < newec; i++) {
                igraph_integer_t tmp = VECTOR(edges)[i];
                VECTOR(edges)[i] = VECTOR(*permutation)[tmp];
            }
        }
        IGRAPH_CHECK(igraph_create(new_graph, &edges, no_of_nodes, directed));
        igraph_vector_int_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&newedges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&add, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&delete, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges * 2);

    IGRAPH_CHECK(igraph_get_edgelist(old_graph, &edges, /* bycol= */ 0));
    /* The sampling method used is analogous to the one in igraph_erdos_renyi_game_gnp(),
     * and assumes that the edge list of the old graph is in order of increasing "codes".
     * Even IGRAPH_EDGEORDER_TO does not guarantee this, therefore we sort explicitly.
     */
    sort_edges(&edges, old_graph);

    RNG_BEGIN();

    if (p_del > 0) {
        last = RNG_GEOM(p_del);
        while (last < no_of_edges) {
            IGRAPH_CHECK(igraph_vector_push_back(&delete, last));
            last += RNG_GEOM(p_del);
            last += 1;
        }
    }
    no_del = igraph_vector_size(&delete);

    if (p_add > 0) {
        last = RNG_GEOM(p_add);
        while (last < no_of_missing) {
            IGRAPH_CHECK(igraph_vector_push_back(&add, last));
            last += RNG_GEOM(p_add);
            last += 1;
        }
    }
    no_add = igraph_vector_size(&add);

    RNG_END();

    /* Now we are merging the original edges, the edges that are removed,
       and the new edges. We have the following pointers:
       - p_a: the next edge to add
       - p_d: the next edge to delete
       - p_e: the next original edge
       - next_e: the code of the next edge in 'edges'
       - next_a: the code of the next edge to add
       - next_d: the code of the next edge to delete */

#define CODEE() (CODE(VECTOR(edges)[2*p_e], VECTOR(edges)[2*p_e+1]))

    /* First we (re)code the edges to delete */

    for (i = 0; i < no_del; i++) {
        igraph_integer_t td = VECTOR(delete)[i];
        igraph_integer_t from = VECTOR(edges)[2 * td];
        igraph_integer_t to = VECTOR(edges)[2 * td + 1];
        VECTOR(delete)[i] = CODE(from, to);
    }

    IGRAPH_CHECK(igraph_vector_int_reserve(&newedges,
                                       (no_of_edges - no_del + no_add) * 2));

    /* Now we can do the merge. Additional edges are tricky, because
       the code must be shifted by the edges in the original graph. */

#define UPD_E() \
    { if (p_e < no_of_edges) { next_e=CODEE(); } else { next_e = IGRAPH_INFINITY; } }
#define UPD_A() \
    { if (p_a < no_add) { \
            next_a = VECTOR(add)[p_a] + p_e; } else { next_a = IGRAPH_INFINITY; } }
#define UPD_D() \
    { if (p_d < no_del) { \
            next_d = VECTOR(delete)[p_d]; } else { next_d = IGRAPH_INFINITY; } }

    UPD_E(); UPD_A(); UPD_D();

    while (next_e != IGRAPH_INFINITY || next_a != IGRAPH_INFINITY || next_d != IGRAPH_INFINITY) {
        IGRAPH_ALLOW_INTERRUPTION();
        if (next_e <= next_a && next_e < next_d) {

            /* keep an edge */
            IGRAPH_CHECK(igraph_vector_int_push_back(&newedges, VECTOR(edges)[2 * p_e]));
            IGRAPH_CHECK(igraph_vector_int_push_back(&newedges, VECTOR(edges)[2 * p_e + 1]));
            p_e ++; UPD_E(); UPD_A()

        } else if (next_e <= next_a && next_e == next_d) {

            /* delete an edge */
            p_e ++; UPD_E(); UPD_A();
            p_d++; UPD_D();

        } else {

            /* add an edge */
            igraph_integer_t to, from;
            IGRAPH_ASSERT(isfinite(next_a));
            if (directed) {
                to = floor(next_a / no_of_nodes);
                from = next_a - ((igraph_real_t)to) * no_of_nodes;
                if (from == to) {
                    to = no_of_nodes - 1;
                }
            } else {
                to = floor((sqrt(8 * next_a + 1) + 1) / 2);
                from = next_a - (((igraph_real_t)to) * (to - 1)) / 2;
            }
            IGRAPH_CHECK(igraph_vector_int_push_back(&newedges, from));
            IGRAPH_CHECK(igraph_vector_int_push_back(&newedges, to));
            p_a++; UPD_A();

        }
    }

    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&add);
    igraph_vector_destroy(&delete);
    IGRAPH_FINALLY_CLEAN(3);

    if (permutation) {
        newec = igraph_vector_int_size(&newedges);
        for (i = 0; i < newec; i++) {
            igraph_integer_t tmp = VECTOR(newedges)[i];
            VECTOR(newedges)[i] = VECTOR(*permutation)[tmp];
        }
    }

    IGRAPH_CHECK(igraph_create(new_graph, &newedges, no_of_nodes, directed));

    igraph_vector_int_destroy(&newedges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

#undef D_CODE
#undef U_CODE
#undef CODE
#undef CODEE
#undef UPD_E
#undef UPD_A
#undef UPD_D

/**
 * \function igraph_correlated_pair_game
 * \brief Generates pairs of correlated random graphs.
 *
 * Sample two random graphs, with given correlation.
 *
 * \param graph1 The first graph will be stored here.
 * \param graph2 The second graph will be stored here.
 * \param n The number of vertices in both graphs.
 * \param corr A scalar in the unit interval, the target Pearson
 *        correlation between the adjacency matrices of the original the
 *        generated graph (the adjacency matrix being used as a vector).
 * \param p A numeric scalar, the probability of an edge between two
 *        vertices, it must in the open (0,1) interval.
 * \param directed Whether to generate directed graphs.
 * \param permutation A permutation to apply to the vertices of the
 *        second graph. It can also be a null pointer, in which case
 *        the vertices will not be permuted.
 * \return Error code
 *
 * \sa \ref igraph_correlated_game() for generating a correlated pair
 * to a given graph.
 */
igraph_error_t igraph_correlated_pair_game(igraph_t *graph1, igraph_t *graph2,
                                igraph_integer_t n, igraph_real_t corr, igraph_real_t p,
                                igraph_bool_t directed,
                                const igraph_vector_int_t *permutation) {

    IGRAPH_CHECK(igraph_erdos_renyi_game_gnp(graph1, n, p, directed, IGRAPH_NO_LOOPS));
    IGRAPH_CHECK(igraph_correlated_game(graph1, graph2, corr, p, permutation));
    return IGRAPH_SUCCESS;
}
