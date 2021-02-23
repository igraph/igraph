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
#include "igraph_random.h"

/**
 * \function igraph_correlated_game
 * \brief Generates a random graph correlated to an existing graph.
 *
 * Sample a new graph by perturbing the adjacency matrix of a
 * given graph and shuffling its vertices.
 *
 * \param old_graph The original graph.
 * \param new_graph The new graph will be stored here.
 * \param corr A scalar in the unit interval, the target Pearson
 *        correlation between the adjacency matrices of the original the
 *        generated graph (the adjacency matrix being used as a vector).
 * \param p A numeric scalar, the probability of an edge between two
 *        vertices, it must in the open (0,1) interval.
 * \param permutation A permutation to apply to the vertices of the
 *        generated graph. It can also be a null pointer, in which case
 *        the vertices will not be permuted.
 * \return Error code
 *
 * \sa \ref igraph_correlated_pair_game() for generating a pair
 * of correlated random graphs in one go.
 */
int igraph_correlated_game(const igraph_t *old_graph, igraph_t *new_graph,
                           igraph_real_t corr, igraph_real_t p,
                           const igraph_vector_t *permutation) {

    int no_of_nodes = igraph_vcount(old_graph);
    int no_of_edges = igraph_ecount(old_graph);
    igraph_bool_t directed = igraph_is_directed(old_graph);
    igraph_real_t no_of_all = directed ? no_of_nodes * (no_of_nodes - 1) :
                              no_of_nodes * (no_of_nodes - 1) / 2;
    igraph_real_t no_of_missing = no_of_all - no_of_edges;
    igraph_real_t q = p + corr * (1 - p);
    igraph_real_t p_del = 1 - q;
    igraph_real_t p_add = ((1 - q) * (p / (1 - p)));
    igraph_vector_t add, delete, edges, newedges;
    igraph_real_t last;
    int p_e = 0, p_a = 0, p_d = 0, no_add, no_del;
    igraph_real_t inf = IGRAPH_INFINITY;
    igraph_real_t next_e, next_a, next_d;
    int i;

    if (corr < -1 || corr > 1) {
        IGRAPH_ERROR("Correlation must be in [-1,1] in correlated "
                     "Erdos-Renyi game", IGRAPH_EINVAL);
    }
    if (p <= 0 || p >= 1) {
        IGRAPH_ERROR("Edge probability must be in (0,1) in correlated "
                     "Erdos-Renyi game", IGRAPH_EINVAL);
    }
    if (permutation) {
        if (igraph_vector_size(permutation) != no_of_nodes) {
            IGRAPH_ERROR("Invalid permutation length in correlated Erdos-Renyi game",
                         IGRAPH_EINVAL);
        }
    }

    /* Special cases */

    if (corr == 0) {
        return igraph_erdos_renyi_game(new_graph, IGRAPH_ERDOS_RENYI_GNP,
                                       no_of_nodes, p, directed,
                                       IGRAPH_NO_LOOPS);
    }
    if (corr == 1) {
        /* We don't copy, because we don't need the attributes.... */
        IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);
        IGRAPH_CHECK(igraph_get_edgelist(old_graph, &edges, /* bycol= */ 0));
        if (permutation) {
            int newec = igraph_vector_size(&edges);
            for (i = 0; i < newec; i++) {
                int tmp = VECTOR(edges)[i];
                VECTOR(edges)[i] = VECTOR(*permutation)[tmp];
            }
        }
        IGRAPH_CHECK(igraph_create(new_graph, &edges, no_of_nodes, directed));
        igraph_vector_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);
        return 0;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&newedges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&add, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&delete, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);

    IGRAPH_CHECK(igraph_get_edgelist(old_graph, &edges, /* bycol= */ 0));

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

    IGRAPH_CHECK(igraph_get_edgelist(old_graph, &edges, /* bycol= */ 0));

    /* Now we are merging the original edges, the edges that are removed,
       and the new edges. We have the following pointers:
       - p_a: the next edge to add
       - p_d: the next edge to delete
       - p_e: the next original edge
       - next_e: the code of the next edge in 'edges'
       - next_a: the code of the next edge to add
       - next_d: the code of the next edge to delete */

#define D_CODE(f,t) (((t)==no_of_nodes-1 ? f : t) * no_of_nodes + (f))
#define U_CODE(f,t) ((t) * ((t)-1) / 2 + (f))
#define CODE(f,t) (directed ? D_CODE(f,t) : U_CODE(f,t))
#define CODEE() (CODE(VECTOR(edges)[2*p_e], VECTOR(edges)[2*p_e+1]))

    /* First we (re)code the edges to delete */

    for (i = 0; i < no_del; i++) {
        int td = VECTOR(delete)[i];
        int from = VECTOR(edges)[2 * td];
        int to = VECTOR(edges)[2 * td + 1];
        VECTOR(delete)[i] = CODE(from, to);
    }

    IGRAPH_CHECK(igraph_vector_reserve(&newedges,
                                       (no_of_edges - no_del + no_add) * 2));

    /* Now we can do the merge. Additional edges are tricky, because
       the code must be shifted by the edges in the original graph. */

#define UPD_E()                             \
    { if (p_e < no_of_edges) { next_e=CODEE(); } else { next_e = inf; } }
#define UPD_A()                             \
{ if (p_a < no_add) { \
            next_a = VECTOR(add)[p_a] + p_e; } else { next_a = inf; } }
#define UPD_D()                             \
{ if (p_d < no_del) { \
            next_d = VECTOR(delete)[p_d]; } else { next_d = inf; } }

    UPD_E(); UPD_A(); UPD_D();

    while (next_e != inf || next_a != inf || next_d != inf) {
        if (next_e <= next_a && next_e < next_d) {

            /* keep an edge */
            IGRAPH_CHECK(igraph_vector_push_back(&newedges, VECTOR(edges)[2 * p_e]));
            IGRAPH_CHECK(igraph_vector_push_back(&newedges, VECTOR(edges)[2 * p_e + 1]));
            p_e ++; UPD_E(); UPD_A()

        } else if (next_e <= next_a && next_e == next_d) {

            /* delete an edge */
            p_e ++; UPD_E(); UPD_A();
            p_d++; UPD_D();

        } else {

            /* add an edge */
            int to, from;
            if (directed) {
                to = (int) floor(next_a / no_of_nodes);
                from = (int) (next_a - ((igraph_real_t)to) * no_of_nodes);
                if (from == to) {
                    to = no_of_nodes - 1;
                }
            } else {
                to = (int) floor((sqrt(8 * next_a + 1) + 1) / 2);
                from = (int) (next_a - (((igraph_real_t)to) * (to - 1)) / 2);
            }
            IGRAPH_CHECK(igraph_vector_push_back(&newedges, from));
            IGRAPH_CHECK(igraph_vector_push_back(&newedges, to));
            p_a++; UPD_A();

        }
    }

    igraph_vector_destroy(&edges);
    igraph_vector_destroy(&add);
    igraph_vector_destroy(&delete);
    IGRAPH_FINALLY_CLEAN(3);

    if (permutation) {
        int newec = igraph_vector_size(&newedges);
        for (i = 0; i < newec; i++) {
            int tmp = VECTOR(newedges)[i];
            VECTOR(newedges)[i] = VECTOR(*permutation)[tmp];
        }
    }

    IGRAPH_CHECK(igraph_create(new_graph, &newedges, no_of_nodes, directed));

    igraph_vector_destroy(&newedges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
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
int igraph_correlated_pair_game(igraph_t *graph1, igraph_t *graph2,
                                igraph_integer_t n, igraph_real_t corr, igraph_real_t p,
                                igraph_bool_t directed,
                                const igraph_vector_t *permutation) {

    IGRAPH_CHECK(igraph_erdos_renyi_game(graph1, IGRAPH_ERDOS_RENYI_GNP, n, p,
                                         directed, IGRAPH_NO_LOOPS));
    IGRAPH_CHECK(igraph_correlated_game(graph1, graph2, corr, p, permutation));
    return 0;
}
