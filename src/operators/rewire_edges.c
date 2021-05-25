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

#include "graph/attributes.h"

static int igraph_i_rewire_edges_no_multiple(igraph_t *graph, igraph_real_t prob,
                                             igraph_bool_t loops,
                                             igraph_vector_t *edges) {

    int no_verts = igraph_vcount(graph);
    int no_edges = igraph_ecount(graph);
    igraph_vector_t eorder, tmp;
    igraph_vector_int_t first, next, prev, marked;
    int i, to_rewire, last_other = -1;

    /* Create our special graph representation */

# define ADD_STUB(vertex, stub) do {                \
        if (VECTOR(first)[(vertex)]) {              \
            VECTOR(prev)[(int) VECTOR(first)[(vertex)]-1]=(stub)+1;   \
        }                               \
        VECTOR(next)[(stub)]=VECTOR(first)[(vertex)];       \
        VECTOR(prev)[(stub)]=0;                 \
        VECTOR(first)[(vertex)]=(stub)+1;               \
    } while (0)

# define DEL_STUB(vertex, stub) do {                    \
        if (VECTOR(next)[(stub)]) {                     \
            VECTOR(prev)[VECTOR(next)[(stub)]-1]=VECTOR(prev)[(stub)];    \
        }                                   \
        if (VECTOR(prev)[(stub)]) {                     \
            VECTOR(next)[VECTOR(prev)[(stub)]-1]=VECTOR(next)[(stub)];    \
        } else {                                \
            VECTOR(first)[(vertex)]=VECTOR(next)[(stub)];         \
        }                                   \
    } while (0)

# define MARK_NEIGHBORS(vertex) do {                \
        int xxx_ =VECTOR(first)[(vertex)];              \
        while (xxx_) {                      \
            int o= (int) VECTOR(*edges)[xxx_ % 2 ? xxx_ : xxx_-2];    \
            VECTOR(marked)[o]=other+1;                \
            xxx_=VECTOR(next)[xxx_-1];                \
        }                               \
    } while (0)

    IGRAPH_CHECK(igraph_vector_int_init(&first, no_verts));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &first);
    IGRAPH_CHECK(igraph_vector_int_init(&next, no_edges * 2));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &next);
    IGRAPH_CHECK(igraph_vector_int_init(&prev, no_edges * 2));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &prev);
    IGRAPH_CHECK(igraph_get_edgelist(graph, edges, /*bycol=*/ 0));
    IGRAPH_VECTOR_INIT_FINALLY(&eorder, no_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, no_edges);
    for (i = 0; i < no_edges; i++) {
        int idx1 = 2 * i, idx2 = idx1 + 1,
            from = (int) VECTOR(*edges)[idx1], to = (int) VECTOR(*edges)[idx2];
        VECTOR(tmp)[i] = from;
        ADD_STUB(from, idx1);
        ADD_STUB(to, idx2);
    }
    IGRAPH_CHECK(igraph_vector_order1(&tmp, &eorder, no_verts));
    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_vector_int_init(&marked, no_verts));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &marked);

    /* Rewire the stubs, part I */

    to_rewire = (int) RNG_GEOM(prob);
    while (to_rewire < no_edges) {
        int stub = (int) (2 * VECTOR(eorder)[to_rewire] + 1);
        int v = (int) VECTOR(*edges)[stub];
        int ostub = stub - 1;
        int other = (int) VECTOR(*edges)[ostub];
        int pot;
        if (last_other != other) {
            MARK_NEIGHBORS(other);
        }
        /* Do the rewiring */
        do {
            if (loops) {
                pot = (int) RNG_INTEGER(0, no_verts - 1);
            } else {
                pot = (int) RNG_INTEGER(0, no_verts - 2);
                pot = pot != other ? pot : no_verts - 1;
            }
        } while (VECTOR(marked)[pot] == other + 1 && pot != v);

        if (pot != v) {
            DEL_STUB(v, stub);
            ADD_STUB(pot, stub);
            VECTOR(marked)[v] = 0;
            VECTOR(marked)[pot] = other + 1;
            VECTOR(*edges)[stub] = pot;
        }

        to_rewire += RNG_GEOM(prob) + 1;
        last_other = other;
    }

    /* Create the new index, from the potentially rewired stubs */

    IGRAPH_VECTOR_INIT_FINALLY(&tmp, no_edges);
    for (i = 0; i < no_edges; i++) {
        VECTOR(tmp)[i] = VECTOR(*edges)[2 * i + 1];
    }
    IGRAPH_CHECK(igraph_vector_order1(&tmp, &eorder, no_verts));
    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);

    /* Rewire the stubs, part II */

    igraph_vector_int_null(&marked);
    last_other = -1;

    to_rewire = (int) RNG_GEOM(prob);
    while (to_rewire < no_edges) {
        int stub = (int) (2 * VECTOR(eorder)[to_rewire]);
        int v = (int) VECTOR(*edges)[stub];
        int ostub = stub + 1;
        int other = (int) VECTOR(*edges)[ostub];
        int pot;
        if (last_other != other) {
            MARK_NEIGHBORS(other);
        }
        /* Do the rewiring */
        do {
            if (loops) {
                pot = (int) RNG_INTEGER(0, no_verts - 1);
            } else {
                pot = (int) RNG_INTEGER(0, no_verts - 2);
                pot = pot != other ? pot : no_verts - 1;
            }
        } while (VECTOR(marked)[pot] == other + 1 && pot != v);
        if (pot != v) {
            DEL_STUB(v, stub);
            ADD_STUB(pot, stub);
            VECTOR(marked)[v] = 0;
            VECTOR(marked)[pot] = other + 1;
            VECTOR(*edges)[stub] = pot;
        }

        to_rewire += RNG_GEOM(prob) + 1;
        last_other = other;
    }

    igraph_vector_int_destroy(&marked);
    igraph_vector_int_destroy(&prev);
    igraph_vector_int_destroy(&next);
    igraph_vector_int_destroy(&first);
    igraph_vector_destroy(&eorder);
    IGRAPH_FINALLY_CLEAN(5);

    return 0;
}

#undef ADD_STUB
#undef DEL_STUB
#undef MARK_NEIGHBORS

/**
 * \function igraph_rewire_edges
 * \brief Rewires the edges of a graph with constant probability.
 *
 * This function rewires the edges of a graph with a constant
 * probability. More precisely each end point of each edge is rewired
 * to a uniformly randomly chosen vertex with constant probability \p
 * prob.
 *
 * </para><para> Note that this function modifies the input \p graph,
 * call \ref igraph_copy() if you want to keep it.
 *
 * \param graph The input graph, this will be rewired, it can be
 *    directed or undirected.
 * \param prob The rewiring probability a constant between zero and
 *    one (inclusive).
 * \param loops Boolean, whether loop edges are allowed in the new
 *    graph, or not.
 * \param multiple Boolean, whether multiple edges are allowed in the
 *    new graph.
 * \return Error code.
 *
 * \sa \ref igraph_watts_strogatz_game() uses this function for the
 * rewiring.
 *
 * Time complexity: O(|V|+|E|).
 */
int igraph_rewire_edges(igraph_t *graph, igraph_real_t prob,
                        igraph_bool_t loops, igraph_bool_t multiple) {

    igraph_t newgraph;
    long int no_of_edges = igraph_ecount(graph);
    long int no_of_nodes = igraph_vcount(graph);
    long int endpoints = no_of_edges * 2;
    long int to_rewire;
    igraph_vector_t edges;

    if (prob < 0 || prob > 1) {
        IGRAPH_ERROR("Rewiring probability should be between zero and one",
                     IGRAPH_EINVAL);
    }

    if (prob == 0) {
        /* This is easy, just leave things as they are */
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, endpoints);

    RNG_BEGIN();

    if (prob != 0 && no_of_edges > 0) {
        if (multiple) {
            /* If multiple edges are allowed, then there is an easy and fast
            method. Each endpoint of an edge is rewired with probability p,
             so the "skips" between the really rewired endpoints follow a
             geometric distribution. */
            IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));
            to_rewire = (long int) RNG_GEOM(prob);
            while (to_rewire < endpoints) {
                if (loops) {
                    VECTOR(edges)[to_rewire] = RNG_INTEGER(0, no_of_nodes - 1);
                } else {
                    long int opos = to_rewire % 2 ? to_rewire - 1 : to_rewire + 1;
                    long int nei = (long int) VECTOR(edges)[opos];
                    long int r = RNG_INTEGER(0, no_of_nodes - 2);
                    VECTOR(edges)[ to_rewire ] = (r != nei ? r : no_of_nodes - 1);
                }
                to_rewire += RNG_GEOM(prob) + 1;
            }

        } else {
            IGRAPH_CHECK(igraph_i_rewire_edges_no_multiple(graph, prob, loops,
                         &edges));
        }
    }

    RNG_END();

    IGRAPH_CHECK(igraph_create(&newgraph, &edges, (igraph_integer_t) no_of_nodes,
                               igraph_is_directed(graph)));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_FINALLY(igraph_destroy, &newgraph);
    IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
    IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, 1, 1, 1);
    IGRAPH_FINALLY_CLEAN(1);
    igraph_destroy(graph);
    *graph = newgraph;

    return 0;
}

/**
 * \function igraph_rewire_directed_edges
 * \brief Rewires the chosen endpoint of directed edges.
 *
 * This function rewires either the start or end of directed edges in a graph
 * with a constant probability. Correspondingly, either the in-degree sequence
 * or the out-degree sequence of the graph will be preserved.
 *
 * </para><para> Note that this function modifies the input \p graph,
 * call \ref igraph_copy() if you want to keep it.
 *
 * \param graph The input graph, this will be rewired, it can be
 *    directed or undirected. If it is directed, \ref igraph_rewire_edges()
 *    will be called.
 * \param prob The rewiring probability, a constant between zero and
 *    one (inclusive).
 * \param loops Boolean, whether loop edges are allowed in the new
 *    graph, or not.
 * \param mode The endpoints of directed edges to rewire. It is ignored for
 *    undirected graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          rewire the end of each directed edge
 *        \cli IGRAPH_IN
 *          rewire the start of each directed edge
 *        \cli IGRAPH_ALL
 *          rewire both endpoints of each edge
 *        \endclist
 * \return Error code.
 *
 * \sa \ref igraph_rewire_edges(), \ref igraph_rewire()
 *
 * Time complexity: O(|E|).
 */
int igraph_rewire_directed_edges(igraph_t *graph, igraph_real_t prob,
                                 igraph_bool_t loops, igraph_neimode_t mode) {

    if (prob < 0 || prob > 1) {
        IGRAPH_ERROR("Rewiring probability should be between zero and one",
                     IGRAPH_EINVAL);
    }

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
    }

    if (prob == 0) {
        return IGRAPH_SUCCESS;
    }

    if (igraph_is_directed(graph) && mode != IGRAPH_ALL) {
        igraph_t newgraph;
        long int no_of_edges = igraph_ecount(graph);
        long int no_of_nodes = igraph_vcount(graph);
        long int to_rewire;
        long int offset = 0;
        igraph_vector_t edges;

        IGRAPH_VECTOR_INIT_FINALLY(&edges, 2 * no_of_edges);

        switch (mode) {
        case IGRAPH_IN:
            offset = 0;
            break;
        case IGRAPH_OUT:
            offset = 1;
            break;
        case IGRAPH_ALL:
            break; /* suppress compiler warning */
        }

        IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));

        RNG_BEGIN();

        to_rewire = RNG_GEOM(prob);
        while (to_rewire < no_of_edges) {
            if (loops) {
                VECTOR(edges)[2 * to_rewire + offset] = RNG_INTEGER(0, no_of_nodes - 1);
            } else {
                long int nei = (long int) VECTOR(edges)[2 * to_rewire + (1 - offset)];
                long int r = RNG_INTEGER(0, no_of_nodes - 2);
                VECTOR(edges)[2 * to_rewire + offset] = (r != nei ? r : no_of_nodes - 1);
            }
            to_rewire += RNG_GEOM(prob) + 1;
        }

        RNG_END();

        IGRAPH_CHECK(igraph_create(&newgraph, &edges, (igraph_integer_t) no_of_nodes,
                                   igraph_is_directed(graph)));
        igraph_vector_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);

        IGRAPH_FINALLY(igraph_destroy, &newgraph);
        IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
        IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, 1, 1, 1);
        IGRAPH_FINALLY_CLEAN(1);
        igraph_destroy(graph);
        *graph = newgraph;

    } else {
        IGRAPH_CHECK(igraph_rewire_edges(graph, prob, loops, /* multiple = */ 0));
    }

    return 0;
}
