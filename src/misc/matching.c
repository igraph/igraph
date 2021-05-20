/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2012  Tamas Nepusz <ntamas@gmail.com>

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

#include "igraph_matching.h"

#include "igraph_adjlist.h"
#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_structural.h"

#include <math.h>

/* #define MATCHING_DEBUG */

#ifdef _MSC_VER
/* MSVC does not support variadic macros */
#include <stdarg.h>
static void debug(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
#ifdef MATCHING_DEBUG
    vfprintf(stderr, fmt, args);
#endif
    va_end(args);
}
#else
#ifdef MATCHING_DEBUG
    #define debug(...) fprintf(stderr, __VA_ARGS__)
#else
    #define debug(...)
#endif
#endif

/**
 * \function igraph_is_matching
 * Checks whether the given matching is valid for the given graph.
 *
 * This function checks a matching vector and verifies whether its length
 * matches the number of vertices in the given graph, its values are between
 * -1 (inclusive) and the number of vertices (exclusive), and whether there
 * exists a corresponding edge in the graph for every matched vertex pair.
 * For bipartite graphs, it also verifies whether the matched vertices are
 * in different parts of the graph.
 *
 * \param graph The input graph. It can be directed but the edge directions
 *              will be ignored.
 * \param types If the graph is bipartite and you are interested in bipartite
 *              matchings only, pass the vertex types here. If the graph is
 *              non-bipartite, simply pass \c NULL.
 * \param matching The matching itself. It must be a vector where element i
 *                 contains the ID of the vertex that vertex i is matched to,
 *                 or -1 if vertex i is unmatched.
 * \param result Pointer to a boolean variable, the result will be returned
 *               here.
 *
 * \sa \ref igraph_is_maximal_matching() if you are also interested in whether
 *     the matching is maximal (i.e. non-extendable).
 *
 * Time complexity: O(|V|+|E|) where |V| is the number of vertices and
 * |E| is the number of edges.
 *
 * \example examples/simple/igraph_maximum_bipartite_matching.c
 */
int igraph_is_matching(const igraph_t* graph,
                       const igraph_vector_bool_t* types, const igraph_vector_long_t* matching,
                       igraph_bool_t* result) {
    long int i, j, no_of_nodes = igraph_vcount(graph);
    igraph_bool_t conn;

    /* Checking match vector length */
    if (igraph_vector_long_size(matching) != no_of_nodes) {
        *result = 0; return IGRAPH_SUCCESS;
    }

    for (i = 0; i < no_of_nodes; i++) {
        j = VECTOR(*matching)[i];

        /* Checking range of each element in the match vector */
        if (j < -1 || j >= no_of_nodes) {
            *result = 0; return IGRAPH_SUCCESS;
        }
        /* When i is unmatched, we're done */
        if (j == -1) {
            continue;
        }
        /* Matches must be mutual */
        if (VECTOR(*matching)[j] != i) {
            *result = 0; return IGRAPH_SUCCESS;
        }
        /* Matched vertices must be connected */
        IGRAPH_CHECK(igraph_are_connected(graph, (igraph_integer_t) i,
                                          (igraph_integer_t) j, &conn));
        if (!conn) {
            /* Try the other direction -- for directed graphs */
            IGRAPH_CHECK(igraph_are_connected(graph, (igraph_integer_t) j,
                                              (igraph_integer_t) i, &conn));
            if (!conn) {
                *result = 0; return IGRAPH_SUCCESS;
            }
        }
    }

    if (types != 0) {
        /* Matched vertices must be of different types */
        for (i = 0; i < no_of_nodes; i++) {
            j = VECTOR(*matching)[i];
            if (j == -1) {
                continue;
            }
            if (VECTOR(*types)[i] == VECTOR(*types)[j]) {
                *result = 0; return IGRAPH_SUCCESS;
            }
        }
    }

    *result = 1;
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_maximal_matching
 * Checks whether a matching in a graph is maximal.
 *
 * A matching is maximal if and only if there exists no unmatched vertex in a
 * graph such that one of its neighbors is also unmatched.
 *
 * \param graph The input graph. It can be directed but the edge directions
 *              will be ignored.
 * \param types If the graph is bipartite and you are interested in bipartite
 *              matchings only, pass the vertex types here. If the graph is
 *              non-bipartite, simply pass \c NULL.
 * \param matching The matching itself. It must be a vector where element i
 *                 contains the ID of the vertex that vertex i is matched to,
 *                 or -1 if vertex i is unmatched.
 * \param result Pointer to a boolean variable, the result will be returned
 *               here.
 *
 * \sa \ref igraph_is_matching() if you are only interested in whether a
 *     matching vector is valid for a given graph.
 *
 * Time complexity: O(|V|+|E|) where |V| is the number of vertices and
 * |E| is the number of edges.
 *
 * \example examples/simple/igraph_maximum_bipartite_matching.c
 */
int igraph_is_maximal_matching(const igraph_t* graph,
                               const igraph_vector_bool_t* types, const igraph_vector_long_t* matching,
                               igraph_bool_t* result) {
    long int i, j, n, no_of_nodes = igraph_vcount(graph);
    igraph_vector_t neis;
    igraph_bool_t valid;

    IGRAPH_CHECK(igraph_is_matching(graph, types, matching, &valid));
    if (!valid) {
        *result = 0; return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

    valid = 1;
    for (i = 0; i < no_of_nodes; i++) {
        j = VECTOR(*matching)[i];
        if (j != -1) {
            continue;
        }

        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) i,
                                      IGRAPH_ALL));
        n = igraph_vector_size(&neis);
        for (j = 0; j < n; j++) {
            if (VECTOR(*matching)[(long int)VECTOR(neis)[j]] == -1) {
                if (types == 0 ||
                    VECTOR(*types)[i] != VECTOR(*types)[(long int)VECTOR(neis)[j]]) {
                    valid = 0; break;
                }
            }
        }
    }

    igraph_vector_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(1);

    *result = valid;
    return IGRAPH_SUCCESS;
}

static int igraph_i_maximum_bipartite_matching_unweighted(
        const igraph_t* graph,
        const igraph_vector_bool_t* types, igraph_integer_t* matching_size,
        igraph_vector_long_t* matching);
static int igraph_i_maximum_bipartite_matching_weighted(
        const igraph_t* graph,
        const igraph_vector_bool_t* types, igraph_integer_t* matching_size,
        igraph_real_t* matching_weight, igraph_vector_long_t* matching,
        const igraph_vector_t* weights, igraph_real_t eps);

#define MATCHED(v) (VECTOR(match)[v] != -1)
#define UNMATCHED(v) (!MATCHED(v))

/**
 * \function igraph_maximum_bipartite_matching
 * Calculates a maximum matching in a bipartite graph.
 *
 * A matching in a bipartite graph is a partial assignment of vertices
 * of the first kind to vertices of the second kind such that each vertex of
 * the first kind is matched to at most one vertex of the second kind and
 * vice versa, and matched vertices must be connected by an edge in the graph.
 * The size (or cardinality) of a matching is the number of edges.
 * A matching is a maximum matching if there exists no other matching with
 * larger cardinality. For weighted graphs, a maximum matching is a matching
 * whose edges have the largest possible total weight among all possible
 * matchings.
 *
 * </para><para>
 * Maximum matchings in bipartite graphs are found by the push-relabel algorithm
 * with greedy initialization and a global relabeling after every n/2 steps where
 * n is the number of vertices in the graph.
 *
 * </para><para>
 * References: Cherkassky BV, Goldberg AV, Martin P, Setubal JC and Stolfi J:
 * Augment or push: A computational study of bipartite matching and
 * unit-capacity flow algorithms. ACM Journal of Experimental Algorithmics 3,
 * 1998.
 *
 * </para><para>
 * Kaya K, Langguth J, Manne F and Ucar B: Experiments on push-relabel-based
 * maximum cardinality matching algorithms for bipartite graphs. Technical
 * Report TR/PA/11/33 of the Centre Europeen de Recherche et de Formation
 * Avancee en Calcul Scientifique, 2011.
 *
 * \param graph The input graph. It can be directed but the edge directions
 *              will be ignored.
 * \param types Boolean vector giving the vertex types of the graph.
 * \param matching_size The size of the matching (i.e. the number of matched
 *                      vertex pairs will be returned here). It may be \c NULL
 *                      if you don't need this.
 * \param matching_weight The weight of the matching if the edges are weighted,
 *                        or the size of the matching again if the edges are
 *                        unweighted. It may be \c NULL if you don't need this.
 * \param matching The matching itself. It must be a vector where element i
 *                 contains the ID of the vertex that vertex i is matched to,
 *                 or -1 if vertex i is unmatched.
 * \param weights A null pointer (=no edge weights), or a vector giving the
 *                weights of the edges. Note that the algorithm is stable
 *                only for integer weights.
 * \param eps A small real number used in equality tests in the weighted
 *            bipartite matching algorithm. Two real numbers are considered
 *            equal in the algorithm if their difference is smaller than
 *            \c eps. This is required to avoid the accumulation of numerical
 *            errors. It is advised to pass a value derived from the
 *            \c DBL_EPSILON constant in \c float.h here. If you are
 *            running the algorithm with no \c weights vector, this argument
 *            is ignored.
 * \return Error code.
 *
 * Time complexity: O(sqrt(|V|) |E|) for unweighted graphs (according to the
 * technical report referenced above), O(|V||E|) for weighted graphs.
 *
 * \example examples/simple/igraph_maximum_bipartite_matching.c
 */
int igraph_maximum_bipartite_matching(const igraph_t* graph,
                                      const igraph_vector_bool_t* types, igraph_integer_t* matching_size,
                                      igraph_real_t* matching_weight, igraph_vector_long_t* matching,
                                      const igraph_vector_t* weights, igraph_real_t eps) {

    /* Sanity checks */
    if (igraph_vector_bool_size(types) < igraph_vcount(graph)) {
        IGRAPH_ERROR("types vector too short", IGRAPH_EINVAL);
    }
    if (weights && igraph_vector_size(weights) < igraph_ecount(graph)) {
        IGRAPH_ERROR("weights vector too short", IGRAPH_EINVAL);
    }

    if (weights == 0) {
        IGRAPH_CHECK(igraph_i_maximum_bipartite_matching_unweighted(graph, types,
                     matching_size, matching));
        if (matching_weight != 0) {
            *matching_weight = *matching_size;
        }
        return IGRAPH_SUCCESS;
    } else {
        IGRAPH_CHECK(igraph_i_maximum_bipartite_matching_weighted(graph, types,
                     matching_size, matching_weight, matching, weights, eps));
        return IGRAPH_SUCCESS;
    }
}

static int igraph_i_maximum_bipartite_matching_unweighted_relabel(
        const igraph_t* graph,
        const igraph_vector_bool_t* types, igraph_vector_t* labels,
        igraph_vector_long_t* matching, igraph_bool_t smaller_set);

/**
 * Finding maximum bipartite matchings on bipartite graphs using the
 * push-relabel algorithm.
 *
 * The implementation follows the pseudocode in Algorithm 1 of the
 * following paper:
 *
 * Kaya K, Langguth J, Manne F and Ucar B: Experiments on push-relabel-based
 * maximum cardinality matching algorithms for bipartite graphs. Technical
 * Report TR/PA/11/33 of CERFACS (Centre Européen de Recherche et de Formation
 * Avancée en Calcul Scientifique).
 * http://www.cerfacs.fr/algor/reports/2011/TR_PA_11_33.pdf
 */
static int igraph_i_maximum_bipartite_matching_unweighted(
        const igraph_t* graph,
        const igraph_vector_bool_t* types, igraph_integer_t* matching_size,
        igraph_vector_long_t* matching) {
    long int i, j, k, n, no_of_nodes = igraph_vcount(graph);
    long int num_matched;             /* number of matched vertex pairs */
    igraph_vector_long_t match;       /* will store the matching */
    igraph_vector_t labels;           /* will store the labels */
    igraph_vector_t neis;             /* used to retrieve the neighbors of a node */
    igraph_dqueue_long_t q;           /* a FIFO for push ordering */
    igraph_bool_t smaller_set;        /* denotes which part of the bipartite graph is smaller */
    long int label_changed = 0;       /* Counter to decide when to run a global relabeling */
    long int relabeling_freq = no_of_nodes / 2;

    /* We will use:
     * - FIFO push ordering
     * - global relabeling frequency: n/2 steps where n is the number of nodes
     * - simple greedy matching for initialization
     */

    /* (1) Initialize data structures */
    IGRAPH_CHECK(igraph_vector_long_init(&match, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &match);
    IGRAPH_VECTOR_INIT_FINALLY(&labels, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_dqueue_long_init(&q, 0));
    IGRAPH_FINALLY(igraph_dqueue_long_destroy, &q);

    /* (2) Initially, every node is unmatched */
    igraph_vector_long_fill(&match, -1);

    /* (3) Find an initial matching in a greedy manner.
     *     At the same time, find which side of the graph is smaller. */
    num_matched = 0; j = 0;
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(*types)[i]) {
            j++;
        }
        if (MATCHED(i)) {
            continue;
        }
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) i,
                                      IGRAPH_ALL));
        n = igraph_vector_size(&neis);
        for (j = 0; j < n; j++) {
            k = (long int) VECTOR(neis)[j];
            if (VECTOR(*types)[k] == VECTOR(*types)[i]) {
                IGRAPH_ERROR("Graph is not bipartite with supplied types vector", IGRAPH_EINVAL);
            }
            if (UNMATCHED(k)) {
                /* We match vertex i to vertex VECTOR(neis)[j] */
                VECTOR(match)[k] = i;
                VECTOR(match)[i] = k;
                num_matched++;
                break;
            }
        }
    }
    smaller_set = (j <= no_of_nodes / 2);

    /* (4) Set the initial labeling -- lines 1 and 2 in the tech report */
    IGRAPH_CHECK(igraph_i_maximum_bipartite_matching_unweighted_relabel(
                     graph, types, &labels, &match, smaller_set));

    /* (5) Fill the push queue with the unmatched nodes from the smaller set. */
    for (i = 0; i < no_of_nodes; i++) {
        if (UNMATCHED(i) && VECTOR(*types)[i] == smaller_set) {
            IGRAPH_CHECK(igraph_dqueue_long_push(&q, i));
        }
    }

    /* (6) Main loop from the referenced tech report -- lines 4--13 */
    label_changed = 0;
    while (!igraph_dqueue_long_empty(&q)) {
        long int v = igraph_dqueue_long_pop(&q);             /* Line 13 */
        long int u = -1, label_u = 2 * no_of_nodes;
        long int w;

        if (label_changed >= relabeling_freq) {
            /* Run global relabeling */
            IGRAPH_CHECK(igraph_i_maximum_bipartite_matching_unweighted_relabel(
                             graph, types, &labels, &match, smaller_set));
            label_changed = 0;
        }

        debug("Considering vertex %ld\n", v);

        /* Line 5: find row u among the neighbors of v s.t. label(u) is minimal */
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) v,
                                      IGRAPH_ALL));
        n = igraph_vector_size(&neis);
        for (i = 0; i < n; i++) {
            if (VECTOR(labels)[(long int)VECTOR(neis)[i]] < label_u) {
                u = (long int) VECTOR(neis)[i];
                label_u = (long int) VECTOR(labels)[u];
                label_changed++;
            }
        }

        debug("  Neighbor with smallest label: %ld (label=%ld)\n", u, label_u);

        if (label_u < no_of_nodes) {                         /* Line 6 */
            VECTOR(labels)[v] = VECTOR(labels)[u] + 1;         /* Line 7 */
            if (MATCHED(u)) {                                  /* Line 8 */
                w = VECTOR(match)[u];
                debug("  Vertex %ld is matched to %ld, performing a double push\n", u, w);
                if (w != v) {
                    VECTOR(match)[u] = -1; VECTOR(match)[w] = -1;  /* Line 9 */
                    IGRAPH_CHECK(igraph_dqueue_long_push(&q, w));  /* Line 10 */
                    debug("  Unmatching & activating vertex %ld\n", w);
                    num_matched--;
                }
            }
            VECTOR(match)[u] = v; VECTOR(match)[v] = u;      /* Line 11 */
            num_matched++;
            VECTOR(labels)[u] += 2;                          /* Line 12 */
            label_changed++;
        }
    }

    /* Fill the output parameters */
    if (matching != 0) {
        IGRAPH_CHECK(igraph_vector_long_update(matching, &match));
    }
    if (matching_size != 0) {
        *matching_size = (igraph_integer_t) num_matched;
    }

    /* Release everything */
    igraph_dqueue_long_destroy(&q);
    igraph_vector_destroy(&neis);
    igraph_vector_destroy(&labels);
    igraph_vector_long_destroy(&match);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

static int igraph_i_maximum_bipartite_matching_unweighted_relabel(
        const igraph_t* graph,
        const igraph_vector_bool_t* types, igraph_vector_t* labels,
        igraph_vector_long_t* match, igraph_bool_t smaller_set) {
    long int i, j, n, no_of_nodes = igraph_vcount(graph), matched_to;
    igraph_dqueue_long_t q;
    igraph_vector_t neis;

    debug("Running global relabeling.\n");

    /* Set all the labels to no_of_nodes first */
    igraph_vector_fill(labels, no_of_nodes);

    /* Allocate vector for neighbors */
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

    /* Create a FIFO for the BFS and initialize it with the unmatched rows
     * (i.e. members of the larger set) */
    IGRAPH_CHECK(igraph_dqueue_long_init(&q, 0));
    IGRAPH_FINALLY(igraph_dqueue_long_destroy, &q);
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(*types)[i] != smaller_set && VECTOR(*match)[i] == -1) {
            IGRAPH_CHECK(igraph_dqueue_long_push(&q, i));
            VECTOR(*labels)[i] = 0;
        }
    }

    /* Run the BFS */
    while (!igraph_dqueue_long_empty(&q)) {
        long int v = igraph_dqueue_long_pop(&q);
        long int w;

        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) v,
                                      IGRAPH_ALL));

        n = igraph_vector_size(&neis);
        for (j = 0; j < n; j++) {
            w = (long int) VECTOR(neis)[j];
            if (VECTOR(*labels)[w] == no_of_nodes) {
                VECTOR(*labels)[w] = VECTOR(*labels)[v] + 1;
                matched_to = VECTOR(*match)[w];
                if (matched_to != -1 && VECTOR(*labels)[matched_to] == no_of_nodes) {
                    IGRAPH_CHECK(igraph_dqueue_long_push(&q, matched_to));
                    VECTOR(*labels)[matched_to] = VECTOR(*labels)[w] + 1;
                }
            }
        }
    }

    igraph_dqueue_long_destroy(&q);
    igraph_vector_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * Finding maximum bipartite matchings on bipartite graphs using the
 * Hungarian algorithm (a.k.a. Kuhn-Munkres algorithm).
 *
 * The algorithm uses a maximum cardinality matching on a subset of
 * tight edges as a starting point. This is achieved by
 * \c igraph_i_maximum_bipartite_matching_unweighted on the restricted
 * graph.
 *
 * The algorithm works reliably only if the weights are integers. The
 * \c eps parameter should specity a very small number; if the slack on
 * an edge falls below \c eps, it will be considered tight. If all your
 * weights are integers, you can safely set \c eps to zero.
 */
static int igraph_i_maximum_bipartite_matching_weighted(
        const igraph_t* graph,
        const igraph_vector_bool_t* types, igraph_integer_t* matching_size,
        igraph_real_t* matching_weight, igraph_vector_long_t* matching,
        const igraph_vector_t* weights, igraph_real_t eps) {
    long int i, j, k, n, no_of_nodes, no_of_edges;
    igraph_integer_t u, v, w, msize;
    igraph_t newgraph;
    igraph_vector_long_t match;       /* will store the matching */
    igraph_vector_t slack;            /* will store the slack on each edge */
    igraph_vector_t parent;           /* parent vertices during a BFS */
    igraph_vector_t vec1, vec2;       /* general temporary vectors */
    igraph_vector_t labels;           /* will store the labels */
    igraph_dqueue_long_t q;           /* a FIFO for BST */
    igraph_bool_t smaller_set_type;   /* denotes which part of the bipartite graph is smaller */
    igraph_vector_t smaller_set;      /* stores the vertex IDs of the smaller set */
    igraph_vector_t larger_set;       /* stores the vertex IDs of the larger set */
    long int smaller_set_size;        /* size of the smaller set */
    long int larger_set_size;         /* size of the larger set */
    igraph_real_t dual;               /* solution of the dual problem */
    igraph_adjlist_t tight_phantom_edges; /* adjacency list to manage tight phantom edges */
    igraph_integer_t alternating_path_endpoint;
    igraph_vector_int_t* neis;
    igraph_vector_int_t *neis2;
    igraph_inclist_t inclist;         /* incidence list of the original graph */

    /* The Hungarian algorithm is originally for complete bipartite graphs.
     * For non-complete bipartite graphs, a phantom edge of weight zero must be
     * added between every pair of non-connected vertices. We don't do this
     * explicitly of course. See the comments below about how phantom edges
     * are taken into account. */

    no_of_nodes = igraph_vcount(graph);
    no_of_edges = igraph_ecount(graph);
    if (eps < 0) {
        IGRAPH_WARNING("negative epsilon given, clamping to zero");
        eps = 0;
    }

    /* (1) Initialize data structures */
    IGRAPH_CHECK(igraph_vector_long_init(&match, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &match);
    IGRAPH_CHECK(igraph_vector_init(&slack, no_of_edges));
    IGRAPH_FINALLY(igraph_vector_destroy, &slack);
    IGRAPH_VECTOR_INIT_FINALLY(&vec1, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&vec2, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&labels, no_of_nodes);
    IGRAPH_CHECK(igraph_dqueue_long_init(&q, 0));
    IGRAPH_FINALLY(igraph_dqueue_long_destroy, &q);
    IGRAPH_VECTOR_INIT_FINALLY(&parent, no_of_nodes);
    IGRAPH_CHECK(igraph_adjlist_init_empty(&tight_phantom_edges,
                                           (igraph_integer_t) no_of_nodes));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &tight_phantom_edges);
    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    IGRAPH_VECTOR_INIT_FINALLY(&smaller_set, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&larger_set, 0);

    /* (2) Find which set is the smaller one */
    j = 0;
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(*types)[i] == 0) {
            j++;
        }
    }
    smaller_set_type = (j > no_of_nodes / 2);
    smaller_set_size = smaller_set_type ? (no_of_nodes - j) : j;
    larger_set_size = no_of_nodes - smaller_set_size;
    IGRAPH_CHECK(igraph_vector_reserve(&smaller_set, smaller_set_size));
    IGRAPH_CHECK(igraph_vector_reserve(&larger_set, larger_set_size));
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(*types)[i] == smaller_set_type) {
            IGRAPH_CHECK(igraph_vector_push_back(&smaller_set, i));
        } else {
            IGRAPH_CHECK(igraph_vector_push_back(&larger_set, i));
        }
    }

    /* (3) Calculate the initial labeling and the set of tight edges. Use the
     *     smaller set only. Here we can assume that there are no phantom edges
     *     among the tight ones. */
    dual = 0;
    for (i = 0; i < no_of_nodes; i++) {
        igraph_real_t max_weight = 0;

        if (VECTOR(*types)[i] != smaller_set_type) {
            VECTOR(labels)[i] = 0;
            continue;
        }

        neis = igraph_inclist_get(&inclist, i);
        n = igraph_vector_int_size(neis);
        for (j = 0, k = 0; j < n; j++) {
            k = (long int) VECTOR(*neis)[j];
            u = IGRAPH_OTHER(graph, k, i);
            if (VECTOR(*types)[u] == VECTOR(*types)[i]) {
                IGRAPH_ERROR("Graph is not bipartite with supplied types vector", IGRAPH_EINVAL);
            }
            if (VECTOR(*weights)[k] > max_weight) {
                max_weight = VECTOR(*weights)[k];
            }
        }

        VECTOR(labels)[i] = max_weight;
        dual += max_weight;
    }

    igraph_vector_clear(&vec1);
    IGRAPH_CHECK(igraph_get_edgelist(graph, &vec2, 0));
#define IS_TIGHT(i) (VECTOR(slack)[i] <= eps)
    for (i = 0, j = 0; i < no_of_edges; i++, j += 2) {
        u = (igraph_integer_t) VECTOR(vec2)[j];
        v = (igraph_integer_t) VECTOR(vec2)[j + 1];
        VECTOR(slack)[i] = VECTOR(labels)[u] + VECTOR(labels)[v] - VECTOR(*weights)[i];
        if (IS_TIGHT(i)) {
            IGRAPH_CHECK(igraph_vector_push_back(&vec1, u));
            IGRAPH_CHECK(igraph_vector_push_back(&vec1, v));
        }
    }
    igraph_vector_clear(&vec2);

    /* (4) Construct a temporary graph on which the initial maximum matching
     *     will be calculated (only on the subset of tight edges) */
    IGRAPH_CHECK(igraph_create(&newgraph, &vec1,
                               (igraph_integer_t) no_of_nodes, 0));
    IGRAPH_FINALLY(igraph_destroy, &newgraph);
    IGRAPH_CHECK(igraph_maximum_bipartite_matching(&newgraph, types, &msize, 0, &match, 0, 0));
    igraph_destroy(&newgraph);
    IGRAPH_FINALLY_CLEAN(1);

    /* (5) Main loop until the matching becomes maximal */
    while (msize < smaller_set_size) {
        igraph_real_t min_slack, min_slack_2;
        igraph_integer_t min_slack_u, min_slack_v;

        /* mark min_slack_u as unused; it is actually used when debugging, but
         * gcc complains when we are not debugging */
        IGRAPH_UNUSED(min_slack_u);

        /* (7) Fill the push queue with the unmatched nodes from the smaller set. */
        igraph_vector_clear(&vec1);
        igraph_vector_clear(&vec2);
        igraph_vector_fill(&parent, -1);
        for (j = 0; j < smaller_set_size; j++) {
            i = VECTOR(smaller_set)[j];
            if (UNMATCHED(i)) {
                IGRAPH_CHECK(igraph_dqueue_long_push(&q, i));
                VECTOR(parent)[i] = i;
                IGRAPH_CHECK(igraph_vector_push_back(&vec1, i));
            }
        }

#ifdef MATCHING_DEBUG
        debug("Matching:");
        igraph_vector_long_print(&match);
        debug("Unmatched vertices are marked by non-negative numbers:\n");
        igraph_vector_print(&parent);
        debug("Labeling:");
        igraph_vector_print(&labels);
        debug("Slacks:");
        igraph_vector_print(&slack);
#endif

        /* (8) Run the BFS */
        alternating_path_endpoint = -1;
        while (!igraph_dqueue_long_empty(&q)) {
            v = (int) igraph_dqueue_long_pop(&q);

            debug("Considering vertex %ld\n", (long int)v);

            /* v is always in the smaller set. Find the neighbors of v, which
             * are all in the larger set. Find the pairs of these nodes in
             * the smaller set and push them to the queue. Mark the traversed
             * nodes as seen.
             *
             * Here we have to be careful as there are two types of incident
             * edges on v: real edges and phantom ones. Real edges are
             * given by igraph_inclist_get. Phantom edges are not given so we
             * (ab)use an adjacency list data structure that lists the
             * vertices connected to v by phantom edges only. */
            neis = igraph_inclist_get(&inclist, v);
            n = igraph_vector_int_size(neis);
            for (i = 0; i < n; i++) {
                j = (long int) VECTOR(*neis)[i];
                /* We only care about tight edges */
                if (!IS_TIGHT(j)) {
                    continue;
                }
                /* Have we seen the other endpoint already? */
                u = IGRAPH_OTHER(graph, j, v);
                if (VECTOR(parent)[u] >= 0) {
                    continue;
                }
                debug("  Reached vertex %ld via edge %ld\n", (long)u, (long)j);
                VECTOR(parent)[u] = v;
                IGRAPH_CHECK(igraph_vector_push_back(&vec2, u));
                w = (int) VECTOR(match)[u];
                if (w == -1) {
                    /* u is unmatched and it is in the larger set. Therefore, we
                     * could improve the matching by following the parents back
                     * from u to the root.
                     */
                    alternating_path_endpoint = u;
                    break;  /* since we don't need any more endpoints that come from v */
                } else {
                    IGRAPH_CHECK(igraph_dqueue_long_push(&q, w));
                    VECTOR(parent)[w] = u;
                }
                IGRAPH_CHECK(igraph_vector_push_back(&vec1, w));
            }

            /* Now do the same with the phantom edges */
            neis2 = igraph_adjlist_get(&tight_phantom_edges, v);
            n = igraph_vector_int_size(neis2);
            for (i = 0; i < n; i++) {
                u = (igraph_integer_t) VECTOR(*neis2)[i];
                /* Have we seen u already? */
                if (VECTOR(parent)[u] >= 0) {
                    continue;
                }
                /* Check if the edge is really tight; it might have happened that the
                 * edge became non-tight in the meanwhile. We do not remove these from
                 * tight_phantom_edges at the moment, so we check them once again here.
                 */
                if (fabs(VECTOR(labels)[(long int)v] + VECTOR(labels)[(long int)u]) > eps) {
                    continue;
                }
                debug("  Reached vertex %ld via tight phantom edge\n", (long)u);
                VECTOR(parent)[u] = v;
                IGRAPH_CHECK(igraph_vector_push_back(&vec2, u));
                w = (int) VECTOR(match)[u];
                if (w == -1) {
                    /* u is unmatched and it is in the larger set. Therefore, we
                     * could improve the matching by following the parents back
                     * from u to the root.
                     */
                    alternating_path_endpoint = u;
                    break;  /* since we don't need any more endpoints that come from v */
                } else {
                    IGRAPH_CHECK(igraph_dqueue_long_push(&q, w));
                    VECTOR(parent)[w] = u;
                }
                IGRAPH_CHECK(igraph_vector_push_back(&vec1, w));
            }
        }

        /* Okay; did we have an alternating path? */
        if (alternating_path_endpoint != -1) {
#ifdef MATCHING_DEBUG
            debug("BFS parent tree:");
            igraph_vector_print(&parent);
#endif
            /* Increase the size of the matching with the alternating path. */
            v = alternating_path_endpoint;
            u = (igraph_integer_t) VECTOR(parent)[v];
            debug("Extending matching with alternating path ending in %ld.\n", (long int)v);

            while (u != v) {
                w = (int) VECTOR(match)[v];
                if (w != -1) {
                    VECTOR(match)[w] = -1;
                }
                VECTOR(match)[v] = u;

                VECTOR(match)[v] = u;
                w = (int) VECTOR(match)[u];
                if (w != -1) {
                    VECTOR(match)[w] = -1;
                }
                VECTOR(match)[u] = v;

                v = (igraph_integer_t) VECTOR(parent)[u];
                u = (igraph_integer_t) VECTOR(parent)[v];
            }

            msize++;

#ifdef MATCHING_DEBUG
            debug("New matching after update:");
            igraph_vector_long_print(&match);
            debug("Matching size is now: %ld\n", (long)msize);
#endif
            continue;
        }

#ifdef MATCHING_DEBUG
        debug("Vertices reachable from unmatched ones via tight edges:\n");
        igraph_vector_print(&vec1);
        igraph_vector_print(&vec2);
#endif

        /* At this point, vec1 contains the nodes in the smaller set (A)
         * reachable from unmatched nodes in A via tight edges only, while vec2
         * contains the nodes in the larger set (B) reachable from unmatched
         * nodes in A via tight edges only. Also, parent[i] >= 0 if node i
         * is reachable */

        /* Check the edges between reachable nodes in A and unreachable
         * nodes in B, and find the minimum slack on them.
         *
         * Since the weights are positive, we do no harm if we first
         * assume that there are no "real" edges between the two sets
         * mentioned above and determine an upper bound for min_slack
         * based on this. */
        min_slack = IGRAPH_INFINITY;
        min_slack_u = min_slack_v = 0;
        n = igraph_vector_size(&vec1);
        for (j = 0; j < larger_set_size; j++) {
            i = VECTOR(larger_set)[j];
            if (VECTOR(labels)[i] < min_slack) {
                min_slack = VECTOR(labels)[i];
                min_slack_v = (igraph_integer_t) i;
            }
        }
        min_slack_2 = IGRAPH_INFINITY;
        for (i = 0; i < n; i++) {
            u = (igraph_integer_t) VECTOR(vec1)[i];
            /* u is surely from the smaller set, but we are interested in it
             * only if it is reachable from an unmatched vertex */
            if (VECTOR(parent)[u] < 0) {
                continue;
            }
            if (VECTOR(labels)[u] < min_slack_2) {
                min_slack_2 = VECTOR(labels)[u];
                min_slack_u = u;
            }
        }
        min_slack += min_slack_2;
        debug("Starting approximation for min_slack = %.4f (based on vertex pair %ld--%ld)\n",
              min_slack, (long int)min_slack_u, (long int)min_slack_v);

        n = igraph_vector_size(&vec1);
        for (i = 0; i < n; i++) {
            u = (igraph_integer_t) VECTOR(vec1)[i];
            /* u is a reachable node in A; get its incident edges.
             *
             * There are two types of incident edges: 1) real edges,
             * 2) phantom edges. Phantom edges were treated earlier
             * when we determined the initial value for min_slack. */
            debug("Trying to expand along vertex %ld\n", (long int)u);
            neis = igraph_inclist_get(&inclist, u);
            k = igraph_vector_int_size(neis);
            for (j = 0; j < k; j++) {
                /* v is the vertex sitting at the other end of an edge incident
                 * on u; check whether it was reached */
                v = IGRAPH_OTHER(graph, VECTOR(*neis)[j], u);
                debug("  Edge %ld -- %ld (ID=%ld)\n", (long int)u, (long int)v, (long int)VECTOR(*neis)[j]);
                if (VECTOR(parent)[v] >= 0) {
                    /* v was reached, so we are not interested in it */
                    debug("    %ld was reached, so we are not interested in it\n", (long int)v);
                    continue;
                }
                /* v is the ID of the edge from now on */
                v = (igraph_integer_t) VECTOR(*neis)[j];
                if (VECTOR(slack)[v] < min_slack) {
                    min_slack = VECTOR(slack)[v];
                    min_slack_u = u;
                    min_slack_v = IGRAPH_OTHER(graph, v, u);
                }
                debug("    Slack of this edge: %.4f, min slack is now: %.4f\n",
                      VECTOR(slack)[v], min_slack);
            }
        }
        debug("Minimum slack: %.4f on edge %d--%d\n", min_slack, (int)min_slack_u, (int)min_slack_v);

        if (min_slack > 0) {
            /* Decrease the label of reachable nodes in A by min_slack.
             * Also update the dual solution */
            n = igraph_vector_size(&vec1);
            for (i = 0; i < n; i++) {
                u = (igraph_integer_t) VECTOR(vec1)[i];
                VECTOR(labels)[u] -= min_slack;
                neis = igraph_inclist_get(&inclist, u);
                k = igraph_vector_int_size(neis);
                for (j = 0; j < k; j++) {
                    debug("  Decreasing slack of edge %ld (%ld--%ld) by %.4f\n",
                          (long)VECTOR(*neis)[j], (long)u,
                          (long)IGRAPH_OTHER(graph, VECTOR(*neis)[j], u), min_slack);
                    VECTOR(slack)[(long int)VECTOR(*neis)[j]] -= min_slack;
                }
                dual -= min_slack;
            }

            /* Increase the label of reachable nodes in B by min_slack.
             * Also update the dual solution */
            n = igraph_vector_size(&vec2);
            for (i = 0; i < n; i++) {
                u = (igraph_integer_t) VECTOR(vec2)[i];
                VECTOR(labels)[u] += min_slack;
                neis = igraph_inclist_get(&inclist, u);
                k = igraph_vector_int_size(neis);
                for (j = 0; j < k; j++) {
                    debug("  Increasing slack of edge %ld (%ld--%ld) by %.4f\n",
                          (long)VECTOR(*neis)[j], (long)u,
                          (long)IGRAPH_OTHER(graph, (long)VECTOR(*neis)[j], u), min_slack);
                    VECTOR(slack)[(long int)VECTOR(*neis)[j]] += min_slack;
                }
                dual += min_slack;
            }
        }

        /* Update the set of tight phantom edges.
         * Note that we must do it even if min_slack is zero; the reason is that
         * it can happen that min_slack is zero in the first step if there are
         * isolated nodes in the input graph.
         *
         * TODO: this is O(n^2) here. Can we do it faster? */
        for (i = 0; i < smaller_set_size; i++) {
            u = VECTOR(smaller_set)[i];
            for (j = 0; j < larger_set_size; j++) {
                v = VECTOR(larger_set)[j];
                if (VECTOR(labels)[(long int)u] + VECTOR(labels)[(long int)v] <= eps) {
                    /* Tight phantom edge found. Note that we don't have to check whether
                     * u and v are connected; if they were, then the slack of this edge
                     * would be negative. */
                    neis2 = igraph_adjlist_get(&tight_phantom_edges, u);
                    if (!igraph_vector_int_binsearch(neis2, v, &k)) {
                        debug("New tight phantom edge: %ld -- %ld\n", (long)u, (long)v);
                        IGRAPH_CHECK(igraph_vector_int_insert(neis2, k, v));
                    }
                }
            }
        }

#ifdef MATCHING_DEBUG
        debug("New labels:");
        igraph_vector_print(&labels);
        debug("Slacks after updating with min_slack:");
        igraph_vector_print(&slack);
#endif
    }

    /* Cleanup: remove phantom edges from the matching */
    for (i = 0; i < smaller_set_size; i++) {
        u = VECTOR(smaller_set)[i];
        v = VECTOR(match)[u];
        if (v != -1) {
            neis2 = igraph_adjlist_get(&tight_phantom_edges, u);
            if (igraph_vector_int_binsearch(neis2, v, 0)) {
                VECTOR(match)[u] = VECTOR(match)[v] = -1;
                msize--;
            }
        }
    }

    /* Fill the output parameters */
    if (matching != 0) {
        IGRAPH_CHECK(igraph_vector_long_update(matching, &match));
    }
    if (matching_size != 0) {
        *matching_size = msize;
    }
    if (matching_weight != 0) {
        *matching_weight = 0;
        for (i = 0; i < no_of_edges; i++) {
            if (IS_TIGHT(i)) {
                IGRAPH_CHECK(igraph_edge(graph, (igraph_integer_t) i, &u, &v));
                if (VECTOR(match)[u] == v) {
                    *matching_weight += VECTOR(*weights)[i];
                }
            }
        }
    }

    /* Release everything */
#undef IS_TIGHT
    igraph_vector_destroy(&larger_set);
    igraph_vector_destroy(&smaller_set);
    igraph_inclist_destroy(&inclist);
    igraph_adjlist_destroy(&tight_phantom_edges);
    igraph_vector_destroy(&parent);
    igraph_dqueue_long_destroy(&q);
    igraph_vector_destroy(&labels);
    igraph_vector_destroy(&vec1);
    igraph_vector_destroy(&vec2);
    igraph_vector_destroy(&slack);
    igraph_vector_long_destroy(&match);
    IGRAPH_FINALLY_CLEAN(11);

    return IGRAPH_SUCCESS;
}

int igraph_maximum_matching(const igraph_t* graph, igraph_integer_t* matching_size,
                            igraph_real_t* matching_weight, igraph_vector_long_t* matching,
                            const igraph_vector_t* weights) {
    IGRAPH_UNUSED(graph);
    IGRAPH_UNUSED(matching_size);
    IGRAPH_UNUSED(matching_weight);
    IGRAPH_UNUSED(matching);
    IGRAPH_UNUSED(weights);
    IGRAPH_ERROR("maximum matching on general graphs not implemented yet",
                 IGRAPH_UNIMPLEMENTED);
}

#ifdef MATCHING_DEBUG
    #undef MATCHING_DEBUG
#endif
