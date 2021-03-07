/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "igraph_flow.h"

#include "igraph_adjlist.h"
#include "igraph_components.h"
#include "igraph_conversion.h"
#include "igraph_constants.h"
#include "igraph_constructors.h"
#include "igraph_dqueue.h"
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_progress.h"
#include "igraph_operators.h"
#include "igraph_structural.h"
#include "igraph_topology.h"

#include "core/buckets.h"
#include "core/cutheap.h"
#include "core/interruption.h"
#include "core/math.h"

#include "config.h"

/*
 * Some general remarks about the functions in this file.
 *
 * The following measures can be calculated:
 * ( 1) s-t maximum flow value, directed graph
 * ( 2) s-t maximum flow value, undirected graph
 * ( 3) s-t maximum flow, directed graph
 * ( 4) s-t maximum flow, undirected graph
 * ( 5) s-t minimum cut value, directed graph
 * ( 6) s-t minimum cut value, undirected graph
 * ( 7) minimum cut value, directed graph
 * ( 8) minimum cut value, undirected graph
 * ( 9) s-t minimum cut, directed graph
 * (10) s-t minimum cut, undirected graph
 * (11) minimum cut, directed graph
 * (12) minimum cut, undirected graph
 * (13) s-t edge connectivity, directed graph
 * (14) s-t edge connectivity, undirected graph
 * (15) edge connectivity, directed graph
 * (16) edge connectivity, undirected graph
 * (17) s-t vertex connectivity, directed graph
 * (18) s-t vertex connectivity, undirected graph
 * (19) vertex connectivity, directed graph
 * (20) vertex connectivity, undirected graph
 * (21) s-t number of edge disjoint paths, directed graph
 * (22) s-t number of edge disjoint paths, undirected graph
 * (23) s-t number of vertex disjoint paths, directed graph
 * (24) s-t number of vertex disjoint paths, undirected graph
 * (25) graph adhesion, directed graph
 * (26) graph adhesion, undirected graph
 * (27) graph cohesion, directed graph
 * (28) graph cohesion, undirected graph
 *
 * This is how they are calculated:
 * ( 1) igraph_maxflow_value, calls igraph_maxflow.
 * ( 2) igraph_maxflow_value, calls igraph_maxflow, this calls
 *      igraph_i_maxflow_undirected. This transforms the graph into a
 *      directed graph, including two mutual edges instead of every
 *      undirected edge, then igraph_maxflow is called again with the
 *      directed graph.
 * ( 3) igraph_maxflow, does the push-relabel algorithm, optionally
 *      calculates the cut, the partitions and the flow itself.
 * ( 4) igraph_maxflow calls igraph_i_maxflow_undirected, this converts
 *      the undirected graph into a directed one, adding two mutual edges
 *      for each undirected edge, then igraph_maxflow is called again,
 *      with the directed graph. After igraph_maxflow returns, we need
 *      to edit the flow (and the cut) to make it sense for the
 *      original graph.
 * ( 5) igraph_st_mincut_value, we just call igraph_maxflow_value
 * ( 6) igraph_st_mincut_value, we just call igraph_maxflow_value
 * ( 7) igraph_mincut_value, we call igraph_maxflow_value (|V|-1)*2
 *      times, from vertex 0 to all other vertices and from all other
 *      vertices to vertex 0
 * ( 8) We call igraph_i_mincut_value_undirected, that calls
 *      igraph_i_mincut_undirected with partition=partition2=cut=NULL
 *      The Stoer-Wagner algorithm is used.
 * ( 9) igraph_st_mincut, just calls igraph_maxflow.
 * (10) igraph_st_mincut, just calls igraph_maxflow.
 * (11) igraph_mincut, calls igraph_i_mincut_directed, which runs
 *      the maximum flow algorithm 2(|V|-1) times, from vertex zero to
 *      and from all other vertices and stores the smallest cut.
 * (12) igraph_mincut, igraph_i_mincut_undirected is called,
 *      this is the Stoer-Wagner algorithm
 * (13) We just call igraph_maxflow_value, back to (1)
 * (14) We just call igraph_maxflow_value, back to (2)
 * (15) We just call igraph_mincut_value (possibly after some basic
 *      checks). Back to (7)
 * (16) We just call igraph_mincut_value (possibly after some basic
 *      checks). Back to (8).
 * (17) We call igraph_i_st_vertex_connectivity_directed.
 *      That creates a new graph with 2*|V| vertices and smartly chosen
 *      edges, so that the s-t edge connectivity of this graph is the
 *      same as the s-t vertex connectivity of the original graph.
 *      So finally it calls igraph_maxflow_value, go to (1)
 * (18) We call igraph_i_st_vertex_connectivity_undirected.
 *      We convert the graph to a directed one,
 *      IGRAPH_TO_DIRECTED_MUTUAL method. Then we call
 *      igraph_i_st_vertex_connectivity_directed, see (17).
 * (19) We call igraph_i_vertex_connectivity_directed.
 *      That calls igraph_st_vertex_connectivity for all pairs of
 *      vertices. Back to (17).
 * (20) We call igraph_i_vertex_connectivity_undirected.
 *      That converts the graph into a directed one
 *      (IGRAPH_TO_DIRECTED_MUTUAL) and calls the directed version,
 *      igraph_i_vertex_connectivity_directed, see (19).
 * (21) igraph_edge_disjoint_paths, we just call igraph_maxflow_value, (1).
 * (22) igraph_edge_disjoint_paths, we just call igraph_maxflow_value, (2).
 * (23) igraph_vertex_disjoint_paths, if there is a connection between
 *      the two vertices, then we remove that (or all of them if there
 *      are many), as this could mess up vertex connectivity
 *      calculation. The we call
 *      igraph_i_st_vertex_connectivity_directed, see (19).
 * (24) igraph_vertex_disjoint_paths, if there is a connection between
 *      the two vertices, then we remove that (or all of them if there
 *      are many), as this could mess up vertex connectivity
 *      calculation. The we call
 *      igraph_i_st_vertex_connectivity_undirected, see (20).
 * (25) We just call igraph_edge_connectivity, see (15).
 * (26) We just call igraph_edge_connectivity, see (16).
 * (27) We just call igraph_vertex_connectivity, see (19).
 * (28) We just call igraph_vertex_connectivity, see (20).
 */

/*
 * This is an internal function that calculates the maximum flow value
 * on undirected graphs, either for an s-t vertex pair or for the
 * graph (i.e. all vertex pairs).
 *
 * It does it by converting the undirected graph to a corresponding
 * directed graph, including reciprocal directed edges instead of each
 * undirected edge.
 */

static int igraph_i_maxflow_undirected(const igraph_t *graph,
                                       igraph_real_t *value,
                                       igraph_vector_t *flow,
                                       igraph_vector_t *cut,
                                       igraph_vector_t *partition,
                                       igraph_vector_t *partition2,
                                       igraph_integer_t source,
                                       igraph_integer_t target,
                                       const igraph_vector_t *capacity,
                                       igraph_maxflow_stats_t *stats) {
    igraph_integer_t no_of_edges = (igraph_integer_t) igraph_ecount(graph);
    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph);
    igraph_vector_t edges;
    igraph_vector_t newcapacity;
    igraph_t newgraph;
    long int i;

    /* We need to convert this to directed by hand, since we need to be
       sure that the edge ids will be handled properly to build the new
       capacity vector. */

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&newcapacity, no_of_edges * 2);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 4));
    IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));
    IGRAPH_CHECK(igraph_vector_resize(&edges, no_of_edges * 4));
    for (i = 0; i < no_of_edges; i++) {
        VECTOR(edges)[no_of_edges * 2 + i * 2] = VECTOR(edges)[i * 2 + 1];
        VECTOR(edges)[no_of_edges * 2 + i * 2 + 1] = VECTOR(edges)[i * 2];
        VECTOR(newcapacity)[i] = VECTOR(newcapacity)[no_of_edges + i] =
                                     capacity ? VECTOR(*capacity)[i] : 1.0;
    }

    IGRAPH_CHECK(igraph_create(&newgraph, &edges, no_of_nodes, IGRAPH_DIRECTED));
    IGRAPH_FINALLY(igraph_destroy, &newgraph);

    IGRAPH_CHECK(igraph_maxflow(&newgraph, value, flow, cut, partition,
                                partition2, source, target, &newcapacity, stats));

    if (cut) {
        long int i, cs = igraph_vector_size(cut);
        for (i = 0; i < cs; i++) {
            if (VECTOR(*cut)[i] >= no_of_edges) {
                VECTOR(*cut)[i] -= no_of_edges;
            }
        }
    }

    /* The flow has one non-zero value for each real-nonreal edge pair,
       by definition, we convert it to a positive-negative vector. If
       for an edge the flow is negative that means that it is going
       from the bigger vertex id to the smaller one. For positive
       values the direction is the opposite. */
    if (flow) {
        long int i;
        for (i = 0; i < no_of_edges; i++) {
            VECTOR(*flow)[i] -= VECTOR(*flow)[i + no_of_edges];
        }
        IGRAPH_CHECK(igraph_vector_resize(flow, no_of_edges));
    }

    igraph_destroy(&newgraph);
    igraph_vector_destroy(&edges);
    igraph_vector_destroy(&newcapacity);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

#define FIRST(i)       (VECTOR(*first)[(i)])
#define LAST(i)        (VECTOR(*first)[(i)+1])
#define CURRENT(i)     (VECTOR(*current)[(i)])
#define RESCAP(i)      (VECTOR(*rescap)[(i)])
#define REV(i)         (VECTOR(*rev)[(i)])
#define HEAD(i)        (VECTOR(*to)[(i)])
#define EXCESS(i)      (VECTOR(*excess)[(i)])
#define DIST(i)        (VECTOR(*distance)[(i)])
#define DISCHARGE(v)   (igraph_i_mf_discharge((v), &current, &first, &rescap, \
                        &to, &distance, &excess,        \
                        no_of_nodes, source, target,    \
                        &buckets, &ibuckets,        \
                        &rev, stats, &npushsince,       \
                        &nrelabelsince))
#define PUSH(v,e,n)    (igraph_i_mf_push((v), (e), (n), current, rescap,      \
                        excess, target, source, buckets,     \
                        ibuckets, distance, rev, stats,      \
                        npushsince))
#define RELABEL(v)     (igraph_i_mf_relabel((v), no_of_nodes, distance,       \
                        first, rescap, to, current,       \
                        stats, nrelabelsince))
#define GAP(b)         (igraph_i_mf_gap((b), stats, buckets, ibuckets,        \
                                        no_of_nodes, distance))
#define BFS()          (igraph_i_mf_bfs(&bfsq, source, target, no_of_nodes,   \
                                        &buckets, &ibuckets, &distance,       \
                                        &first, &current, &to, &excess,       \
                                        &rescap, &rev))

static void igraph_i_mf_gap(long int b, igraph_maxflow_stats_t *stats,
                            igraph_buckets_t *buckets, igraph_dbuckets_t *ibuckets,
                            long int no_of_nodes,
                            igraph_vector_long_t *distance) {

    IGRAPH_UNUSED(buckets);

    long int bo;
    (stats->nogap)++;
    for (bo = b + 1; bo <= no_of_nodes; bo++) {
        while (!igraph_dbuckets_empty_bucket(ibuckets, bo)) {
            long int n = igraph_dbuckets_pop(ibuckets, bo);
            (stats->nogapnodes)++;
            DIST(n) = no_of_nodes;
        }
    }
}

static void igraph_i_mf_relabel(long int v, long int no_of_nodes,
                                igraph_vector_long_t *distance,
                                igraph_vector_long_t *first,
                                igraph_vector_t *rescap, igraph_vector_long_t *to,
                                igraph_vector_long_t *current,
                                igraph_maxflow_stats_t *stats, int *nrelabelsince) {

    long int min = no_of_nodes;
    long int k, l, min_edge = 0;
    (stats->norelabel)++; (*nrelabelsince)++;
    DIST(v) = no_of_nodes;
    for (k = FIRST(v), l = LAST(v); k < l; k++) {
        if (RESCAP(k) > 0 && DIST(HEAD(k)) < min) {
            min = DIST(HEAD(k));
            min_edge = k;
        }
    }
    min++;
    if (min < no_of_nodes) {
        DIST(v) = min;
        CURRENT(v) = min_edge;
    }
}

static void igraph_i_mf_push(long int v, long int e, long int n,
                             igraph_vector_long_t *current,
                             igraph_vector_t *rescap, igraph_vector_t *excess,
                             long int target, long int source,
                             igraph_buckets_t *buckets, igraph_dbuckets_t *ibuckets,
                             igraph_vector_long_t *distance,
                             igraph_vector_long_t *rev, igraph_maxflow_stats_t *stats,
                             int *npushsince) {


    IGRAPH_UNUSED(current);
    IGRAPH_UNUSED(source);

    igraph_real_t delta =
        RESCAP(e) < EXCESS(v) ? RESCAP(e) : EXCESS(v);
    (stats->nopush)++; (*npushsince)++;
    if (EXCESS(n) == 0 && n != target) {
        igraph_dbuckets_delete(ibuckets, DIST(n), n);
        igraph_buckets_add(buckets, (long int) DIST(n), n);
    }
    RESCAP(e) -= delta;
    RESCAP(REV(e)) += delta;
    EXCESS(n) += delta;
    EXCESS(v) -= delta;
}

static void igraph_i_mf_discharge(long int v,
                                  igraph_vector_long_t *current,
                                  igraph_vector_long_t *first,
                                  igraph_vector_t *rescap,
                                  igraph_vector_long_t *to,
                                  igraph_vector_long_t *distance,
                                  igraph_vector_t *excess,
                                  long int no_of_nodes, long int source,
                                  long int target, igraph_buckets_t *buckets,
                                  igraph_dbuckets_t *ibuckets,
                                  igraph_vector_long_t *rev,
                                  igraph_maxflow_stats_t *stats,
                                  int *npushsince, int *nrelabelsince) {
    do {
        long int i;
        long int start = (long int) CURRENT(v);
        long int stop = (long int) LAST(v);
        for (i = start; i < stop; i++) {
            if (RESCAP(i) > 0) {
                long int nei = HEAD(i);
                if (DIST(v) == DIST(nei) + 1) {
                    PUSH((v), i, nei);
                    if (EXCESS(v) == 0) {
                        break;
                    }
                }
            }
        }
        if (i == stop) {
            long int origdist = DIST(v);
            RELABEL(v);
            if (igraph_buckets_empty_bucket(buckets, origdist) &&
                igraph_dbuckets_empty_bucket(ibuckets, origdist)) {
                GAP(origdist);
            }
            if (DIST(v) == no_of_nodes) {
                break;
            }
        } else {
            CURRENT(v) = i;
            igraph_dbuckets_add(ibuckets, DIST(v), v);
            break;
        }
    } while (1);
}

static void igraph_i_mf_bfs(igraph_dqueue_long_t *bfsq,
                            long int source, long int target,
                            long int no_of_nodes, igraph_buckets_t *buckets,
                            igraph_dbuckets_t *ibuckets,
                            igraph_vector_long_t *distance,
                            igraph_vector_long_t *first, igraph_vector_long_t *current,
                            igraph_vector_long_t *to, igraph_vector_t *excess,
                            igraph_vector_t *rescap, igraph_vector_long_t *rev) {

    long int k, l;

    IGRAPH_UNUSED(source);

    igraph_buckets_clear(buckets);
    igraph_dbuckets_clear(ibuckets);
    igraph_vector_long_fill(distance, no_of_nodes);
    DIST(target) = 0;

    igraph_dqueue_long_push(bfsq, target);
    while (!igraph_dqueue_long_empty(bfsq)) {
        long int node = igraph_dqueue_long_pop(bfsq);
        long int ndist = DIST(node) + 1;
        for (k = FIRST(node), l = LAST(node); k < l; k++) {
            if (RESCAP(REV(k)) > 0) {
                long int nei = HEAD(k);
                if (DIST(nei) == no_of_nodes) {
                    DIST(nei) = ndist;
                    CURRENT(nei) = FIRST(nei);
                    if (EXCESS(nei) > 0) {
                        igraph_buckets_add(buckets, ndist, nei);
                    } else {
                        igraph_dbuckets_add(ibuckets, ndist, nei);
                    }
                    igraph_dqueue_long_push(bfsq, nei);
                }
            }
        }
    }
}

/**
 * \function igraph_maxflow
 * Maximum network flow between a pair of vertices
 *
 * </para><para>This function implements the Goldberg-Tarjan algorithm for
 * calculating value of the maximum flow in a directed or undirected
 * graph. The algorithm was given in Andrew V. Goldberg, Robert
 * E. Tarjan: A New Approach to the Maximum-Flow Problem, Journal of
 * the ACM, 35(4), 921-940, 1988. </para>
 *
 * <para> The input of the function is a graph, a vector
 * of real numbers giving the capacity of the edges and two vertices
 * of the graph, the source and the target. A flow is a function
 * assigning positive real numbers to the edges and satisfying two
 * requirements: (1) the flow value is less than the capacity of the
 * edge and (2) at each vertex except the source and the target, the
 * incoming flow (i.e. the sum of the flow on the incoming edges) is
 * the same as the outgoing flow (i.e. the sum of the flow on the
 * outgoing edges). The value of the flow is the incoming flow at the
 * target vertex. The maximum flow is the flow with the maximum
 * value.
 *
 * \param graph The input graph, either directed or undirected.
 * \param value Pointer to a real number, the value of the maximum
 *        will be placed here, unless it is a null pointer.
 * \param flow If not a null pointer, then it must be a pointer to an
 *        initialized vector. The vector will be resized, and the flow
 *        on each edge will be placed in it, in the order of the edge
 *        ids. For undirected graphs this argument is bit trickier,
 *        since for these the flow direction is not predetermined by
 *        the edge direction. For these graphs the elements of the
 *        \p flow vector can be negative, this means that the flow
 *        goes from the bigger vertex id to the smaller one. Positive
 *        values mean that the flow goes from the smaller vertex id to
 *        the bigger one.
 * \param cut A null pointer or a pointer to an initialized vector.
 *        If not a null pointer, then the minimum cut corresponding to
 *        the maximum flow is stored here, i.e. all edge ids that are
 *        part of the minimum cut are stored in the vector.
 * \param partition A null pointer or a pointer to an initialized
 *        vector. If not a null pointer, then the first partition of
 *        the minimum cut that corresponds to the maximum flow will be
 *        placed here. The first partition is always the one that
 *        contains the source vertex.
 * \param partition2 A null pointer or a pointer to an initialized
 *        vector. If not a null pointer, then the second partition of
 *        the minimum cut that corresponds to the maximum flow will be
 *        placed here. The second partition is always the one that
 *        contains the target vertex.
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \param capacity Vector containing the capacity of the edges. If NULL, then
 *        every edge is considered to have capacity 1.0.
 * \param stats Counts of the number of different operations
 *        preformed by the algorithm are stored here.
 * \return Error code.
 *
 * Time complexity: O(|V|^3). In practice it is much faster, but i
 * cannot prove a better lower bound for the data structure i've
 * used. In fact, this implementation runs much faster than the
 * \c hi_pr implementation discussed in
 * B. V. Cherkassky and A. V. Goldberg: On implementing the
 * push-relabel method for the maximum flow problem, (Algorithmica,
 * 19:390--410, 1997) on all the graph classes i've tried.
 *
 * \sa \ref igraph_mincut_value(), \ref igraph_edge_connectivity(),
 * \ref igraph_vertex_connectivity() for
 * properties based on the maximum flow.
 *
 * \example examples/simple/flow.c
 * \example examples/simple/flow2.c
 */

int igraph_maxflow(const igraph_t *graph, igraph_real_t *value,
                   igraph_vector_t *flow, igraph_vector_t *cut,
                   igraph_vector_t *partition, igraph_vector_t *partition2,
                   igraph_integer_t source, igraph_integer_t target,
                   const igraph_vector_t *capacity,
                   igraph_maxflow_stats_t *stats) {

    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph);
    igraph_integer_t no_of_orig_edges = (igraph_integer_t) igraph_ecount(graph);
    igraph_integer_t no_of_edges = 2 * no_of_orig_edges;

    igraph_vector_t rescap, excess;
    igraph_vector_long_t from, to, rev, distance;
    igraph_vector_t edges, rank;
    igraph_vector_long_t current, first;
    igraph_buckets_t buckets;
    igraph_dbuckets_t ibuckets;

    igraph_dqueue_long_t bfsq;

    long int i, j, idx;
    int npushsince = 0, nrelabelsince = 0;

    igraph_maxflow_stats_t local_stats;   /* used if the user passed a null pointer for stats */

    if (stats == 0) {
        stats = &local_stats;
    }

    if (!igraph_is_directed(graph)) {
        IGRAPH_CHECK(igraph_i_maxflow_undirected(graph, value, flow, cut,
                     partition, partition2, source,
                     target, capacity, stats));
        return 0;
    }

    if (capacity && igraph_vector_size(capacity) != no_of_orig_edges) {
        IGRAPH_ERROR("Invalid capacity vector", IGRAPH_EINVAL);
    }
    if (source < 0 || source >= no_of_nodes || target < 0 || target >= no_of_nodes) {
        IGRAPH_ERROR("Invalid source or target vertex", IGRAPH_EINVAL);
    }

    stats->nopush = stats->norelabel = stats->nogap = stats->nogapnodes =
                                           stats->nobfs = 0;

    /*
     * The data structure:
     * - First of all, we consider every edge twice, first the edge
     *   itself, but also its opposite.
     * - (from, to) contain all edges (original + opposite), ordered by
     *   the id of the source vertex. During the algorithm we just need
     *   'to', so from is destroyed soon. We only need it in the
     *   beginning, to create the 'first' pointers.
     * - 'first' is a pointer vector for 'to', first[i] points to the
     *   first neighbor of vertex i and first[i+1]-1 is the last
     *   neighbor of vertex i. (Unless vertex i is isolate, in which
     *   case first[i]==first[i+1]).
     * - 'rev' contains a mapping from an edge to its opposite pair
     * - 'rescap' contains the residual capacities of the edges, this is
     *   initially equal to the capacity of the edges for the original
     *   edges and it is zero for the opposite edges.
     * - 'excess' contains the excess flow for the vertices. I.e. the flow
     *   that is coming in, but it is not going out.
     * - 'current' stores the next neighboring vertex to check, for every
     *   vertex, when excess flow is being pushed to neighbors.
     * - 'distance' stores the distance of the vertices from the source.
     * - 'rank' and 'edges' are only needed temporarily, for ordering and
     *   storing the edges.
     * - we use an igraph_buckets_t data structure ('buckets') to find
     *   the vertices with the highest 'distance' values quickly.
     *   This always contains the vertices that have a positive excess
     *   flow.
     */
#undef FIRST
#undef LAST
#undef CURRENT
#undef RESCAP
#undef REV
#undef HEAD
#undef EXCESS
#undef DIST
#define FIRST(i)       (VECTOR(first)[(i)])
#define LAST(i)        (VECTOR(first)[(i)+1])
#define CURRENT(i)     (VECTOR(current)[(i)])
#define RESCAP(i)      (VECTOR(rescap)[(i)])
#define REV(i)         (VECTOR(rev)[(i)])
#define HEAD(i)        (VECTOR(to)[(i)])
#define EXCESS(i)      (VECTOR(excess)[(i)])
#define DIST(i)        (VECTOR(distance)[(i)])

    igraph_dqueue_long_init(&bfsq,             no_of_nodes);
    IGRAPH_FINALLY(igraph_dqueue_long_destroy, &bfsq);
    IGRAPH_VECTOR_LONG_INIT_FINALLY(&to,       no_of_edges);
    IGRAPH_VECTOR_LONG_INIT_FINALLY(&rev,      no_of_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&rescap,        no_of_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&excess,        no_of_nodes);
    IGRAPH_VECTOR_LONG_INIT_FINALLY(&distance, no_of_nodes);
    IGRAPH_VECTOR_LONG_INIT_FINALLY(&first,    no_of_nodes + 1);

    IGRAPH_VECTOR_INIT_FINALLY(&rank,          no_of_edges);
    IGRAPH_VECTOR_LONG_INIT_FINALLY(&from,     no_of_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&edges,         no_of_edges);

    /* Create the basic data structure */
    IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));
    IGRAPH_CHECK(igraph_vector_rank(&edges, &rank, no_of_nodes));

    for (i = 0; i < no_of_edges; i += 2) {
        long int pos = (long int) VECTOR(rank)[i];
        long int pos2 = (long int) VECTOR(rank)[i + 1];
        VECTOR(from)[pos] = VECTOR(edges)[i];
        VECTOR(to)[pos]   = VECTOR(edges)[i + 1];
        VECTOR(from)[pos2] = VECTOR(edges)[i + 1];
        VECTOR(to)[pos2]   = VECTOR(edges)[i];
        VECTOR(rev)[pos] = pos2;
        VECTOR(rev)[pos2] = pos;
        VECTOR(rescap)[pos] = capacity ? VECTOR(*capacity)[i / 2] : 1.0;
        VECTOR(rescap)[pos2] = 0.0;
    }

    /* The first pointers. This is a but trickier, than one would
       think, because of the possible isolate vertices. */

    idx = -1;
    for (i = 0; i <= VECTOR(from)[0]; i++) {
        idx++; VECTOR(first)[idx] = 0;
    }
    for (i = 1; i < no_of_edges; i++) {
        long int n = (long int) (VECTOR(from)[i] -
                                 VECTOR(from)[ (long int) VECTOR(first)[idx] ]);
        for (j = 0; j < n; j++) {
            idx++; VECTOR(first)[idx] = i;
        }
    }
    idx++;
    while (idx < no_of_nodes + 1) {
        VECTOR(first)[idx++] = no_of_edges;
    }

    igraph_vector_long_destroy(&from);
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(2);

    if (!flow) {
        igraph_vector_destroy(&rank);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* And the current pointers, initially the same as the first */
    IGRAPH_VECTOR_LONG_INIT_FINALLY(&current, no_of_nodes);
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(current)[i] = VECTOR(first)[i];
    }

    /* OK, the graph is set up, initialization */

    IGRAPH_CHECK(igraph_buckets_init(&buckets, no_of_nodes + 1, no_of_nodes));
    IGRAPH_FINALLY(igraph_buckets_destroy, &buckets);
    IGRAPH_CHECK(igraph_dbuckets_init(&ibuckets, no_of_nodes + 1, no_of_nodes));
    IGRAPH_FINALLY(igraph_dbuckets_destroy, &ibuckets);

    /* Send as much flow as possible from the source to its neighbors */
    for (i = FIRST(source), j = LAST(source); i < j; i++) {
        if (HEAD(i) != source) {
            igraph_real_t delta = RESCAP(i);
            RESCAP(i) = 0;
            RESCAP(REV(i)) += delta;
            EXCESS(HEAD(i)) += delta;
        }
    }

    BFS();
    (stats->nobfs)++;

    while (!igraph_buckets_empty(&buckets)) {
        long int vertex = igraph_buckets_popmax(&buckets);
        DISCHARGE(vertex);
        if (npushsince > no_of_nodes / 2 && nrelabelsince > no_of_nodes) {
            (stats->nobfs)++;
            BFS();
            npushsince = nrelabelsince = 0;
        }
    }

    /* Store the result */
    if (value) {
        *value = EXCESS(target);
    }

    /* If we also need the minimum cut */
    if (cut || partition || partition2) {
        /* We need to find all vertices from which the target is reachable
           in the residual graph. We do a breadth-first search, going
           backwards. */
        igraph_dqueue_t Q;
        igraph_vector_bool_t added;
        long int marked = 0;

        IGRAPH_CHECK(igraph_vector_bool_init(&added, no_of_nodes));
        IGRAPH_FINALLY(igraph_vector_bool_destroy, &added);

        IGRAPH_CHECK(igraph_dqueue_init(&Q, 100));
        IGRAPH_FINALLY(igraph_dqueue_destroy, &Q);

        igraph_dqueue_push(&Q, target);
        VECTOR(added)[(long int)target] = 1;
        marked++;
        while (!igraph_dqueue_empty(&Q)) {
            long int actnode = (long int) igraph_dqueue_pop(&Q);
            for (i = FIRST(actnode), j = LAST(actnode); i < j; i++) {
                long int nei = HEAD(i);
                if (!VECTOR(added)[nei] && RESCAP(REV(i)) > 0.0) {
                    VECTOR(added)[nei] = 1;
                    marked++;
                    IGRAPH_CHECK(igraph_dqueue_push(&Q, nei));
                }
            }
        }
        igraph_dqueue_destroy(&Q);
        IGRAPH_FINALLY_CLEAN(1);

        /* Now we marked each vertex that is on one side of the cut,
           check the crossing edges */

        if (cut) {
            igraph_vector_clear(cut);
            for (i = 0; i < no_of_orig_edges; i++) {
                long int f = IGRAPH_FROM(graph, i);
                long int t = IGRAPH_TO(graph, i);
                if (!VECTOR(added)[f] && VECTOR(added)[t]) {
                    IGRAPH_CHECK(igraph_vector_push_back(cut, i));
                }
            }
        }

        if (partition2) {
            long int x = 0;
            IGRAPH_CHECK(igraph_vector_resize(partition2, marked));
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(added)[i]) {
                    VECTOR(*partition2)[x++] = i;
                }
            }
        }

        if (partition) {
            long int x = 0;
            IGRAPH_CHECK(igraph_vector_resize(partition,
                                              no_of_nodes - marked));
            for (i = 0; i < no_of_nodes; i++) {
                if (!VECTOR(added)[i]) {
                    VECTOR(*partition)[x++] = i;
                }
            }
        }

        igraph_vector_bool_destroy(&added);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (flow) {
        /* Initialize the backward distances, with a breadth-first search
           from the source */
        igraph_dqueue_t Q;
        igraph_vector_int_t added;
        long int j, k, l;
        igraph_t flow_graph;
        igraph_vector_t flow_edges;
        igraph_bool_t dag;

        IGRAPH_CHECK(igraph_vector_int_init(&added, no_of_nodes));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &added);
        IGRAPH_CHECK(igraph_dqueue_init(&Q, 100));
        IGRAPH_FINALLY(igraph_dqueue_destroy, &Q);

        igraph_dqueue_push(&Q, source);
        igraph_dqueue_push(&Q, 0);
        VECTOR(added)[(long int)source] = 1;
        while (!igraph_dqueue_empty(&Q)) {
            long int actnode = (long int) igraph_dqueue_pop(&Q);
            long int actdist = (long int) igraph_dqueue_pop(&Q);
            DIST(actnode) = actdist;

            for (i = FIRST(actnode), j = LAST(actnode); i < j; i++) {
                long int nei = HEAD(i);
                if (!VECTOR(added)[nei] && RESCAP(REV(i)) > 0.0) {
                    VECTOR(added)[nei] = 1;
                    IGRAPH_CHECK(igraph_dqueue_push(&Q, nei));
                    IGRAPH_CHECK(igraph_dqueue_push(&Q, actdist + 1));
                }
            }
        } /* !igraph_dqueue_empty(&Q) */

        igraph_vector_int_destroy(&added);
        igraph_dqueue_destroy(&Q);
        IGRAPH_FINALLY_CLEAN(2);

        /* Reinitialize the buckets */
        igraph_buckets_clear(&buckets);
        for (i = 0; i < no_of_nodes; i++) {
            if (EXCESS(i) > 0.0 && i != source && i != target) {
                igraph_buckets_add(&buckets, (long int) DIST(i), i);
            }
        }

        /* Now we return the flow to the source */
        while (!igraph_buckets_empty(&buckets)) {
            long int vertex = igraph_buckets_popmax(&buckets);

            /* DISCHARGE(vertex) comes here */
            do {
                for (i = (long int) CURRENT(vertex), j = LAST(vertex); i < j; i++) {
                    if (RESCAP(i) > 0) {
                        long int nei = HEAD(i);

                        if (DIST(vertex) == DIST(nei) + 1) {
                            igraph_real_t delta =
                                RESCAP(i) < EXCESS(vertex) ? RESCAP(i) : EXCESS(vertex);
                            RESCAP(i) -= delta;
                            RESCAP(REV(i)) += delta;

                            if (nei != source && EXCESS(nei) == 0.0 &&
                                DIST(nei) != no_of_nodes) {
                                igraph_buckets_add(&buckets, (long int) DIST(nei), nei);
                            }

                            EXCESS(nei) += delta;
                            EXCESS(vertex) -= delta;

                            if (EXCESS(vertex) == 0) {
                                break;
                            }

                        }
                    }
                }

                if (i == j) {

                    /* RELABEL(vertex) comes here */
                    igraph_real_t min;
                    long int min_edge = 0;
                    DIST(vertex) = min = no_of_nodes;
                    for (k = FIRST(vertex), l = LAST(vertex); k < l; k++) {
                        if (RESCAP(k) > 0) {
                            if (DIST(HEAD(k)) < min) {
                                min = DIST(HEAD(k));
                                min_edge = k;
                            }
                        }
                    }

                    min++;

                    if (min < no_of_nodes) {
                        DIST(vertex) = min;
                        CURRENT(vertex) = min_edge;
                        /* Vertex is still active */
                        igraph_buckets_add(&buckets, (long int) DIST(vertex), vertex);
                    }

                    /* TODO: gap heuristics here ??? */

                } else {
                    CURRENT(vertex) = FIRST(vertex);
                }

                break;

            } while (1);
        }

        /* We need to eliminate flow cycles now. Before that we check that
           there is a cycle in the flow graph.

           First we do a couple of DFSes from the source vertex to the
           target and factor out the paths we find. If there is no more
           path to the target, then all remaining flow must be in flow
           cycles, so we don't need it at all.

           Some details. 'stack' contains the whole path of the DFS, both
           the vertices and the edges, they are alternating in the stack.
           'current' helps finding the next outgoing edge of a vertex
           quickly, the next edge of 'v' is FIRST(v)+CURRENT(v). If this
           is LAST(v), then there are no more edges to try.

           The 'added' vector contains 0 if the vertex was not visited
           before, 1 if it is currently in 'stack', and 2 if it is not in
           'stack', but it was visited before. */

        IGRAPH_VECTOR_INIT_FINALLY(&flow_edges, 0);
        for (i = 0, j = 0; i < no_of_edges; i += 2, j++) {
            long int pos = (long int) VECTOR(rank)[i];
            if ((capacity ? VECTOR(*capacity)[j] : 1.0) > RESCAP(pos)) {
                IGRAPH_CHECK(igraph_vector_push_back(&flow_edges,
                                                     IGRAPH_FROM(graph, j)));
                IGRAPH_CHECK(igraph_vector_push_back(&flow_edges,
                                                     IGRAPH_TO(graph, j)));
            }
        }
        IGRAPH_CHECK(igraph_create(&flow_graph, &flow_edges, no_of_nodes,
                                   IGRAPH_DIRECTED));
        igraph_vector_destroy(&flow_edges);
        IGRAPH_FINALLY_CLEAN(1);
        IGRAPH_FINALLY(igraph_destroy, &flow_graph);
        IGRAPH_CHECK(igraph_is_dag(&flow_graph, &dag));
        igraph_destroy(&flow_graph);
        IGRAPH_FINALLY_CLEAN(1);

        if (!dag) {
            igraph_vector_long_t stack;
            igraph_vector_t mycap;

            IGRAPH_CHECK(igraph_vector_long_init(&stack, 0));
            IGRAPH_FINALLY(igraph_vector_long_destroy, &stack);
            IGRAPH_CHECK(igraph_vector_int_init(&added, no_of_nodes));
            IGRAPH_FINALLY(igraph_vector_int_destroy, &added);
            IGRAPH_VECTOR_INIT_FINALLY(&mycap, no_of_edges);

#define MYCAP(i)      (VECTOR(mycap)[(i)])

            for (i = 0; i < no_of_edges; i += 2) {
                long int pos = (long int) VECTOR(rank)[i];
                long int pos2 = (long int) VECTOR(rank)[i + 1];
                MYCAP(pos) = (capacity ? VECTOR(*capacity)[i / 2] : 1.0) - RESCAP(pos);
                MYCAP(pos2) = 0.0;
            }

            do {
                igraph_vector_long_null(&current);
                igraph_vector_long_clear(&stack);
                igraph_vector_int_null(&added);

                IGRAPH_CHECK(igraph_vector_long_push_back(&stack, -1));
                IGRAPH_CHECK(igraph_vector_long_push_back(&stack, source));
                VECTOR(added)[(long int)source] = 1;
                while (!igraph_vector_long_empty(&stack) &&
                       igraph_vector_long_tail(&stack) != target) {
                    long int actnode = igraph_vector_long_tail(&stack);
                    long int edge = FIRST(actnode) + (long int) CURRENT(actnode);
                    long int nei;
                    while (edge < LAST(actnode) && MYCAP(edge) == 0.0) {
                        edge++;
                    }
                    nei = edge < LAST(actnode) ? HEAD(edge) : -1;

                    if (edge < LAST(actnode) && !VECTOR(added)[nei]) {
                        /* Go forward along next edge, if the vertex was not
                           visited before */
                        IGRAPH_CHECK(igraph_vector_long_push_back(&stack, edge));
                        IGRAPH_CHECK(igraph_vector_long_push_back(&stack, nei));
                        VECTOR(added)[nei] = 1;
                        CURRENT(actnode) += 1;
                    } else if (edge < LAST(actnode) && VECTOR(added)[nei] == 1) {
                        /* We found a flow cycle, factor it out. Go back in stack
                           until we find 'nei' again, determine the flow along the
                           cycle. */
                        igraph_real_t thisflow = MYCAP(edge);
                        long int idx;
                        for (idx = igraph_vector_long_size(&stack) - 2;
                             idx >= 0 && VECTOR(stack)[idx + 1] != nei; idx -= 2) {
                            long int e = VECTOR(stack)[idx];
                            igraph_real_t rcap = e >= 0 ? MYCAP(e) : MYCAP(edge);
                            if (rcap < thisflow) {
                                thisflow = rcap;
                            }
                        }
                        MYCAP(edge) -= thisflow; RESCAP(edge) += thisflow;
                        for (idx = igraph_vector_long_size(&stack) - 2;
                             idx >= 0 && VECTOR(stack)[idx + 1] != nei; idx -= 2) {
                            long int e = VECTOR(stack)[idx];
                            if (e >= 0) {
                                MYCAP(e) -= thisflow;
                                RESCAP(e) += thisflow;
                            }
                        }
                        CURRENT(actnode) += 1;
                    } else if (edge < LAST(actnode)) { /* && VECTOR(added)[nei]==2 */
                        /* The next edge leads to a vertex that was visited before,
                           but it is currently not in 'stack' */
                        CURRENT(actnode) += 1;
                    } else {
                        /* Go backward, take out the node and the edge that leads to it */
                        igraph_vector_long_pop_back(&stack);
                        igraph_vector_long_pop_back(&stack);
                        VECTOR(added)[actnode] = 2;
                    }
                }

                /* If non-empty, then it contains a path from source to target
                   in the residual graph. We factor out this path from the flow. */
                if (!igraph_vector_long_empty(&stack)) {
                    long int pl = igraph_vector_long_size(&stack);
                    igraph_real_t thisflow = EXCESS(target);
                    for (i = 2; i < pl; i += 2) {
                        long int edge = VECTOR(stack)[i];
                        igraph_real_t rcap = MYCAP(edge);
                        if (rcap < thisflow) {
                            thisflow = rcap;
                        }
                    }
                    for (i = 2; i < pl; i += 2) {
                        long int edge = VECTOR(stack)[i];
                        MYCAP(edge) -= thisflow;
                    }
                }

            } while (!igraph_vector_long_empty(&stack));

            igraph_vector_destroy(&mycap);
            igraph_vector_int_destroy(&added);
            igraph_vector_long_destroy(&stack);
            IGRAPH_FINALLY_CLEAN(3);
        }

        /* ----------------------------------------------------------- */

        IGRAPH_CHECK(igraph_vector_resize(flow, no_of_orig_edges));
        for (i = 0, j = 0; i < no_of_edges; i += 2, j++) {
            long int pos = (long int) VECTOR(rank)[i];
            VECTOR(*flow)[j] = (capacity ? VECTOR(*capacity)[j] : 1.0) -
                               RESCAP(pos);
        }

        igraph_vector_destroy(&rank);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_dbuckets_destroy(&ibuckets);
    igraph_buckets_destroy(&buckets);
    igraph_vector_long_destroy(&current);
    igraph_vector_long_destroy(&first);
    igraph_vector_long_destroy(&distance);
    igraph_vector_destroy(&excess);
    igraph_vector_destroy(&rescap);
    igraph_vector_long_destroy(&rev);
    igraph_vector_long_destroy(&to);
    igraph_dqueue_long_destroy(&bfsq);
    IGRAPH_FINALLY_CLEAN(10);

    return 0;
}

/**
 * \function igraph_maxflow_value
 * \brief Maximum flow in a network with the push/relabel algorithm
 *
 * </para><para>This function implements the Goldberg-Tarjan algorithm for
 * calculating value of the maximum flow in a directed or undirected
 * graph. The algorithm was given in Andrew V. Goldberg, Robert
 * E. Tarjan: A New Approach to the Maximum-Flow Problem, Journal of
 * the ACM, 35(4), 921-940, 1988. </para>
 *
 * <para> The input of the function is a graph, a vector
 * of real numbers giving the capacity of the edges and two vertices
 * of the graph, the source and the target. A flow is a function
 * assigning positive real numbers to the edges and satisfying two
 * requirements: (1) the flow value is less than the capacity of the
 * edge and (2) at each vertex except the source and the target, the
 * incoming flow (i.e. the sum of the flow on the incoming edges) is
 * the same as the outgoing flow (i.e. the sum of the flow on the
 * outgoing edges). The value of the flow is the incoming flow at the
 * target vertex. The maximum flow is the flow with the maximum
 * value. </para>
 *
 * <para> According to a theorem by Ford and Fulkerson
 * (L. R. Ford Jr. and D. R. Fulkerson. Maximal flow through a
 * network. Canadian J. Math., 8:399-404, 1956.) the maximum flow
 * between two vertices is the same as the
 * minimum cut between them (also called the minimum s-t cut). So \ref
 * igraph_st_mincut_value() gives the same result in all cases as \c
 * igraph_maxflow_value().</para>
 *
 * <para> Note that the value of the maximum flow is the same as the
 * minimum cut in the graph.
 * \param graph The input graph, either directed or undirected.
 * \param value Pointer to a real number, the result will be placed here.
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \param capacity Vector containing the capacity of the edges. If NULL, then
 *        every edge is considered to have capacity 1.0.
 * \param stats Counts of the number of different operations
 *        preformed by the algorithm are stored here.
 * \return Error code.
 *
 * Time complexity: O(|V|^3).
 *
 * \sa \ref igraph_maxflow() to calculate the actual flow.
 * \ref igraph_mincut_value(), \ref igraph_edge_connectivity(),
 * \ref igraph_vertex_connectivity() for
 * properties based on the maximum flow.
 */

int igraph_maxflow_value(const igraph_t *graph, igraph_real_t *value,
                         igraph_integer_t source, igraph_integer_t target,
                         const igraph_vector_t *capacity,
                         igraph_maxflow_stats_t *stats) {

    return igraph_maxflow(graph, value, /*flow=*/ 0, /*cut=*/ 0,
                          /*partition=*/ 0, /*partition1=*/ 0,
                          source, target, capacity, stats);
}

/**
 * \function igraph_st_mincut_value
 * \brief The minimum s-t cut in a graph
 *
 * </para><para> The minimum s-t cut in a weighted (=valued) graph is the
 * total minimum edge weight needed to remove from the graph to
 * eliminate all paths from a given vertex (\c source) to
 * another vertex (\c target). Directed paths are considered in
 * directed graphs, and undirected paths in undirected graphs.  </para>
 *
 * <para> The minimum s-t cut between two vertices is known to be same
 * as the maximum flow between these two vertices. So this function
 * calls \ref igraph_maxflow_value() to do the calculation.
 * \param graph The input graph.
 * \param value Pointer to a real variable, the result will be stored
 *        here.
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \param capacity Pointer to the capacity vector, it should contain
 *        non-negative numbers and its length should be the same the
 *        the number of edges in the graph. It can be a null pointer, then
 *        every edge has unit capacity.
 * \return Error code.
 *
 * Time complexity: O(|V|^3), see also the discussion for \ref
 * igraph_maxflow_value(), |V| is the number of vertices.
 */

int igraph_st_mincut_value(const igraph_t *graph, igraph_real_t *value,
                           igraph_integer_t source, igraph_integer_t target,
                           const igraph_vector_t *capacity) {

    if (source == target) {
        IGRAPH_ERROR("source and target vertices are the same", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_maxflow_value(graph, value, source, target, capacity, 0));

    return 0;
}

/**
 * \function igraph_st_mincut
 * Minimum cut between a source and a target vertex
 *
 * Finds the edge set that has the smallest total capacity among all
 * edge sets that disconnect the source and target vertices.
 *
 * </para><para>The calculation is performed using maximum flow
 * techniques, by calling \ref igraph_maxflow().
 * \param graph The input graph.
 * \param value Pointer to a real variable, the value of the cut is
 *        stored here.
 * \param cut Pointer to a real vector, the edge ids that are included
 *        in the cut are stored here. This argument is ignored if it
 *        is a null pointer.
 * \param partition Pointer to a real vector, the vertex ids of the
 *        vertices in the first partition of the cut are stored
 *        here. The first partition is always the one that contains the
 *        source vertex. This argument is ignored if it is a null pointer.
 * \param partition2 Pointer to a real vector, the vertex ids of the
 *        vertices in the second partition of the cut are stored here.
 *        The second partition is always the one that contains the
 *        target vertex. This argument is ignored if it is a null pointer.
 * \param source Integer, the id of the source vertex.
 * \param target Integer, the id of the target vertex.
 * \param capacity Vector containing the capacity of the edges. If a
 *        null pointer, then every edge is considered to have capacity
 *        1.0.
 * \return Error code.
 *
 * \sa \ref igraph_maxflow().
 *
 * Time complexity: see \ref igraph_maxflow().
 */

int igraph_st_mincut(const igraph_t *graph, igraph_real_t *value,
                     igraph_vector_t *cut, igraph_vector_t *partition,
                     igraph_vector_t *partition2,
                     igraph_integer_t source, igraph_integer_t target,
                     const igraph_vector_t *capacity) {

    return igraph_maxflow(graph, value, /*flow=*/ 0,
                          cut, partition, partition2,
                          source, target, capacity, 0);
}

/* This is a flow-based version, but there is a better one
   for undirected graphs */

/* int igraph_i_mincut_value_undirected(const igraph_t *graph, */
/*                   igraph_real_t *res, */
/*                   const igraph_vector_t *capacity) { */

/*   long int no_of_edges=igraph_ecount(graph); */
/*   long int no_of_nodes=igraph_vcount(graph); */
/*   igraph_vector_t edges; */
/*   igraph_vector_t newcapacity; */
/*   igraph_t newgraph; */
/*   long int i; */

/*   /\* We need to convert this to directed by hand, since we need to be */
/*      sure that the edge ids will be handled properly to build the new */
/*      capacity vector. *\/ */

/*   IGRAPH_VECTOR_INIT_FINALLY(&edges, 0); */
/*   IGRAPH_VECTOR_INIT_FINALLY(&newcapacity, no_of_edges*2); */
/*   IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges*4)); */
/*   IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0)); */
/*   IGRAPH_CHECK(igraph_vector_resize(&edges, no_of_edges*4)); */
/*   for (i=0; i<no_of_edges; i++) { */
/*     VECTOR(edges)[no_of_edges*2+i*2] = VECTOR(edges)[i*2+1]; */
/*     VECTOR(edges)[no_of_edges*2+i*2+1] = VECTOR(edges)[i*2]; */
/*     VECTOR(newcapacity)[i] = VECTOR(newcapacity)[no_of_edges+i] =  */
/*       capacity ? VECTOR(*capacity)[i] : 1.0 ; */
/*   } */

/*   IGRAPH_CHECK(igraph_create(&newgraph, &edges, no_of_nodes, IGRAPH_DIRECTED)); */
/*   IGRAPH_FINALLY(igraph_destroy, &newgraph); */

/*   IGRAPH_CHECK(igraph_mincut_value(&newgraph, res, &newcapacity)); */

/*   igraph_destroy(&newgraph); */
/*   igraph_vector_destroy(&edges); */
/*   igraph_vector_destroy(&newcapacity); */
/*   IGRAPH_FINALLY_CLEAN(3); */

/*   return 0; */
/* } */

/*
 * This is the Stoer-Wagner algorithm, it works for calculating the
 * minimum cut for undirected graphs, for the whole graph.
 * I.e. this is basically the edge-connectivity of the graph.
 * It can also calculate the cut itself, not just the cut value.
 */

static int igraph_i_mincut_undirected(const igraph_t *graph,
                                      igraph_real_t *res,
                                      igraph_vector_t *partition,
                                      igraph_vector_t *partition2,
                                      igraph_vector_t *cut,
                                      const igraph_vector_t *capacity) {

    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph);
    igraph_integer_t no_of_edges = (igraph_integer_t) igraph_ecount(graph);

    igraph_i_cutheap_t heap;
    igraph_real_t mincut = IGRAPH_INFINITY; /* infinity */
    long int i;

    igraph_adjlist_t adjlist;
    igraph_inclist_t inclist;

    igraph_vector_t mergehist;
    igraph_bool_t calc_cut = partition || partition2 || cut;
    long int act_step = 0, mincut_step = 0;

    if (capacity && igraph_vector_size(capacity) != no_of_edges) {
        IGRAPH_ERROR("Invalid capacity vector size", IGRAPH_EINVAL);
    }

    /* Check if the graph is connected at all */
    {
        igraph_vector_t memb, csize;
        igraph_integer_t no;
        IGRAPH_VECTOR_INIT_FINALLY(&memb, 0);
        IGRAPH_VECTOR_INIT_FINALLY(&csize, 0);
        IGRAPH_CHECK(igraph_clusters(graph, &memb, &csize, &no,
                                     /*mode=*/ IGRAPH_WEAK));
        if (no != 1) {
            if (res) {
                *res = 0;
            }
            if (cut) {
                igraph_vector_clear(cut);
            }
            if (partition) {
                int j = 0;
                IGRAPH_CHECK(igraph_vector_resize(partition,
                                                  (long int) VECTOR(csize)[0]));
                for (i = 0; i < no_of_nodes; i++) {
                    if (VECTOR(memb)[i] == 0) {
                        VECTOR(*partition)[j++] = i;
                    }
                }
            }
            if (partition2) {
                int j = 0;
                IGRAPH_CHECK(igraph_vector_resize(partition2, no_of_nodes -
                                                  (long int) VECTOR(csize)[0]));
                for (i = 0; i < no_of_nodes; i++) {
                    if (VECTOR(memb)[i] != 0) {
                        VECTOR(*partition2)[j++] = i;
                    }
                }
            }
        }
        igraph_vector_destroy(&csize);
        igraph_vector_destroy(&memb);
        IGRAPH_FINALLY_CLEAN(2);

        if (no != 1) {
            return 0;
        }
    }

    if (calc_cut) {
        IGRAPH_VECTOR_INIT_FINALLY(&mergehist, 0);
        IGRAPH_CHECK(igraph_vector_reserve(&mergehist, no_of_nodes * 2));
    }

    IGRAPH_CHECK(igraph_i_cutheap_init(&heap, no_of_nodes));
    IGRAPH_FINALLY(igraph_i_cutheap_destroy, &heap);

    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    while (igraph_i_cutheap_size(&heap) >= 2) {

        long int last;
        igraph_real_t acut;
        long int a, n;

        igraph_vector_int_t *edges, *edges2;
        igraph_vector_int_t *neis, *neis2;

        do {
            a = igraph_i_cutheap_popmax(&heap);

            /* update the weights of the active vertices connected to a */
            edges = igraph_inclist_get(&inclist, a);
            neis = igraph_adjlist_get(&adjlist, a);
            n = igraph_vector_int_size(edges);
            for (i = 0; i < n; i++) {
                igraph_integer_t edge = (igraph_integer_t) VECTOR(*edges)[i];
                igraph_integer_t to = (igraph_integer_t) VECTOR(*neis)[i];
                igraph_real_t weight = capacity ? VECTOR(*capacity)[(long int)edge] : 1.0;
                igraph_i_cutheap_update(&heap, to, weight);
            }

        } while (igraph_i_cutheap_active_size(&heap) > 1);

        /* Now, there is only one active vertex left,
           calculate the cut of the phase */
        acut = igraph_i_cutheap_maxvalue(&heap);
        last = igraph_i_cutheap_popmax(&heap);

        if (acut < mincut) {
            mincut = acut;
            mincut_step = act_step;
        }

        if (mincut == 0) {
            break;
        }

        /* And contract the last and the remaining vertex (a and last) */
        /* Before actually doing that, make some notes */
        act_step++;
        if (calc_cut) {
            IGRAPH_CHECK(igraph_vector_push_back(&mergehist, a));
            IGRAPH_CHECK(igraph_vector_push_back(&mergehist, last));
        }
        /* First remove the a--last edge if there is one, a is still the
           last deactivated vertex */
        edges = igraph_inclist_get(&inclist, a);
        neis = igraph_adjlist_get(&adjlist, a);
        n = igraph_vector_int_size(edges);
        for (i = 0; i < n; ) {
            if (VECTOR(*neis)[i] == last) {
                VECTOR(*neis)[i] = VECTOR(*neis)[n - 1];
                VECTOR(*edges)[i] = VECTOR(*edges)[n - 1];
                igraph_vector_int_pop_back(neis);
                igraph_vector_int_pop_back(edges);
                n--;
            } else {
                i++;
            }
        }

        edges = igraph_inclist_get(&inclist, last);
        neis = igraph_adjlist_get(&adjlist, last);
        n = igraph_vector_int_size(edges);
        for (i = 0; i < n; ) {
            if (VECTOR(*neis)[i] == a) {
                VECTOR(*neis)[i] = VECTOR(*neis)[n - 1];
                VECTOR(*edges)[i] = VECTOR(*edges)[n - 1];
                igraph_vector_int_pop_back(neis);
                igraph_vector_int_pop_back(edges);
                n--;
            } else {
                i++;
            }
        }

        /* Now rewrite the edge lists of last's neighbors */
        neis = igraph_adjlist_get(&adjlist, last);
        n = igraph_vector_int_size(neis);
        for (i = 0; i < n; i++) {
            igraph_integer_t nei = (igraph_integer_t) VECTOR(*neis)[i];
            long int n2, j;
            neis2 = igraph_adjlist_get(&adjlist, nei);
            n2 = igraph_vector_int_size(neis2);
            for (j = 0; j < n2; j++) {
                if (VECTOR(*neis2)[j] == last) {
                    VECTOR(*neis2)[j] = a;
                }
            }
        }

        /* And append the lists of last to the lists of a */
        edges = igraph_inclist_get(&inclist, a);
        neis = igraph_adjlist_get(&adjlist, a);
        edges2 = igraph_inclist_get(&inclist, last);
        neis2 = igraph_adjlist_get(&adjlist, last);
        IGRAPH_CHECK(igraph_vector_int_append(edges, edges2));
        IGRAPH_CHECK(igraph_vector_int_append(neis, neis2));
        igraph_vector_int_clear(edges2); /* TODO: free it */
        igraph_vector_int_clear(neis2);  /* TODO: free it */

        /* Remove the deleted vertex from the heap entirely */
        igraph_i_cutheap_reset_undefine(&heap, last);
    }

    *res = mincut;

    igraph_inclist_destroy(&inclist);
    igraph_adjlist_destroy(&adjlist);
    igraph_i_cutheap_destroy(&heap);
    IGRAPH_FINALLY_CLEAN(3);

    if (calc_cut) {
        long int bignode = (long int) VECTOR(mergehist)[2 * mincut_step + 1];
        long int i, idx;
        long int size = 1;
        char *mark;
        mark = IGRAPH_CALLOC(no_of_nodes, char);
        if (!mark) {
            IGRAPH_ERROR("Not enough memory for minimum cut", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, mark);

        /* first count the vertices in the partition */
        mark[bignode] = 1;
        for (i = mincut_step - 1; i >= 0; i--) {
            if ( mark[ (long int) VECTOR(mergehist)[2 * i] ] ) {
                size++;
                mark [ (long int) VECTOR(mergehist)[2 * i + 1] ] = 1;
            }
        }

        /* now store them, if requested */
        if (partition) {
            IGRAPH_CHECK(igraph_vector_resize(partition, size));
            idx = 0;
            VECTOR(*partition)[idx++] = bignode;
            for (i = mincut_step - 1; i >= 0; i--) {
                if (mark[ (long int) VECTOR(mergehist)[2 * i] ]) {
                    VECTOR(*partition)[idx++] = VECTOR(mergehist)[2 * i + 1];
                }
            }
        }

        /* The other partition too? */
        if (partition2) {
            IGRAPH_CHECK(igraph_vector_resize(partition2, no_of_nodes - size));
            idx = 0;
            for (i = 0; i < no_of_nodes; i++) {
                if (!mark[i]) {
                    VECTOR(*partition2)[idx++] = i;
                }
            }
        }

        /* The edges in the cut are also requested? */
        /* We want as few memory allocated for 'cut' as possible,
           so we first collect the edges in mergehist, we don't
           need that anymore. Then we copy it to 'cut';  */
        if (cut) {
            igraph_integer_t from, to;
            igraph_vector_clear(&mergehist);
            for (i = 0; i < no_of_edges; i++) {
                igraph_edge(graph, (igraph_integer_t) i, &from, &to);
                if ((mark[(long int)from] && !mark[(long int)to]) ||
                    (mark[(long int)to] && !mark[(long int)from])) {
                    IGRAPH_CHECK(igraph_vector_push_back(&mergehist, i));
                }
            }
            igraph_vector_clear(cut);
            IGRAPH_CHECK(igraph_vector_append(cut, &mergehist));
        }

        igraph_free(mark);
        igraph_vector_destroy(&mergehist);
        IGRAPH_FINALLY_CLEAN(2);
    }

    return 0;
}

static int igraph_i_mincut_directed(const igraph_t *graph,
                                    igraph_real_t *value,
                                    igraph_vector_t *partition,
                                    igraph_vector_t *partition2,
                                    igraph_vector_t *cut,
                                    const igraph_vector_t *capacity) {
    long int i;
    long int no_of_nodes = igraph_vcount(graph);
    igraph_real_t flow;
    igraph_real_t minmaxflow = IGRAPH_INFINITY;
    igraph_vector_t mypartition, mypartition2, mycut;
    igraph_vector_t *ppartition = 0, *ppartition2 = 0, *pcut = 0;
    igraph_vector_t bestpartition, bestpartition2, bestcut;

    if (partition) {
        IGRAPH_VECTOR_INIT_FINALLY(&bestpartition, 0);
    }
    if (partition2) {
        IGRAPH_VECTOR_INIT_FINALLY(&bestpartition2, 0);
    }
    if (cut) {
        IGRAPH_VECTOR_INIT_FINALLY(&bestcut, 0);
    }

    if (partition) {
        IGRAPH_VECTOR_INIT_FINALLY(&mypartition, 0);
        ppartition = &mypartition;
    }
    if (partition2) {
        IGRAPH_VECTOR_INIT_FINALLY(&mypartition2, 0);
        ppartition2 = &mypartition2;
    }
    if (cut) {
        IGRAPH_VECTOR_INIT_FINALLY(&mycut, 0);
        pcut = &mycut;
    }

    for (i = 1; i < no_of_nodes; i++) {
        IGRAPH_CHECK(igraph_maxflow(graph, /*value=*/ &flow, /*flow=*/ 0,
                                    pcut, ppartition, ppartition2, /*source=*/ 0,
                                    /*target=*/ (igraph_integer_t) i, capacity, 0));
        if (flow < minmaxflow) {
            minmaxflow = flow;
            if (cut) {
                IGRAPH_CHECK(igraph_vector_update(&bestcut, &mycut));
            }
            if (partition) {
                IGRAPH_CHECK(igraph_vector_update(&bestpartition, &mypartition));
            }
            if (partition2) {
                IGRAPH_CHECK(igraph_vector_update(&bestpartition2, &mypartition2));
            }

            if (minmaxflow == 0) {
                break;
            }
        }
        IGRAPH_CHECK(igraph_maxflow(graph, /*value=*/ &flow, /*flow=*/ 0,
                                    pcut, ppartition, ppartition2,
                                    /*source=*/ (igraph_integer_t) i,
                                    /*target=*/ 0, capacity, 0));
        if (flow < minmaxflow) {
            minmaxflow = flow;
            if (cut) {
                IGRAPH_CHECK(igraph_vector_update(&bestcut, &mycut));
            }
            if (partition) {
                IGRAPH_CHECK(igraph_vector_update(&bestpartition, &mypartition));
            }
            if (partition2) {
                IGRAPH_CHECK(igraph_vector_update(&bestpartition2, &mypartition2));
            }

            if (minmaxflow == 0) {
                break;
            }
        }
    }

    if (value) {
        *value = minmaxflow;
    }

    if (cut) {
        igraph_vector_destroy(&mycut);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (partition) {
        igraph_vector_destroy(&mypartition);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (partition2) {
        igraph_vector_destroy(&mypartition2);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (cut) {
        IGRAPH_CHECK(igraph_vector_update(cut, &bestcut));
        igraph_vector_destroy(&bestcut);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (partition2) {
        IGRAPH_CHECK(igraph_vector_update(partition2, &bestpartition2));
        igraph_vector_destroy(&bestpartition2);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (partition) {
        IGRAPH_CHECK(igraph_vector_update(partition, &bestpartition));
        igraph_vector_destroy(&bestpartition);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/**
 * \function igraph_mincut
 * \brief Calculates the minimum cut in a graph.
 *
 * This function calculates the minimum cut in a graph.
 * The minimum cut is the minimum set of edges which needs to be
 * removed to disconnect the graph. The minimum is calculated using
 * the weights (\p capacity) of the edges, so the cut with the minimum
 * total capacity is calculated.
 *
 * </para><para> For directed graphs an implementation based on
 * calculating 2|V|-2 maximum flows is used.
 * For undirected graphs we use the Stoer-Wagner
 * algorithm, as described in M. Stoer and F. Wagner: A simple min-cut
 * algorithm, Journal of the ACM, 44 585-591, 1997.
 *
 * </para><para>
 * The first implementation of the actual cut calculation for
 * undirected graphs was made by Gregory Benison, thanks Greg.
 * \param graph The input graph.
 * \param value Pointer to a float, the value of the cut will be
 *    stored here.
 * \param partition Pointer to an initialized vector, the ids
 *    of the vertices in the first partition after separating the
 *    graph will be stored here. The vector will be resized as
 *    needed. This argument is ignored if it is a NULL pointer.
 * \param partition2 Pointer to an initialized vector the ids
 *    of the vertices in the second partition will be stored here.
 *    The vector will be resized as needed. This argument is ignored
 *    if it is a NULL pointer.
 * \param cut Pointer to an initialized vector, the ids of the edges
 *    in the cut will be stored here. This argument is ignored if it
 *    is a NULL pointer.
 * \param capacity A numeric vector giving the capacities of the
 *    edges. If a null pointer then all edges have unit capacity.
 * \return Error code.
 *
 * \sa \ref igraph_mincut_value(), a simpler interface for calculating
 * the value of the cut only.
 *
 * Time complexity: for directed graphs it is O(|V|^4), but see the
 * remarks at \ref igraph_maxflow(). For undirected graphs it is
 * O(|V||E|+|V|^2 log|V|). |V| and |E| are the number of vertices and
 * edges respectively.
 *
 * \example examples/simple/igraph_mincut.c
 */

int igraph_mincut(const igraph_t *graph,
                  igraph_real_t *value,
                  igraph_vector_t *partition,
                  igraph_vector_t *partition2,
                  igraph_vector_t *cut,
                  const igraph_vector_t *capacity) {

    if (igraph_is_directed(graph)) {
        if (partition || partition2 || cut) {
            igraph_i_mincut_directed(graph, value, partition, partition2, cut,
                                     capacity);
        } else {
            return igraph_mincut_value(graph, value, capacity);
        }
    } else {
        IGRAPH_CHECK(igraph_i_mincut_undirected(graph, value, partition,
                                                partition2, cut, capacity));
        return IGRAPH_SUCCESS;
    }

    return 0;
}


static int igraph_i_mincut_value_undirected(const igraph_t *graph,
                                            igraph_real_t *res,
                                            const igraph_vector_t *capacity) {
    return igraph_i_mincut_undirected(graph, res, 0, 0, 0, capacity);
}

/**
 * \function igraph_mincut_value
 * \brief The minimum edge cut in a graph
 *
 * </para><para> The minimum edge cut in a graph is the total minimum
 * weight of the edges needed to remove from the graph to make the
 * graph \em not strongly connected. (If the original graph is not
 * strongly connected then this is zero.) Note that in undirected
 * graphs strong connectedness is the same as weak connectedness. </para>
 *
 * <para> The minimum cut can be calculated with maximum flow
 * techniques, although the current implementation does this only for
 * directed graphs and a separate non-flow based implementation is
 * used for undirected graphs. See Mechthild Stoer and Frank Wagner: A
 * simple min-cut algorithm, Journal of the ACM 44 585--591, 1997.
 * For directed graphs
 * the maximum flow is calculated between a fixed vertex and all the
 * other vertices in the graph and this is done in both
 * directions. Then the minimum is taken to get the minimum cut.
 *
 * \param graph The input graph.
 * \param res Pointer to a real variable, the result will be stored
 *    here.
 * \param capacity Pointer to the capacity vector, it should contain
 *    the same number of non-negative numbers as the number of edges in
 *    the graph. If a null pointer then all edges will have unit capacity.
 * \return Error code.
 *
 * \sa \ref igraph_mincut(), \ref igraph_maxflow_value(), \ref
 * igraph_st_mincut_value().
 *
 * Time complexity: O(log(|V|)*|V|^2) for undirected graphs and
 * O(|V|^4) for directed graphs, but see also the discussion at the
 * documentation of \ref igraph_maxflow_value().
 */

int igraph_mincut_value(const igraph_t *graph, igraph_real_t *res,
                        const igraph_vector_t *capacity) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_real_t minmaxflow, flow;
    long int i;

    minmaxflow = IGRAPH_INFINITY;

    if (!igraph_is_directed(graph)) {
        IGRAPH_CHECK(igraph_i_mincut_value_undirected(graph, res, capacity));
        return 0;
    }

    for (i = 1; i < no_of_nodes; i++) {
        IGRAPH_CHECK(igraph_maxflow_value(graph, &flow, 0, (igraph_integer_t) i,
                                          capacity, 0));
        if (flow < minmaxflow) {
            minmaxflow = flow;
            if (flow == 0) {
                break;
            }
        }
        IGRAPH_CHECK(igraph_maxflow_value(graph, &flow, (igraph_integer_t) i, 0,
                                          capacity, 0));
        if (flow < minmaxflow) {
            minmaxflow = flow;
            if (flow == 0) {
                break;
            }
        }
    }

    if (res) {
        *res = minmaxflow;
    }

    return 0;
}

static int igraph_i_st_vertex_connectivity_directed(const igraph_t *graph,
                                                    igraph_integer_t *res,
                                                    igraph_integer_t source,
                                                    igraph_integer_t target,
                                                    igraph_vconn_nei_t neighbors) {

    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph);
    igraph_integer_t no_of_edges = (igraph_integer_t) igraph_ecount(graph);
    igraph_vector_t edges;
    igraph_real_t real_res;
    igraph_t newgraph;
    long int i;
    igraph_bool_t conn1;

    if (source < 0 || source >= no_of_nodes || target < 0 || target >= no_of_nodes) {
        IGRAPH_ERROR("Invalid source or target vertex", IGRAPH_EINVAL);
    }

    switch (neighbors) {
    case IGRAPH_VCONN_NEI_ERROR:
        IGRAPH_CHECK(igraph_are_connected(graph, source, target, &conn1));
        if (conn1) {
            IGRAPH_ERROR("vertices connected", IGRAPH_EINVAL);
        }
        break;
    case IGRAPH_VCONN_NEI_NEGATIVE:
        IGRAPH_CHECK(igraph_are_connected(graph, source, target, &conn1));
        if (conn1) {
            *res = -1;
            return 0;
        }
        break;
    case IGRAPH_VCONN_NEI_NUMBER_OF_NODES:
        IGRAPH_CHECK(igraph_are_connected(graph, source, target, &conn1));
        if (conn1) {
            *res = no_of_nodes;
            return 0;
        }
        break;
    case IGRAPH_VCONN_NEI_IGNORE:
        break;
    default:
        IGRAPH_ERROR("Unknown `igraph_vconn_nei_t'", IGRAPH_EINVAL);
        break;
    }

    /* Create the new graph */

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, 2 * (no_of_edges + no_of_nodes)));
    IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));
    IGRAPH_CHECK(igraph_vector_resize(&edges, 2 * (no_of_edges + no_of_nodes)));

    for (i = 0; i < 2 * no_of_edges; i += 2) {
        igraph_integer_t to = (igraph_integer_t) VECTOR(edges)[i + 1];
        if (to != source && to != target) {
            VECTOR(edges)[i + 1] = no_of_nodes + to;
        }
    }

    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(edges)[ 2 * (no_of_edges + i)   ] = no_of_nodes + i;
        VECTOR(edges)[ 2 * (no_of_edges + i) + 1 ] = i;
    }

    IGRAPH_CHECK(igraph_create(&newgraph, &edges, 2 * no_of_nodes,
                               igraph_is_directed(graph)));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_destroy, &newgraph);

    /* Do the maximum flow */

    IGRAPH_CHECK(igraph_maxflow_value(&newgraph, &real_res,
                                      source, target, 0, 0));
    *res = (igraph_integer_t)real_res;

    igraph_destroy(&newgraph);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_st_vertex_connectivity_undirected(const igraph_t *graph,
                                                      igraph_integer_t *res,
                                                      igraph_integer_t source,
                                                      igraph_integer_t target,
                                                      igraph_vconn_nei_t neighbors) {

    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph);
    igraph_t newgraph;
    igraph_bool_t conn;

    if (source < 0 || source >= no_of_nodes || target < 0 || target >= no_of_nodes) {
        IGRAPH_ERROR("Invalid source or target vertex", IGRAPH_EINVAL);
    }

    switch (neighbors) {
    case IGRAPH_VCONN_NEI_ERROR:
        IGRAPH_CHECK(igraph_are_connected(graph, source, target, &conn));
        if (conn) {
            IGRAPH_ERROR("vertices connected", IGRAPH_EINVAL);
        }
        break;
    case IGRAPH_VCONN_NEI_NEGATIVE:
        IGRAPH_CHECK(igraph_are_connected(graph, source, target, &conn));
        if (conn) {
            *res = -1;
            return 0;
        }
        break;
    case IGRAPH_VCONN_NEI_NUMBER_OF_NODES:
        IGRAPH_CHECK(igraph_are_connected(graph, source, target, &conn));
        if (conn) {
            *res = no_of_nodes;
            return 0;
        }
        break;
    case IGRAPH_VCONN_NEI_IGNORE:
        break;
    default:
        IGRAPH_ERROR("Unknown `igraph_vconn_nei_t'", IGRAPH_EINVAL);
        break;
    }

    IGRAPH_CHECK(igraph_copy(&newgraph, graph));
    IGRAPH_FINALLY(igraph_destroy, &newgraph);
    IGRAPH_CHECK(igraph_to_directed(&newgraph, IGRAPH_TO_DIRECTED_MUTUAL));

    IGRAPH_CHECK(igraph_i_st_vertex_connectivity_directed(&newgraph, res,
                 source, target,
                 IGRAPH_VCONN_NEI_IGNORE));

    igraph_destroy(&newgraph);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_st_vertex_connectivity
 * \brief The vertex connectivity of a pair of vertices
 *
 * </para><para>The vertex connectivity of two vertices (\c source and
 * \c target) is the minimum number of vertices that have to be
 * deleted to eliminate all paths from \c source to \c
 * target. Directed paths are considered in directed graphs.</para>
 *
 * <para>The vertex connectivity of a pair is the same as the number
 * of different (i.e. node-independent) paths from source to
 * target.</para>
 *
 * <para>The current implementation uses maximum flow calculations to
 * obtain the result.
 * \param graph The input graph.
 * \param res Pointer to an integer, the result will be stored here.
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \param neighbors A constant giving what to do if the two vertices
 *     are connected. Possible values:
 *     \c IGRAPH_VCONN_NEI_ERROR, stop with an error message,
 *     \c IGRAPH_VCONN_NEGATIVE, return -1.
 *     \c IGRAPH_VCONN_NUMBER_OF_NODES, return the number of nodes.
 *     \c IGRAPH_VCONN_IGNORE, ignore the fact that the two vertices
 *        are connected and calculate the number of vertices needed
 *        to eliminate all paths except for the trivial (direct) paths
 *        between \p source and \p vertex. TODO: what about neighbors?
 * \return Error code.
 *
 * Time complexity: O(|V|^3), but see the discussion at \ref
 * igraph_maxflow_value().
 *
 * \sa \ref igraph_vertex_connectivity(),
 * \ref igraph_edge_connectivity(),
 * \ref igraph_maxflow_value().
 */

int igraph_st_vertex_connectivity(const igraph_t *graph,
                                  igraph_integer_t *res,
                                  igraph_integer_t source,
                                  igraph_integer_t target,
                                  igraph_vconn_nei_t neighbors) {

    if (source == target) {
        IGRAPH_ERROR("source and target vertices are the same", IGRAPH_EINVAL);
    }

    if (igraph_is_directed(graph)) {
        IGRAPH_CHECK(igraph_i_st_vertex_connectivity_directed(graph, res,
                     source, target,
                     neighbors));
    } else {
        IGRAPH_CHECK(igraph_i_st_vertex_connectivity_undirected(graph, res,
                     source, target,
                     neighbors));
    }

    return 0;
}

static int igraph_i_vertex_connectivity_directed(const igraph_t *graph,
                                                 igraph_integer_t *res) {

    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph);
    long int i, j;
    igraph_integer_t minconn = no_of_nodes - 1, conn = 0;

    for (i = 0; i < no_of_nodes; i++) {
        for (j = 0; j < no_of_nodes; j++) {
            if (i == j) {
                continue;
            }

            IGRAPH_ALLOW_INTERRUPTION();

            IGRAPH_CHECK(igraph_st_vertex_connectivity(graph, &conn,
                         (igraph_integer_t) i,
                         (igraph_integer_t) j,
                         IGRAPH_VCONN_NEI_NUMBER_OF_NODES));
            if (conn < minconn) {
                minconn = conn;
                if (conn == 0) {
                    break;
                }
            }
        }
        if (conn == 0) {
            break;
        }
    }

    if (res) {
        *res = minconn;
    }

    return 0;
}

static int igraph_i_vertex_connectivity_undirected(const igraph_t *graph,
                                                   igraph_integer_t *res) {
    igraph_t newgraph;

    IGRAPH_CHECK(igraph_copy(&newgraph, graph));
    IGRAPH_FINALLY(igraph_destroy, &newgraph);
    IGRAPH_CHECK(igraph_to_directed(&newgraph, IGRAPH_TO_DIRECTED_MUTUAL));

    IGRAPH_CHECK(igraph_i_vertex_connectivity_directed(&newgraph, res));

    igraph_destroy(&newgraph);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/* Use that vertex.connectivity(G) <= edge.connectivity(G) <= min(degree(G)) */
static int igraph_i_connectivity_checks(const igraph_t *graph,
                                        igraph_integer_t *res,
                                        igraph_bool_t *found) {
    igraph_bool_t conn;
    *found = 0;

    if (igraph_vcount(graph) == 0) {
        *res = 0;
        *found = 1;
        return 0;
    }

    IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_STRONG));
    if (!conn) {
        *res = 0;
        *found = 1;
    } else {
        igraph_vector_t degree;
        IGRAPH_VECTOR_INIT_FINALLY(&degree, 0);
        if (!igraph_is_directed(graph)) {
            IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
                                       IGRAPH_OUT, IGRAPH_LOOPS));
            if (igraph_vector_min(&degree) == 1) {
                *res = 1;
                *found = 1;
            }
        } else {
            /* directed, check both in- & out-degree */
            IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
                                       IGRAPH_OUT, IGRAPH_LOOPS));
            if (igraph_vector_min(&degree) == 1) {
                *res = 1;
                *found = 1;
            } else {
                IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
                                           IGRAPH_IN, IGRAPH_LOOPS));
                if (igraph_vector_min(&degree) == 1) {
                    *res = 1;
                    *found = 1;
                }
            }
        }
        igraph_vector_destroy(&degree);
        IGRAPH_FINALLY_CLEAN(1);
    }
    return 0;
}

/**
 * \function igraph_vertex_connectivity
 * The vertex connectivity of a graph
 *
 * </para><para> The vertex connectivity of a graph is the minimum
 * vertex connectivity along each pairs of vertices in the graph.
 * </para>
 * <para> The vertex connectivity of a graph is the same as group
 * cohesion as defined in Douglas R. White and Frank Harary: The
 * cohesiveness of blocks in social networks: node connectivity and
 * conditional density, Sociological Methodology 31:305--359, 2001.
 * \param graph The input graph.
 * \param res Pointer to an integer, the result will be stored here.
 * \param checks Logical constant. Whether to check that the graph is
 *    connected and also the degree of the vertices. If the graph is
 *    not (strongly) connected then the connectivity is obviously zero. Otherwise
 *    if the minimum degree is one then the vertex connectivity is also
 *    one. It is a good idea to perform these checks, as they can be
 *    done quickly compared to the connectivity calculation itself.
 *    They were suggested by Peter McMahan, thanks Peter.
 * \return Error code.
 *
 * Time complexity: O(|V|^5).
 *
 * \sa \ref igraph_st_vertex_connectivity(), \ref igraph_maxflow_value(),
 * and \ref igraph_edge_connectivity().
 */

int igraph_vertex_connectivity(const igraph_t *graph, igraph_integer_t *res,
                               igraph_bool_t checks) {

    igraph_bool_t ret = 0;

    if (checks) {
        IGRAPH_CHECK(igraph_i_connectivity_checks(graph, res, &ret));
    }

    /* Are we done yet? */
    if (!ret) {
        if (igraph_is_directed(graph)) {
            IGRAPH_CHECK(igraph_i_vertex_connectivity_directed(graph, res));
        } else {
            IGRAPH_CHECK(igraph_i_vertex_connectivity_undirected(graph, res));
        }
    }

    return 0;
}

/**
 * \function igraph_st_edge_connectivity
 * \brief Edge connectivity of a pair of vertices
 *
 * </para><para> The edge connectivity of two vertices (\c source and
 * \c target) in a graph is the minimum number of edges that
 * have to be deleted from the graph to eliminate all paths from \c
 * source to \c target.</para>
 *
 * <para>This function uses the maximum flow algorithm to calculate
 * the edge connectivity.
 * \param graph The input graph, it has to be directed.
 * \param res Pointer to an integer, the result will be stored here.
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \return Error code.
 *
 * Time complexity: O(|V|^3).
 *
 * \sa \ref igraph_maxflow_value(), \ref igraph_edge_connectivity(),
 * \ref igraph_st_vertex_connectivity(), \ref
 * igraph_vertex_connectivity().
 */

int igraph_st_edge_connectivity(const igraph_t *graph, igraph_integer_t *res,
                                igraph_integer_t source,
                                igraph_integer_t target) {
    igraph_real_t flow;

    if (source == target) {
        IGRAPH_ERROR("source and target vertices are the same", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_maxflow_value(graph, &flow, source, target, 0, 0));
    *res = (igraph_integer_t) flow;

    return 0;
}


/**
 * \function igraph_edge_connectivity
 * \brief The minimum edge connectivity in a graph.
 *
 * </para><para> This is the minimum of the edge connectivity over all
 * pairs of vertices in the graph. </para>
 *
 * <para>
 * The edge connectivity of a graph is the same as group adhesion as
 * defined in Douglas R. White and Frank Harary: The cohesiveness of
 * blocks in social networks: node connectivity and conditional
 * density, Sociological Methodology 31:305--359, 2001.
 * \param graph The input graph.
 * \param res Pointer to an integer, the result will be stored here.
 * \param checks Logical constant. Whether to check that the graph is
 *    connected and also the degree of the vertices. If the graph is
 *    not (strongly) connected then the connectivity is obviously zero. Otherwise
 *    if the minimum degree is one then the edge connectivity is also
 *    one. It is a good idea to perform these checks, as they can be
 *    done quickly compared to the connectivity calculation itself.
 *    They were suggested by Peter McMahan, thanks Peter.
 * \return Error code.
 *
 * Time complexity: O(log(|V|)*|V|^2) for undirected graphs and
 * O(|V|^4) for directed graphs, but see also the discussion at the
 * documentation of \ref igraph_maxflow_value().
 *
 * \sa \ref igraph_st_edge_connectivity(), \ref igraph_maxflow_value(),
 * \ref igraph_vertex_connectivity().
 */

int igraph_edge_connectivity(const igraph_t *graph, igraph_integer_t *res,
                             igraph_bool_t checks) {
    igraph_bool_t ret = 0;
    igraph_integer_t number_of_nodes = igraph_vcount(graph);

    /* igraph_mincut_value returns infinity for the singleton graph,
     * which cannot be cast to an integer. We catch this case early
     * and postulate the edge-connectivity of this graph to be 0.
     * This is consistent with what other software packages return. */
    if (number_of_nodes <= 1) {
        *res = 0;
        return 0;
    }

    /* Use that vertex.connectivity(G) <= edge.connectivity(G) <= min(degree(G)) */
    if (checks) {
        IGRAPH_CHECK(igraph_i_connectivity_checks(graph, res, &ret));
    }

    if (!ret) {
        igraph_real_t real_res;
        IGRAPH_CHECK(igraph_mincut_value(graph, &real_res, 0));
        *res = (igraph_integer_t)real_res;
    }

    return 0;
}

/**
 * \function igraph_edge_disjoint_paths
 * \brief The maximum number of edge-disjoint paths between two vertices.
 *
 * </para><para> A set of paths between two vertices is called
 * edge-disjoint if they do not share any edges. The maximum number of
 * edge-disjoint paths are calculated by this function using maximum
 * flow techniques. Directed paths are considered in directed
 * graphs. </para>
 *
 * <para> Note that the number of disjoint paths is the same as the
 * edge connectivity of the two vertices using uniform edge weights.
 * \param graph The input graph, can be directed or undirected.
 * \param res Pointer to an integer variable, the result will be
 *        stored here.
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \return Error code.
 *
 * Time complexity: O(|V|^3), but see the discussion at \ref
 * igraph_maxflow_value().
 *
 * \sa \ref igraph_vertex_disjoint_paths(), \ref
 * igraph_st_edge_connectivity(), \ref igraph_maxflow_value().
 */

int igraph_edge_disjoint_paths(const igraph_t *graph, igraph_integer_t *res,
                               igraph_integer_t source,
                               igraph_integer_t target) {

    igraph_real_t flow;

    if (source == target) {
        IGRAPH_ERROR("Not implemented for source=target", IGRAPH_UNIMPLEMENTED);
    }

    IGRAPH_CHECK(igraph_maxflow_value(graph, &flow, source, target, 0, 0));

    *res = (igraph_integer_t) flow;

    return 0;
}

/**
 * \function igraph_vertex_disjoint_paths
 * \brief Maximum number of vertex-disjoint paths between two vertices.
 *
 * </para><para> A set of paths between two vertices is called
 * vertex-disjoint if they share no vertices. The calculation is
 * performed by using maximum flow techniques. </para>
 *
 * <para> Note that the number of vertex-disjoint paths is the same as
 * the vertex connectivity of the two vertices in most cases (if the
 * two vertices are not connected by an edge).
 * \param graph The input graph.
 * \param res Pointer to an integer variable, the result will be
 *        stored here.
 * \param source The id of the source vertex.
 * \param target The id of the target vertex.
 * \return Error code.
 *
 * Time complexity: O(|V|^3).
 *
 * \sa \ref igraph_edge_disjoint_paths(), \ref
 * igraph_vertex_connectivity(), \ref igraph_maxflow_value().
 */

int igraph_vertex_disjoint_paths(const igraph_t *graph, igraph_integer_t *res,
                                 igraph_integer_t source,
                                 igraph_integer_t target) {

    igraph_bool_t conn;

    if (source == target) {
        IGRAPH_ERROR("The source==target case is not implemented",
                     IGRAPH_UNIMPLEMENTED);
    }

    igraph_are_connected(graph, source, target, &conn);
    if (conn) {
        /* We need to remove every (possibly directed) edge between source
           and target and calculate the disjoint paths on the new
           graph. Finally we add 1 for the removed connection(s).  */
        igraph_es_t es;
        igraph_vector_t v;
        igraph_t newgraph;
        IGRAPH_VECTOR_INIT_FINALLY(&v, 2);
        VECTOR(v)[0] = source;
        VECTOR(v)[1] = target;
        IGRAPH_CHECK(igraph_es_multipairs(&es, &v, IGRAPH_DIRECTED));
        IGRAPH_FINALLY(igraph_es_destroy, &es);

        IGRAPH_CHECK(igraph_copy(&newgraph, graph));
        IGRAPH_FINALLY(igraph_destroy, &newgraph);
        IGRAPH_CHECK(igraph_delete_edges(&newgraph, es));

        if (igraph_is_directed(graph)) {
            IGRAPH_CHECK(igraph_i_st_vertex_connectivity_directed(&newgraph, res,
                         source, target,
                         IGRAPH_VCONN_NEI_IGNORE));
        } else {
            IGRAPH_CHECK(igraph_i_st_vertex_connectivity_undirected(&newgraph, res,
                         source, target,
                         IGRAPH_VCONN_NEI_IGNORE));
        }

        if (res) {
            *res += 1;
        }

        IGRAPH_FINALLY_CLEAN(3);
        igraph_destroy(&newgraph);
        igraph_es_destroy(&es);
        igraph_vector_destroy(&v);
    }

    /* These do nothing if the two vertices are connected,
       so it is safe to call them. */

    if (igraph_is_directed(graph)) {
        IGRAPH_CHECK(igraph_i_st_vertex_connectivity_directed(graph, res,
                     source, target,
                     IGRAPH_VCONN_NEI_IGNORE));
    } else {
        IGRAPH_CHECK(igraph_i_st_vertex_connectivity_undirected(graph, res,
                     source, target,
                     IGRAPH_VCONN_NEI_IGNORE));
    }

    return 0;
}

/**
 * \function igraph_adhesion
 * \brief Graph adhesion, this is (almost) the same as edge connectivity.
 *
 * </para><para> This quantity is defined by White and Harary in
 * The cohesiveness of blocks in social networks: node connectivity and
 * conditional density, (Sociological Methodology 31:305--359, 2001)
 * and basically it is the edge connectivity of the graph
 * with uniform edge weights.
 * \param graph The input graph, either directed or undirected.
 * \param res Pointer to an integer, the result will be stored here.
 * \param checks Logical constant. Whether to check that the graph is
 *    connected and also the degree of the vertices. If the graph is
 *    not (strongly) connected then the adhesion is obviously zero. Otherwise
 *    if the minimum degree is one then the adhesion is also
 *    one. It is a good idea to perform these checks, as they can be
 *    done quickly compared to the edge connectivity calculation itself.
 *    They were suggested by Peter McMahan, thanks Peter.
* \return Error code.
 *
 * Time complexity: O(log(|V|)*|V|^2) for undirected graphs and
 * O(|V|^4) for directed graphs, but see also the discussion at the
 * documentation of \ref igraph_maxflow_value().
 *
 * \sa \ref igraph_cohesion(), \ref igraph_maxflow_value(), \ref
 * igraph_edge_connectivity(), \ref igraph_mincut_value().
 */

int igraph_adhesion(const igraph_t *graph, igraph_integer_t *res,
                    igraph_bool_t checks) {
    return igraph_edge_connectivity(graph, res, checks);
}

/**
 * \function igraph_cohesion
 * \brief Graph cohesion, this is the same as vertex connectivity.
 *
 * </para><para> This quantity was defined by White and Harary in <quote>The
 * cohesiveness of blocks in social networks: node connectivity and
 * conditional density</quote>, (Sociological Methodology 31:305--359, 2001)
 * and it is the same as the vertex connectivity of a
 * graph.
 * \param graph The input graph.
 * \param res Pointer to an integer variable, the result will be
 *        stored here.
 * \param checks Logical constant. Whether to check that the graph is
 *    connected and also the degree of the vertices. If the graph is
 *    not (strongly) connected then the cohesion is obviously zero. Otherwise
 *    if the minimum degree is one then the cohesion is also
 *    one. It is a good idea to perform these checks, as they can be
 *    done quickly compared to the vertex connectivity calculation itself.
 *    They were suggested by Peter McMahan, thanks Peter.
 * \return Error code.
 *
 * Time complexity: O(|V|^4), |V| is the number of vertices. In
 * practice it is more like O(|V|^2), see \ref igraph_maxflow_value().
 *
 * \sa \ref igraph_vertex_connectivity(), \ref igraph_adhesion(),
 * \ref igraph_maxflow_value().
 */

int igraph_cohesion(const igraph_t *graph, igraph_integer_t *res,
                    igraph_bool_t checks) {

    IGRAPH_CHECK(igraph_vertex_connectivity(graph, res, checks));
    return 0;
}

/**
 * \function igraph_gomory_hu_tree
 * \brief Gomory-Hu tree of a graph.
 *
 * </para><para>
 * The Gomory-Hu tree is a concise representation of the value of all the
 * maximum flows (or minimum cuts) in a graph. The vertices of the tree
 * correspond exactly to the vertices of the original graph in the same order.
 * Edges of the Gomory-Hu tree are annotated by flow values.  The value of
 * the maximum flow (or minimum cut) between an arbitrary (u,v) vertex
 * pair in the original graph is then given by the minimum flow value (i.e.
 * edge annotation) along the shortest path between u and v in the
 * Gomory-Hu tree.
 *
 * </para><para>This implementation uses Gusfield's algorithm to construct the
 * Gomory-Hu tree. See the following paper for more details:
 *
 * </para><para>
 * Gusfield D: Very simple methods for all pairs network flow analysis. SIAM J
 * Comput 19(1):143-155, 1990.
 *
 * \param graph The input graph.
 * \param tree  Pointer to an uninitialized graph; the result will be
 *              stored here.
 * \param flows Pointer to an uninitialized vector; the flow values
 *              corresponding to each edge in the Gomory-Hu tree will
 *              be returned here. You may pass a NULL pointer here if you are
 *              not interested in the flow values.
 * \param capacity Vector containing the capacity of the edges. If NULL, then
 *        every edge is considered to have capacity 1.0.
 * \return Error code.
 *
 * Time complexity: O(|V|^4) since it performs a max-flow calculation
 * between vertex zero and every other vertex and max-flow is
 * O(|V|^3).
 *
 * \sa \ref igraph_maxflow()
 */
int igraph_gomory_hu_tree(const igraph_t *graph, igraph_t *tree,
                          igraph_vector_t *flows, const igraph_vector_t *capacity) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t source, target, mid, i, n;
    igraph_vector_t neighbors;
    igraph_vector_t flow_values;
    igraph_vector_t partition;
    igraph_vector_t partition2;
    igraph_real_t flow_value;

    if (igraph_is_directed(graph)) {
        IGRAPH_ERROR("Gomory-Hu tree can only be calculated for undirected graphs",
                     IGRAPH_EINVAL);
    }

    /* Allocate memory */
    IGRAPH_VECTOR_INIT_FINALLY(&neighbors, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&flow_values, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&partition, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&partition2, 0);

    /* Initialize the tree: every edge points to node 0 */
    /* Actually, this is done implicitly since both 'neighbors' and 'flow_values' are
     * initialized to zero already */

    /* For each source vertex except vertex zero... */
    for (source = 1; source < no_of_nodes; source++) {
        IGRAPH_ALLOW_INTERRUPTION();
        IGRAPH_PROGRESS("Gomory-Hu tree", (100.0 * (source - 1)) / (no_of_nodes - 1), 0);

        /* Find its current neighbor in the tree */
        target = VECTOR(neighbors)[(long int)source];

        /* Find the maximum flow between source and target */
        IGRAPH_CHECK(igraph_maxflow(graph, &flow_value, 0, 0, &partition, &partition2,
                                    source, target, capacity, 0));

        /* Store the maximum flow and determine which side each node is on */
        VECTOR(flow_values)[(long int)source] = flow_value;

        /* Update the tree */
        /* igraph_maxflow() guarantees that the source vertex will be in &partition
         * and not in &partition2 */
        n = igraph_vector_size(&partition);
        for (i = 0; i < n; i++) {
            mid = VECTOR(partition)[i];
            if (mid > source && VECTOR(neighbors)[(long int)mid] == target) {
                VECTOR(neighbors)[(long int)mid] = source;
            }
        }
    }

    IGRAPH_PROGRESS("Gomory-Hu tree", 100.0, 0);

    /* Re-use the 'partition' vector as an edge list now */
    IGRAPH_CHECK(igraph_vector_resize(&partition, 2 * (no_of_nodes - 1)));
    for (i = 1, mid = 0; i < no_of_nodes; i++, mid += 2) {
        VECTOR(partition)[(long int)mid]   = i;
        VECTOR(partition)[(long int)mid + 1] = VECTOR(neighbors)[(long int)i];
    }

    /* Create the tree graph; we use igraph_subgraph_edges here to keep the
     * graph and vertex attributes */
    IGRAPH_CHECK(igraph_subgraph_edges(graph, tree, igraph_ess_none(), 0));
    IGRAPH_CHECK(igraph_add_edges(tree, &partition, 0));

    /* Free the allocated memory */
    igraph_vector_destroy(&partition2);
    igraph_vector_destroy(&partition);
    igraph_vector_destroy(&neighbors);
    IGRAPH_FINALLY_CLEAN(3);

    /* Return the flow values to the caller */
    if (flows != 0) {
        IGRAPH_CHECK(igraph_vector_update(flows, &flow_values));
        if (no_of_nodes > 0) {
            igraph_vector_remove(flows, 0);
        }
    }

    /* Free the remaining allocated memory */
    igraph_vector_destroy(&flow_values);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
