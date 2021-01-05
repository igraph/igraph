/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_structural.h"
#include "igraph_transitivity.h"
#include "igraph_paths.h"
#include "core/math.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_centrality.h"
#include "igraph_components.h"
#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_dqueue.h"
#include "igraph_attributes.h"
#include "igraph_neighborhood.h"
#include "igraph_stack.h"
#include "igraph_topology.h"
#include "igraph_qsort.h"
#include "igraph_error.h"
#include "config.h"

#include <string.h>
#include <limits.h>

#include "core/fixed_vectorlist.h"
#include "core/indheap.h"
#include "core/interruption.h"
#include "properties/properties_internal.h"

/**
 * \section about_structural
 *
 * <para>These functions usually calculate some structural property
 * of a graph, like its diameter, the degree of the nodes, etc.</para>
 */


/**
 * \function igraph_path_length_hist
 * Create a histogram of all shortest path lengths.
 *
 * This function calculates a histogram, by calculating the
 * shortest path length between each pair of vertices. For directed
 * graphs both directions might be considered and then every pair of vertices
 * appears twice in the histogram.
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the result is stored
 *     here. The first (i.e. zeroth) element contains the number of
 *     shortest paths of length 1, etc. The supplied vector is resized
 *     as needed.
 * \param unconnected Pointer to a real number, the number of
 *     pairs for which the second vertex is not reachable from the
 *     first is stored here.
 * \param directed Whether to consider directed paths in a directed
 *     graph (if not zero). This argument is ignored for undirected
 *     graphs.
 * \return Error code.
 *
 * Time complexity: O(|V||E|), the number of vertices times the number
 * of edges.
 *
 * \sa \ref igraph_average_path_length() and \ref igraph_shortest_paths()
 */

int igraph_path_length_hist(const igraph_t *graph, igraph_vector_t *res,
                            igraph_real_t *unconnected, igraph_bool_t directed) {

    long int no_of_nodes = igraph_vcount(graph);
    long int i, j, n;
    igraph_vector_long_t already_added;
    long int nodes_reached;

    igraph_dqueue_t q = IGRAPH_DQUEUE_NULL;
    igraph_vector_int_t *neis;
    igraph_neimode_t dirmode;
    igraph_adjlist_t allneis;
    igraph_real_t unconn = 0;
    long int ressize;

    if (directed) {
        dirmode = IGRAPH_OUT;
    } else {
        dirmode = IGRAPH_ALL;
    }

    IGRAPH_CHECK(igraph_vector_long_init(&already_added, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &already_added);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, dirmode));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

    IGRAPH_CHECK(igraph_vector_resize(res, 0));
    ressize = 0;

    for (i = 0; i < no_of_nodes; i++) {
        nodes_reached = 1;      /* itself */
        IGRAPH_CHECK(igraph_dqueue_push(&q, i));
        IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
        VECTOR(already_added)[i] = i + 1;

        IGRAPH_PROGRESS("Path-hist: ", 100.0 * i / no_of_nodes, NULL);

        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);
            long int actdist = (long int) igraph_dqueue_pop(&q);

            neis = igraph_adjlist_get(&allneis, actnode);
            n = igraph_vector_int_size(neis);
            for (j = 0; j < n; j++) {
                long int neighbor = (long int) VECTOR(*neis)[j];
                if (VECTOR(already_added)[neighbor] == i + 1) {
                    continue;
                }
                VECTOR(already_added)[neighbor] = i + 1;
                nodes_reached++;
                if (actdist + 1 > ressize) {
                    IGRAPH_CHECK(igraph_vector_resize(res, actdist + 1));
                    for (; ressize < actdist + 1; ressize++) {
                        VECTOR(*res)[ressize] = 0;
                    }
                }
                VECTOR(*res)[actdist] += 1;

                IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
            }
        } /* while !igraph_dqueue_empty */

        unconn += (no_of_nodes - nodes_reached);

    } /* for i<no_of_nodes */

    IGRAPH_PROGRESS("Path-hist: ", 100.0, NULL);

    /* count every pair only once for an undirected graph */
    if (!directed || !igraph_is_directed(graph)) {
        for (i = 0; i < ressize; i++) {
            VECTOR(*res)[i] /= 2;
        }
        unconn /= 2;
    }

    igraph_vector_long_destroy(&already_added);
    igraph_dqueue_destroy(&q);
    igraph_adjlist_destroy(&allneis);
    IGRAPH_FINALLY_CLEAN(3);

    if (unconnected) {
        *unconnected = unconn;
    }

    return 0;
}

/**
 * \ingroup structural
 * \function igraph_subcomponent
 * \brief The vertices in the same component as a given vertex.
 *
 * \param graph The graph object.
 * \param res The result, vector with the ids of the vertices in the
 *        same component.
 * \param vertex The id of the vertex of which the component is
 *        searched.
 * \param mode Type of the component for directed graphs, possible
 *        values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the set of vertices reachable \em from the
 *          \p vertex,
 *        \cli IGRAPH_IN
 *          the set of vertices from which the
 *          \p vertex is reachable.
 *        \cli IGRAPH_ALL
 *          the graph is considered as an
 *          undirected graph. Note that this is \em not the same
 *          as the union of the previous two.
 *        \endclist
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *          not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           \p vertex is an invalid vertex id
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument passed.
 *        \endclist
 *
 * Time complexity: O(|V|+|E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph.
 *
 * \sa \ref igraph_induced_subgraph() if you want a graph object consisting only
 * a given set of vertices and the edges between them.
 */

int igraph_subcomponent(const igraph_t *graph, igraph_vector_t *res, igraph_real_t vertex,
                        igraph_neimode_t mode) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_t q = IGRAPH_DQUEUE_NULL;
    char *already_added;
    long int i, vsize;
    igraph_vector_t tmp = IGRAPH_VECTOR_NULL;

    if (!IGRAPH_FINITE(vertex) || vertex < 0 || vertex >= no_of_nodes) {
        IGRAPH_ERROR("subcomponent failed", IGRAPH_EINVVID);
    }
    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("invalid mode argument", IGRAPH_EINVMODE);
    }

    already_added = igraph_Calloc(no_of_nodes, char);
    if (already_added == 0) {
        IGRAPH_ERROR("subcomponent failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, already_added);

    igraph_vector_clear(res);

    IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_dqueue_push(&q, vertex));
    IGRAPH_CHECK(igraph_vector_push_back(res, vertex));
    already_added[(long int)vertex] = 1;

    while (!igraph_dqueue_empty(&q)) {
        long int actnode = (long int) igraph_dqueue_pop(&q);

        IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_CHECK(igraph_neighbors(graph, &tmp, (igraph_integer_t) actnode,
                                      mode));
        vsize = igraph_vector_size(&tmp);
        for (i = 0; i < vsize; i++) {
            long int neighbor = (long int) VECTOR(tmp)[i];

            if (already_added[neighbor]) {
                continue;
            }
            already_added[neighbor] = 1;
            IGRAPH_CHECK(igraph_vector_push_back(res, neighbor));
            IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
        }
    }

    igraph_dqueue_destroy(&q);
    igraph_vector_destroy(&tmp);
    igraph_Free(already_added);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

/**
 * Subgraph creation, old version: it copies the graph and then deletes
 * unneeded vertices.
 */
int igraph_i_subgraph_copy_and_delete(const igraph_t *graph, igraph_t *res,
                                      const igraph_vs_t vids,
                                      igraph_vector_t *map,
                                      igraph_vector_t *invmap) {
    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t delete = IGRAPH_VECTOR_NULL;
    char *remain;
    long int i;
    igraph_vit_t vit;

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    IGRAPH_VECTOR_INIT_FINALLY(&delete, 0);
    remain = igraph_Calloc(no_of_nodes, char);
    if (remain == 0) {
        IGRAPH_ERROR("subgraph failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, remain);
    IGRAPH_CHECK(igraph_vector_reserve(&delete, no_of_nodes - IGRAPH_VIT_SIZE(vit)));

    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        remain[ (long int) IGRAPH_VIT_GET(vit) ] = 1;
    }

    for (i = 0; i < no_of_nodes; i++) {

        IGRAPH_ALLOW_INTERRUPTION();

        if (remain[i] == 0) {
            IGRAPH_CHECK(igraph_vector_push_back(&delete, i));
        }
    }

    igraph_Free(remain);
    IGRAPH_FINALLY_CLEAN(1);

    /* must set res->attr to 0 before calling igraph_copy */
    res->attr = 0;         /* Why is this needed? TODO */
    IGRAPH_CHECK(igraph_copy(res, graph));
    IGRAPH_FINALLY(igraph_destroy, res);
    IGRAPH_CHECK(igraph_delete_vertices_idx(res, igraph_vss_vector(&delete),
                                            map, invmap));

    igraph_vector_destroy(&delete);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(3);
    return 0;
}

/**
 * Subgraph creation, new version: creates the new graph instead of
 * copying the old one.
 */
int igraph_i_subgraph_create_from_scratch(const igraph_t *graph,
        igraph_t *res,
        const igraph_vs_t vids,
        igraph_vector_t *map,
        igraph_vector_t *invmap) {
    igraph_bool_t directed = igraph_is_directed(graph);
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_new_nodes = 0;
    long int i, j, n;
    long int to;
    igraph_integer_t eid;
    igraph_vector_t vids_old2new, vids_new2old;
    igraph_vector_t eids_new2old;
    igraph_vector_t nei_edges;
    igraph_vector_t new_edges;
    igraph_vit_t vit;
    igraph_vector_t *my_vids_old2new = &vids_old2new,
                     *my_vids_new2old = &vids_new2old;

    /* The order of initialization is important here, they will be destroyed in the
     * opposite order */
    IGRAPH_VECTOR_INIT_FINALLY(&eids_new2old, 0);
    if (invmap) {
        my_vids_new2old = invmap;
        igraph_vector_clear(my_vids_new2old);
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(&vids_new2old, 0);
    }
    IGRAPH_VECTOR_INIT_FINALLY(&new_edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&nei_edges, 0);
    if (map) {
        my_vids_old2new = map;
        IGRAPH_CHECK(igraph_vector_resize(map, no_of_nodes));
        igraph_vector_null(map);
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(&vids_old2new, no_of_nodes);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    /* Calculate the mapping from the old node IDs to the new ones. The other
     * igraph_simplify implementation in igraph_i_simplify_copy_and_delete
     * ensures that the order of vertex IDs is kept during remapping (i.e.
     * if the old ID of vertex A is less than the old ID of vertex B, then
     * the same will also be true for the new IDs). To ensure compatibility
     * with the other implementation, we have to fetch the vertex IDs into
     * a vector first and then sort it. We temporarily use new_edges for that.
     */
    IGRAPH_CHECK(igraph_vit_as_vector(&vit, &nei_edges));
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_vector_sort(&nei_edges);
    n = igraph_vector_size(&nei_edges);
    for (i = 0; i < n; i++) {
        long int vid = (long int) VECTOR(nei_edges)[i];
        if (VECTOR(*my_vids_old2new)[vid] == 0) {
            IGRAPH_CHECK(igraph_vector_push_back(my_vids_new2old, vid));
            no_of_new_nodes++;
            VECTOR(*my_vids_old2new)[vid] = no_of_new_nodes;
        }
    }

    /* Create the new edge list */
    for (i = 0; i < no_of_new_nodes; i++) {
        long int old_vid = (long int) VECTOR(*my_vids_new2old)[i];
        long int new_vid = i;
        igraph_bool_t skip_loop_edge;

        IGRAPH_CHECK(igraph_incident(graph, &nei_edges, old_vid, IGRAPH_OUT));
        n = igraph_vector_size(&nei_edges);

        if (directed) {
            /* directed graph; this is easier */
            for (j = 0; j < n; j++) {
                eid = (igraph_integer_t) VECTOR(nei_edges)[j];

                to = (long int) VECTOR(*my_vids_old2new)[ (long int)IGRAPH_TO(graph, eid) ];
                if (!to) {
                    continue;
                }

                IGRAPH_CHECK(igraph_vector_push_back(&new_edges, new_vid));
                IGRAPH_CHECK(igraph_vector_push_back(&new_edges, to - 1));
                IGRAPH_CHECK(igraph_vector_push_back(&eids_new2old, eid));
            }
        } else {
            /* undirected graph. We need to be careful with loop edges as each
             * loop edge will appear twice. We use a boolean flag to skip every
             * second loop edge */
            skip_loop_edge = 0;
            for (j = 0; j < n; j++) {
                eid = (igraph_integer_t) VECTOR(nei_edges)[j];

                if (IGRAPH_FROM(graph, eid) != old_vid) {
                    /* avoid processing edges twice */
                    continue;
                }

                to = (long int) VECTOR(*my_vids_old2new)[ (long int)IGRAPH_TO(graph, eid) ];
                if (!to) {
                    continue;
                }
                to -= 1;

                if (new_vid == to) {
                    /* this is a loop edge; check whether we need to skip it */
                    skip_loop_edge = !skip_loop_edge;
                    if (skip_loop_edge) {
                        continue;
                    }
                }

                IGRAPH_CHECK(igraph_vector_push_back(&new_edges, new_vid));
                IGRAPH_CHECK(igraph_vector_push_back(&new_edges, to));
                IGRAPH_CHECK(igraph_vector_push_back(&eids_new2old, eid));
            }
        }
    }

    /* Get rid of some vectors that are not needed anymore */
    if (!map) {
        igraph_vector_destroy(&vids_old2new);
        IGRAPH_FINALLY_CLEAN(1);
    }
    igraph_vector_destroy(&nei_edges);
    IGRAPH_FINALLY_CLEAN(1);

    /* Create the new graph */
    IGRAPH_CHECK(igraph_create(res, &new_edges, (igraph_integer_t)
                               no_of_new_nodes, directed));
    IGRAPH_I_ATTRIBUTE_DESTROY(res);

    /* Now we can also get rid of the new_edges vector */
    igraph_vector_destroy(&new_edges);
    IGRAPH_FINALLY_CLEAN(1);

    /* Make sure that the newly created graph is destroyed if something happens from
     * now on */
    IGRAPH_FINALLY(igraph_destroy, res);

    /* Copy the graph attributes */
    IGRAPH_CHECK(igraph_i_attribute_copy(res, graph,
                                         /* ga = */ 1, /* va = */ 0, /* ea = */ 0));

    /* Copy the vertex attributes */
    IGRAPH_CHECK(igraph_i_attribute_permute_vertices(graph, res,
                 my_vids_new2old));

    /* Copy the edge attributes */
    IGRAPH_CHECK(igraph_i_attribute_permute_edges(graph, res, &eids_new2old));

    if (!invmap) {
        igraph_vector_destroy(my_vids_new2old);
        IGRAPH_FINALLY_CLEAN(1);
    }
    igraph_vector_destroy(&eids_new2old);
    IGRAPH_FINALLY_CLEAN(2);   /* 1 + 1 since we don't need to destroy res */

    return 0;
}

/**
 * \ingroup structural
 * \function igraph_induced_subgraph
 * \brief Creates a subgraph induced by the specified vertices.
 *
 * </para><para>
 * This function collects the specified vertices and all edges between
 * them to a new graph.
 * As the vertex ids in a graph always start with zero, this function
 * very likely needs to reassign ids to the vertices.
 * \param graph The graph object.
 * \param res The subgraph, another graph object will be stored here,
 *        do \em not initialize this object before calling this
 *        function, and call \ref igraph_destroy() on it if you don't need
 *        it any more.
 * \param vids A vertex selector describing which vertices to keep.
 * \param impl This parameter selects which implementation should we
 *        use when constructing the new graph. Basically there are two
 *        possibilities: \c IGRAPH_SUBGRAPH_COPY_AND_DELETE copies the
 *        existing graph and deletes the vertices that are not needed
 *        in the new graph, while \c IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH
 *        constructs the new graph from scratch without copying the old
 *        one. The latter is more efficient if you are extracting a
 *        relatively small subpart of a very large graph, while the
 *        former is better if you want to extract a subgraph whose size
 *        is comparable to the size of the whole graph. There is a third
 *        possibility: \c IGRAPH_SUBGRAPH_AUTO will select one of the
 *        two methods automatically based on the ratio of the number
 *        of vertices in the new and the old graph.
 *
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data.
 *         \c IGRAPH_EINVVID, invalid vertex id in
 *         \p vids.
 *
 * Time complexity: O(|V|+|E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the original graph.
 *
 * \sa \ref igraph_delete_vertices() to delete the specified set of
 * vertices from a graph, the opposite of this function.
 */
int igraph_induced_subgraph(const igraph_t *graph, igraph_t *res,
                            const igraph_vs_t vids, igraph_subgraph_implementation_t impl) {
    return igraph_induced_subgraph_map(graph, res, vids, impl, /* map= */ 0,
                                       /* invmap= */ 0);
}

int igraph_i_induced_subgraph_suggest_implementation(
    const igraph_t *graph, const igraph_vs_t vids,
    igraph_subgraph_implementation_t *result) {
    double ratio;
    igraph_integer_t num_vs;

    if (igraph_vs_is_all(&vids)) {
        ratio = 1.0;
    } else {
        IGRAPH_CHECK(igraph_vs_size(graph, &vids, &num_vs));
        ratio = (igraph_real_t) num_vs / igraph_vcount(graph);
    }

    /* TODO: needs benchmarking; threshold was chosen totally arbitrarily */
    if (ratio > 0.5) {
        *result = IGRAPH_SUBGRAPH_COPY_AND_DELETE;
    } else {
        *result = IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH;
    }

    return 0;
}

int igraph_induced_subgraph_map(const igraph_t *graph, igraph_t *res,
                                const igraph_vs_t vids,
                                igraph_subgraph_implementation_t impl,
                                igraph_vector_t *map,
                                igraph_vector_t *invmap) {

    if (impl == IGRAPH_SUBGRAPH_AUTO) {
        IGRAPH_CHECK(igraph_i_induced_subgraph_suggest_implementation(graph, vids, &impl));
    }

    switch (impl) {
    case IGRAPH_SUBGRAPH_COPY_AND_DELETE:
        return igraph_i_subgraph_copy_and_delete(graph, res, vids, map, invmap);

    case IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH:
        return igraph_i_subgraph_create_from_scratch(graph, res, vids, map,
                invmap);

    default:
        IGRAPH_ERROR("unknown subgraph implementation type", IGRAPH_EINVAL);
    }
}

/**
 * \ingroup structural
 * \function igraph_subgraph_edges
 * \brief Creates a subgraph with the specified edges and their endpoints.
 *
 * </para><para>
 * This function collects the specified edges and their endpoints to a new
 * graph.
 * As the vertex ids in a graph always start with zero, this function
 * very likely needs to reassign ids to the vertices.
 * \param graph The graph object.
 * \param res The subgraph, another graph object will be stored here,
 *        do \em not initialize this object before calling this
 *        function, and call \ref igraph_destroy() on it if you don't need
 *        it any more.
 * \param eids An edge selector describing which edges to keep.
 * \param delete_vertices Whether to delete the vertices not incident on any
 *        of the specified edges as well. If \c FALSE, the number of vertices
 *        in the result graph will always be equal to the number of vertices
 *        in the input graph.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data.
 *         \c IGRAPH_EINVEID, invalid edge id in
 *         \p eids.
 *
 * Time complexity: O(|V|+|E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the original graph.
 *
 * \sa \ref igraph_delete_edges() to delete the specified set of
 * edges from a graph, the opposite of this function.
 */

int igraph_subgraph_edges(const igraph_t *graph, igraph_t *res,
                          const igraph_es_t eids, igraph_bool_t delete_vertices) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_vector_t delete = IGRAPH_VECTOR_NULL;
    char *vremain, *eremain;
    long int i;
    igraph_eit_t eit;

    IGRAPH_CHECK(igraph_eit_create(graph, eids, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    IGRAPH_VECTOR_INIT_FINALLY(&delete, 0);
    vremain = igraph_Calloc(no_of_nodes, char);
    if (vremain == 0) {
        IGRAPH_ERROR("subgraph_edges failed", IGRAPH_ENOMEM);
    }
    eremain = igraph_Calloc(no_of_edges, char);
    if (eremain == 0) {
        IGRAPH_ERROR("subgraph_edges failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, vremain);
    IGRAPH_FINALLY(igraph_free, eremain);
    IGRAPH_CHECK(igraph_vector_reserve(&delete, no_of_edges - IGRAPH_EIT_SIZE(eit)));

    /* Collect the vertex and edge IDs that will remain */
    for (IGRAPH_EIT_RESET(eit); !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
        igraph_integer_t from, to;
        long int eid = (long int) IGRAPH_EIT_GET(eit);
        IGRAPH_CHECK(igraph_edge(graph, (igraph_integer_t) eid, &from, &to));
        eremain[eid] = vremain[(long int)from] = vremain[(long int)to] = 1;
    }

    /* Collect the edge IDs to be deleted */
    for (i = 0; i < no_of_edges; i++) {
        IGRAPH_ALLOW_INTERRUPTION();
        if (eremain[i] == 0) {
            IGRAPH_CHECK(igraph_vector_push_back(&delete, i));
        }
    }

    igraph_Free(eremain);
    IGRAPH_FINALLY_CLEAN(1);

    /* Delete the unnecessary edges */
    /* must set res->attr to 0 before calling igraph_copy */
    res->attr = 0;         /* Why is this needed? TODO */
    IGRAPH_CHECK(igraph_copy(res, graph));
    IGRAPH_FINALLY(igraph_destroy, res);
    IGRAPH_CHECK(igraph_delete_edges(res, igraph_ess_vector(&delete)));

    if (delete_vertices) {
        /* Collect the vertex IDs to be deleted */
        igraph_vector_clear(&delete);
        for (i = 0; i < no_of_nodes; i++) {
            IGRAPH_ALLOW_INTERRUPTION();
            if (vremain[i] == 0) {
                IGRAPH_CHECK(igraph_vector_push_back(&delete, i));
            }
        }
    }

    igraph_Free(vremain);
    IGRAPH_FINALLY_CLEAN(1);

    /* Delete the unnecessary vertices */
    if (delete_vertices) {
        IGRAPH_CHECK(igraph_delete_vertices(res, igraph_vss_vector(&delete)));
    }

    igraph_vector_destroy(&delete);
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(3);
    return 0;
}

void igraph_i_simplify_free(igraph_vector_ptr_t *p);

void igraph_i_simplify_free(igraph_vector_ptr_t *p) {
    long int i, n = igraph_vector_ptr_size(p);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(*p)[i];
        if (v) {
            igraph_vector_destroy(v);
        }
    }
    igraph_vector_ptr_destroy(p);
}

/**
 * \ingroup structural
 * \function igraph_simplify
 * \brief Removes loop and/or multiple edges from the graph.
 *
 * \param graph The graph object.
 * \param multiple Logical, if true, multiple edges will be removed.
 * \param loops Logical, if true, loops (self edges) will be removed.
 * \param edge_comb What to do with the edge attributes. See the igraph
 *        manual section about attributes for details.
 * \return Error code:
 *    \c IGRAPH_ENOMEM if we are out of memory.
 *
 * Time complexity: O(|V|+|E|).
 *
 * \example examples/simple/igraph_simplify.c
 */

int igraph_simplify(igraph_t *graph, igraph_bool_t multiple,
                    igraph_bool_t loops,
                    const igraph_attribute_combination_t *edge_comb) {

    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    long int edge;
    igraph_bool_t attr = edge_comb && igraph_has_attribute_table();
    long int from, to, pfrom = -1, pto = -2;
    igraph_t res;
    igraph_es_t es;
    igraph_eit_t eit;
    igraph_vector_t mergeinto;
    long int actedge;

    if (!multiple && !loops)
        /* nothing to do */
    {
        return IGRAPH_SUCCESS;
    }

    if (!multiple) {
        /* removing loop edges only, this is simple. No need to combine anything
         * and the whole process can be done in-place */
        IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
        IGRAPH_CHECK(igraph_es_all(&es, IGRAPH_EDGEORDER_ID));
        IGRAPH_FINALLY(igraph_es_destroy, &es);
        IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
        IGRAPH_FINALLY(igraph_eit_destroy, &eit);

        while (!IGRAPH_EIT_END(eit)) {
            edge = IGRAPH_EIT_GET(eit);
            from = IGRAPH_FROM(graph, edge);
            to = IGRAPH_TO(graph, edge);
            if (from == to) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, edge));
            }
            IGRAPH_EIT_NEXT(eit);
        }

        igraph_eit_destroy(&eit);
        igraph_es_destroy(&es);
        IGRAPH_FINALLY_CLEAN(2);

        if (igraph_vector_size(&edges) > 0) {
            IGRAPH_CHECK(igraph_delete_edges(graph, igraph_ess_vector(&edges)));
        }

        igraph_vector_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);

        return IGRAPH_SUCCESS;
    }

    if (attr) {
        IGRAPH_VECTOR_INIT_FINALLY(&mergeinto, no_of_edges);
    }
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 2));

    IGRAPH_CHECK(igraph_es_all(&es, IGRAPH_EDGEORDER_FROM));
    IGRAPH_FINALLY(igraph_es_destroy, &es);
    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    for (actedge = -1; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
        edge = IGRAPH_EIT_GET(eit);
        from = IGRAPH_FROM(graph, edge);
        to = IGRAPH_TO(graph, edge);

        if (loops && from == to) {
            /* Loop edge to be removed */
            if (attr) {
                VECTOR(mergeinto)[edge] = -1;
            }
        } else if (multiple && from == pfrom && to == pto) {
            /* Multiple edge to be contracted */
            if (attr) {
                VECTOR(mergeinto)[edge] = actedge;
            }
        } else {
            /* Edge to be kept */
            igraph_vector_push_back(&edges, from);
            igraph_vector_push_back(&edges, to);
            if (attr) {
                actedge++;
                VECTOR(mergeinto)[edge] = actedge;
            }
        }
        pfrom = from; pto = to;
    }

    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(&res, &edges, (igraph_integer_t) no_of_nodes,
                               igraph_is_directed(graph)));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_FINALLY(igraph_destroy, &res);

    IGRAPH_I_ATTRIBUTE_DESTROY(&res);
    IGRAPH_I_ATTRIBUTE_COPY(&res, graph, /*graph=*/ 1,
                            /*vertex=*/ 1, /*edge=*/ 0);

    if (attr) {
        igraph_fixed_vectorlist_t vl;
        IGRAPH_CHECK(igraph_fixed_vectorlist_convert(&vl, &mergeinto,
                     actedge + 1));
        IGRAPH_FINALLY(igraph_fixed_vectorlist_destroy, &vl);

        IGRAPH_CHECK(igraph_i_attribute_combine_edges(graph, &res, &vl.v,
                     edge_comb));

        igraph_fixed_vectorlist_destroy(&vl);
        igraph_vector_destroy(&mergeinto);
        IGRAPH_FINALLY_CLEAN(2);
    }

    IGRAPH_FINALLY_CLEAN(1);
    igraph_destroy(graph);
    *graph = res;

    return 0;
}

/**
 * \ingroup structural
 * \function igraph_reciprocity
 * \brief Calculates the reciprocity of a directed graph.
 *
 * </para><para>
 * The measure of reciprocity defines the proportion of mutual
 * connections, in a directed graph. It is most commonly defined as
 * the probability that the opposite counterpart of a directed edge is
 * also included in the graph. In adjacency matrix notation:
 * <code>sum(i, j, (A.*A')ij) / sum(i, j, Aij)</code>, where
 * <code>A.*A'</code> is the element-wise product of matrix
 * <code>A</code> and its transpose. This measure is
 * calculated if the \p mode argument is \c
 * IGRAPH_RECIPROCITY_DEFAULT.
 *
 * </para><para>
 * Prior to igraph version 0.6, another measure was implemented,
 * defined as the probability of mutual connection between a vertex
 * pair if we know that there is a (possibly non-mutual) connection
 * between them. In other words, (unordered) vertex pairs are
 * classified into three groups: (1) disconnected, (2)
 * non-reciprocally connected, (3) reciprocally connected.
 * The result is the size of group (3), divided by the sum of group
 * sizes (2)+(3). This measure is calculated if \p mode is \c
 * IGRAPH_RECIPROCITY_RATIO.
 *
 * \param graph The graph object.
 * \param res Pointer to an \c igraph_real_t which will contain the result.
 * \param ignore_loops Whether to ignore loop edges.
 * \param mode Type of reciprocity to calculate, possible values are
 *    \c IGRAPH_RECIPROCITY_DEFAULT and \c IGRAPH_RECIPROCITY_RATIO,
 *    please see their description above.
 * \return Error code:
 *         \c IGRAPH_EINVAL: graph has no edges
 *         \c IGRAPH_ENOMEM: not enough memory for
 *         temporary data.
 *
 * Time complexity: O(|V|+|E|), |V| is the number of vertices,
 * |E| is the number of edges.
 *
 * \example examples/simple/igraph_reciprocity.c
 */

int igraph_reciprocity(const igraph_t *graph, igraph_real_t *res,
                       igraph_bool_t ignore_loops,
                       igraph_reciprocity_t mode) {

    igraph_integer_t nonrec = 0, rec = 0, loops = 0;
    igraph_vector_t inneis, outneis;
    long int i;
    long int no_of_nodes = igraph_vcount(graph);

    if (mode != IGRAPH_RECIPROCITY_DEFAULT &&
        mode != IGRAPH_RECIPROCITY_RATIO) {
        IGRAPH_ERROR("Invalid reciprocity type", IGRAPH_EINVAL);
    }

    /* THIS IS AN EXIT HERE !!!!!!!!!!!!!! */
    if (!igraph_is_directed(graph)) {
        *res = 1.0;
        return 0;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&inneis, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&outneis, 0);

    for (i = 0; i < no_of_nodes; i++) {
        long int ip, op;
        igraph_neighbors(graph, &inneis, (igraph_integer_t) i, IGRAPH_IN);
        igraph_neighbors(graph, &outneis, (igraph_integer_t) i, IGRAPH_OUT);

        ip = op = 0;
        while (ip < igraph_vector_size(&inneis) &&
               op < igraph_vector_size(&outneis)) {
            if (VECTOR(inneis)[ip] < VECTOR(outneis)[op]) {
                nonrec += 1;
                ip++;
            } else if (VECTOR(inneis)[ip] > VECTOR(outneis)[op]) {
                nonrec += 1;
                op++;
            } else {

                /* loop edge? */
                if (VECTOR(inneis)[ip] == i) {
                    loops += 1;
                    if (!ignore_loops) {
                        rec += 1;
                    }
                } else {
                    rec += 1;
                }

                ip++;
                op++;
            }
        }
        nonrec += (igraph_vector_size(&inneis) - ip) +
                  (igraph_vector_size(&outneis) - op);
    }

    if (mode == IGRAPH_RECIPROCITY_DEFAULT) {
        if (ignore_loops) {
            *res = (igraph_real_t) rec / (igraph_ecount(graph) - loops);
        } else {
            *res = (igraph_real_t) rec / (igraph_ecount(graph));
        }
    } else if (mode == IGRAPH_RECIPROCITY_RATIO) {
        *res = (igraph_real_t) rec / (rec + nonrec);
    }

    igraph_vector_destroy(&inneis);
    igraph_vector_destroy(&outneis);
    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

/**
 * \function igraph_constraint
 * \brief Burt's constraint scores.
 *
 * </para><para>
 * This function calculates Burt's constraint scores for the given
 * vertices, also known as structural holes.
 *
 * </para><para>
 * Burt's constraint is higher if ego has less, or mutually stronger
 * related (i.e. more redundant) contacts. Burt's measure of
 * constraint, C[i], of vertex i's ego network V[i], is defined for
 * directed and valued graphs,
 * <blockquote><para>
 * C[i] = sum( sum( (p[i,q] p[q,j])^2, q in V[i], q != i,j ), j in
 * V[], j != i)
 * </para></blockquote>
 * for a graph of order (i.e. number of vertices) N, where proportional
 * tie strengths are defined as
 * <blockquote><para>
 * p[i,j]=(a[i,j]+a[j,i]) / sum(a[i,k]+a[k,i], k in V[i], k != i),
 * </para></blockquote>
 * a[i,j] are elements of A and
 * the latter being the graph adjacency matrix. For isolated vertices,
 * constraint is undefined.
 *
 * </para><para>
 * Burt, R.S. (2004). Structural holes and good ideas. American
 * Journal of Sociology 110, 349-399.
 *
 * </para><para>
 * The first R version of this function was contributed by Jeroen
 * Bruggeman.
 * \param graph A graph object.
 * \param res Pointer to an initialized vector, the result will be
 *        stored here. The vector will be resized to have the
 *        appropriate size for holding the result.
 * \param vids Vertex selector containing the vertices for which the
 *        constraint should be calculated.
 * \param weights Vector giving the weights of the edges. If it is
 *        \c NULL then each edge is supposed to have the same weight.
 * \return Error code.
 *
 * Time complexity: O(|V|+E|+n*d^2), n is the number of vertices for
 * which the constraint is calculated and d is the average degree, |V|
 * is the number of vertices, |E| the number of edges in the
 * graph. If the weights argument is \c NULL then the time complexity
 * is O(|V|+n*d^2).
 */

int igraph_constraint(const igraph_t *graph, igraph_vector_t *res,
                      igraph_vs_t vids, const igraph_vector_t *weights) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_vit_t vit;
    long int nodes_to_calc;
    long int a, b, c, i, j, q, vsize, vsize2;
    igraph_integer_t edge, from, to, edge2;

    igraph_vector_t contrib;
    igraph_vector_t degree;
    igraph_vector_t ineis_in, ineis_out, jneis_in, jneis_out;

    if (weights != 0 && igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid length of weight vector", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&contrib, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&ineis_in, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&ineis_out, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&jneis_in, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&jneis_out, 0);

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    nodes_to_calc = IGRAPH_VIT_SIZE(vit);

    if (weights == 0) {
        IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
                                   IGRAPH_ALL, IGRAPH_NO_LOOPS));
    } else {
        for (a = 0; a < no_of_edges; a++) {
            igraph_edge(graph, (igraph_integer_t) a, &from, &to);
            if (from != to) {
                VECTOR(degree)[(long int) from] += VECTOR(*weights)[a];
                VECTOR(degree)[(long int) to  ] += VECTOR(*weights)[a];
            }
        }
    }

    IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
    igraph_vector_null(res);

    for (a = 0; a < nodes_to_calc; a++, IGRAPH_VIT_NEXT(vit)) {
        i = IGRAPH_VIT_GET(vit);

        /* get neighbors of i */
        IGRAPH_CHECK(igraph_incident(graph, &ineis_in, (igraph_integer_t) i,
                                     IGRAPH_IN));
        IGRAPH_CHECK(igraph_incident(graph, &ineis_out, (igraph_integer_t) i,
                                     IGRAPH_OUT));

        /* NaN for isolates */
        if (igraph_vector_size(&ineis_in) == 0 &&
            igraph_vector_size(&ineis_out) == 0) {
            VECTOR(*res)[a] = IGRAPH_NAN;
        }

        /* zero their contribution */
        vsize = igraph_vector_size(&ineis_in);
        for (b = 0; b < vsize; b++) {
            edge = (igraph_integer_t) VECTOR(ineis_in)[b];
            j = (long int) IGRAPH_OTHER(graph, edge, i);
            VECTOR(contrib)[j] = 0.0;
        }
        vsize = igraph_vector_size(&ineis_out);
        for (b = 0; b < vsize; b++) {
            edge = (igraph_integer_t) VECTOR(ineis_out)[b];
            j = (long int) IGRAPH_OTHER(graph, edge, i);
            VECTOR(contrib)[j] = 0.0;
        }

        /* add the direct contributions, in-neighbors and out-neighbors */
        vsize = igraph_vector_size(&ineis_in);
        for (b = 0; b < vsize; b++) {
            edge = (igraph_integer_t) VECTOR(ineis_in)[b];
            j = (long int) IGRAPH_OTHER(graph, edge, i);
            if (i != j) {     /* excluding loops */
                if (weights) {
                    VECTOR(contrib)[j] +=
                        VECTOR(*weights)[(long int)edge] / VECTOR(degree)[i];
                } else {
                    VECTOR(contrib)[j] += 1.0 / VECTOR(degree)[i];
                }
            }
        }
        if (igraph_is_directed(graph)) {
            vsize = igraph_vector_size(&ineis_out);
            for (b = 0; b < vsize; b++) {
                edge = (igraph_integer_t) VECTOR(ineis_out)[b];
                j = (long int) IGRAPH_OTHER(graph, edge, i);
                if (i != j) {
                    if (weights) {
                        VECTOR(contrib)[j] +=
                            VECTOR(*weights)[(long int)edge] / VECTOR(degree)[i];
                    } else {
                        VECTOR(contrib)[j] += 1.0 / VECTOR(degree)[i];
                    }
                }
            }
        }

        /* add the indirect contributions, in-in, in-out, out-in, out-out */
        vsize = igraph_vector_size(&ineis_in);
        for (b = 0; b < vsize; b++) {
            edge = (igraph_integer_t) VECTOR(ineis_in)[b];
            j = (long int) IGRAPH_OTHER(graph, edge, i);
            if (i == j) {
                continue;
            }
            IGRAPH_CHECK(igraph_incident(graph, &jneis_in, (igraph_integer_t) j,
                                         IGRAPH_IN));
            IGRAPH_CHECK(igraph_incident(graph, &jneis_out, (igraph_integer_t) j,
                                         IGRAPH_OUT));
            vsize2 = igraph_vector_size(&jneis_in);
            for (c = 0; c < vsize2; c++) {
                edge2 = (igraph_integer_t) VECTOR(jneis_in)[c];
                q = (long int) IGRAPH_OTHER(graph, edge2, j);
                if (j != q) {
                    if (weights) {
                        VECTOR(contrib)[q] +=
                            VECTOR(*weights)[(long int)edge] *
                            VECTOR(*weights)[(long int)edge2] /
                            VECTOR(degree)[i] / VECTOR(degree)[j];
                    } else {
                        VECTOR(contrib)[q] += 1 / VECTOR(degree)[i] / VECTOR(degree)[j];
                    }
                }
            }
            if (igraph_is_directed(graph)) {
                vsize2 = igraph_vector_size(&jneis_out);
                for (c = 0; c < vsize2; c++) {
                    edge2 = (igraph_integer_t) VECTOR(jneis_out)[c];
                    q = (long int) IGRAPH_OTHER(graph, edge2, j);
                    if (j != q) {
                        if (weights) {
                            VECTOR(contrib)[q] +=
                                VECTOR(*weights)[(long int)edge] *
                                VECTOR(*weights)[(long int)edge2] /
                                VECTOR(degree)[i] / VECTOR(degree)[j];
                        } else {
                            VECTOR(contrib)[q] += 1 / VECTOR(degree)[i] / VECTOR(degree)[j];
                        }
                    }
                }
            }
        }
        if (igraph_is_directed(graph)) {
            vsize = igraph_vector_size(&ineis_out);
            for (b = 0; b < vsize; b++) {
                edge = (igraph_integer_t) VECTOR(ineis_out)[b];
                j = (long int) IGRAPH_OTHER(graph, edge, i);
                if (i == j) {
                    continue;
                }
                IGRAPH_CHECK(igraph_incident(graph, &jneis_in, (igraph_integer_t) j,
                                             IGRAPH_IN));
                IGRAPH_CHECK(igraph_incident(graph, &jneis_out, (igraph_integer_t) j,
                                             IGRAPH_OUT));
                vsize2 = igraph_vector_size(&jneis_in);
                for (c = 0; c < vsize2; c++) {
                    edge2 = (igraph_integer_t) VECTOR(jneis_in)[c];
                    q = (long int) IGRAPH_OTHER(graph, edge2, j);
                    if (j != q) {
                        if (weights) {
                            VECTOR(contrib)[q] +=
                                VECTOR(*weights)[(long int)edge] *
                                VECTOR(*weights)[(long int)edge2] /
                                VECTOR(degree)[i] / VECTOR(degree)[j];
                        } else {
                            VECTOR(contrib)[q] += 1 / VECTOR(degree)[i] / VECTOR(degree)[j];
                        }
                    }
                }
                vsize2 = igraph_vector_size(&jneis_out);
                for (c = 0; c < vsize2; c++) {
                    edge2 = (igraph_integer_t) VECTOR(jneis_out)[c];
                    q = (long int) IGRAPH_OTHER(graph, edge2, j);
                    if (j != q) {
                        if (weights) {
                            VECTOR(contrib)[q] +=
                                VECTOR(*weights)[(long int)edge] *
                                VECTOR(*weights)[(long int)edge2] /
                                VECTOR(degree)[i] / VECTOR(degree)[j];
                        } else {
                            VECTOR(contrib)[q] += 1 / VECTOR(degree)[i] / VECTOR(degree)[j];
                        }
                    }
                }
            }
        }

        /* squared sum of the contributions */
        vsize = igraph_vector_size(&ineis_in);
        for (b = 0; b < vsize; b++) {
            edge = (igraph_integer_t) VECTOR(ineis_in)[b];
            j = (long int) IGRAPH_OTHER(graph, edge, i);
            if (i == j) {
                continue;
            }
            VECTOR(*res)[a] += VECTOR(contrib)[j] * VECTOR(contrib)[j];
            VECTOR(contrib)[j] = 0.0;
        }
        if (igraph_is_directed(graph)) {
            vsize =  igraph_vector_size(&ineis_out);
            for (b = 0; b < vsize; b++) {
                edge = (igraph_integer_t) VECTOR(ineis_out)[b];
                j = (long int) IGRAPH_OTHER(graph, edge, i);
                if (i == j) {
                    continue;
                }
                VECTOR(*res)[a] += VECTOR(contrib)[j] * VECTOR(contrib)[j];
                VECTOR(contrib)[j] = 0.0;
            }
        }
    }

    igraph_vit_destroy(&vit);
    igraph_vector_destroy(&jneis_out);
    igraph_vector_destroy(&jneis_in);
    igraph_vector_destroy(&ineis_out);
    igraph_vector_destroy(&ineis_in);
    igraph_vector_destroy(&degree);
    igraph_vector_destroy(&contrib);
    IGRAPH_FINALLY_CLEAN(7);

    return 0;
}

/**
 * \function igraph_maxdegree
 * \brief Calculate the maximum degree in a graph (or set of vertices).
 *
 * </para><para>
 * The largest in-, out- or total degree of the specified vertices is
 * calculated.
 * \param graph The input graph.
 * \param res Pointer to an integer (\c igraph_integer_t), the result
 *        will be stored here.
 * \param vids Vector giving the vertex IDs for which the maximum degree will
 *        be calculated.
 * \param mode Defines the type of the degree.
 *        \c IGRAPH_OUT, out-degree,
 *        \c IGRAPH_IN, in-degree,
 *        \c IGRAPH_ALL, total degree (sum of the
 *        in- and out-degree).
 *        This parameter is ignored for undirected graphs.
 * \param loops Boolean, gives whether the self-loops should be
 *        counted.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex id.
 *         \c IGRAPH_EINVMODE: invalid mode argument.
 *
 * Time complexity: O(v) if
 * loops is
 * TRUE, and
 * O(v*d)
 * otherwise. v is the number
 * vertices for which the degree will be calculated, and
 * d is their (average) degree.
 */

int igraph_maxdegree(const igraph_t *graph, igraph_integer_t *res,
                     igraph_vs_t vids, igraph_neimode_t mode,
                     igraph_bool_t loops) {

    igraph_vector_t tmp;

    IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);

    igraph_degree(graph, &tmp, vids, mode, loops);
    *res = (igraph_integer_t) igraph_vector_max(&tmp);

    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_density
 * Calculate the density of a graph.
 *
 * </para><para>The density of a graph is simply the ratio number of
 * edges and the number of possible edges. Note that density is
 * ill-defined for graphs with multiple and/or loop edges, so consider
 * calling \ref igraph_simplify() on the graph if you know that it
 * contains multiple or loop edges.
 * \param graph The input graph object.
 * \param res Pointer to a real number, the result will be stored
 *   here.
 * \param loops Logical constant, whether to include loops in the
 *   calculation. If this constant is TRUE then
 *   loop edges are thought to be possible in the graph (this does not
 *   necessarily mean that the graph really contains any loops). If
 *   this is FALSE then the result is only correct if the graph does not
 *   contain loops.
 * \return Error code.
 *
 * Time complexity: O(1).
 */

int igraph_density(const igraph_t *graph, igraph_real_t *res,
                   igraph_bool_t loops) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_real_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);

    if (no_of_nodes == 0) {
        *res = IGRAPH_NAN;
        return 0;
    }

    if (!loops) {
        if (no_of_nodes == 1) {
            *res = IGRAPH_NAN;
        } else if (directed) {
            *res = no_of_edges / no_of_nodes / (no_of_nodes - 1);
        } else {
            *res = no_of_edges / no_of_nodes * 2.0 / (no_of_nodes - 1);
        }
    } else {
        if (directed) {
            *res = no_of_edges / no_of_nodes / no_of_nodes;
        } else {
            *res = no_of_edges / no_of_nodes * 2.0 / (no_of_nodes + 1);
        }
    }

    return 0;
}

/**
 * \function igraph_neighborhood_size
 * \brief Calculates the size of the neighborhood of a given vertex.
 *
 * The neighborhood of a given order of a vertex includes all vertices
 * which are closer to the vertex than the order. I.e., order 0 is
 * always the vertex itself, order 1 is the vertex plus its immediate
 * neighbors, order 2 is order 1 plus the immediate neighbors of the
 * vertices in order 1, etc.
 *
 * </para><para> This function calculates the size of the neighborhood
 * of the given order for the given vertices.
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the result will be
 *    stored here. It will be resized as needed.
 * \param vids The vertices for which the calculation is performed.
 * \param order Integer giving the order of the neighborhood.
 * \param mode Specifies how to use the direction of the edges if a
 *   directed graph is analyzed. For \c IGRAPH_OUT only the outgoing
 *   edges are followed, so all vertices reachable from the source
 *   vertex in at most \c order steps are counted. For \c IGRAPH_IN
 *   all vertices from which the source vertex is reachable in at most
 *   \c order steps are counted. \c IGRAPH_ALL ignores the direction
 *   of the edges. This argument is ignored for undirected graphs.
 * \param mindist The minimum distance to include a vertex in the counting.
 *   If this is one, then the starting vertex is not counted. If this is
 *   two, then its neighbors are not counted, either, etc.
 * \return Error code.
 *
 * \sa \ref igraph_neighborhood() for calculating the actual neighborhood,
 * \ref igraph_neighborhood_graphs() for creating separate graphs from
 * the neighborhoods.
 *
 * Time complexity: O(n*d*o), where n is the number vertices for which
 * the calculation is performed, d is the average degree, o is the order.
 */

int igraph_neighborhood_size(const igraph_t *graph, igraph_vector_t *res,
                             igraph_vs_t vids, igraph_integer_t order,
                             igraph_neimode_t mode,
                             igraph_integer_t mindist) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_t q;
    igraph_vit_t vit;
    long int i, j;
    long int *added;
    igraph_vector_t neis;

    if (order < 0) {
        IGRAPH_ERROR("Negative order in neighborhood size", IGRAPH_EINVAL);
    }

    if (mindist < 0 || mindist > order) {
        IGRAPH_ERROR("Minimum distance should be between zero and order",
                     IGRAPH_EINVAL);
    }

    added = igraph_Calloc(no_of_nodes, long int);
    if (added == 0) {
        IGRAPH_ERROR("Cannot calculate neighborhood size", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));

    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        long int node = IGRAPH_VIT_GET(vit);
        long int size = mindist == 0 ? 1 : 0;
        added[node] = i + 1;
        igraph_dqueue_clear(&q);
        if (order > 0) {
            igraph_dqueue_push(&q, node);
            igraph_dqueue_push(&q, 0);
        }

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);
            long int actdist = (long int) igraph_dqueue_pop(&q);
            long int n;
            igraph_neighbors(graph, &neis, (igraph_integer_t) actnode, mode);
            n = igraph_vector_size(&neis);

            if (actdist < order - 1) {
                /* we add them to the q */
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
                        if (actdist + 1 >= mindist) {
                            size++;
                        }
                    }
                }
            } else {
                /* we just count them, but don't add them */
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (actdist + 1 >= mindist) {
                            size++;
                        }
                    }
                }
            }

        } /* while q not empty */

        VECTOR(*res)[i] = size;
    } /* for VIT, i */

    igraph_vector_destroy(&neis);
    igraph_vit_destroy(&vit);
    igraph_dqueue_destroy(&q);
    igraph_Free(added);
    IGRAPH_FINALLY_CLEAN(4);

    return 0;
}

/**
 * \function igraph_neighborhood
 * Calculate the neighborhood of vertices.
 *
 * The neighborhood of a given order of a vertex includes all vertices
 * which are closer to the vertex than the order. I.e., order 0 is
 * always the vertex itself, order 1 is the vertex plus its immediate
 * neighbors, order 2 is order 1 plus the immediate neighbors of the
 * vertices in order 1, etc.
 *
 * </para><para> This function calculates the vertices within the
 * neighborhood of the specified vertices.
 * \param graph The input graph.
 * \param res An initialized pointer vector. Note that the objects
 *    (pointers) in the vector will \em not be freed, but the pointer
 *    vector will be resized as needed. The result of the calculation
 *    will be stored here in \ref igraph_vector_t objects.
 * \param vids The vertices for which the calculation is performed.
 * \param order Integer giving the order of the neighborhood.
 * \param mode Specifies how to use the direction of the edges if a
 *   directed graph is analyzed. For \c IGRAPH_OUT only the outgoing
 *   edges are followed, so all vertices reachable from the source
 *   vertex in at most \p order steps are included. For \c IGRAPH_IN
 *   all vertices from which the source vertex is reachable in at most
 *   \p order steps are included. \c IGRAPH_ALL ignores the direction
 *   of the edges. This argument is ignored for undirected graphs.
 * \param mindist The minimum distance to include a vertex in the counting.
 *   If this is one, then the starting vertex is not counted. If this is
 *   two, then its neighbors are not counted, either, etc.
 * \return Error code.
 *
 * \sa \ref igraph_neighborhood_size() to calculate the size of the
 * neighborhood, \ref igraph_neighborhood_graphs() for creating
 * graphs from the neighborhoods.
 *
 * Time complexity: O(n*d*o), n is the number of vertices for which
 * the calculation is performed, d is the average degree, o is the
 * order.
 */

int igraph_neighborhood(const igraph_t *graph, igraph_vector_ptr_t *res,
                        igraph_vs_t vids, igraph_integer_t order,
                        igraph_neimode_t mode, igraph_integer_t mindist) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_t q;
    igraph_vit_t vit;
    long int i, j;
    long int *added;
    igraph_vector_t neis;
    igraph_vector_t tmp;
    igraph_vector_t *newv;

    if (order < 0) {
        IGRAPH_ERROR("Negative order in neighborhood size", IGRAPH_EINVAL);
    }

    if (mindist < 0 || mindist > order) {
        IGRAPH_ERROR("Minimum distance should be between zero and order",
                     IGRAPH_EINVAL);
    }

    added = igraph_Calloc(no_of_nodes, long int);
    if (added == 0) {
        IGRAPH_ERROR("Cannot calculate neighborhood size", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
    IGRAPH_CHECK(igraph_vector_ptr_resize(res, IGRAPH_VIT_SIZE(vit)));

    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        long int node = IGRAPH_VIT_GET(vit);
        added[node] = i + 1;
        igraph_vector_clear(&tmp);
        if (mindist == 0) {
            IGRAPH_CHECK(igraph_vector_push_back(&tmp, node));
        }
        if (order > 0) {
            igraph_dqueue_push(&q, node);
            igraph_dqueue_push(&q, 0);
        }

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);
            long int actdist = (long int) igraph_dqueue_pop(&q);
            long int n;
            igraph_neighbors(graph, &neis, (igraph_integer_t) actnode, mode);
            n = igraph_vector_size(&neis);

            if (actdist < order - 1) {
                /* we add them to the q */
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_push_back(&tmp, nei));
                        }
                    }
                }
            } else {
                /* we just count them but don't add them to q */
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_push_back(&tmp, nei));
                        }
                    }
                }
            }

        } /* while q not empty */

        newv = igraph_Calloc(1, igraph_vector_t);
        if (newv == 0) {
            IGRAPH_ERROR("Cannot calculate neighborhood", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, newv);
        IGRAPH_CHECK(igraph_vector_copy(newv, &tmp));
        VECTOR(*res)[i] = newv;
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_destroy(&tmp);
    igraph_vector_destroy(&neis);
    igraph_vit_destroy(&vit);
    igraph_dqueue_destroy(&q);
    igraph_Free(added);
    IGRAPH_FINALLY_CLEAN(5);

    return 0;
}

/**
 * \function igraph_neighborhood_graphs
 * Create graphs from the neighborhood(s) of some vertex/vertices.
 *
 * The neighborhood of a given order of a vertex includes all vertices
 * which are closer to the vertex than the order. Ie. order 0 is
 * always the vertex itself, order 1 is the vertex plus its immediate
 * neighbors, order 2 is order 1 plus the immediate neighbors of the
 * vertices in order 1, etc.
 *
 * </para><para> This function finds every vertex in the neighborhood
 * of a given parameter vertex and creates a graph from these
 * vertices.
 *
 * </para><para> The first version of this function was written by
 * Vincent Matossian, thanks Vincent.
 * \param graph The input graph.
 * \param res Pointer to a pointer vector, the result will be stored
 *   here, ie. \p res will contain pointers to \c igraph_t
 *   objects. It will be resized if needed but note that the
 *   objects in the pointer vector will not be freed.
 * \param vids The vertices for which the calculation is performed.
 * \param order Integer giving the order of the neighborhood.
 * \param mode Specifies how to use the direction of the edges if a
 *   directed graph is analyzed. For \c IGRAPH_OUT only the outgoing
 *   edges are followed, so all vertices reachable from the source
 *   vertex in at most \p order steps are counted. For \c IGRAPH_IN
 *   all vertices from which the source vertex is reachable in at most
 *   \p order steps are counted. \c IGRAPH_ALL ignores the direction
 *   of the edges. This argument is ignored for undirected graphs.
 * \param mindist The minimum distance to include a vertex in the counting.
 *   If this is one, then the starting vertex is not counted. If this is
 *   two, then its neighbors are not counted, either, etc.
 * \return Error code.
 *
 * \sa \ref igraph_neighborhood_size() for calculating the neighborhood
 * sizes only, \ref igraph_neighborhood() for calculating the
 * neighborhoods (but not creating graphs).
 *
 * Time complexity: O(n*(|V|+|E|)), where n is the number vertices for
 * which the calculation is performed, |V| and |E| are the number of
 * vertices and edges in the original input graph.
 */

int igraph_neighborhood_graphs(const igraph_t *graph, igraph_vector_ptr_t *res,
                               igraph_vs_t vids, igraph_integer_t order,
                               igraph_neimode_t mode,
                               igraph_integer_t mindist) {
    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_t q;
    igraph_vit_t vit;
    long int i, j;
    long int *added;
    igraph_vector_t neis;
    igraph_vector_t tmp;
    igraph_t *newg;

    if (order < 0) {
        IGRAPH_ERROR("Negative order in neighborhood size", IGRAPH_EINVAL);
    }

    if (mindist < 0 || mindist > order) {
        IGRAPH_ERROR("Minimum distance should be between zero and order",
                     IGRAPH_EINVAL);
    }

    added = igraph_Calloc(no_of_nodes, long int);
    if (added == 0) {
        IGRAPH_ERROR("Cannot calculate neighborhood size", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
    IGRAPH_CHECK(igraph_vector_ptr_resize(res, IGRAPH_VIT_SIZE(vit)));

    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        long int node = IGRAPH_VIT_GET(vit);
        added[node] = i + 1;
        igraph_vector_clear(&tmp);
        if (mindist == 0) {
            IGRAPH_CHECK(igraph_vector_push_back(&tmp, node));
        }
        if (order > 0) {
            igraph_dqueue_push(&q, node);
            igraph_dqueue_push(&q, 0);
        }

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);
            long int actdist = (long int) igraph_dqueue_pop(&q);
            long int n;
            igraph_neighbors(graph, &neis, (igraph_integer_t) actnode, mode);
            n = igraph_vector_size(&neis);

            if (actdist < order - 1) {
                /* we add them to the q */
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_push_back(&tmp, nei));
                        }
                    }
                }
            } else {
                /* we just count them but don't add them to q */
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_push_back(&tmp, nei));
                        }
                    }
                }
            }

        } /* while q not empty */

        newg = igraph_Calloc(1, igraph_t);
        if (newg == 0) {
            IGRAPH_ERROR("Cannot create neighborhood graph", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, newg);
        if (igraph_vector_size(&tmp) < no_of_nodes) {
            IGRAPH_CHECK(igraph_induced_subgraph(graph, newg,
                                                 igraph_vss_vector(&tmp),
                                                 IGRAPH_SUBGRAPH_AUTO));
        } else {
            IGRAPH_CHECK(igraph_copy(newg, graph));
        }
        VECTOR(*res)[i] = newg;
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_destroy(&tmp);
    igraph_vector_destroy(&neis);
    igraph_vit_destroy(&vit);
    igraph_dqueue_destroy(&q);
    igraph_Free(added);
    IGRAPH_FINALLY_CLEAN(5);

    return 0;
}

/**
 * \function igraph_topological_sorting
 * \brief Calculate a possible topological sorting of the graph.
 *
 * </para><para>
 * A topological sorting of a directed acyclic graph (DAG) is a linear ordering
 * of its vertices where each vertex comes before all nodes to which it has
 * edges. Every DAG has at least one topological sort, and may have many.
 * This function returns one possible topological sort among them. If the
 * graph is not acyclic (it has at least one cycle), an error is raised.
 *
 * \param graph The input graph.
 * \param res Pointer to a vector, the result will be stored here.
 *   It will be resized if needed.
 * \param mode Specifies how to use the direction of the edges.
 *   For \c IGRAPH_OUT, the sorting order ensures that each vertex comes
 *   before all vertices to which it has edges, so vertices with no incoming
 *   edges go first. For \c IGRAPH_IN, it is quite the opposite: each
 *   vertex comes before all vertices from which it receives edges. Vertices
 *   with no outgoing edges go first.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), where |V| and |E| are the number of
 * vertices and edges in the original input graph.
 *
 * \sa \ref igraph_is_dag() if you are only interested in whether a given
 *     graph is a DAG or not, or \ref igraph_feedback_arc_set() to find a
 *     set of edges whose removal makes the graph acyclic.
 *
 * \example examples/simple/igraph_topological_sorting.c
 */
int igraph_topological_sorting(const igraph_t* graph, igraph_vector_t *res,
                               igraph_neimode_t mode) {
    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t degrees, neis;
    igraph_dqueue_t sources;
    igraph_neimode_t deg_mode;
    long int node, i, j;

    if (mode == IGRAPH_ALL || !igraph_is_directed(graph)) {
        IGRAPH_ERROR("Topological sorting does not make sense for undirected graphs", IGRAPH_EINVAL);
    } else if (mode == IGRAPH_OUT) {
        deg_mode = IGRAPH_IN;
    } else if (mode == IGRAPH_IN) {
        deg_mode = IGRAPH_OUT;
    } else {
        IGRAPH_ERROR("Invalid mode", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&degrees, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_dqueue_init(&sources, 0));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &sources);
    IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), deg_mode, 0));

    igraph_vector_clear(res);

    /* Do we have nodes with no incoming vertices? */
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(degrees)[i] == 0) {
            IGRAPH_CHECK(igraph_dqueue_push(&sources, i));
        }
    }

    /* Take all nodes with no incoming vertices and remove them */
    while (!igraph_dqueue_empty(&sources)) {
        igraph_real_t tmp = igraph_dqueue_pop(&sources); node = (long) tmp;
        /* Add the node to the result vector */
        igraph_vector_push_back(res, node);
        /* Exclude the node from further source searches */
        VECTOR(degrees)[node] = -1;
        /* Get the neighbors and decrease their degrees by one */
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) node, mode));
        j = igraph_vector_size(&neis);
        for (i = 0; i < j; i++) {
            VECTOR(degrees)[(long)VECTOR(neis)[i]]--;
            if (VECTOR(degrees)[(long)VECTOR(neis)[i]] == 0) {
                IGRAPH_CHECK(igraph_dqueue_push(&sources, VECTOR(neis)[i]));
            }
        }
    }

    if (igraph_vector_size(res) < no_of_nodes) {
        IGRAPH_ERROR("The graph has cycles; topological sorting is only possible in acyclic graphs", IGRAPH_EINVAL);
    }

    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&neis);
    igraph_dqueue_destroy(&sources);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

/**
 * \function igraph_is_dag
 * Checks whether a graph is a directed acyclic graph (DAG) or not.
 *
 * </para><para>
 * A directed acyclic graph (DAG) is a directed graph with no cycles.
 *
 * \param graph The input graph.
 * \param res Pointer to a boolean constant, the result
 *     is stored here.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), where |V| and |E| are the number of
 * vertices and edges in the original input graph.
 *
 * \sa \ref igraph_topological_sorting() to get a possible topological
 *     sorting of a DAG.
 */
int igraph_is_dag(const igraph_t* graph, igraph_bool_t *res) {
    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t degrees, neis;
    igraph_dqueue_t sources;
    long int node, i, j, nei, vertices_left;

    if (!igraph_is_directed(graph)) {
        *res = 0;
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&degrees, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_dqueue_init(&sources, 0));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &sources);
    IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), IGRAPH_OUT, 1));

    vertices_left = no_of_nodes;

    /* Do we have nodes with no incoming edges? */
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(degrees)[i] == 0) {
            IGRAPH_CHECK(igraph_dqueue_push(&sources, i));
        }
    }

    /* Take all nodes with no incoming edges and remove them */
    while (!igraph_dqueue_empty(&sources)) {
        igraph_real_t tmp = igraph_dqueue_pop(&sources); node = (long) tmp;
        /* Exclude the node from further source searches */
        VECTOR(degrees)[node] = -1;
        vertices_left--;
        /* Get the neighbors and decrease their degrees by one */
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) node,
                                      IGRAPH_IN));
        j = igraph_vector_size(&neis);
        for (i = 0; i < j; i++) {
            nei = (long)VECTOR(neis)[i];
            if (nei == node) {
                continue;
            }
            VECTOR(degrees)[nei]--;
            if (VECTOR(degrees)[nei] == 0) {
                IGRAPH_CHECK(igraph_dqueue_push(&sources, nei));
            }
        }
    }

    *res = (vertices_left == 0);
    if (vertices_left < 0) {
        IGRAPH_WARNING("vertices_left < 0 in igraph_is_dag, possible bug");
    }

    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&neis);
    igraph_dqueue_destroy(&sources);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_simple
 * \brief Decides whether the input graph is a simple graph.
 *
 * </para><para>
 * A graph is a simple graph if it does not contain loop edges and
 * multiple edges.
 *
 * \param graph The input graph.
 * \param res Pointer to a boolean constant, the result
 *     is stored here.
 * \return Error code.
 *
 * \sa \ref igraph_is_loop() and \ref igraph_is_multiple() to
 * find the loops and multiple edges, \ref igraph_simplify() to
 * get rid of them, or \ref igraph_has_multiple() to decide whether
 * there is at least one multiple edge.
 *
 * Time complexity: O(|V|+|E|).
 */

int igraph_is_simple(const igraph_t *graph, igraph_bool_t *res) {
    long int vc = igraph_vcount(graph);
    long int ec = igraph_ecount(graph);

    if (vc == 0 || ec == 0) {
        *res = 1;
    } else {
        igraph_vector_t neis;
        long int i, j, n;
        igraph_bool_t found = 0;
        IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
        for (i = 0; i < vc; i++) {
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) i, IGRAPH_OUT));
            n = igraph_vector_size(&neis);
            for (j = 0; j < n; j++) {
                if (VECTOR(neis)[j] == i) {
                    found = 1; break;
                }
                if (j > 0 && VECTOR(neis)[j - 1] == VECTOR(neis)[j]) {
                    found = 1; break;
                }
            }
        }
        *res = !found;
        igraph_vector_destroy(&neis);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/**
 * \function igraph_has_loop
 * \brief Returns whether the graph has at least one loop edge.
 *
 * </para><para>
 * A loop edge is an edge from a vertex to itself.
 * \param graph The input graph.
 * \param res Pointer to an initialized boolean vector for storing the result.
 *
 * \sa \ref igraph_simplify() to get rid of loop edges.
 *
 * Time complexity: O(e), the number of edges to check.
 *
 * \example examples/simple/igraph_has_loop.c
 */

int igraph_has_loop(const igraph_t *graph, igraph_bool_t *res) {
    long int i, m = igraph_ecount(graph);

    *res = 0;

    for (i = 0; i < m; i++) {
        if (IGRAPH_FROM(graph, i) == IGRAPH_TO(graph, i)) {
            *res = 1;
            break;
        }
    }

    return 0;
}

/**
 * \function igraph_is_loop
 * \brief Find the loop edges in a graph.
 *
 * </para><para>
 * A loop edge is an edge from a vertex to itself.
 * \param graph The input graph.
 * \param res Pointer to an initialized boolean vector for storing the result,
 *         it will be resized as needed.
 * \param es The edges to check, for all edges supply \ref igraph_ess_all() here.
 * \return Error code.
 *
 * \sa \ref igraph_simplify() to get rid of loop edges.
 *
 * Time complexity: O(e), the number of edges to check.
 *
 * \example examples/simple/igraph_is_loop.c
 */

int igraph_is_loop(const igraph_t *graph, igraph_vector_bool_t *res,
                   igraph_es_t es) {
    igraph_eit_t eit;
    long int i;
    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    IGRAPH_CHECK(igraph_vector_bool_resize(res, IGRAPH_EIT_SIZE(eit)));

    for (i = 0; !IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
        long int e = IGRAPH_EIT_GET(eit);
        VECTOR(*res)[i] = (IGRAPH_FROM(graph, e) == IGRAPH_TO(graph, e)) ? 1 : 0;
    }

    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_has_multiple
 * \brief Check whether the graph has at least one multiple edge.
 *
 * </para><para>
 * An edge is a multiple edge if there is another
 * edge with the same head and tail vertices in the graph.
 *
 * \param graph The input graph.
 * \param res Pointer to a boolean variable, the result will be stored here.
 * \return Error code.
 *
 * \sa \ref igraph_count_multiple(), \ref igraph_is_multiple() and \ref igraph_simplify().
 *
 * Time complexity: O(e*d), e is the number of edges to check and d is the
 * average degree (out-degree in directed graphs) of the vertices at the
 * tail of the edges.
 *
 * \example examples/simple/igraph_has_multiple.c
 */

int igraph_has_multiple(const igraph_t *graph, igraph_bool_t *res) {
    long int vc = igraph_vcount(graph);
    long int ec = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);

    if (vc == 0 || ec == 0) {
        *res = 0;
    } else {
        igraph_vector_t neis;
        long int i, j, n;
        igraph_bool_t found = 0;
        IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
        for (i = 0; i < vc && !found; i++) {
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) i,
                                          IGRAPH_OUT));
            n = igraph_vector_size(&neis);
            for (j = 1; j < n; j++) {
                if (VECTOR(neis)[j - 1] == VECTOR(neis)[j]) {
                    /* If the graph is undirected, loop edges appear twice in the neighbor
                     * list, so check the next item as well */
                    if (directed) {
                        /* Directed, so this is a real multiple edge */
                        found = 1; break;
                    } else if (VECTOR(neis)[j - 1] != i) {
                        /* Undirected, but not a loop edge */
                        found = 1; break;
                    } else if (j < n - 1 && VECTOR(neis)[j] == VECTOR(neis)[j + 1]) {
                        /* Undirected, loop edge, multiple times */
                        found = 1; break;
                    }
                }
            }
        }
        *res = found;
        igraph_vector_destroy(&neis);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/**
 * \function igraph_is_multiple
 * \brief Find the multiple edges in a graph.
 *
 * </para><para>
 * An edge is a multiple edge if there is another
 * edge with the same head and tail vertices in the graph.
 *
 * </para><para>
 * Note that this function returns true only for the second or more
 * appearances of the multiple edges.
 * \param graph The input graph.
 * \param res Pointer to a boolean vector, the result will be stored
 *        here. It will be resized as needed.
 * \param es The edges to check. Supply \ref igraph_ess_all() if you want
 *        to check all edges.
 * \return Error code.
 *
 * \sa \ref igraph_count_multiple(), \ref igraph_has_multiple() and \ref igraph_simplify().
 *
 * Time complexity: O(e*d), e is the number of edges to check and d is the
 * average degree (out-degree in directed graphs) of the vertices at the
 * tail of the edges.
 *
 * \example examples/simple/igraph_is_multiple.c
 */

int igraph_is_multiple(const igraph_t *graph, igraph_vector_bool_t *res,
                       igraph_es_t es) {
    igraph_eit_t eit;
    long int i;
    igraph_lazy_inclist_t inclist;

    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    IGRAPH_CHECK(igraph_vector_bool_resize(res, IGRAPH_EIT_SIZE(eit)));

    for (i = 0; !IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
        long int e = IGRAPH_EIT_GET(eit);
        long int from = IGRAPH_FROM(graph, e);
        long int to = IGRAPH_TO(graph, e);
        igraph_vector_t *neis = igraph_lazy_inclist_get(&inclist,
                                (igraph_integer_t) from);
        long int j, n = igraph_vector_size(neis);
        VECTOR(*res)[i] = 0;
        for (j = 0; j < n; j++) {
            long int e2 = (long int) VECTOR(*neis)[j];
            long int to2 = IGRAPH_OTHER(graph, e2, from);
            if (to2 == to && e2 < e) {
                VECTOR(*res)[i] = 1;
            }
        }
    }

    igraph_lazy_inclist_destroy(&inclist);
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}


/**
 * \function igraph_count_multiple
 * \brief Count the number of appearances of the edges in a graph.
 *
 * </para><para>
 * If the graph has no multiple edges then the result vector will be
 * filled with ones.
 * (An edge is a multiple edge if there is another
 * edge with the same head and tail vertices in the graph.)
 *
 * </para><para>
 * \param graph The input graph.
 * \param res Pointer to a vector, the result will be stored
 *        here. It will be resized as needed.
 * \param es The edges to check. Supply \ref igraph_ess_all() if you want
 *        to check all edges.
 * \return Error code.
 *
 * \sa \ref igraph_is_multiple() and \ref igraph_simplify().
 *
 * Time complexity: O(E d), E is the number of edges to check and d is the
 * average degree (out-degree in directed graphs) of the vertices at the
 * tail of the edges.
 */

int igraph_count_multiple(const igraph_t *graph, igraph_vector_t *res, igraph_es_t es) {
    igraph_eit_t eit;
    long int i;
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_lazy_inclist_t inclist;

    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_EIT_SIZE(eit)));

    for (i = 0; !IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
        long int e = IGRAPH_EIT_GET(eit);
        long int from = IGRAPH_FROM(graph, e);
        long int to = IGRAPH_TO(graph, e);
        igraph_vector_t *neis = igraph_lazy_inclist_get(&inclist,
                                (igraph_integer_t) from);
        long int j, n = igraph_vector_size(neis);
        VECTOR(*res)[i] = 0;
        for (j = 0; j < n; j++) {
            long int e2 = (long int) VECTOR(*neis)[j];
            long int to2 = IGRAPH_OTHER(graph, e2, from);
            if (to2 == to) {
                VECTOR(*res)[i] += 1;
            }
        }
        /* for loop edges, divide the result by two */
        if (!directed && to == from) {
            VECTOR(*res)[i] /= 2;
        }
    }

    igraph_lazy_inclist_destroy(&inclist);
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_girth
 * \brief The girth of a graph is the length of the shortest cycle in it.
 *
 * </para><para>
 * The current implementation works for undirected graphs only,
 * directed graphs are treated as undirected graphs. Self-loops and
 * multiple edges are ignored.
 *
 * </para><para>
 * For graphs that contain no cycles, and only for such graphs, 
 * zero is returned. Note that in some applications, it is customary 
 * to define the girth of acyclic graphs to be infinity. However, infinity
 * is not representable as an \c igraph_integer_t, therefore zero is used
 * for this case.
 *
 * </para><para>
 * This implementation is based on Alon Itai and Michael Rodeh:
 * Finding a minimum circuit in a graph
 * \emb Proceedings of the ninth annual ACM symposium on Theory of
 * computing \eme, 1-10, 1977. The first implementation of this
 * function was done by Keith Briggs, thanks Keith.
 * \param graph The input graph.
 * \param girth Pointer to an integer, if not \c NULL then the result
 *     will be stored here.
 * \param circle Pointer to an initialized vector, the vertex ids in
 *     the shortest circle will be stored here. If \c NULL then it is
 *     ignored.
 * \return Error code.
 *
 * Time complexity: O((|V|+|E|)^2), |V| is the number of vertices, |E|
 * is the number of edges in the general case. If the graph has no
 * cycles at all then the function needs O(|V|+|E|) time to realize
 * this and then it stops.
 *
 * \example examples/simple/igraph_girth.c
 */

int igraph_girth(const igraph_t *graph, igraph_integer_t *girth,
                 igraph_vector_t *circle) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_t q;
    igraph_lazy_adjlist_t adjlist;
    long int mincirc = LONG_MAX, minvertex = 0;
    long int node;
    igraph_bool_t triangle = 0;
    igraph_vector_t *neis;
    igraph_vector_long_t level;
    long int stoplevel = no_of_nodes + 1;
    igraph_bool_t anycircle = 0;
    long int t1 = 0, t2 = 0;

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, IGRAPH_ALL,
                                          IGRAPH_SIMPLIFY));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vector_long_init(&level, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &level);

    for (node = 0; !triangle && node < no_of_nodes; node++) {

        /* Are there circles in this graph at all? */
        if (node == 1 && anycircle == 0) {
            igraph_bool_t conn;
            IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_WEAK));
            if (conn) {
                /* No, there are none */
                break;
            }
        }

        anycircle = 0;
        igraph_dqueue_clear(&q);
        igraph_vector_long_null(&level);
        IGRAPH_CHECK(igraph_dqueue_push(&q, node));
        VECTOR(level)[node] = 1;

        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);
            long int actlevel = VECTOR(level)[actnode];
            long int i, n;

            if (actlevel >= stoplevel) {
                break;
            }

            neis = igraph_lazy_adjlist_get(&adjlist, (igraph_integer_t) actnode);
            n = igraph_vector_size(neis);
            for (i = 0; i < n; i++) {
                long int nei = (long int) VECTOR(*neis)[i];
                long int neilevel = VECTOR(level)[nei];
                if (neilevel != 0) {
                    if (neilevel == actlevel - 1) {
                        continue;
                    } else {
                        /* found circle */
                        stoplevel = neilevel;
                        anycircle = 1;
                        if (actlevel < mincirc) {
                            /* Is it a minimum circle? */
                            mincirc = actlevel + neilevel - 1;
                            minvertex = node;
                            t1 = actnode; t2 = nei;
                            if (neilevel == 2) {
                                /* Is it a triangle? */
                                triangle = 1;
                            }
                        }
                        if (neilevel == actlevel) {
                            break;
                        }
                    }
                } else {
                    igraph_dqueue_push(&q, nei);
                    VECTOR(level)[nei] = actlevel + 1;
                }
            }

        } /* while q !empty */
    } /* node */

    if (girth) {
        if (mincirc == LONG_MAX) {
            *girth = mincirc = 0;
        } else {
            *girth = (igraph_integer_t) mincirc;
        }
    }

    /* Store the actual circle, if needed */
    if (circle) {
        IGRAPH_CHECK(igraph_vector_resize(circle, mincirc));
        if (mincirc != 0) {
            long int i, n, idx = 0;
            igraph_dqueue_clear(&q);
            igraph_vector_long_null(&level); /* used for father pointers */
#define FATHER(x) (VECTOR(level)[(x)])
            IGRAPH_CHECK(igraph_dqueue_push(&q, minvertex));
            FATHER(minvertex) = minvertex;
            while (FATHER(t1) == 0 || FATHER(t2) == 0) {
                long int actnode = (long int) igraph_dqueue_pop(&q);
                neis = igraph_lazy_adjlist_get(&adjlist, (igraph_integer_t) actnode);
                n = igraph_vector_size(neis);
                for (i = 0; i < n; i++) {
                    long int nei = (long int) VECTOR(*neis)[i];
                    if (FATHER(nei) == 0) {
                        FATHER(nei) = actnode + 1;
                        igraph_dqueue_push(&q, nei);
                    }
                }
            }  /* while q !empty */
            /* Ok, now use FATHER to create the path */
            while (t1 != minvertex) {
                VECTOR(*circle)[idx++] = t1;
                t1 = FATHER(t1) - 1;
            }
            VECTOR(*circle)[idx] = minvertex;
            idx = mincirc - 1;
            while (t2 != minvertex) {
                VECTOR(*circle)[idx--] = t2;
                t2 = FATHER(t2) - 1;
            }
        } /* anycircle */
    } /* circle */
#undef FATHER

    igraph_vector_long_destroy(&level);
    igraph_dqueue_destroy(&q);
    igraph_lazy_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

static int igraph_i_linegraph_undirected(const igraph_t *graph, igraph_t *linegraph);

static int igraph_i_linegraph_directed(const igraph_t *graph, igraph_t *linegraph);

/* Note to self: tried using adjacency lists instead of igraph_incident queries,
 * with minimal performance improvements on a graph with 70K vertices and 360K
 * edges. (1.09s instead of 1.10s). I think it's not worth the fuss. */
static int igraph_i_linegraph_undirected(const igraph_t *graph, igraph_t *linegraph) {
    long int no_of_edges = igraph_ecount(graph);
    long int i, j, n;
    igraph_vector_t adjedges, adjedges2;
    igraph_vector_t edges;
    long int prev = -1;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&adjedges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&adjedges2, 0);

    for (i = 0; i < no_of_edges; i++) {
        long int from = IGRAPH_FROM(graph, i);
        long int to = IGRAPH_TO(graph, i);

        IGRAPH_ALLOW_INTERRUPTION();

        if (from != prev) {
            IGRAPH_CHECK(igraph_incident(graph, &adjedges, (igraph_integer_t) from,
                                         IGRAPH_ALL));
        }
        n = igraph_vector_size(&adjedges);
        for (j = 0; j < n; j++) {
            long int e = (long int) VECTOR(adjedges)[j];
            if (e < i) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, e));
            }
        }

        IGRAPH_CHECK(igraph_incident(graph, &adjedges2, (igraph_integer_t) to,
                                     IGRAPH_ALL));
        n = igraph_vector_size(&adjedges2);
        for (j = 0; j < n; j++) {
            long int e = (long int) VECTOR(adjedges2)[j];
            if (e < i) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, e));
            }
        }

        prev = from;
    }

    igraph_vector_destroy(&adjedges);
    igraph_vector_destroy(&adjedges2);
    IGRAPH_FINALLY_CLEAN(2);

    igraph_create(linegraph, &edges, (igraph_integer_t) no_of_edges,
                  igraph_is_directed(graph));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_linegraph_directed(const igraph_t *graph, igraph_t *linegraph) {
    long int no_of_edges = igraph_ecount(graph);
    long int i, j, n;
    igraph_vector_t adjedges;
    igraph_vector_t edges;
    long int prev = -1;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&adjedges, 0);

    for (i = 0; i < no_of_edges; i++) {
        long int from = IGRAPH_FROM(graph, i);

        IGRAPH_ALLOW_INTERRUPTION();

        if (from != prev) {
            IGRAPH_CHECK(igraph_incident(graph, &adjedges, (igraph_integer_t) from,
                                         IGRAPH_IN));
        }
        n = igraph_vector_size(&adjedges);
        for (j = 0; j < n; j++) {
            long int e = (long int) VECTOR(adjedges)[j];
            IGRAPH_CHECK(igraph_vector_push_back(&edges, e));
            IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
        }

        prev = from;
    }

    igraph_vector_destroy(&adjedges);
    IGRAPH_FINALLY_CLEAN(1);
    igraph_create(linegraph, &edges, (igraph_integer_t) no_of_edges, igraph_is_directed(graph));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_linegraph
 * \brief Create the line graph of a graph.
 *
 * The line graph L(G) of a G undirected graph is defined as follows.
 * L(G) has one vertex for each edge in G and two vertices in L(G) are connected
 * by an edge if their corresponding edges share an end point.
 *
 * </para><para>
 * The line graph L(G) of a G directed graph is slightly different,
 * L(G) has one vertex for each edge in G and two vertices in L(G) are connected
 * by a directed edge if the target of the first vertex's corresponding edge
 * is the same as the source of the second vertex's corresponding edge.
 *
 * </para><para>
 * Edge \em i  in the original graph will correspond to vertex \em i
 * in the line graph.
 *
 * </para><para>
 * The first version of this function was contributed by Vincent Matossian,
 * thanks.
 * \param graph The input graph, may be directed or undirected.
 * \param linegraph Pointer to an uninitialized graph object, the
 *        result is stored here.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), the number of edges plus the number of vertices.
 */

int igraph_linegraph(const igraph_t *graph, igraph_t *linegraph) {

    if (igraph_is_directed(graph)) {
        return igraph_i_linegraph_directed(graph, linegraph);
    } else {
        return igraph_i_linegraph_undirected(graph, linegraph);
    }
}

/**
 * \function igraph_add_edge
 * \brief Adds a single edge to a graph.
 *
 * </para><para>
 * For directed graphs the edge points from \p from to \p to.
 *
 * </para><para>
 * Note that if you want to add many edges to a big graph, then it is
 * inefficient to add them one by one, it is better to collect them into
 * a vector and add all of them via a single \ref igraph_add_edges() call.
 * \param igraph The graph.
 * \param from The id of the first vertex of the edge.
 * \param to The id of the second vertex of the edge.
 * \return Error code.
 *
 * \sa \ref igraph_add_edges() to add many edges, \ref
 * igraph_delete_edges() to remove edges and \ref
 * igraph_add_vertices() to add vertices.
 *
 * Time complexity: O(|V|+|E|), the number of edges plus the number of
 * vertices.
 */

int igraph_add_edge(igraph_t *graph, igraph_integer_t from, igraph_integer_t to) {

    igraph_vector_t edges;
    int ret;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 2);

    VECTOR(edges)[0] = from;
    VECTOR(edges)[1] = to;
    IGRAPH_CHECK(ret = igraph_add_edges(graph, &edges, 0));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return ret;
}

/*
 * \example examples/simple/graph_convergence_degree.c
 */

int igraph_convergence_degree(const igraph_t *graph, igraph_vector_t *result,
                              igraph_vector_t *ins, igraph_vector_t *outs) {
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    long int i, j, k, n;
    long int *geodist;
    igraph_vector_int_t *eids;
    igraph_vector_t *ins_p, *outs_p, ins_v, outs_v;
    igraph_dqueue_t q;
    igraph_inclist_t inclist;
    igraph_bool_t directed = igraph_is_directed(graph);

    if (result != 0) {
        IGRAPH_CHECK(igraph_vector_resize(result, no_of_edges));
    }
    IGRAPH_CHECK(igraph_dqueue_init(&q, 100));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &q);

    if (ins == 0) {
        ins_p = &ins_v;
        IGRAPH_VECTOR_INIT_FINALLY(ins_p, no_of_edges);
    } else {
        ins_p = ins;
        IGRAPH_CHECK(igraph_vector_resize(ins_p, no_of_edges));
        igraph_vector_null(ins_p);
    }

    if (outs == 0) {
        outs_p = &outs_v;
        IGRAPH_VECTOR_INIT_FINALLY(outs_p, no_of_edges);
    } else {
        outs_p = outs;
        IGRAPH_CHECK(igraph_vector_resize(outs_p, no_of_edges));
        igraph_vector_null(outs_p);
    }

    geodist = igraph_Calloc(no_of_nodes, long int);
    if (geodist == 0) {
        IGRAPH_ERROR("Cannot calculate convergence degrees", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, geodist);

    /* Collect shortest paths originating from/to every node to correctly
     * determine input field sizes */
    for (k = 0; k < (directed ? 2 : 1); k++) {
        igraph_neimode_t neimode = (k == 0) ? IGRAPH_OUT : IGRAPH_IN;
        igraph_real_t *vec;
        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, neimode));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
        vec = (k == 0) ? VECTOR(*ins_p) : VECTOR(*outs_p);
        for (i = 0; i < no_of_nodes; i++) {
            igraph_dqueue_clear(&q);
            memset(geodist, 0, sizeof(long int) * (size_t) no_of_nodes);
            geodist[i] = 1;
            IGRAPH_CHECK(igraph_dqueue_push(&q, i));
            IGRAPH_CHECK(igraph_dqueue_push(&q, 0.0));
            while (!igraph_dqueue_empty(&q)) {
                long int actnode = (long int) igraph_dqueue_pop(&q);
                long int actdist = (long int) igraph_dqueue_pop(&q);
                IGRAPH_ALLOW_INTERRUPTION();
                eids = igraph_inclist_get(&inclist, actnode);
                n = igraph_vector_int_size(eids);
                for (j = 0; j < n; j++) {
                    long int neighbor = IGRAPH_OTHER(graph, VECTOR(*eids)[j], actnode);
                    if (geodist[neighbor] != 0) {
                        /* we've already seen this node, another shortest path? */
                        if (geodist[neighbor] - 1 == actdist + 1) {
                            /* Since this edge is in the BFS tree rooted at i, we must
                             * increase either the size of the infield or the outfield */
                            if (!directed) {
                                if (actnode < neighbor) {
                                    VECTOR(*ins_p)[(long int)VECTOR(*eids)[j]] += 1;
                                } else {
                                    VECTOR(*outs_p)[(long int)VECTOR(*eids)[j]] += 1;
                                }
                            } else {
                                vec[(long int)VECTOR(*eids)[j]] += 1;
                            }
                        } else if (geodist[neighbor] - 1 < actdist + 1) {
                            continue;
                        }
                    } else {
                        /* we haven't seen this node yet */
                        IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
                        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
                        /* Since this edge is in the BFS tree rooted at i, we must
                         * increase either the size of the infield or the outfield */
                        if (!directed) {
                            if (actnode < neighbor) {
                                VECTOR(*ins_p)[(long int)VECTOR(*eids)[j]] += 1;
                            } else {
                                VECTOR(*outs_p)[(long int)VECTOR(*eids)[j]] += 1;
                            }
                        } else {
                            vec[(long int)VECTOR(*eids)[j]] += 1;
                        }
                        geodist[neighbor] = actdist + 2;
                    }
                }
            }
        }

        igraph_inclist_destroy(&inclist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (result != 0) {
        for (i = 0; i < no_of_edges; i++)
            VECTOR(*result)[i] = (VECTOR(*ins_p)[i] - VECTOR(*outs_p)[i]) /
                                 (VECTOR(*ins_p)[i] + VECTOR(*outs_p)[i]);
        if (!directed) {
            for (i = 0; i < no_of_edges; i++)
                if (VECTOR(*result)[i] < 0) {
                    VECTOR(*result)[i] = -VECTOR(*result)[i];
                }
        }
    }

    if (ins == 0) {
        igraph_vector_destroy(ins_p);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (outs == 0) {
        igraph_vector_destroy(outs_p);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_Free(geodist);
    igraph_dqueue_destroy(&q);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

/**
 * \function igraph_unfold_tree
 * Unfolding a graph into a tree, by possibly multiplicating its vertices.
 *
 * A graph is converted into a tree (or forest, if it is unconnected),
 * by performing a breadth-first search on it, and replicating
 * vertices that were found a second, third, etc. time.
 * \param graph The input graph, it can be either directed or
 *   undirected.
 * \param tree Pointer to an uninitialized graph object, the result is
 *   stored here.
 * \param mode For directed graphs; whether to follow paths along edge
 *    directions (\c IGRAPH_OUT), or the opposite (\c IGRAPH_IN), or
 *    ignore edge directions completely (\c IGRAPH_ALL). It is ignored
 *    for undirected graphs.
 * \param roots A numeric vector giving the root vertex, or vertices
 *   (if the graph is not connected), to start from.
 * \param vertex_index Pointer to an initialized vector, or a null
 *   pointer. If not a null pointer, then a mapping from the vertices
 *   in the new graph to the ones in the original is created here.
 * \return Error code.
 *
 * Time complexity: O(n+m), linear in the number vertices and edges.
 *
 */

int igraph_unfold_tree(const igraph_t *graph, igraph_t *tree,
                       igraph_neimode_t mode, const igraph_vector_t *roots,
                       igraph_vector_t *vertex_index) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    long int no_of_roots = igraph_vector_size(roots);
    long int tree_vertex_count = no_of_nodes;

    igraph_vector_t edges;
    igraph_vector_bool_t seen_vertices;
    igraph_vector_bool_t seen_edges;

    igraph_dqueue_t Q;
    igraph_vector_t neis;

    long int i, n, r, v_ptr = no_of_nodes;

    /* TODO: handle not-connected graphs, multiple root vertices */

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    igraph_vector_reserve(&edges, no_of_edges * 2);
    IGRAPH_DQUEUE_INIT_FINALLY(&Q, 100);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&seen_vertices, no_of_nodes);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&seen_edges, no_of_edges);

    if (vertex_index) {
        IGRAPH_CHECK(igraph_vector_resize(vertex_index, no_of_nodes));
        for (i = 0; i < no_of_nodes; i++) {
            VECTOR(*vertex_index)[i] = i;
        }
    }

    for (r = 0; r < no_of_roots; r++) {

        long int root = (long int) VECTOR(*roots)[r];
        VECTOR(seen_vertices)[root] = 1;
        igraph_dqueue_push(&Q, root);

        while (!igraph_dqueue_empty(&Q)) {
            long int actnode = (long int) igraph_dqueue_pop(&Q);

            IGRAPH_CHECK(igraph_incident(graph, &neis, (igraph_integer_t) actnode, mode));
            n = igraph_vector_size(&neis);
            for (i = 0; i < n; i++) {

                long int edge = (long int) VECTOR(neis)[i];
                long int from = IGRAPH_FROM(graph, edge);
                long int to = IGRAPH_TO(graph, edge);
                long int nei = IGRAPH_OTHER(graph, edge, actnode);

                if (! VECTOR(seen_edges)[edge]) {

                    VECTOR(seen_edges)[edge] = 1;

                    if (! VECTOR(seen_vertices)[nei]) {

                        igraph_vector_push_back(&edges, from);
                        igraph_vector_push_back(&edges, to);

                        VECTOR(seen_vertices)[nei] = 1;
                        IGRAPH_CHECK(igraph_dqueue_push(&Q, nei));

                    } else {

                        tree_vertex_count++;
                        if (vertex_index) {
                            IGRAPH_CHECK(igraph_vector_push_back(vertex_index, nei));
                        }

                        if (from == nei) {
                            igraph_vector_push_back(&edges, v_ptr++);
                            igraph_vector_push_back(&edges, to);
                        } else {
                            igraph_vector_push_back(&edges, from);
                            igraph_vector_push_back(&edges, v_ptr++);
                        }
                    }
                }

            } /* for i<n */

        } /* ! igraph_dqueue_empty(&Q) */

    } /* r < igraph_vector_size(roots) */

    igraph_vector_bool_destroy(&seen_edges);
    igraph_vector_bool_destroy(&seen_vertices);
    igraph_vector_destroy(&neis);
    igraph_dqueue_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(4);

    IGRAPH_CHECK(igraph_create(tree, &edges, tree_vertex_count, igraph_is_directed(graph)));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_is_mutual
 * Check whether the edges of a directed graph are mutual.
 *
 * An (A,B) edge is mutual if the graph contains the (B,A) edge, too.
 * </para>
 *
 * <para>An undirected graph only has mutual edges, by definition.
 * </para>
 *
 * <para>Edge multiplicity is not considered here, e.g. if there are two
 * (A,B) edges and one (B,A) edge, then all three are considered to be
 * mutual.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the result is stored
 *        here.
 * \param es The sequence of edges to check. Supply
 *        <code>igraph_ess_all()</code> for all edges, see \ref
 *        igraph_ess_all().
 * \return Error code.
 *
 * Time complexity: O(n log(d)), n is the number of edges supplied, d
 * is the maximum in-degree of the vertices that are targets of the
 * supplied edges. An upper limit of the time complexity is O(n log(|E|)),
 * |E| is the number of edges in the graph.
 */

int igraph_is_mutual(igraph_t *graph, igraph_vector_bool_t *res, igraph_es_t es) {

    igraph_eit_t eit;
    igraph_lazy_adjlist_t adjlist;
    long int i;

    /* How many edges do we have? */
    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);
    IGRAPH_CHECK(igraph_vector_bool_resize(res, IGRAPH_EIT_SIZE(eit)));

    /* An undirected graph has mutual edges by definition,
       res is already properly resized */
    if (! igraph_is_directed(graph)) {
        igraph_vector_bool_fill(res, 1);
        igraph_eit_destroy(&eit);
        IGRAPH_FINALLY_CLEAN(1);
        return 0;
    }

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, IGRAPH_OUT, IGRAPH_DONT_SIMPLIFY));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);

    for (i = 0; ! IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
        long int edge = IGRAPH_EIT_GET(eit);
        long int from = IGRAPH_FROM(graph, edge);
        long int to = IGRAPH_TO(graph, edge);

        /* Check whether there is a to->from edge, search for from in the
           out-list of to. We don't search an empty vector, because
           vector_binsearch seems to have a bug with this. */
        igraph_vector_t *neis = igraph_lazy_adjlist_get(&adjlist,
                                (igraph_integer_t) to);
        if (igraph_vector_empty(neis)) {
            VECTOR(*res)[i] = 0;
        } else {
            VECTOR(*res)[i] = igraph_vector_binsearch2(neis, from);
        }
    }

    igraph_lazy_adjlist_destroy(&adjlist);
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

static int igraph_i_avg_nearest_neighbor_degree_weighted(const igraph_t *graph,
        igraph_vs_t vids,
        igraph_neimode_t mode,
        igraph_neimode_t neighbor_degree_mode,
        igraph_vector_t *knn,
        igraph_vector_t *knnk,
        const igraph_vector_t *weights);

static int igraph_i_avg_nearest_neighbor_degree_weighted(const igraph_t *graph,
        igraph_vs_t vids,
        igraph_neimode_t mode,
        igraph_neimode_t neighbor_degree_mode,
        igraph_vector_t *knn,
        igraph_vector_t *knnk,
        const igraph_vector_t *weights) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t neis, edge_neis;
    long int i, j, no_vids;
    igraph_vit_t vit;
    igraph_vector_t my_knn_v, *my_knn = knn;
    igraph_vector_t strength, deg;
    igraph_integer_t maxdeg;
    igraph_vector_t deghist;
    igraph_real_t mynan = IGRAPH_NAN;

    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid weight vector size", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    no_vids = IGRAPH_VIT_SIZE(vit);

    if (!knn) {
        IGRAPH_VECTOR_INIT_FINALLY(&my_knn_v, no_vids);
        my_knn = &my_knn_v;
    } else {
        IGRAPH_CHECK(igraph_vector_resize(knn, no_vids));
    }

    /* Get degree of neighbours */
    IGRAPH_VECTOR_INIT_FINALLY(&deg, no_of_nodes);
    IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(),
                               neighbor_degree_mode, IGRAPH_LOOPS));
    IGRAPH_VECTOR_INIT_FINALLY(&strength, no_of_nodes);

    /* Get strength of all nodes */
    IGRAPH_CHECK(igraph_strength(graph, &strength, igraph_vss_all(),
                                 mode, IGRAPH_LOOPS, weights));

    /* Get maximum degree for initialization */
    IGRAPH_CHECK(igraph_maxdegree(graph, &maxdeg, igraph_vss_all(),
                                  mode, IGRAPH_LOOPS));
    IGRAPH_VECTOR_INIT_FINALLY(&neis, (long int)maxdeg);
    IGRAPH_VECTOR_INIT_FINALLY(&edge_neis, (long int)maxdeg);
    igraph_vector_resize(&neis, 0);
    igraph_vector_resize(&edge_neis, 0);

    if (knnk) {
        IGRAPH_CHECK(igraph_vector_resize(knnk, (long int)maxdeg));
        igraph_vector_null(knnk);
        IGRAPH_VECTOR_INIT_FINALLY(&deghist, (long int)maxdeg);
    }

    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_real_t sum = 0.0;
        long int v = IGRAPH_VIT_GET(vit);
        long int nv;
        igraph_real_t str = VECTOR(strength)[v];
        /* Get neighbours and incident edges */
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) v, mode));
        IGRAPH_CHECK(igraph_incident(graph, &edge_neis, (igraph_integer_t) v, mode));
        nv = igraph_vector_size(&neis);
        for (j = 0; j < nv; j++) {
            long int nei = (long int) VECTOR(neis)[j];
            long int e = (long int) VECTOR(edge_neis)[j];
            double w = VECTOR(*weights)[e];
            sum += w * VECTOR(deg)[nei];
        }
        if (str != 0.0) {
            VECTOR(*my_knn)[i] = sum / str;
        } else {
            VECTOR(*my_knn)[i] = mynan;
        }
        if (knnk && nv > 0) {
            VECTOR(*knnk)[nv - 1] += VECTOR(*my_knn)[i];
            VECTOR(deghist)[nv - 1] += 1;
        }
    }

    igraph_vector_destroy(&edge_neis);
    igraph_vector_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(2);

    if (knnk) {
        for (i = 0; i < maxdeg; i++) {
            igraph_real_t dh = VECTOR(deghist)[i];
            if (dh != 0) {
                VECTOR(*knnk)[i] /= dh;
            } else {
                VECTOR(*knnk)[i] = mynan;
            }
        }

        igraph_vector_destroy(&deghist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_destroy(&strength);
    igraph_vector_destroy(&deg);
    IGRAPH_FINALLY_CLEAN(2);

    if (!knn) {
        igraph_vector_destroy(&my_knn_v);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_avg_nearest_neighbor_degree
 * Average neighbor degree.
 *
 * Calculates the average degree of the neighbors for each vertex (\p knn), and
 * optionally, the same quantity as a function of the vertex degree (\p knnk).
 *
 * </para><para>
 * For isolated vertices \p knn is set to NaN.
 * The same is done in \p knnk for vertex degrees that
 * don't appear in the graph.
 *
 * </para><para>
 * The weighted version computes a weighted average of the neighbor degrees as
 *
 * <code>k_nn_u = 1/s_u sum_v w_uv k_v</code>,
 *
 * where <code>s_u = sum_v w_uv</code> is the sum of the incident edge weights
 * of vertex \c u, i.e. its strength.
 * The sum runs over the neighbors \c v of vertex \c u
 * as indicated by \p mode. <code>w_uv</code> denotes the weighted adjacency matrix
 * and <code>k_v</code> is the neighbors' degree, specified by \p neighbor_degree_mode.
 *
 * </para><para>
 * Reference:
 * A. Barrat, M. Barthélemy, R. Pastor-Satorras, and A. Vespignani,
 * The architecture of complex weighted networks,
 * Proc. Natl. Acad. Sci. USA 101, 3747 (2004).
 * https://dx.doi.org/10.1073/pnas.0400087101
 *
 * \param graph The input graph. It may be directed.
 * \param vids The vertices for which the calculation is performed.
 * \param mode The type of neighbors to consider in directed graphs.
 *   \c IGRAPH_OUT considers out-neighbors, \c IGRAPH_IN in-neighbors
 *   and \c IGRAPH_ALL ignores edge directions.
 * \param neighbor_degree_mode The type of degree to average in directed graphs.
 *   \c IGRAPH_OUT averages out-degrees, \c IGRAPH_IN averages in-degrees
 *   and \c IGRAPH_ALL ignores edge directions for the degree calculation.
 * \param vids The vertices for which the calculation is performed.
 * \param knn Pointer to an initialized vector, the result will be
 *   stored here. It will be resized as needed. Supply a \c NULL pointer
 *   here, if you only want to calculate \c knnk.
 * \param knnk Pointer to an initialized vector, the average
 *   neighbor degree as a function of the vertex degree is stored
 *   here. The first (zeroth) element is for degree one vertices,
 *   etc. Supply a \c NULL pointer here if you don't want to calculate
 *   this.
 * \param weights Optional edge weights. Supply a null pointer here
 *   for the non-weighted version.
 *
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 *
 * \example examples/simple/igraph_knn.c
 */

int igraph_avg_nearest_neighbor_degree(const igraph_t *graph,
                                       igraph_vs_t vids,
                                       igraph_neimode_t mode,
                                       igraph_neimode_t neighbor_degree_mode,
                                       igraph_vector_t *knn,
                                       igraph_vector_t *knnk,
                                       const igraph_vector_t *weights) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t neis;
    long int i, j, no_vids;
    igraph_vit_t vit;
    igraph_vector_t my_knn_v, *my_knn = knn;
    igraph_vector_t deg;
    igraph_integer_t maxdeg;
    igraph_vector_t deghist;
    igraph_real_t mynan = IGRAPH_NAN;
    igraph_bool_t simple;

    IGRAPH_CHECK(igraph_is_simple(graph, &simple));
    if (!simple) {
        IGRAPH_ERROR("Average nearest neighbor degree works only with "
                     "simple graphs", IGRAPH_EINVAL);
    }

    if (weights) {
        return igraph_i_avg_nearest_neighbor_degree_weighted(graph, vids,
                mode, neighbor_degree_mode, knn, knnk, weights);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    no_vids = IGRAPH_VIT_SIZE(vit);

    if (!knn) {
        IGRAPH_VECTOR_INIT_FINALLY(&my_knn_v, no_vids);
        my_knn = &my_knn_v;
    } else {
        IGRAPH_CHECK(igraph_vector_resize(knn, no_vids));
    }

    IGRAPH_VECTOR_INIT_FINALLY(&deg, no_of_nodes);
    IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(),
                               neighbor_degree_mode, IGRAPH_LOOPS));
    igraph_maxdegree(graph, &maxdeg, igraph_vss_all(), mode, IGRAPH_LOOPS);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, maxdeg);
    igraph_vector_resize(&neis, 0);

    if (knnk) {
        IGRAPH_CHECK(igraph_vector_resize(knnk, (long int)maxdeg));
        igraph_vector_null(knnk);
        IGRAPH_VECTOR_INIT_FINALLY(&deghist, (long int)maxdeg);
    }

    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_real_t sum = 0.0;
        long int v = IGRAPH_VIT_GET(vit);
        long int nv;
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) v, mode));
        nv = igraph_vector_size(&neis);
        for (j = 0; j < nv; j++) {
            long int nei = (long int) VECTOR(neis)[j];
            sum += VECTOR(deg)[nei];
        }
        if (nv != 0) {
            VECTOR(*my_knn)[i] = sum / nv;
        } else {
            VECTOR(*my_knn)[i] = mynan;
        }
        if (knnk && nv > 0) {
            VECTOR(*knnk)[nv - 1] += VECTOR(*my_knn)[i];
            VECTOR(deghist)[nv - 1] += 1;
        }
    }

    if (knnk) {
        for (i = 0; i < maxdeg; i++) {
            long int dh = (long int) VECTOR(deghist)[i];
            if (dh != 0) {
                VECTOR(*knnk)[i] /= dh;
            } else {
                VECTOR(*knnk)[i] = mynan;
            }
        }
        igraph_vector_destroy(&deghist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_destroy(&neis);
    igraph_vector_destroy(&deg);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(3);

    if (!knn) {
        igraph_vector_destroy(&my_knn_v);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/**
 * \function igraph_strength
 * Strength of the vertices, weighted vertex degree in other words.
 *
 * In a weighted network the strength of a vertex is the sum of the
 * weights of all incident edges. In a non-weighted network this is
 * exactly the vertex degree.
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the result is stored
 *   here. It will be resized as needed.
 * \param vids The vertices for which the calculation is performed.
 * \param mode Gives whether to count only outgoing (\c IGRAPH_OUT),
 *   incoming (\c IGRAPH_IN) edges or both (\c IGRAPH_ALL).
 * \param loops A logical scalar, whether to count loop edges as well.
 * \param weights A vector giving the edge weights. If this is a NULL
 *   pointer, then \ref igraph_degree() is called to perform the
 *   calculation.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number vertices and
 * edges.
 *
 * \sa \ref igraph_degree() for the traditional, non-weighted version.
 */

int igraph_strength(const igraph_t *graph, igraph_vector_t *res,
                    const igraph_vs_t vids, igraph_neimode_t mode,
                    igraph_bool_t loops, const igraph_vector_t *weights) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_vit_t vit;
    long int no_vids;
    igraph_vector_t neis;
    long int i;

    if (!weights) {
        return igraph_degree(graph, res, vids, mode, loops);
    }

    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    no_vids = IGRAPH_VIT_SIZE(vit);

    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&neis, no_of_nodes));
    IGRAPH_CHECK(igraph_vector_resize(res, no_vids));
    igraph_vector_null(res);

    if (loops) {
        for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
            long int vid = IGRAPH_VIT_GET(vit);
            long int j, n;
            IGRAPH_CHECK(igraph_incident(graph, &neis, (igraph_integer_t) vid, mode));
            n = igraph_vector_size(&neis);
            for (j = 0; j < n; j++) {
                long int edge = (long int) VECTOR(neis)[j];
                VECTOR(*res)[i] += VECTOR(*weights)[edge];
            }
        }
    } else {
        for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
            long int vid = IGRAPH_VIT_GET(vit);
            long int j, n;
            IGRAPH_CHECK(igraph_incident(graph, &neis, (igraph_integer_t) vid, mode));
            n = igraph_vector_size(&neis);
            for (j = 0; j < n; j++) {
                long int edge = (long int) VECTOR(neis)[j];
                long int from = IGRAPH_FROM(graph, edge);
                long int to = IGRAPH_TO(graph, edge);
                if (from != to) {
                    VECTOR(*res)[i] += VECTOR(*weights)[edge];
                }
            }
        }
    }

    igraph_vit_destroy(&vit);
    igraph_vector_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}


/**
 * \function igraph_sort_vertex_ids_by_degree
 * \brief Calculate a list of vertex ids sorted by degree of the corresponding vertex.
 *
 * The list of vertex ids is returned in a vector that is sorted
 * in ascending or descending order of vertex degree.
 *
 * \param graph The input graph.
 * \param outvids Pointer to an initialized vector that will be
 *        resized and will contain the ordered vertex ids.
 * \param vids Input vertex selector of vertex ids to include in
 *        calculation.
 * \param mode Defines the type of the degree.
 *        \c IGRAPH_OUT, out-degree,
 *        \c IGRAPH_IN, in-degree,
 *        \c IGRAPH_ALL, total degree (sum of the
 *        in- and out-degree).
 *        This parameter is ignored for undirected graphs.
 * \param loops Boolean, gives whether the self-loops should be
 *        counted.
 * \param order Specifies whether the ordering should be ascending
 *        (\c IGRAPH_ASCENDING) or descending (\c IGRAPH_DESCENDING).
 * \param only_indices If true, then return a sorted list of indices
 *        into a vector corresponding to \c vids, rather than a list
 *        of vertex ids. This parameter is ignored if \c vids is set
 *        to all vertices via igraph_vs_all() or igraph_vss_all(),
 *        because in this case the indices and vertex ids are the
 *        same.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex id.
 *         \c IGRAPH_EINVMODE: invalid mode argument.
 *
 */

int igraph_sort_vertex_ids_by_degree(const igraph_t *graph,
                                     igraph_vector_t *outvids,
                                     igraph_vs_t vids,
                                     igraph_neimode_t mode,
                                     igraph_bool_t loops,
                                     igraph_order_t order,
                                     igraph_bool_t only_indices) {
    long int i;
    igraph_vector_t degrees, vs_vec;
    IGRAPH_VECTOR_INIT_FINALLY(&degrees, 0);
    IGRAPH_CHECK(igraph_degree(graph, &degrees, vids, mode, loops));
    IGRAPH_CHECK((int) igraph_vector_qsort_ind(&degrees, outvids,
                 order == IGRAPH_DESCENDING));
    if (only_indices || igraph_vs_is_all(&vids) ) {
        igraph_vector_destroy(&degrees);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(&vs_vec, 0);
        IGRAPH_CHECK(igraph_vs_as_vector(graph, vids, &vs_vec));
        for (i = 0; i < igraph_vector_size(outvids); i++) {
            VECTOR(*outvids)[i] = VECTOR(vs_vec)[(long int)VECTOR(*outvids)[i]];
        }
        igraph_vector_destroy(&vs_vec);
        igraph_vector_destroy(&degrees);
        IGRAPH_FINALLY_CLEAN(2);
    }
    return 0;
}

/**
 * \function igraph_contract_vertices
 * Replace multiple vertices with a single one.
 *
 * This function creates a new graph, by merging several
 * vertices into one. The vertices in the new graph correspond
 * to sets of vertices in the input graph.
 * \param graph The input graph, it can be directed or
 *        undirected.
 * \param mapping A vector giving the mapping. For each
 *        vertex in the original graph, it should contain
 *        its id in the new graph.
 * \param vertex_comb What to do with the vertex attributes.
 *        See the igraph manual section about attributes for
 *        details.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number
 * or vertices plus edges.
 */

int igraph_contract_vertices(igraph_t *graph,
                             const igraph_vector_t *mapping,
                             const igraph_attribute_combination_t *vertex_comb) {
    igraph_vector_t edges;
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_bool_t vattr = vertex_comb && igraph_has_attribute_table();
    igraph_t res;
    long int e, last = -1;
    long int no_new_vertices;

    if (igraph_vector_size(mapping) != no_of_nodes) {
        IGRAPH_ERROR("Invalid mapping vector length",
                     IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 2));

    if (no_of_nodes > 0) {
        last = (long int) igraph_vector_max(mapping);
    }

    for (e = 0; e < no_of_edges; e++) {
        long int from = IGRAPH_FROM(graph, e);
        long int to = IGRAPH_TO(graph, e);

        long int nfrom = (long int) VECTOR(*mapping)[from];
        long int nto = (long int) VECTOR(*mapping)[to];

        igraph_vector_push_back(&edges, nfrom);
        igraph_vector_push_back(&edges, nto);

        if (nfrom > last) {
            last = nfrom;
        }
        if (nto   > last) {
            last = nto;
        }
    }

    no_new_vertices = last + 1;

    IGRAPH_CHECK(igraph_create(&res, &edges, (igraph_integer_t) no_new_vertices,
                               igraph_is_directed(graph)));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_FINALLY(igraph_destroy, &res);

    IGRAPH_I_ATTRIBUTE_DESTROY(&res);
    IGRAPH_I_ATTRIBUTE_COPY(&res, graph, /*graph=*/ 1,
                            /*vertex=*/ 0, /*edge=*/ 1);

    if (vattr) {
        long int i;
        igraph_vector_ptr_t merges;
        igraph_vector_t sizes;
        igraph_vector_t *vecs;

        vecs = igraph_Calloc(no_new_vertices, igraph_vector_t);
        if (!vecs) {
            IGRAPH_ERROR("Cannot combine attributes while contracting"
                         " vertices", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, vecs);
        IGRAPH_CHECK(igraph_vector_ptr_init(&merges, no_new_vertices));
        IGRAPH_FINALLY(igraph_i_simplify_free, &merges);
        IGRAPH_VECTOR_INIT_FINALLY(&sizes, no_new_vertices);

        for (i = 0; i < no_of_nodes; i++) {
            long int to = (long int) VECTOR(*mapping)[i];
            VECTOR(sizes)[to] += 1;
        }
        for (i = 0; i < no_new_vertices; i++) {
            igraph_vector_t *v = &vecs[i];
            IGRAPH_CHECK(igraph_vector_init(v, (long int) VECTOR(sizes)[i]));
            igraph_vector_clear(v);
            VECTOR(merges)[i] = v;
        }
        for (i = 0; i < no_of_nodes; i++) {
            long int to = (long int) VECTOR(*mapping)[i];
            igraph_vector_t *v = &vecs[to];
            igraph_vector_push_back(v, i);
        }

        IGRAPH_CHECK(igraph_i_attribute_combine_vertices(graph, &res,
                     &merges,
                     vertex_comb));

        igraph_vector_destroy(&sizes);
        igraph_i_simplify_free(&merges);
        igraph_free(vecs);
        IGRAPH_FINALLY_CLEAN(3);
    }

    IGRAPH_FINALLY_CLEAN(1);
    igraph_destroy(graph);
    *graph = res;

    return 0;
}

/* Create the transitive closure of a tree graph.
   This is fairly simple, we just collect all ancestors of a vertex
   using a depth-first search.
 */

int igraph_transitive_closure_dag(const igraph_t *graph,
                                  igraph_t *closure) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t deg;
    igraph_vector_t new_edges;
    igraph_vector_t ancestors;
    long int root;
    igraph_vector_t neighbors;
    igraph_stack_t path;
    igraph_vector_bool_t done;

    if (!igraph_is_directed(graph)) {
        IGRAPH_ERROR("Tree transitive closure of a directed graph",
                     IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&new_edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&deg, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&ancestors, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&neighbors, 0);
    IGRAPH_CHECK(igraph_stack_init(&path, 0));
    IGRAPH_FINALLY(igraph_stack_destroy, &path);
    IGRAPH_CHECK(igraph_vector_bool_init(&done, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &done);

    IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(),
                               IGRAPH_OUT, IGRAPH_LOOPS));

#define STAR (-1)

    for (root = 0; root < no_of_nodes; root++) {
        if (VECTOR(deg)[root] != 0) {
            continue;
        }
        IGRAPH_CHECK(igraph_stack_push(&path, root));

        while (!igraph_stack_empty(&path)) {
            long int node = (long int) igraph_stack_top(&path);
            if (node == STAR) {
                /* Leaving a node */
                long int j, n;
                igraph_stack_pop(&path);
                node = (long int) igraph_stack_pop(&path);
                if (!VECTOR(done)[node]) {
                    igraph_vector_pop_back(&ancestors);
                    VECTOR(done)[node] = 1;
                }
                n = igraph_vector_size(&ancestors);
                for (j = 0; j < n; j++) {
                    IGRAPH_CHECK(igraph_vector_push_back(&new_edges, node));
                    IGRAPH_CHECK(igraph_vector_push_back(&new_edges,
                                                         VECTOR(ancestors)[j]));
                }
            } else {
                /* Getting into a node */
                long int n, j;
                if (!VECTOR(done)[node]) {
                    IGRAPH_CHECK(igraph_vector_push_back(&ancestors, node));
                }
                IGRAPH_CHECK(igraph_neighbors(graph, &neighbors,
                                              (igraph_integer_t) node, IGRAPH_IN));
                n = igraph_vector_size(&neighbors);
                IGRAPH_CHECK(igraph_stack_push(&path, STAR));
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neighbors)[j];
                    IGRAPH_CHECK(igraph_stack_push(&path, nei));
                }
            }
        }
    }

#undef STAR

    igraph_vector_bool_destroy(&done);
    igraph_stack_destroy(&path);
    igraph_vector_destroy(&neighbors);
    igraph_vector_destroy(&ancestors);
    igraph_vector_destroy(&deg);
    IGRAPH_FINALLY_CLEAN(5);

    IGRAPH_CHECK(igraph_create(closure, &new_edges, (igraph_integer_t)no_of_nodes,
                               IGRAPH_DIRECTED));

    igraph_vector_destroy(&new_edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_diversity
 * Structural diversity index of the vertices
 *
 * This measure was defined in Nathan Eagle, Michael Macy and Rob
 * Claxton: Network Diversity and Economic Development, Science 328,
 * 1029--1031, 2010.
 *
 * </para><para>
 * It is simply the (normalized) Shannon entropy of the
 * incident edges' weights. D(i)=H(i)/log(k[i]), and
 * H(i) = -sum(p[i,j] log(p[i,j]), j=1..k[i]),
 * where p[i,j]=w[i,j]/sum(w[i,l], l=1..k[i]),  k[i] is the (total)
 * degree of vertex i, and w[i,j] is the weight of the edge(s) between
 * vertex i and j.
 * \param graph The input graph, edge directions are ignored.
 * \param weights The edge weights, in the order of the edge ids, must
 *    have appropriate length.
 * \param res An initialized vector, the results are stored here.
 * \param vids Vector with the vertex ids for which to calculate the
 *    measure.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear.
 *
 */

int igraph_diversity(igraph_t *graph, const igraph_vector_t *weights,
                     igraph_vector_t *res, const igraph_vs_t vids) {

    int no_of_nodes = igraph_vcount(graph);
    int no_of_edges = igraph_ecount(graph);
    igraph_vector_t incident;
    igraph_vit_t vit;
    igraph_real_t s, ent, w;
    int i, j, k;

    if (!weights) {
        IGRAPH_ERROR("Edge weights must be given", IGRAPH_EINVAL);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid edge weight vector length", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&incident, 10);

    if (igraph_vs_is_all(&vids)) {
        IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
        for (i = 0; i < no_of_nodes; i++) {
            s = ent = 0.0;
            IGRAPH_CHECK(igraph_incident(graph, &incident, i, /*mode=*/ IGRAPH_ALL));
            for (j = 0, k = (int) igraph_vector_size(&incident); j < k; j++) {
                w = VECTOR(*weights)[(long int)VECTOR(incident)[j]];
                s += w;
                ent += (w * log(w));
            }
            VECTOR(*res)[i] = (log(s) - ent / s) / log(k);
        }
    } else {
        IGRAPH_CHECK(igraph_vector_resize(res, 0));
        IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
        IGRAPH_FINALLY(igraph_vit_destroy, &vit);

        for (IGRAPH_VIT_RESET(vit), i = 0;
             !IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), i++) {
            long int v = IGRAPH_VIT_GET(vit);
            s = ent = 0.0;
            IGRAPH_CHECK(igraph_incident(graph, &incident, (igraph_integer_t) v,
                                         /*mode=*/ IGRAPH_ALL));
            for (j = 0, k = (int) igraph_vector_size(&incident); j < k; j++) {
                w = VECTOR(*weights)[(long int)VECTOR(incident)[j]];
                s += w;
                ent += (w * log(w));
            }
            IGRAPH_CHECK(igraph_vector_push_back(res, (log(s) - ent / s) / log(k)));
        }

        igraph_vit_destroy(&vit);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_destroy(&incident);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}


/* igraph_is_tree -- check if a graph is a tree */

/* count the number of vertices reachable from the root */
static int igraph_i_is_tree_visitor(igraph_integer_t root, const igraph_adjlist_t *al, igraph_integer_t *visited_count) {
    igraph_stack_int_t stack;
    igraph_vector_bool_t visited;
    long i;

    IGRAPH_CHECK(igraph_vector_bool_init(&visited, igraph_adjlist_size(al)));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &visited);

    IGRAPH_CHECK(igraph_stack_int_init(&stack, 0));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &stack);

    *visited_count = 0;

    /* push the root into the stack */
    IGRAPH_CHECK(igraph_stack_int_push(&stack, root));

    while (! igraph_stack_int_empty(&stack)) {
        igraph_integer_t u;
        igraph_vector_int_t *neighbors;
        long ncount;

        /* take a vertex from the stack, mark it as visited */
        u = igraph_stack_int_pop(&stack);
        if (IGRAPH_LIKELY(! VECTOR(visited)[u])) {
            VECTOR(visited)[u] = 1;
            *visited_count += 1;
        }

        /* register all its yet-unvisited neighbours for future processing */
        neighbors = igraph_adjlist_get(al, u);
        ncount = igraph_vector_int_size(neighbors);
        for (i = 0; i < ncount; ++i) {
            igraph_integer_t v = VECTOR(*neighbors)[i];
            if (! VECTOR(visited)[v]) {
                IGRAPH_CHECK(igraph_stack_int_push(&stack, v));
            }
        }
    }

    igraph_stack_int_destroy(&stack);
    igraph_vector_bool_destroy(&visited);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup structural
 * \function igraph_is_tree
 * \brief Decides whether the graph is a tree.
 *
 * An undirected graph is a tree if it is connected and has no cycles.
 * </para><para>
 *
 * In the directed case, a possible additional requirement is that all
 * edges are oriented away from a root (out-tree or arborescence) or all edges
 * are oriented towards a root (in-tree or anti-arborescence).
 * This test can be controlled using the \p mode parameter.
 * </para><para>
 *
 * By convention, the null graph (i.e. the graph with no vertices) is considered not to be a tree.
 *
 * \param graph The graph object to analyze.
 * \param res Pointer to a logical variable, the result will be stored
 *        here.
 * \param root If not \c NULL, the root node will be stored here. When \p mode
 *        is \c IGRAPH_ALL or the graph is undirected, any vertex can be the root
 *        and \p root is set to 0 (the first vertex). When \p mode is \c IGRAPH_OUT
 *        or \c IGRAPH_IN, the root is set to the vertex with zero in- or out-degree,
 *        respectively.
 * \param mode For a directed graph this specifies whether to test for an
 *        out-tree, an in-tree or ignore edge directions. The respective
 *        possible values are:
 *        \c IGRAPH_OUT, \c IGRAPH_IN, \c IGRAPH_ALL. This argument is
 *        ignored for undirected graphs.
 * \return Error code:
 *        \c IGRAPH_EINVAL: invalid mode argument.
 *
 * Time complexity: At most O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa igraph_is_weakly_connected()
 *
 * \example examples/simple/igraph_tree.c
 */

int igraph_is_tree(const igraph_t *graph, igraph_bool_t *res, igraph_integer_t *root, igraph_neimode_t mode) {
    igraph_adjlist_t al;
    igraph_integer_t iroot = 0;
    igraph_integer_t visited_count;
    igraph_integer_t vcount, ecount;

    vcount = igraph_vcount(graph);
    ecount = igraph_ecount(graph);

    /* A tree must have precisely vcount-1 edges. */
    /* By convention, the zero-vertex graph will not be considered a tree. */
    if (ecount != vcount - 1) {
        *res = 0;
        return IGRAPH_SUCCESS;
    }

    /* The single-vertex graph is a tree, provided it has no edges (checked in the previous if (..)) */
    if (vcount == 1) {
        *res = 1;
        if (root) {
            *root = 0;
        }
        return IGRAPH_SUCCESS;
    }

    /* For higher vertex counts we cannot short-circuit due to the possibility
     * of loops or multi-edges even when the edge count is correct. */

    /* Ignore mode for undirected graphs. */
    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    IGRAPH_CHECK(igraph_adjlist_init(graph, &al, mode));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);

    /* The main algorithm:
     * We find a root and check that all other vertices are reachable from it.
     * We have already checked the number of edges, so with the additional
     * reachability condition we can verify if the graph is a tree.
     *
     * For directed graphs, the root is the node with no incoming/outgoing
     * connections, depending on 'mode'. For undirected, it is arbitrary, so
     * we choose 0.
     */

    *res = 1; /* assume success */

    switch (mode) {
    case IGRAPH_ALL:
        iroot = 0;
        break;

    case IGRAPH_IN:
    case IGRAPH_OUT: {
        igraph_vector_t degree;
        igraph_integer_t i;

        IGRAPH_CHECK(igraph_vector_init(&degree, 0));
        IGRAPH_FINALLY(igraph_vector_destroy, &degree);

        IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), mode == IGRAPH_IN ? IGRAPH_OUT : IGRAPH_IN, /* loops = */ 1));

        for (i = 0; i < vcount; ++i)
            if (VECTOR(degree)[i] == 0) {
                break;
            }

        /* if no suitable root is found, the graph is not a tree */
        if (i == vcount) {
            *res = 0;
        } else {
            iroot = i;
        }

        igraph_vector_destroy(&degree);
        IGRAPH_FINALLY_CLEAN(1);
    }

    break;
    default:
        IGRAPH_ERROR("Invalid mode", IGRAPH_EINVMODE);
    }

    /* if no suitable root was found, skip visiting vertices */
    if (*res) {
        IGRAPH_CHECK(igraph_i_is_tree_visitor(iroot, &al, &visited_count));
        *res = visited_count == vcount;
    }

    if (root) {
        *root = iroot;
    }

    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
