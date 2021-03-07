/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2021 The igraph development team

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

#include "igraph_operators.h"

#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_memory.h"

#include "graph/attributes.h"
#include "core/interruption.h"

/**
 * Subgraph creation, old version: it copies the graph and then deletes
 * unneeded vertices.
 */
static int igraph_i_subgraph_copy_and_delete(const igraph_t *graph, igraph_t *res,
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
    remain = IGRAPH_CALLOC(no_of_nodes, char);
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

    IGRAPH_FREE(remain);
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
static int igraph_i_subgraph_create_from_scratch(const igraph_t *graph,
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

static int igraph_i_induced_subgraph_suggest_implementation(
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
    vremain = IGRAPH_CALLOC(no_of_nodes, char);
    if (vremain == 0) {
        IGRAPH_ERROR("subgraph_edges failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, vremain);
    eremain = IGRAPH_CALLOC(no_of_edges, char);
    if (eremain == 0) {
        IGRAPH_ERROR("subgraph_edges failed", IGRAPH_ENOMEM);
    }
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

    IGRAPH_FREE(eremain);
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

    IGRAPH_FREE(vremain);
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
