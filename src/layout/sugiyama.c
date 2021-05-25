/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_layout.h"
#include "igraph_centrality.h"
#include "igraph_components.h"
#include "igraph_constants.h"
#include "igraph_constructors.h"
#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_structural.h"
#include "igraph_types.h"

#include "internal/glpk_support.h"
#include "misc/feedback_arc_set.h"

#include "config.h"

#include <limits.h>

/* #define SUGIYAMA_DEBUG */

#ifdef _MSC_VER
/* MSVC does not support variadic macros */
#include <stdarg.h>
static void debug(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
#ifdef SUGIYAMA_DEBUG
    vfprintf(stderr, fmt, args);
#endif
    va_end(args);
}
#else
#ifdef SUGIYAMA_DEBUG
    #define debug(...) fprintf(stderr, __VA_ARGS__)
#else
    #define debug(...)
#endif
#endif

/* MSVC uses __forceinline instead of inline */
#ifdef _MSC_VER
    #define INLINE __forceinline
#else
    #define INLINE inline
#endif

/*
 * Implementation of the Sugiyama layout algorithm as described in:
 *
 * [1] K. Sugiyama, S. Tagawa and M. Toda, "Methods for Visual Understanding of
 * Hierarchical Systems". IEEE Transactions on Systems, Man and Cybernetics
 * 11(2):109-125, 1981.
 *
 * The layering (if not given in advance) is calculated by ... TODO
 *
 * [2] TODO
 *
 * The X coordinates of nodes within a layer are calculated using the method of
 * Brandes & Köpf:
 *
 * [3] U. Brandes and B. Köpf, "Fast and Simple Horizontal Coordinate
 * Assignment".  In: Lecture Notes in Computer Science 2265:31-44, 2002.
 *
 * Layer compaction is done according to:
 *
 * [4] N.S. Nikolov and A. Tarassov, "Graph layering by promotion of nodes".
 * Journal of Discrete Applied Mathematics, special issue: IV ALIO/EURO
 * workshop on applied combinatorial optimization, 154(5).
 *
 * The steps of the algorithm are as follows:
 *
 *   1. Cycle removal by finding an approximately minimal feedback arc set
 *      and reversing the direction of edges in the set.  Algorithms for
 *      finding minimal feedback arc sets are as follows:
 *
 *        - Find a cycle and find its minimum weight edge. Decrease the weight
 *          of all the edges by w. Remove those edges whose weight became zero.
 *          Repeat until there are no cycles. Re-introduce removed edges in
 *          decreasing order of weights, ensuring that no cycles are created.
 *
 *        - Order the vertices somehow and remove edges which point backwards
 *          in the ordering. Eades et al proposed the following procedure:
 *
 *            1. Iteratively remove sinks and prepend them to a vertex sequence
 *               s2.
 *
 *            2. Iteratively remove sources and append them to a vertex sequence
 *               s1.
 *
 *            3. Choose a vertex u s.t. the difference between the number of
 *               rightward arcs and the number of leftward arcs is the largest,
 *               remove u and append it to s1. Goto step 1 if there are still
 *               more vertices.
 *
 *            4. Concatenate s1 with s2.
 *
 *          This algorithm is known to produce feedback arc sets at most the
 *          size of m/2 - n/6, where m is the number of edges. Further
 *          improvements are possible in step 3 which bring down the size of
 *          the set to at most m/4 for cubic directed graphs, see Eades (1995).
 *
 *        - For undirected graphs, find a maximum weight spanning tree and
 *          remove all the edges not in the spanning tree. For directed graphs,
 *          find minimal cuts iteratively and remove edges pointing from A to
 *          B or from B to A in the cut, depending on which one is smaller. Yes,
 *          this is time-consuming.
 *
 *   2. Assigning vertices to layers according to [2].
 *
 *   3. Extracting weakly connected components. The remaining steps are
 *      executed for each component.
 *
 *   4. Compacting the layering using the method of [4]. TODO
 *      Steps 2-4 are performed only when no layering is given in advance.
 *
 *   5. Adding dummy nodes to ensure that each edge spans at most one layer
 *      only.
 *
 *   6. Finding an optimal ordering of vertices within a layer using the
 *      Sugiyama framework [1].
 *
 *   7. Assigning horizontal coordinates to each vertex using [3].
 *
 *   8. ???
 *
 *   9. Profit!
 */

/**
 * Data structure to store a layering of the graph.
 */
typedef struct {
    igraph_vector_ptr_t layers;
} igraph_i_layering_t;

/**
 * Initializes a layering.
 */
static int igraph_i_layering_init(igraph_i_layering_t* layering,
                                  const igraph_vector_t* membership) {
    long int i, n, num_layers;

    if (igraph_vector_size(membership) == 0) {
        num_layers = 0;
    } else {
        num_layers = (long int) igraph_vector_max(membership) + 1;
    }

    IGRAPH_CHECK(igraph_vector_ptr_init(&layering->layers, num_layers));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &layering->layers);

    for (i = 0; i < num_layers; i++) {
        igraph_vector_t* vec = IGRAPH_CALLOC(1, igraph_vector_t);
        IGRAPH_VECTOR_INIT_FINALLY(vec, 0);
        VECTOR(layering->layers)[i] = vec;
        IGRAPH_FINALLY_CLEAN(1);
    }
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&layering->layers, igraph_vector_destroy);

    n = igraph_vector_size(membership);
    for (i = 0; i < n; i++) {
        long int l = (long int) VECTOR(*membership)[i];
        igraph_vector_t* vec = VECTOR(layering->layers)[l];
        IGRAPH_CHECK(igraph_vector_push_back(vec, i));
    }

    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * Destroys a layering.
 */
static void igraph_i_layering_destroy(igraph_i_layering_t* layering) {
    igraph_vector_ptr_destroy_all(&layering->layers);
}

/**
 * Returns the number of layers in a layering.
 */
static int igraph_i_layering_num_layers(const igraph_i_layering_t* layering) {
    return (int) igraph_vector_ptr_size(&layering->layers);
}

/**
 * Returns the list of vertices in a given layer
 */
static igraph_vector_t* igraph_i_layering_get(const igraph_i_layering_t* layering,
                                       long int index) {
    return (igraph_vector_t*)VECTOR(layering->layers)[index];
}


/**
 * Forward declarations
 */

static int igraph_i_layout_sugiyama_place_nodes_vertically(const igraph_t* graph,
        const igraph_vector_t* weights, igraph_vector_t* membership);
static int igraph_i_layout_sugiyama_order_nodes_horizontally(const igraph_t* graph,
        igraph_matrix_t* layout, const igraph_i_layering_t* layering,
        long int maxiter);
static int igraph_i_layout_sugiyama_place_nodes_horizontally(const igraph_t* graph,
        igraph_matrix_t* layout, const igraph_i_layering_t* layering,
        igraph_real_t hgap, igraph_integer_t no_of_real_nodes);

/**
 * Calculated the median of four numbers (not necessarily sorted).
 */
static INLINE igraph_real_t igraph_i_median_4(igraph_real_t x1,
        igraph_real_t x2, igraph_real_t x3, igraph_real_t x4) {
    igraph_real_t arr[4] = { x1, x2, x3, x4 };
    igraph_vector_t vec;
    igraph_vector_view(&vec, arr, 4);
    igraph_vector_sort(&vec);
    return (arr[1] + arr[2]) / 2.0;
}


/**
 * \ingroup layout
 * \function igraph_layout_sugiyama
 * \brief Sugiyama layout algorithm for layered directed acyclic graphs.
 *
 * </para><para>
 * This layout algorithm is designed for directed acyclic graphs where each
 * vertex is assigned to a layer. Layers are indexed from zero, and vertices
 * of the same layer will be placed on the same horizontal line. The X coordinates
 * of vertices within each layer are decided by the heuristic proposed by
 * Sugiyama et al to minimize edge crossings.
 *
 * </para><para>
 * You can also try to lay out undirected graphs, graphs containing cycles, or
 * graphs without an a priori layered assignment with this algorithm. igraph
 * will try to eliminate cycles and assign vertices to layers, but there is no
 * guarantee on the quality of the layout in such cases.
 *
 * </para><para>
 * The Sugiyama layout may introduce "bends" on the edges in order to obtain a
 * visually more pleasing layout. This is achieved by adding dummy nodes to
 * edges spanning more than one layer. The resulting layout assigns coordinates
 * not only to the nodes of the original graph but also to the dummy nodes.
 * The layout algorithm will also return the extended graph with the dummy nodes.
 * An edge in the original graph may either be mapped to a single edge in the
 * extended graph or a \em path that starts and ends in the original
 * source and target vertex and passes through multiple dummy vertices. In
 * such cases, the user may also request the mapping of the edges of the extended
 * graph back to the edges of the original graph.
 *
 * </para><para>
 * For more details, see K. Sugiyama, S. Tagawa and M. Toda, "Methods for Visual
 * Understanding of Hierarchical Systems". IEEE Transactions on Systems, Man and
 * Cybernetics 11(2):109-125, 1981.
 *
 * \param graph Pointer to an initialized graph object.
 * \param res   Pointer to an initialized matrix object. This will contain
 *              the result and will be resized as needed. The first |V| rows
 *              of the layout will contain the coordinates of the original graph,
 *              the remaining rows contain the positions of the dummy nodes.
 *              Therefore, you can use the result both with \p graph or with
 *              \p extended_graph.
 * \param extended_graph Pointer to an uninitialized graph object or \c NULL.
 *                       The extended graph with the added dummy nodes will be
 *                       returned here. In this graph, each edge points downwards
 *                       to lower layers, spans exactly one layer and the first
 *                       |V| vertices coincide with the vertices of the
 *                       original graph.
 * \param extd_to_orig_eids Pointer to a vector or \c NULL. If not \c NULL, the
 *                          mapping from the edge IDs of the extended graph back
 *                          to the edge IDs of the original graph will be stored
 *                          here.
 * \param layers  The layer index for each vertex or \c NULL if the layers should
 *                be determined automatically by igraph.
 * \param hgap  The preferred minimum horizontal gap between vertices in the same
 *              layer.
 * \param vgap  The distance between layers.
 * \param maxiter Maximum number of iterations in the crossing minimization stage.
 *                100 is a reasonable default; if you feel that you have too
 *                many edge crossings, increase this.
 * \param weights Weights of the edges. These are used only if the graph contains
 *                cycles; igraph will tend to reverse edges with smaller
 *                weights when breaking the cycles.
 */
int igraph_layout_sugiyama(const igraph_t *graph, igraph_matrix_t *res,
                           igraph_t *extd_graph, igraph_vector_t *extd_to_orig_eids,
                           const igraph_vector_t* layers, igraph_real_t hgap, igraph_real_t vgap,
                           long int maxiter, const igraph_vector_t *weights) {
    long int i, j, k, l, m, nei;
    long int no_of_nodes = (long int)igraph_vcount(graph);
    long int comp_idx;
    long int next_extd_vertex_id = no_of_nodes;
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_integer_t no_of_components;  /* number of components of the original graph */
    igraph_vector_t membership;         /* components of the original graph */
    igraph_vector_t extd_edgelist;   /* edge list of the extended graph */
    igraph_vector_t layers_own;  /* layer indices after having eliminated empty layers */
    igraph_real_t dx = 0, dx2 = 0; /* displacement of the current component on the X axis */
    igraph_vector_t layer_to_y; /* mapping from layer indices to final Y coordinates */

    if (layers && igraph_vector_size(layers) != no_of_nodes) {
        IGRAPH_ERROR("layer vector too short or too long", IGRAPH_EINVAL);
    }

    if (extd_graph != 0) {
        IGRAPH_VECTOR_INIT_FINALLY(&extd_edgelist, 0);
        if (extd_to_orig_eids != 0) {
            igraph_vector_clear(extd_to_orig_eids);
        }
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));
    IGRAPH_VECTOR_INIT_FINALLY(&membership, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&layer_to_y, 0);

    /* 1. Find a feedback arc set if we don't have a layering yet. If we do have
     *    a layering, we can leave all the edges as is as they will be re-oriented
     *    to point downwards only anyway. */
    if (layers == 0) {
        IGRAPH_VECTOR_INIT_FINALLY(&layers_own, no_of_nodes);
        IGRAPH_CHECK(igraph_i_layout_sugiyama_place_nodes_vertically(
                         graph, weights, &layers_own));
    } else {
        IGRAPH_CHECK(igraph_vector_copy(&layers_own, layers));
        IGRAPH_FINALLY(igraph_vector_destroy, &layers_own);
    }

    /* Normalize layering, eliminate empty layers */
    if (no_of_nodes > 0) {
        igraph_vector_t inds;
        IGRAPH_VECTOR_INIT_FINALLY(&inds, 0);
        IGRAPH_CHECK((int) igraph_vector_qsort_ind(&layers_own, &inds, 0));
        j = -1; dx = VECTOR(layers_own)[(long int)VECTOR(inds)[0]] - 1;
        for (i = 0; i < no_of_nodes; i++) {
            k = (long int)VECTOR(inds)[i];
            if (VECTOR(layers_own)[k] > dx) {
                /* New layer starts here */
                dx = VECTOR(layers_own)[k];
                j++;
                IGRAPH_CHECK(igraph_vector_push_back(&layer_to_y, dx * vgap));
            }
            VECTOR(layers_own)[k] = j;
        }
        igraph_vector_destroy(&inds);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* 2. Find the connected components. */
    IGRAPH_CHECK(igraph_clusters(graph, &membership, 0, &no_of_components,
                                 IGRAPH_WEAK));

    /* 3. For each component... */
    dx = 0;
    for (comp_idx = 0; comp_idx < no_of_components; comp_idx++) {
        /* Extract the edges of the comp_idx'th component and add dummy nodes for edges
         * spanning more than one layer. */
        long int component_size, next_new_vertex_id;
        igraph_vector_t old2new_vertex_ids;
        igraph_vector_t new2old_vertex_ids;
        igraph_vector_t new_layers;
        igraph_vector_t edgelist;
        igraph_vector_t neis;

        IGRAPH_VECTOR_INIT_FINALLY(&edgelist, 0);
        IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
        IGRAPH_VECTOR_INIT_FINALLY(&new2old_vertex_ids, no_of_nodes);
        IGRAPH_VECTOR_INIT_FINALLY(&old2new_vertex_ids, no_of_nodes);
        IGRAPH_VECTOR_INIT_FINALLY(&new_layers, 0);

        igraph_vector_fill(&old2new_vertex_ids, -1);

        /* Construct a mapping from the old vertex ids to the new ones */
        for (i = 0, next_new_vertex_id = 0; i < no_of_nodes; i++) {
            if (VECTOR(membership)[i] == comp_idx) {
                IGRAPH_CHECK(igraph_vector_push_back(&new_layers, VECTOR(layers_own)[i]));
                VECTOR(new2old_vertex_ids)[next_new_vertex_id] = i;
                VECTOR(old2new_vertex_ids)[i] = next_new_vertex_id;
                next_new_vertex_id++;
            }
        }
        component_size = next_new_vertex_id;

        /* Construct a proper layering of the component in new_graph where each edge
         * points downwards and spans exactly one layer. */
        for (i = 0; i < no_of_nodes; i++) {
            if (VECTOR(membership)[i] != comp_idx) {
                continue;
            }

            /* Okay, this vertex is in the component we are considering.
             * Add the neighbors of this vertex, excluding loops */
            IGRAPH_CHECK(igraph_incident(graph, &neis, (igraph_integer_t) i,
                                         IGRAPH_OUT));
            j = igraph_vector_size(&neis);
            for (k = 0; k < j; k++) {
                long int eid = (long int) VECTOR(neis)[k];
                if (directed) {
                    nei = IGRAPH_TO(graph, eid);
                } else {
                    nei = IGRAPH_OTHER(graph, eid, i);
                    if (nei < i) { /* to avoid considering edges twice */
                        continue;
                    }
                }
                if (VECTOR(layers_own)[i] == VECTOR(layers_own)[nei]) {
                    /* Edge goes within the same layer, we don't need this in the
                     * layered graph, but we need it in the extended graph */
                    if (extd_graph != 0) {
                        IGRAPH_CHECK(igraph_vector_push_back(&extd_edgelist, i));
                        IGRAPH_CHECK(igraph_vector_push_back(&extd_edgelist, nei));
                        if (extd_to_orig_eids != 0) {
                            IGRAPH_CHECK(igraph_vector_push_back(extd_to_orig_eids, eid));
                        }
                    }
                } else if (VECTOR(layers_own)[i] > VECTOR(layers_own)[nei]) {
                    /* Edge goes upwards, we have to flip it */
                    IGRAPH_CHECK(igraph_vector_push_back(&edgelist,
                                                         VECTOR(old2new_vertex_ids)[nei]));
                    for (l = (long int) VECTOR(layers_own)[nei] + 1;
                         l < VECTOR(layers_own)[i]; l++) {
                        IGRAPH_CHECK(igraph_vector_push_back(&new_layers, l));
                        IGRAPH_CHECK(igraph_vector_push_back(&edgelist, next_new_vertex_id));
                        IGRAPH_CHECK(igraph_vector_push_back(&edgelist, next_new_vertex_id++));
                    }
                    IGRAPH_CHECK(igraph_vector_push_back(&edgelist,
                                                         VECTOR(old2new_vertex_ids)[i]));
                    /* Also add the edge to the extended graph if needed, but this time
                     * with the proper orientation */
                    if (extd_graph != 0) {
                        IGRAPH_CHECK(igraph_vector_push_back(&extd_edgelist, i));
                        next_extd_vertex_id += VECTOR(layers_own)[i] - VECTOR(layers_own)[nei] - 1;
                        for (l = (long int) VECTOR(layers_own)[i] - 1, m = 1;
                             l > VECTOR(layers_own)[nei]; l--, m++) {
                            IGRAPH_CHECK(igraph_vector_push_back(&extd_edgelist, next_extd_vertex_id - m));
                            IGRAPH_CHECK(igraph_vector_push_back(&extd_edgelist, next_extd_vertex_id - m));
                            if (extd_to_orig_eids != 0) {
                                IGRAPH_CHECK(igraph_vector_push_back(extd_to_orig_eids, eid));
                            }
                        }
                        IGRAPH_CHECK(igraph_vector_push_back(&extd_edgelist, nei));
                        if (extd_to_orig_eids != 0) {
                            IGRAPH_CHECK(igraph_vector_push_back(extd_to_orig_eids, eid));
                        }
                    }
                } else {
                    /* Edge goes downwards */
                    IGRAPH_CHECK(igraph_vector_push_back(&edgelist,
                                                         VECTOR(old2new_vertex_ids)[i]));
                    for (l = (long int) VECTOR(layers_own)[i] + 1;
                         l < VECTOR(layers_own)[nei]; l++) {
                        IGRAPH_CHECK(igraph_vector_push_back(&new_layers, l));
                        IGRAPH_CHECK(igraph_vector_push_back(&edgelist, next_new_vertex_id));
                        IGRAPH_CHECK(igraph_vector_push_back(&edgelist, next_new_vertex_id++));
                    }
                    IGRAPH_CHECK(igraph_vector_push_back(&edgelist,
                                                         VECTOR(old2new_vertex_ids)[nei]));
                    /* Also add the edge to the extended graph */
                    if (extd_graph != 0) {
                        IGRAPH_CHECK(igraph_vector_push_back(&extd_edgelist, i));
                        for (l = (long int) VECTOR(layers_own)[i] + 1;
                             l < VECTOR(layers_own)[nei]; l++) {
                            IGRAPH_CHECK(igraph_vector_push_back(&extd_edgelist, next_extd_vertex_id));
                            IGRAPH_CHECK(igraph_vector_push_back(&extd_edgelist, next_extd_vertex_id++));
                            if (extd_to_orig_eids != 0) {
                                IGRAPH_CHECK(igraph_vector_push_back(extd_to_orig_eids, eid));
                            }
                        }
                        IGRAPH_CHECK(igraph_vector_push_back(&extd_edgelist, nei));
                        if (extd_to_orig_eids != 0) {
                            IGRAPH_CHECK(igraph_vector_push_back(extd_to_orig_eids, eid));
                        }
                    }
                }
            }
        }

        /* At this point, we have the subgraph with the dummy nodes and
         * edges, so we can run Sugiyama's algorithm on it. */
        {
            igraph_matrix_t layout;
            igraph_i_layering_t layering;
            igraph_t subgraph;

            IGRAPH_CHECK(igraph_matrix_init(&layout, next_new_vertex_id, 2));
            IGRAPH_FINALLY(igraph_matrix_destroy, &layout);
            IGRAPH_CHECK(igraph_create(&subgraph, &edgelist, (igraph_integer_t)
                                       next_new_vertex_id, 1));
            IGRAPH_FINALLY(igraph_destroy, &subgraph);

            /*
            igraph_vector_print(&edgelist);
            igraph_vector_print(&new_layers);
            */

            /* Assign the vertical coordinates */
            for (i = 0; i < next_new_vertex_id; i++) {
                MATRIX(layout, i, 1) = VECTOR(new_layers)[i];
            }

            /* Create a layering */
            IGRAPH_CHECK(igraph_i_layering_init(&layering, &new_layers));
            IGRAPH_FINALLY(igraph_i_layering_destroy, &layering);

            /* Find the order in which the nodes within a layer should be placed */
            IGRAPH_CHECK(igraph_i_layout_sugiyama_order_nodes_horizontally(&subgraph, &layout,
                         &layering, maxiter));

            /* Assign the horizontal coordinates. This is according to the algorithm
             * of Brandes & Köpf */
            IGRAPH_CHECK(igraph_i_layout_sugiyama_place_nodes_horizontally(&subgraph, &layout,
                         &layering, hgap, (igraph_integer_t) component_size));

            /* Re-assign rows into the result matrix, and at the same time, */
            /* adjust dx so that the next component does not overlap this one */
            j = next_new_vertex_id - component_size;
            k = igraph_matrix_nrow(res);
            IGRAPH_CHECK(igraph_matrix_add_rows(res, j));
            dx2 = dx;
            for (i = 0; i < component_size; i++) {
                l = (long int)VECTOR(new2old_vertex_ids)[i];
                MATRIX(*res, l, 0) = MATRIX(layout, i, 0) + dx;
                MATRIX(*res, l, 1) = VECTOR(layer_to_y)[(long)MATRIX(layout, i, 1)];
                if (dx2 < MATRIX(*res, l, 0)) {
                    dx2 = MATRIX(*res, l, 0);
                }
            }
            for (i = component_size; i < next_new_vertex_id; i++) {
                MATRIX(*res, k, 0) = MATRIX(layout, i, 0) + dx;
                MATRIX(*res, k, 1) = VECTOR(layer_to_y)[(long)MATRIX(layout, i, 1)];
                if (dx2 < MATRIX(*res, k, 0)) {
                    dx2 = MATRIX(*res, k, 0);
                }
                k++;
            }
            dx = dx2 + hgap;

            igraph_destroy(&subgraph);
            igraph_i_layering_destroy(&layering);
            igraph_matrix_destroy(&layout);
            IGRAPH_FINALLY_CLEAN(3);
        }

        igraph_vector_destroy(&new_layers);
        igraph_vector_destroy(&old2new_vertex_ids);
        igraph_vector_destroy(&new2old_vertex_ids);
        igraph_vector_destroy(&edgelist);
        igraph_vector_destroy(&neis);
        IGRAPH_FINALLY_CLEAN(5);
    }

    igraph_vector_destroy(&layers_own);
    igraph_vector_destroy(&layer_to_y);
    igraph_vector_destroy(&membership);
    IGRAPH_FINALLY_CLEAN(3);

    if (extd_graph != 0) {
        IGRAPH_CHECK(igraph_create(extd_graph, &extd_edgelist, (igraph_integer_t)
                                   next_extd_vertex_id, igraph_is_directed(graph)));
        igraph_vector_destroy(&extd_edgelist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static int igraph_i_layout_sugiyama_place_nodes_vertically(const igraph_t* graph,
        const igraph_vector_t* weights, igraph_vector_t* membership) {
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    IGRAPH_CHECK(igraph_vector_resize(membership, no_of_nodes));

    if (no_of_edges == 0) {
        igraph_vector_fill(membership, 0);
        return IGRAPH_SUCCESS;
    }

#ifdef HAVE_GLPK
    if (igraph_is_directed(graph) && no_of_nodes <= 1000) {
        /* Network simplex algorithm of Gansner et al, using the original linear
         * programming formulation */
        long int i, j;
        igraph_vector_t outdegs, indegs, feedback_edges;
        glp_prob *ip;
        glp_smcp parm;

        /* Allocate storage and create the problem */
        ip = glp_create_prob();
        IGRAPH_FINALLY(glp_delete_prob, ip);
        IGRAPH_VECTOR_INIT_FINALLY(&feedback_edges, 0);
        IGRAPH_VECTOR_INIT_FINALLY(&outdegs, no_of_nodes);
        IGRAPH_VECTOR_INIT_FINALLY(&indegs, no_of_nodes);

        /* Find an approximate feedback edge set */
        IGRAPH_CHECK(igraph_i_feedback_arc_set_eades(graph, &feedback_edges, weights, 0));
        igraph_vector_sort(&feedback_edges);

        /* Calculate in- and out-strengths for the remaining edges */
        IGRAPH_CHECK(igraph_strength(graph, &indegs, igraph_vss_all(),
                                     IGRAPH_IN, 1, weights));
        IGRAPH_CHECK(igraph_strength(graph, &outdegs, igraph_vss_all(),
                                     IGRAPH_IN, 1, weights));
        j = igraph_vector_size(&feedback_edges);
        for (i = 0; i < j; i++) {
            long int eid = (long int) VECTOR(feedback_edges)[i];
            long int from = IGRAPH_FROM(graph, eid);
            long int to = IGRAPH_TO(graph, eid);
            VECTOR(outdegs)[from] -= weights ? VECTOR(*weights)[eid] : 1;
            VECTOR(indegs)[to] -= weights ? VECTOR(*weights)[eid] : 1;
        }

        /* Configure GLPK */
        glp_term_out(GLP_OFF);
        glp_init_smcp(&parm);
        parm.msg_lev = GLP_MSG_OFF;
        parm.presolve = GLP_OFF;

        /* Set up variables and objective function coefficients */
        glp_set_obj_dir(ip, GLP_MIN);
        glp_add_cols(ip, (int) no_of_nodes);
        IGRAPH_CHECK(igraph_vector_sub(&outdegs, &indegs));
        for (i = 1; i <= no_of_nodes; i++) {
            glp_set_col_kind(ip, (int) i, GLP_IV);
            glp_set_col_bnds(ip, (int) i, GLP_LO, 0.0, 0.0);
            glp_set_obj_coef(ip, (int) i, VECTOR(outdegs)[i - 1]);
        }
        igraph_vector_destroy(&indegs);
        igraph_vector_destroy(&outdegs);
        IGRAPH_FINALLY_CLEAN(2);

        /* Add constraints */
        glp_add_rows(ip, (int) no_of_edges);
        IGRAPH_CHECK(igraph_vector_push_back(&feedback_edges, -1));
        j = 0;
        for (i = 0; i < no_of_edges; i++) {
            int ind[3];
            double val[3] = {0, -1, 1};
            ind[1] = IGRAPH_FROM(graph, i) + 1;
            ind[2] = IGRAPH_TO(graph, i) + 1;

            if (ind[1] == ind[2]) {
                if (VECTOR(feedback_edges)[j] == i) {
                    j++;
                }
                continue;
            }

            if (VECTOR(feedback_edges)[j] == i) {
                /* This is a feedback edge, add it reversed */
                glp_set_row_bnds(ip, (int) i + 1, GLP_UP, -1, -1);
                j++;
            } else {
                glp_set_row_bnds(ip, (int) i + 1, GLP_LO, 1, 1);
            }
            glp_set_mat_row(ip, (int) i + 1, 2, ind, val);
        }

        /* Solve the problem */
        IGRAPH_GLPK_CHECK(glp_simplex(ip, &parm),
                          "Vertical arrangement step using IP failed");

        /* The problem is totally unimodular, therefore the output of the simplex
         * solver can be converted to an integer solution easily */
        for (i = 0; i < no_of_nodes; i++) {
            VECTOR(*membership)[i] = floor(glp_get_col_prim(ip, (int) i + 1));
        }

        glp_delete_prob(ip);
        igraph_vector_destroy(&feedback_edges);
        IGRAPH_FINALLY_CLEAN(2);
    } else if (igraph_is_directed(graph)) {
        IGRAPH_CHECK(igraph_i_feedback_arc_set_eades(graph, 0, weights, membership));
    } else {
        IGRAPH_CHECK(igraph_i_feedback_arc_set_undirected(graph, 0, weights, membership));
    }
#else
    if (igraph_is_directed(graph)) {
        IGRAPH_CHECK(igraph_i_feedback_arc_set_eades(graph, 0, weights, membership));
    } else {
        IGRAPH_CHECK(igraph_i_feedback_arc_set_undirected(graph, 0, weights, membership));
    }
#endif

    return IGRAPH_SUCCESS;
}

static int igraph_i_layout_sugiyama_calculate_barycenters(const igraph_t* graph,
        const igraph_i_layering_t* layering, long int layer_index,
        igraph_neimode_t direction, const igraph_matrix_t* layout,
        igraph_vector_t* barycenters) {
    long int i, j, m, n;
    igraph_vector_t* layer_members = igraph_i_layering_get(layering, layer_index);
    igraph_vector_t neis;

    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

    n = igraph_vector_size(layer_members);
    IGRAPH_CHECK(igraph_vector_resize(barycenters, n));
    igraph_vector_null(barycenters);

    for (i = 0; i < n; i++) {
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t)
                                      VECTOR(*layer_members)[i], direction));
        m = igraph_vector_size(&neis);
        if (m == 0) {
            /* No neighbors in this direction. Just use the current X coordinate */
            VECTOR(*barycenters)[i] = MATRIX(*layout, i, 0);
        } else {
            for (j = 0; j < m; j++) {
                VECTOR(*barycenters)[i] += MATRIX(*layout, (long)VECTOR(neis)[j], 0);
            }
            VECTOR(*barycenters)[i] /= m;
        }
    }

    igraph_vector_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * Given a properly layered graph where each edge points downwards and spans
 * exactly one layer, arranges the nodes in each layer horizontally in a way
 * that strives to minimize edge crossings.
 */
static int igraph_i_layout_sugiyama_order_nodes_horizontally(const igraph_t* graph,
        igraph_matrix_t* layout, const igraph_i_layering_t* layering,
        long int maxiter) {
    long int i, n, nei;
    long int no_of_vertices = igraph_vcount(graph);
    long int no_of_layers = igraph_i_layering_num_layers(layering);
    long int iter, layer_index;
    igraph_vector_t* layer_members;
    igraph_vector_t neis, barycenters, sort_indices;
    igraph_bool_t changed;

    /* The first column of the matrix will serve as the ordering */
    /* Start with a first-seen ordering within each layer */
    {
        long int *xs = IGRAPH_CALLOC(no_of_layers, long int);
        if (xs == 0) {
            IGRAPH_ERROR("cannot order nodes horizontally", IGRAPH_ENOMEM);
        }
        for (i = 0; i < no_of_vertices; i++) {
            MATRIX(*layout, i, 0) = xs[(long int)MATRIX(*layout, i, 1)]++;
        }
        free(xs);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&barycenters, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&sort_indices, 0);

    /* Start the effective part of the Sugiyama algorithm */
    iter = 0; changed = 1;
    while (changed && iter < maxiter) {
        changed = 0;

        /* Phase 1 */

        /* Moving downwards and sorting by upper barycenters */
        for (layer_index = 1; layer_index < no_of_layers; layer_index++) {
            layer_members = igraph_i_layering_get(layering, layer_index);
            n = igraph_vector_size(layer_members);

            igraph_i_layout_sugiyama_calculate_barycenters(graph,
                    layering, layer_index, IGRAPH_IN, layout, &barycenters);

#ifdef SUGIYAMA_DEBUG
            printf("Layer %ld, aligning to upper barycenters\n", layer_index);
            printf("Vertices: "); igraph_vector_print(layer_members);
            printf("Barycenters: "); igraph_vector_print(&barycenters);
#endif
            IGRAPH_CHECK((int) igraph_vector_qsort_ind(&barycenters,
                         &sort_indices, 0));
            for (i = 0; i < n; i++) {
                nei = (long)VECTOR(*layer_members)[(long)VECTOR(sort_indices)[i]];
                VECTOR(barycenters)[i] = nei;
                MATRIX(*layout, nei, 0) = i;
            }
            if (!igraph_vector_all_e(layer_members, &barycenters)) {
                IGRAPH_CHECK(igraph_vector_update(layer_members, &barycenters));
#ifdef SUGIYAMA_DEBUG
                printf("New vertex order: "); igraph_vector_print(layer_members);
#endif
                changed = 1;
            } else {
#ifdef SUGIYAMA_DEBUG
                printf("Order did not change.\n");
#endif
            }
        }

        /* Moving upwards and sorting by lower barycenters */
        for (layer_index = no_of_layers - 2; layer_index >= 0; layer_index--) {
            layer_members = igraph_i_layering_get(layering, layer_index);
            n = igraph_vector_size(layer_members);

            igraph_i_layout_sugiyama_calculate_barycenters(graph,
                    layering, layer_index, IGRAPH_OUT, layout, &barycenters);

#ifdef SUGIYAMA_DEBUG
            printf("Layer %ld, aligning to lower barycenters\n", layer_index);
            printf("Vertices: "); igraph_vector_print(layer_members);
            printf("Barycenters: "); igraph_vector_print(&barycenters);
#endif

            IGRAPH_CHECK((int) igraph_vector_qsort_ind(&barycenters,
                         &sort_indices, 0));
            for (i = 0; i < n; i++) {
                nei = (long)VECTOR(*layer_members)[(long)VECTOR(sort_indices)[i]];
                VECTOR(barycenters)[i] = nei;
                MATRIX(*layout, nei, 0) = i;
            }
            if (!igraph_vector_all_e(layer_members, &barycenters)) {
                IGRAPH_CHECK(igraph_vector_update(layer_members, &barycenters));
#ifdef SUGIYAMA_DEBUG
                printf("New vertex order: "); igraph_vector_print(layer_members);
#endif
                changed = 1;
            } else {
#ifdef SUGIYAMA_DEBUG
                printf("Order did not change.\n");
#endif
            }
        }

#ifdef SUGIYAMA_DEBUG
        printf("==== Finished iteration %ld\n", iter);
#endif

        iter++;
    }

    igraph_vector_destroy(&barycenters);
    igraph_vector_destroy(&neis);
    igraph_vector_destroy(&sort_indices);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

#define IS_DUMMY(v) ((v >= no_of_real_nodes))
#define IS_INNER_SEGMENT(u, v) (IS_DUMMY(u) && IS_DUMMY(v))
#define X_POS(v) (MATRIX(*layout, v, 0))

static int igraph_i_layout_sugiyama_vertical_alignment(const igraph_t* graph,
        const igraph_i_layering_t* layering, const igraph_matrix_t* layout,
        const igraph_vector_bool_t* ignored_edges,
        igraph_bool_t reverse, igraph_bool_t align_right,
        igraph_vector_t* roots, igraph_vector_t* align);
static int igraph_i_layout_sugiyama_horizontal_compaction(const igraph_t* graph,
        const igraph_vector_t* vertex_to_the_left,
        const igraph_vector_t* roots, const igraph_vector_t* align,
        igraph_real_t hgap, igraph_vector_t* xs);
static int igraph_i_layout_sugiyama_horizontal_compaction_place_block(long int v,
        const igraph_vector_t* vertex_to_the_left,
        const igraph_vector_t* roots, const igraph_vector_t* align,
        igraph_vector_t* sinks, igraph_vector_t* shifts,
        igraph_real_t hgap, igraph_vector_t* xs);

static int igraph_i_layout_sugiyama_place_nodes_horizontally(const igraph_t* graph,
        igraph_matrix_t* layout, const igraph_i_layering_t* layering,
        igraph_real_t hgap, igraph_integer_t no_of_real_nodes) {

    long int i, j, k, l, n;
    long int no_of_layers = igraph_i_layering_num_layers(layering);
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_vector_t neis1, neis2;
    igraph_vector_t xs[4];
    igraph_vector_t roots, align;
    igraph_vector_t vertex_to_the_left;
    igraph_vector_bool_t ignored_edges;

    /*
    {
      igraph_vector_t edgelist;
      IGRAPH_VECTOR_INIT_FINALLY(&edgelist, 0);
      IGRAPH_CHECK(igraph_get_edgelist(graph, &edgelist, 0));
      igraph_vector_print(&edgelist);
      igraph_vector_destroy(&edgelist);
      IGRAPH_FINALLY_CLEAN(1);

      for (i = 0; i < no_of_layers; i++) {
        igraph_vector_t* layer = igraph_i_layering_get(layering, i);
        igraph_vector_print(layer);
      }
    }
    */

    IGRAPH_CHECK(igraph_vector_bool_init(&ignored_edges, no_of_edges));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &ignored_edges);

    IGRAPH_VECTOR_INIT_FINALLY(&vertex_to_the_left, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&neis1, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&neis2, 0);

    /* First, find all type 1 conflicts and mark one of the edges participating
     * in the conflict as being ignored. If one of the edges in the conflict
     * is a non-inner segment and the other is an inner segment, we ignore the
     * non-inner segment as we want to keep inner segments vertical.
     */
    for (i = 0; i < no_of_layers - 1; i++) {
        igraph_vector_t* vertices = igraph_i_layering_get(layering, i);
        n = igraph_vector_size(vertices);

        /* Find all the edges from this layer to the next */
        igraph_vector_clear(&neis1);
        for (j = 0; j < n; j++) {
            IGRAPH_CHECK(igraph_neighbors(graph, &neis2, (igraph_integer_t)
                                          VECTOR(*vertices)[j], IGRAPH_OUT));
            IGRAPH_CHECK(igraph_vector_append(&neis1, &neis2));
        }

        /* Consider all pairs of edges and check whether they are in a type 1
         * conflict */
        n = igraph_vector_size(&neis1);
        for (j = 0; j < n; j++) {
            long int u = IGRAPH_FROM(graph, j);
            long int v = IGRAPH_TO(graph, j);
            igraph_bool_t j_inner = IS_INNER_SEGMENT(u, v);
            igraph_bool_t crossing;

            for (k = j + 1; k < n; k++) {
                long int w = IGRAPH_FROM(graph, k);
                long int x = IGRAPH_TO(graph, k);
                if (IS_INNER_SEGMENT(w, x) == j_inner) {
                    continue;
                }
                /* Do the u --> v and w --> x edges cross? */
                crossing = (u == w || v == x);
                if (!crossing) {
                    if (X_POS(u) <= X_POS(w)) {
                        crossing = X_POS(v) >= X_POS(x);
                    } else {
                        crossing = X_POS(v) <= X_POS(x);
                    }
                }
                if (crossing) {
                    if (j_inner) {
                        VECTOR(ignored_edges)[k] = 1;
                    } else {
                        VECTOR(ignored_edges)[j] = 1;
                    }
                }
            }
        }
    }

    igraph_vector_destroy(&neis1);
    igraph_vector_destroy(&neis2);
    IGRAPH_FINALLY_CLEAN(2);

    /*
     * Prepare vertex_to_the_left where the ith element stores
     * the index of the vertex to the left of vertex i, or i itself if the
     * vertex is the leftmost vertex in a layer.
     */
    for (i = 0; i < no_of_layers; i++) {
        igraph_vector_t* vertices = igraph_i_layering_get(layering, i);
        n = igraph_vector_size(vertices);
        if (n == 0) {
            continue;
        }

        k = l = (long int)VECTOR(*vertices)[0];
        VECTOR(vertex_to_the_left)[k] = k;
        for (j = 1; j < n; j++) {
            k = (long int)VECTOR(*vertices)[j];
            VECTOR(vertex_to_the_left)[k] = l;
            l = k;
        }
    }

    /* Type 1 conflicts found, ignored edges chosen, vertex_to_the_left
     * prepared. Run vertical alignment for all four combinations */
    for (i = 0; i < 4; i++) {
        IGRAPH_VECTOR_INIT_FINALLY(&xs[i], no_of_nodes);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&roots, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&align, no_of_nodes);

    for (i = 0; i < 4; i++) {
        IGRAPH_CHECK(igraph_i_layout_sugiyama_vertical_alignment(graph,
                     layering, layout, &ignored_edges,
                     /* reverse = */ (igraph_bool_t) i / 2, /* align_right = */ i % 2,
                     &roots, &align));
        IGRAPH_CHECK(igraph_i_layout_sugiyama_horizontal_compaction(graph,
                     &vertex_to_the_left, &roots, &align, hgap, &xs[i]));
    }

    {
        igraph_real_t width, min_width, mins[4], maxs[4], diff;
        /* Find the alignment with the minimum width */
        min_width = IGRAPH_INFINITY; j = 0;
        for (i = 0; i < 4; i++) {
            mins[i] = igraph_vector_min(&xs[i]);
            maxs[i] = igraph_vector_max(&xs[i]);
            width = maxs[i] - mins[i];
            if (width < min_width) {
                min_width = width;
                j = i;
            }
        }

        /* Leftmost alignments: align them s.t. the min X coordinate is equal to
         * the minimum X coordinate of the alignment with the smallest width.
         * Rightmost alignments: align them s.t. the max X coordinate is equal to
         * the max X coordinate of the alignment with the smallest width.
         */
        for (i = 0; i < 4; i++) {
            if (j == i) {
                continue;
            }
            if (i % 2 == 0) {
                /* Leftmost alignment */
                diff = mins[j] - mins[i];
            } else {
                /* Rightmost alignment */
                diff = maxs[j] - maxs[i];
            }
            igraph_vector_add_constant(&xs[i], diff);
        }
    }

    /* For every vertex, find the median of the X coordinates in the four
     * alignments */
    for (i = 0; i < no_of_nodes; i++) {
        X_POS(i) = igraph_i_median_4(VECTOR(xs[0])[i], VECTOR(xs[1])[i],
                                     VECTOR(xs[2])[i], VECTOR(xs[3])[i]);
    }

    igraph_vector_destroy(&roots);
    igraph_vector_destroy(&align);
    IGRAPH_FINALLY_CLEAN(2);

    for (i = 0; i < 4; i++) {
        igraph_vector_destroy(&xs[i]);
    }
    IGRAPH_FINALLY_CLEAN(4);

    igraph_vector_destroy(&vertex_to_the_left);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_vector_bool_destroy(&ignored_edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static int igraph_i_layout_sugiyama_vertical_alignment(const igraph_t* graph,
        const igraph_i_layering_t* layering, const igraph_matrix_t* layout,
        const igraph_vector_bool_t* ignored_edges,
        igraph_bool_t reverse, igraph_bool_t align_right,
        igraph_vector_t* roots, igraph_vector_t* align) {
    long int i, j, k, n, di, dj, i_limit, j_limit, r;
    long int no_of_layers = igraph_i_layering_num_layers(layering);
    long int no_of_nodes = igraph_vcount(graph);
    igraph_neimode_t neimode = (reverse ? IGRAPH_OUT : IGRAPH_IN);
    igraph_vector_t neis, xs, inds;

    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&xs, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&inds, 0);

    IGRAPH_CHECK(igraph_vector_resize(roots, no_of_nodes));
    IGRAPH_CHECK(igraph_vector_resize(align, no_of_nodes));

    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(*roots)[i] = VECTOR(*align)[i] = i;
    }

    /* When reverse = False, we are aligning "upwards" in the tree, hence we
     * have to loop i from 1 to no_of_layers-1 (inclusive) and use neimode=IGRAPH_IN.
     * When reverse = True, we are aligning "downwards", hence we have to loop
     * i from no_of_layers-2 to 0 (inclusive) and use neimode=IGRAPH_OUT.
     */
    i       = reverse ? (no_of_layers - 2) : 1;
    di      = reverse ? -1 : 1;
    i_limit = reverse ? -1 : no_of_layers;
    for (; i != i_limit; i += di) {
        igraph_vector_t *layer = igraph_i_layering_get(layering, i);

        /* r = 0 in the paper, but C arrays are indexed from 0 */
        r = align_right ? LONG_MAX : -1;

        /* If align_right is 1, we have to process the layer in reverse order */
        j       = align_right ? (igraph_vector_size(layer) - 1) : 0;
        dj      = align_right ? -1 : 1;
        j_limit = align_right ? -1 : igraph_vector_size(layer);
        for (; j != j_limit; j += dj) {
            long int medians[2];
            long int vertex = (long int) VECTOR(*layer)[j];
            long int pos;

            if (VECTOR(*align)[vertex] != vertex)
                /* This vertex is already aligned with some other vertex,
                 * so there's nothing to do */
            {
                continue;
            }

            /* Find the neighbors of vertex j in layer i */
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) vertex,
                                          neimode));

            n = igraph_vector_size(&neis);
            if (n == 0)
                /* No neighbors in this direction, continue */
            {
                continue;
            }
            if (n == 1) {
                /* Just one neighbor; the median is trivial */
                medians[0] = (long int) VECTOR(neis)[0];
                medians[1] = -1;
            } else {
                /* Sort the neighbors by their X coordinates */
                IGRAPH_CHECK(igraph_vector_resize(&xs, n));
                for (k = 0; k < n; k++) {
                    VECTOR(xs)[k] = X_POS((long int)VECTOR(neis)[k]);
                }
                IGRAPH_CHECK((int) igraph_vector_qsort_ind(&xs, &inds, 0));

                if (n % 2 == 1) {
                    /* Odd number of neighbors, so the median is unique */
                    medians[0] = (long int) VECTOR(neis)[(long int)VECTOR(inds)[n / 2]];
                    medians[1] = -1;
                } else {
                    /* Even number of neighbors, so we have two medians. The order
                     * depends on whether we are processing the layer in leftmost
                     * or rightmost fashion. */
                    if (align_right) {
                        medians[0] = (long int) VECTOR(neis)[(long int)VECTOR(inds)[n / 2]];
                        medians[1] = (long int) VECTOR(neis)[(long int)VECTOR(inds)[n / 2 - 1]];
                    } else {
                        medians[0] = (long int) VECTOR(neis)[(long int)VECTOR(inds)[n / 2 - 1]];
                        medians[1] = (long int) VECTOR(neis)[(long int)VECTOR(inds)[n / 2]];
                    }
                }
            }

            /* Try aligning with the medians */
            for (k = 0; k < 2; k++) {
                igraph_integer_t eid;
                if (medians[k] < 0) {
                    continue;
                }
                if (VECTOR(*align)[vertex] != vertex) {
                    /* Vertex already aligned, continue */
                    continue;
                }
                /* Is the edge between medians[k] and vertex ignored
                 * because of a type 1 conflict? */
                IGRAPH_CHECK(igraph_get_eid(graph, &eid, (igraph_integer_t) vertex,
                                            (igraph_integer_t) medians[k], 0, 1));
                if (VECTOR(*ignored_edges)[(long int)eid]) {
                    continue;
                }
                /* Okay, align with the median if possible */
                pos = (long int) X_POS(medians[k]);
                if ((align_right && r > pos) || (!align_right && r < pos)) {
                    VECTOR(*align)[medians[k]] = vertex;
                    VECTOR(*roots)[vertex] = VECTOR(*roots)[medians[k]];
                    VECTOR(*align)[vertex] = VECTOR(*roots)[medians[k]];
                    r = pos;
                }
            }
        }
    }

    igraph_vector_destroy(&inds);
    igraph_vector_destroy(&neis);
    igraph_vector_destroy(&xs);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/*
 * Runs a horizontal compaction given a vertical alignment (in `align`)
 * and the roots (in `roots`). These come out directly from
 * igraph_i_layout_sugiyama_vertical_alignment.
 *
 * Returns the X coordinates for each vertex in `xs`.
 *
 * `graph` is the input graph, `layering` is the layering on which we operate.
 * `hgap` is the preferred horizontal gap between vertices.
 */
static int igraph_i_layout_sugiyama_horizontal_compaction(const igraph_t* graph,
        const igraph_vector_t* vertex_to_the_left,
        const igraph_vector_t* roots, const igraph_vector_t* align,
        igraph_real_t hgap, igraph_vector_t* xs) {
    long int i;
    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t sinks, shifts, old_xs;
    igraph_real_t shift;

    /* Initialization */

    IGRAPH_VECTOR_INIT_FINALLY(&sinks, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&shifts, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&old_xs, no_of_nodes);

    IGRAPH_CHECK(igraph_vector_resize(xs, no_of_nodes));

    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(sinks)[i] = i;
    }
    igraph_vector_fill(&shifts, IGRAPH_INFINITY);
    igraph_vector_fill(xs, -1);

    /* Calculate the coordinates of the vertices relative to their sinks
     * in their own class. At the end of this for loop, xs will contain the
     * relative displacement of a vertex from its sink, while the shifts list
     * will contain the absolute displacement of the sinks.
     * (For the sinks only, of course, the rest is undefined and unused)
     */
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(*roots)[i] == i) {
            IGRAPH_CHECK(
                igraph_i_layout_sugiyama_horizontal_compaction_place_block(i,
                        vertex_to_the_left, roots, align, &sinks, &shifts, hgap, xs)
            );
        }
    }

    /* In "sinks", only those indices `i` matter for which `i` is in `roots`.
     * All the other values will never be touched.
     */

    /* Calculate the absolute coordinates */
    IGRAPH_CHECK(igraph_vector_update(&old_xs, xs));
    for (i = 0; i < no_of_nodes; i++) {
        long int root = (long int) VECTOR(*roots)[i];
        VECTOR(*xs)[i] = VECTOR(old_xs)[root];
        shift = VECTOR(shifts)[(long int)VECTOR(sinks)[root]];
        if (shift < IGRAPH_INFINITY) {
            VECTOR(*xs)[i] += shift;
        }
    }

    igraph_vector_destroy(&sinks);
    igraph_vector_destroy(&shifts);
    igraph_vector_destroy(&old_xs);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

static int igraph_i_layout_sugiyama_horizontal_compaction_place_block(long int v,
        const igraph_vector_t* vertex_to_the_left,
        const igraph_vector_t* roots, const igraph_vector_t* align,
        igraph_vector_t* sinks, igraph_vector_t* shifts,
        igraph_real_t hgap, igraph_vector_t* xs) {
    long int u, w;
    long int u_sink, v_sink;

    if (VECTOR(*xs)[v] >= 0) {
        return IGRAPH_SUCCESS;
    }

    VECTOR(*xs)[v] = 0;

    w = v;
    do {
        /* Check whether vertex w is the leftmost in its own layer */
        u = (long int) VECTOR(*vertex_to_the_left)[w];
        if (u != w) {
            /* Get the root of u (proceeding all the way upwards in the block) */
            u = (long int) VECTOR(*roots)[u];
            /* Place the block of u recursively */
            IGRAPH_CHECK(
                igraph_i_layout_sugiyama_horizontal_compaction_place_block(u,
                        vertex_to_the_left, roots, align, sinks, shifts, hgap, xs)
            );

            u_sink = (long int) VECTOR(*sinks)[u];
            v_sink = (long int) VECTOR(*sinks)[v];
            /* If v is its own sink yet, set its sink to the sink of u */
            if (v_sink == v) {
                VECTOR(*sinks)[v] = v_sink = u_sink;
            }
            /* If v and u have different sinks (i.e. they are in different classes),
             * shift the sink of u so that the two blocks are separated by the
             * preferred gap
             */
            if (v_sink != u_sink) {
                if (VECTOR(*shifts)[u_sink] > VECTOR(*xs)[v] - VECTOR(*xs)[u] - hgap) {
                    VECTOR(*shifts)[u_sink] = VECTOR(*xs)[v] - VECTOR(*xs)[u] - hgap;
                }
            } else {
                /* v and u have the same sink, i.e. they are in the same class. Make sure
                 * that v is separated from u by at least hgap.
                 */
                if (VECTOR(*xs)[v] < VECTOR(*xs)[u] + hgap) {
                    VECTOR(*xs)[v] = VECTOR(*xs)[u] + hgap;
                }
            }
        }

        /* Follow the alignment */
        w = (long int) VECTOR(*align)[w];
    } while (w != v);

    return IGRAPH_SUCCESS;
}

#undef IS_INNER_SEGMENT
#undef IS_DUMMY
#undef X_POS

#ifdef SUGIYAMA_DEBUG
    #undef SUGIYAMA_DEBUG
#endif
