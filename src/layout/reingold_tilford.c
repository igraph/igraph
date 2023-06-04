/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2003-2020  The igraph development team

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

#include "igraph_adjlist.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_paths.h"
#include "igraph_progress.h"
#include "igraph_structural.h"

#include "core/math.h"

static igraph_error_t igraph_i_layout_reingold_tilford_unreachable(
    const igraph_t *graph,
    igraph_neimode_t mode,
    igraph_integer_t real_root,
    igraph_integer_t no_of_nodes,
    igraph_vector_int_t *pnewedges) {

    igraph_integer_t no_of_newedges;
    igraph_vector_bool_t visited;
    igraph_integer_t i, j, n;
    igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;
    igraph_adjlist_t allneis;
    igraph_vector_int_t *neis;

    igraph_vector_int_clear(pnewedges);

    /* traverse from real_root and see what nodes you cannot reach */
    no_of_newedges = 0;
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&visited, no_of_nodes);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, mode, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

    /* start from real_root and go BFS */
    IGRAPH_CHECK(igraph_dqueue_int_push(&q, real_root));
    while (!igraph_dqueue_int_empty(&q)) {
        igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
        neis = igraph_adjlist_get(&allneis, actnode);
        n = igraph_vector_int_size(neis);
        VECTOR(visited)[actnode] = true;
        for (j = 0; j < n; j++) {
            igraph_integer_t neighbor = VECTOR(*neis)[j];
            if (!VECTOR(visited)[neighbor]) {
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
            }
        }
    }

    for (j = 0; j < no_of_nodes; j++) {
        no_of_newedges += VECTOR(visited)[j] ? 0 : 1;
    }

    /* if any nodes are unreachable, add edges between them and real_root */
    if (no_of_newedges != 0) {

        igraph_vector_int_resize(pnewedges, no_of_newedges * 2);
        j = 0;
        for (i = 0; i < no_of_nodes; i++) {
            if (!VECTOR(visited)[i]) {
                if (mode != IGRAPH_IN) {
                    VECTOR(*pnewedges)[2 * j] = real_root;
                    VECTOR(*pnewedges)[2 * j + 1] = i;
                } else {
                    VECTOR(*pnewedges)[2 * j] = i;
                    VECTOR(*pnewedges)[2 * j + 1] = real_root;
                }
                j++;
            }
        }
    }

    igraph_dqueue_int_destroy(&q);
    igraph_adjlist_destroy(&allneis);
    igraph_vector_bool_destroy(&visited);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}


/* Internal structure for Reingold-Tilford layout */
struct igraph_i_reingold_tilford_vertex {
    igraph_integer_t parent;        /* Parent node index */
    igraph_integer_t level;         /* Level of the node */
    igraph_real_t offset;     /* X offset from parent node */
    igraph_integer_t left_contour;  /* Next left node of the contour
              of the subtree rooted at this node */
    igraph_integer_t right_contour; /* Next right node of the contour
              of the subtree rooted at this node */
    igraph_real_t offset_to_left_contour;  /* X offset when following the left contour */
    igraph_real_t offset_to_right_contour;  /* X offset when following the right contour */
    igraph_integer_t left_extreme;  /* Leftmost node on the deepest layer of the subtree rooted at this node */
    igraph_integer_t right_extreme; /* Rightmost node on the deepest layer of the subtree rooted at this node */
    igraph_real_t offset_to_left_extreme;  /* X offset when jumping to the left extreme node */
    igraph_real_t offset_to_right_extreme;  /* X offset when jumping to the right extreme node */
};

static void igraph_i_layout_reingold_tilford_postorder(struct igraph_i_reingold_tilford_vertex *vdata,
                                                      igraph_integer_t node, igraph_integer_t vcount);
static void igraph_i_layout_reingold_tilford_calc_coords(struct igraph_i_reingold_tilford_vertex *vdata,
                                                        igraph_matrix_t *res, igraph_integer_t node,
                                                        igraph_integer_t vcount, igraph_real_t xpos);

/* uncomment the next line for debugging the Reingold-Tilford layout */
/* #define LAYOUT_RT_DEBUG 1 */

static igraph_error_t igraph_i_layout_reingold_tilford(const igraph_t *graph,
                                            igraph_matrix_t *res,
                                            igraph_neimode_t mode,
                                            igraph_integer_t root) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t i, n, j;
    igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;
    igraph_adjlist_t allneis;
    igraph_vector_int_t *neis;
    struct igraph_i_reingold_tilford_vertex *vdata;

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, mode, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

    vdata = IGRAPH_CALLOC(no_of_nodes, struct igraph_i_reingold_tilford_vertex);
    if (vdata == 0) {
        IGRAPH_ERROR("igraph_layout_reingold_tilford failed", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, vdata);

    for (i = 0; i < no_of_nodes; i++) {
        vdata[i].parent = -1;
        vdata[i].level = -1;
        vdata[i].offset = 0.0;
        vdata[i].left_contour = -1;
        vdata[i].right_contour = -1;
        vdata[i].offset_to_left_contour = 0.0;
        vdata[i].offset_to_right_contour = 0.0;
        vdata[i].left_extreme = i;
        vdata[i].right_extreme = i;
        vdata[i].offset_to_left_extreme = 0.0;
        vdata[i].offset_to_right_extreme = 0.0;
    }
    vdata[root].parent = root;
    vdata[root].level = 0;
    MATRIX(*res, root, 1) = 0;

    /* Step 1: assign Y coordinates based on BFS and setup parents vector */
    IGRAPH_CHECK(igraph_dqueue_int_push(&q, root));
    IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
    while (!igraph_dqueue_int_empty(&q)) {
        igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
        igraph_integer_t actdist = igraph_dqueue_int_pop(&q);
        neis = igraph_adjlist_get(&allneis, actnode);
        n = igraph_vector_int_size(neis);

        for (j = 0; j < n; j++) {
            igraph_integer_t neighbor = VECTOR(*neis)[j];
            if (vdata[neighbor].parent >= 0) {
                continue;
            }
            MATRIX(*res, neighbor, 1) = actdist + 1;
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
            vdata[neighbor].parent = actnode;
            vdata[neighbor].level = actdist + 1;
        }
    }

    /* Step 2: postorder tree traversal, determines the appropriate X
     * offsets for every node */
    igraph_i_layout_reingold_tilford_postorder(vdata, root, no_of_nodes);

    /* Step 3: calculate real coordinates based on X offsets */
    igraph_i_layout_reingold_tilford_calc_coords(vdata, res, root, no_of_nodes, vdata[root].offset);

    igraph_dqueue_int_destroy(&q);
    igraph_adjlist_destroy(&allneis);
    igraph_free(vdata);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_PROGRESS("Reingold-Tilford tree layout", 100.0, NULL);

#ifdef LAYOUT_RT_DEBUG
    for (i = 0; i < no_of_nodes; i++) {
        printf(
            "%3" IGRAPH_PRId ": offset = %.2f, contours = [%" IGRAPH_PRId ", %" IGRAPH_PRId "], contour offsets = [%.2f, %.2f]\n",
            i, vdata[i].offset,
            vdata[i].left_contour, vdata[i].right_contour,
            vdata[i].offset_to_left_contour, vdata[i].offset_to_right_contour
        );
        if (vdata[i].left_extreme != i || vdata[i].right_extreme != i) {
            printf(
                "     extrema = [%" IGRAPH_PRId ", %" IGRAPH_PRId "], offsets to extrema = [%.2f, %.2f]\n",
                vdata[i].left_extreme, vdata[i].right_extreme,
                vdata[i].offset_to_left_extreme, vdata[i].offset_to_right_extreme
            );
        }
    }
#endif

    return IGRAPH_SUCCESS;
}

static void igraph_i_layout_reingold_tilford_calc_coords(
        struct igraph_i_reingold_tilford_vertex *vdata,
        igraph_matrix_t *res, igraph_integer_t node,
        igraph_integer_t vcount, igraph_real_t xpos) {

    MATRIX(*res, node, 0) = xpos;
    for (igraph_integer_t i = 0; i < vcount; i++) {
        if (i == node) {
            continue;
        }
        if (vdata[i].parent == node) {
            igraph_i_layout_reingold_tilford_calc_coords(vdata, res, i, vcount,
                    xpos + vdata[i].offset);
        }
    }
}

static void igraph_i_layout_reingold_tilford_postorder(
        struct igraph_i_reingold_tilford_vertex *vdata,
        igraph_integer_t node, igraph_integer_t vcount) {

    igraph_integer_t childcount, leftroot, leftrootidx;
    const igraph_real_t minsep = 1;
    igraph_real_t avg;

#ifdef LAYOUT_RT_DEBUG
    printf("Starting visiting node %" IGRAPH_PRId "\n", node);
#endif

    /* Check whether this node is a leaf node */
    childcount = 0;
    for (igraph_integer_t i = 0; i < vcount; i++) {
        if (i == node) {
            continue;
        }
        if (vdata[i].parent == node) {
            /* Node i is a child, so visit it recursively */
            childcount++;
            igraph_i_layout_reingold_tilford_postorder(vdata, i, vcount);
        }
    }

    if (childcount == 0) {
        return;
    }

    /* Here we can assume that all of the subtrees have been placed and their
     * left and right contours are calculated. Let's place them next to each
     * other as close as we can.
     * We will take each subtree in an arbitrary order. The root of the
     * first one will be placed at offset 0, the next ones will be placed
     * as close to each other as possible. leftroot stores the root of the
     * rightmost subtree of the already placed subtrees - its right contour
     * will be checked against the left contour of the next subtree */
    leftroot = leftrootidx = -1;
    avg = 0.0;
#ifdef LAYOUT_RT_DEBUG
    printf("Visited node %" IGRAPH_PRId " and arranged its subtrees\n", node);
#endif
    for (igraph_integer_t i = 0, j = 0; i < vcount; i++) {
        if (i == node) {
            continue;
        }
        if (vdata[i].parent == node) {
            if (leftroot >= 0) {
                /* Now we will follow the right contour of leftroot and the
                 * left contour of the subtree rooted at i */
                igraph_integer_t lnode, rnode, auxnode;
                igraph_real_t loffset, roffset, rootsep, newoffset;

#ifdef LAYOUT_RT_DEBUG
                printf("  Placing child %" IGRAPH_PRId " on level %" IGRAPH_PRId ", to the right of %" IGRAPH_PRId "\n", i, vdata[i].level, leftroot);
#endif
                lnode = leftroot; rnode = i;
                rootsep = vdata[leftroot].offset + minsep;
                loffset = vdata[leftroot].offset; roffset = loffset + minsep;

                /* Keep on updating the right contour now that we have attached
                 * a new node to the subtree being built */
                vdata[node].right_contour = i;
                vdata[node].offset_to_right_contour = rootsep;

#ifdef LAYOUT_RT_DEBUG
                printf("    Contour: [%" IGRAPH_PRId ", %" IGRAPH_PRId "], offsets: [%lf, %lf], rootsep: %lf\n",
                       lnode, rnode, loffset, roffset, rootsep);
#endif
                while ((lnode >= 0) && (rnode >= 0)) {
                    /* Step to the next level on the right contour of the left subtree */
                    if (vdata[lnode].right_contour >= 0) {
                        loffset += vdata[lnode].offset_to_right_contour;
                        lnode = vdata[lnode].right_contour;
                    } else {
                        /* Left subtree ended there. The left and right contour
                         * of the left subtree will continue to the next step
                         * on the right subtree. */
                        if (vdata[rnode].left_contour >= 0) {
                            auxnode = vdata[node].left_extreme;

                            /* this is the "threading" step that the original
                             * paper is talking about */
                            newoffset = (vdata[node].offset_to_right_extreme - vdata[node].offset_to_left_extreme) + minsep + vdata[rnode].offset_to_left_contour;
                            vdata[auxnode].left_contour = vdata[rnode].left_contour;
                            vdata[auxnode].right_contour = vdata[rnode].left_contour;
                            vdata[auxnode].offset_to_left_contour = vdata[auxnode].offset_to_right_contour = newoffset;

                            /* since we attached a larger subtree to the
                             * already placed left subtree, we need to update
                             * the extrema of the subtree rooted at 'node' */
                            vdata[node].left_extreme = vdata[i].left_extreme;
                            vdata[node].right_extreme = vdata[i].right_extreme;
                            vdata[node].offset_to_left_extreme = vdata[i].offset_to_left_extreme + rootsep;
                            vdata[node].offset_to_right_extreme = vdata[i].offset_to_right_extreme + rootsep;
#ifdef LAYOUT_RT_DEBUG
                            printf("      Left subtree ended earlier, continuing left subtree's left and right contour on right subtree (node %" IGRAPH_PRId " gets connected to node %" IGRAPH_PRId ")\n", auxnode, vdata[rnode].left_contour);
                            printf("      New contour following offset for node %" IGRAPH_PRId " is %lf\n", auxnode, vdata[auxnode].offset_to_left_contour);
#endif
                        } else {
                            /* Both subtrees are ending at the same time; the
                             * left extreme node of the subtree rooted at
                             * 'node' remains the same but the right extreme
                             * will change */
                            vdata[node].right_extreme = vdata[i].right_extreme;
                            vdata[node].offset_to_right_extreme = vdata[i].offset_to_right_extreme + rootsep;
                        }
                        lnode = -1;
                    }
                    /* Step to the next level on the left contour of the right subtree */
                    if (vdata[rnode].left_contour >= 0) {
                        roffset += vdata[rnode].offset_to_left_contour;
                        rnode = vdata[rnode].left_contour;
                    } else {
                        /* Right subtree ended here. The right contour of the right
                         * subtree will continue to the next step on the left subtree.
                         * Note that lnode has already been advanced here */
                        if (lnode >= 0) {
                            auxnode = vdata[i].right_extreme;

                            /* this is the "threading" step that the original
                             * paper is talking about */
                            newoffset = loffset - rootsep - vdata[i].offset_to_right_extreme;
                            vdata[auxnode].left_contour = lnode;
                            vdata[auxnode].right_contour = lnode;
                            vdata[auxnode].offset_to_left_contour = vdata[auxnode].offset_to_right_contour = newoffset;

                            /* no need to update the extrema of the subtree
                             * rooted at 'node' because the right subtree was
                             * smaller */
#ifdef LAYOUT_RT_DEBUG
                            printf("      Right subtree ended earlier, continuing right subtree's left and right contour on left subtree (node %" IGRAPH_PRId " gets connected to node %" IGRAPH_PRId ")\n", auxnode, lnode);
                            printf("      New contour following offset for node %" IGRAPH_PRId " is %lf\n", auxnode, vdata[auxnode].offset_to_left_contour);
#endif
                        }
                        rnode = -1;
                    }
#ifdef LAYOUT_RT_DEBUG
                    printf("    Contour: [%" IGRAPH_PRId ", %" IGRAPH_PRId "], offsets: [%lf, %lf], rootsep: %lf\n",
                           lnode, rnode, loffset, roffset, rootsep);
#endif

                    /* Push subtrees away if necessary */
                    if ((lnode >= 0) && (rnode >= 0) && (roffset - loffset < minsep)) {
#ifdef LAYOUT_RT_DEBUG
                        printf("    Pushing right subtree away by %lf\n", minsep-roffset+loffset);
#endif
                        rootsep += minsep - roffset + loffset;
                        roffset = loffset + minsep;
                        vdata[node].offset_to_right_contour = rootsep;
                    }
                }

#ifdef LAYOUT_RT_DEBUG
                printf("  Offset of subtree with root node %" IGRAPH_PRId " will be %lf\n", i, rootsep);
#endif
                vdata[i].offset = rootsep;
                vdata[node].offset_to_right_contour = rootsep;
                avg = (avg * j) / (j + 1) + rootsep / (j + 1);
                leftrootidx = j;
                leftroot = i;
            } else {
                /* This is the first child of the node being considered,
                 * so we can simply place the subtree on our virtual canvas. */
#ifdef LAYOUT_RT_DEBUG
                printf("  Placing child %" IGRAPH_PRId " on level %" IGRAPH_PRId " as first child\n", i, vdata[i].level);
#endif
                leftrootidx = j;
                leftroot = i;
                vdata[node].left_contour = i;
                vdata[node].right_contour = i;
                vdata[node].offset_to_left_contour = 0.0;
                vdata[node].offset_to_right_contour = 0.0;
                vdata[node].left_extreme = vdata[i].left_extreme;
                vdata[node].right_extreme = vdata[i].right_extreme;
                vdata[node].offset_to_left_extreme = vdata[i].offset_to_left_extreme;
                vdata[node].offset_to_right_extreme = vdata[i].offset_to_right_extreme;
                avg = vdata[i].offset;
            }
            j++;
        }
    }
#ifdef LAYOUT_RT_DEBUG
    printf("Shifting node %" IGRAPH_PRId " to be centered above children. Shift amount: %lf\n", node, avg);
#endif
    vdata[node].offset_to_left_contour -= avg;
    vdata[node].offset_to_right_contour -= avg;
    vdata[node].offset_to_left_extreme -= avg;
    vdata[node].offset_to_right_extreme -= avg;
    for (igraph_integer_t i = 0; i < vcount; i++) {
        if (i == node) {
            continue;
        }
        if (vdata[i].parent == node) {
            vdata[i].offset -= avg;
        }
    }
}

/* This function computes the number of outgoing (or incoming) connections
 * of clusters, represented as a membership vector. It only works with
 * directed graphs. */
igraph_error_t igraph_i_layout_reingold_tilford_cluster_degrees_directed(
        const igraph_t *graph,
        const igraph_vector_int_t *membership,
        igraph_integer_t no_comps,
        igraph_neimode_t mode,
        igraph_vector_int_t *degrees) {

    igraph_eit_t eit;

    if (! igraph_is_directed(graph) || (mode != IGRAPH_OUT && mode != IGRAPH_IN)) {
        IGRAPH_ERROR("Directed graph expected.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_int_resize(degrees, no_comps));
    igraph_vector_int_null(degrees);

    IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    for (; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
        igraph_integer_t eid = IGRAPH_EIT_GET(eit);

        igraph_integer_t from = IGRAPH_FROM(graph, eid);
        igraph_integer_t to   = IGRAPH_TO(graph, eid);

        igraph_integer_t from_cl = VECTOR(*membership)[from];
        igraph_integer_t to_cl   = VECTOR(*membership)[to];

        igraph_integer_t cl = mode == IGRAPH_OUT ? from_cl : to_cl;

        if (from_cl != to_cl) {
            VECTOR(*degrees)[cl] += 1;
        }
    }

    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* Heuristic method to choose "nice" roots for the Reingold-Tilford layout algorithm.
 *
 * The principle is to select a minimal set of roots so that all other vertices
 * will be reachable from them.
 *
 * In the undirected case, one root is chosen from each connected component.
 * In the directed case, one root is chosen from each strongly connected component
 * that has no incoming (or outgoing) edges (depending on 'mode').
 * When more than one root choice is possible, nodes are prioritized based on
 * either lowest eccentricity (if 'use_eccentricity' is true) or based on
 * highest degree (out- or in-degree in directed mode).
 */

/**
 * \function igraph_roots_for_tree_layout
 * \brief Roots suitable for a nice tree layout.
 *
 * This function chooses a root, or a set of roots suitable for visualizing a tree,
 * or a tree-like graph. It is typically used with \ref igraph_layout_reingold_tilford().
 * The principle is to select a minimal set of roots so that all other vertices
 * will be reachable from them.
 *
 * </para><para>
 * In the undirected case, one root is chosen from each connected component.
 * In the directed case, one root is chosen from each strongly connected component
 * that has no incoming (or outgoing) edges (depending on 'mode'). When more than
 * one root choice is possible, vertices are prioritized based on the given \p heuristic.
 *
 * \param graph The graph, typically a tree, but any graph is accepted.
 * \param mode Whether to interpret the input as undirected, a directed out-tree or in-tree.
 * \param roots An initialized integer vector, the roots will be returned here.
 * \param heuristic The heuristic to use for breaking ties when multiple root
 *   choices are possible.
 *          \clist
 *          \cli IGRAPH_ROOT_CHOICE_DEGREE
 *           Choose the vertices with the highest degree (out- or in-degree
 *           in directed mode). This simple heuristic is fast even in large graphs.
 *          \cli IGRAPH_ROOT_CHOICE_ECCENTRICITY
 *           Choose the vertices with the lowest eccentricity. This usually results
 *           in a "wide and shallow" tree layout. While this heuristic produces
 *           high-quality results, it is slow for large graphs: computing the
 *           eccentricities has quadratic complexity in the number of vertices.
 *          \endclist
 * \return Error code.
 *
 * Time complexity: depends on the heuristic.
 */
igraph_error_t igraph_roots_for_tree_layout(
        const igraph_t *graph,
        igraph_neimode_t mode,
        igraph_vector_int_t *roots,
        igraph_root_choice_t heuristic) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t order, membership;
    igraph_integer_t no_comps;
    igraph_integer_t i, j;
    igraph_bool_t use_eccentricity;

    switch (heuristic) {
    case IGRAPH_ROOT_CHOICE_DEGREE:
        use_eccentricity = false; break;
    case IGRAPH_ROOT_CHOICE_ECCENTRICITY:
        use_eccentricity = true; break;
    default:
        IGRAPH_ERROR("Invalid root choice heuristic given.", IGRAPH_EINVAL);
    }

    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    if (no_of_nodes == 0) {
        igraph_vector_int_clear(roots);
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&order, no_of_nodes);
    if (use_eccentricity) {
        /* Sort vertices by decreasing eccentricity. */

        igraph_vector_t ecc;

        IGRAPH_VECTOR_INIT_FINALLY(&ecc, no_of_nodes);
        IGRAPH_CHECK(igraph_eccentricity(graph, &ecc, igraph_vss_all(), mode));
        IGRAPH_CHECK(igraph_vector_qsort_ind(&ecc, &order, IGRAPH_ASCENDING));

        igraph_vector_destroy(&ecc);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        /* Sort vertices by decreasing degree (out- or in-degree in directed case). */

        IGRAPH_CHECK(igraph_sort_vertex_ids_by_degree(graph, &order,
                     igraph_vss_all(), mode, 0, IGRAPH_DESCENDING, 0));
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&membership, no_of_nodes);
    IGRAPH_CHECK(igraph_connected_components(
        graph, &membership, /*csize=*/ NULL,
        &no_comps, mode == IGRAPH_ALL ? IGRAPH_WEAK : IGRAPH_STRONG
    ));

    IGRAPH_CHECK(igraph_vector_int_resize(roots, no_comps));
    igraph_vector_int_fill(roots, -1); /* -1 signifies a not-yet-determined root for a component */

    if (mode != IGRAPH_ALL) {
        /* Directed case:
         *
         * We break the graph into strongly-connected components and find those components
         * which have no incoming (outgoing) edges. The largest out-degree (in-degree)
         * nodes from these components will be chosen as roots. When the graph is a DAG,
         * these will simply be the source (sink) nodes. */

        igraph_vector_int_t cluster_degrees;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&cluster_degrees, no_of_nodes);
        IGRAPH_CHECK(igraph_i_layout_reingold_tilford_cluster_degrees_directed(
                         graph, &membership, no_comps,
                         mode == IGRAPH_OUT ? IGRAPH_IN : IGRAPH_OUT, /* reverse direction */
                         &cluster_degrees));

        /* Iterate through nodes in decreasing out-degree (or in-degree) order
         * and record largest degree node in each strongly-connected component
         * which has no incoming (outgoing) edges. */
        for (i = 0; i < no_of_nodes; ++i) {
            igraph_integer_t v  = VECTOR(order)[i];
            igraph_integer_t cl = VECTOR(membership)[v];
            if (VECTOR(cluster_degrees)[cl] == 0 && VECTOR(*roots)[cl] == -1) {
                VECTOR(*roots)[cl] = v;
            }
        }

        igraph_vector_int_destroy(&cluster_degrees);
        IGRAPH_FINALLY_CLEAN(1);

        /* Remove remaining -1 indices. These correspond to components that
         * did have some incoming edges. */
        for (i=0, j=0; i < no_comps; ++i) {
            if (VECTOR(*roots)[i] == -1) {
                continue;
            }
            VECTOR(*roots)[j++] = VECTOR(*roots)[i];
        }
        igraph_vector_int_resize(roots, j);

    } else {
        /* Undirected case:
         *
         * Select the highest degree node from each component.
         */

        igraph_integer_t no_seen = 0;

        for (i=0; i < no_of_nodes; ++i) {
            igraph_integer_t v  = VECTOR(order)[i];
            igraph_integer_t cl = VECTOR(membership)[v];
            if (VECTOR(*roots)[cl] == -1) {
                no_seen += 1;
                VECTOR(*roots)[cl] = v;
            }
            if (no_seen == no_comps) {
                /* All components have roots now. */
                break;
            }
        }
    }

    igraph_vector_int_destroy(&membership);
    igraph_vector_int_destroy(&order);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_layout_reingold_tilford
 * \brief Reingold-Tilford layout for tree graphs.
 *
 * </para><para>
 * Arranges the nodes in a tree where the given node is used as the root.
 * The tree is directed downwards and the parents are centered above its
 * children. For the exact algorithm, see:
 *
 * </para><para>
 * Reingold, E and Tilford, J: Tidier drawing of trees.
 * IEEE Trans. Softw. Eng., SE-7(2):223--228, 1981.
 * https://doi.org/10.1109/TSE.1981.234519
 *
 * </para><para>
 * If the given graph is not a tree, a breadth-first search is executed
 * first to obtain a possible spanning tree.
 *
 * \param graph The graph object.
 * \param res The result, the coordinates in a matrix. The parameter
 *   should point to an initialized matrix object and will be resized.
 * \param mode Specifies which edges to consider when building the tree.
 *   If it is \c IGRAPH_OUT then only the outgoing, if it is \c IGRAPH_IN
 *   then only the incoming edges of a parent are considered. If it is
 *   \c IGRAPH_ALL then all edges are used (this was the behavior in
 *   igraph 0.5 and before). This parameter also influences how the root
 *   vertices are calculated, if they are not given. See the \p roots parameter.
 * \param roots The index of the root vertex or root vertices. The set of roots
 *   should be specified so that all vertices of the graph are reachable from them.
 *   Simply put, in the undirected case, one root should be given from each
 *   connected component. If \p roots is \c NULL or a pointer to an empty vector,
 *   then the roots will be selected automatically. Currently, automatic root
 *   selection prefers low eccentricity vertices in graphs with fewer than
 *   500 vertices, and high degree vertices (according to \p mode) in larger graphs.
 *   The root selection heuristic may change without notice. To ensure a consistent
 *   output, please specify the roots manually. The \ref igraph_roots_for_tree_layout()
 *   function gives more control over automatic root selection.
 * \param rootlevel This argument can be useful when drawing forests which are
 *   not trees (i.e. they are unconnected and have tree components). It specifies
 *   the level of the root vertices for every tree in the forest. It is only
 *   considered if not a null pointer and the \p roots argument is also given
 *   (and it is not a null pointer of an empty vector).
 * \return Error code.
 *
 * Added in version 0.2.
 *
 * \sa \ref igraph_layout_reingold_tilford_circular(), \ref igraph_roots_for_tree_layout()
 *
 * \example examples/simple/igraph_layout_reingold_tilford.c
 */
igraph_error_t igraph_layout_reingold_tilford(const igraph_t *graph,
                                   igraph_matrix_t *res,
                                   igraph_neimode_t mode,
                                   const igraph_vector_int_t *roots,
                                   const igraph_vector_int_t *rootlevel) {

    const igraph_integer_t no_of_nodes_orig = igraph_vcount(graph);
    igraph_integer_t no_of_nodes = no_of_nodes_orig;
    igraph_integer_t real_root;
    igraph_t extended;
    const igraph_t *pextended = graph;
    igraph_vector_int_t myroots;
    const igraph_vector_int_t *proots = roots;
    igraph_vector_int_t newedges;


    /* TODO: possible speedup could be achieved if we use a table for storing
     * the children of each node in the tree. (Now the implementation uses a
     * single array containing the parent of each node and a node's children
     * are determined by looking for other nodes that have this node as parent)
     */

    /* at various steps it might be necessary to add edges to the graph */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&newedges, 0);

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    if ( (!roots || igraph_vector_int_size(roots) == 0) &&
         rootlevel && igraph_vector_int_size(rootlevel) != 0 ) {
        IGRAPH_WARNING("Reingold-Tilford layout: 'rootlevel' ignored");
    }

    /* ----------------------------------------------------------------------- */
    /* If root vertices are not given, perform automated root selection. */

    if (!roots || igraph_vector_int_size(roots) == 0) {

        IGRAPH_VECTOR_INT_INIT_FINALLY(&myroots, 0);
        igraph_roots_for_tree_layout(graph, mode, &myroots,
                                     no_of_nodes < 500 ? IGRAPH_ROOT_CHOICE_DEGREE : IGRAPH_ROOT_CHOICE_ECCENTRICITY);
        proots = &myroots;

    } else if (rootlevel && igraph_vector_int_size(rootlevel) > 0 &&
               igraph_vector_int_size(roots) > 1) {

        /* ----------------------------------------------------------------------- */
        /* Many roots were given to us, check 'rootlevel' */

        igraph_integer_t plus_levels = 0;

        if (igraph_vector_int_size(roots) != igraph_vector_int_size(rootlevel)) {
            IGRAPH_ERROR("Reingold-Tilford: 'roots' and 'rootlevel' lengths differ",
                         IGRAPH_EINVAL);
        }

        /* count the rootlevels that are not zero */
        for (igraph_integer_t i = 0; i < igraph_vector_int_size(roots); i++) {
            plus_levels += VECTOR(*rootlevel)[i];
        }

        /* make copy of graph, add vertices/edges */
        if (plus_levels != 0) {
            igraph_integer_t edgeptr = 0;

            pextended = &extended;
            IGRAPH_CHECK(igraph_copy(&extended, graph));
            IGRAPH_FINALLY(igraph_destroy, &extended);
            IGRAPH_CHECK(igraph_add_vertices(&extended, plus_levels, 0));

            igraph_vector_int_resize(&newedges, plus_levels * 2);

            for (igraph_integer_t i = 0; i < igraph_vector_int_size(roots); i++) {
                igraph_integer_t rl = VECTOR(*rootlevel)[i];
                igraph_integer_t rn = VECTOR(*roots)[i];
                igraph_integer_t j;

                /* zero-level roots don't get anything special */
                if (rl == 0) {
                    continue;
                }

                /* for each nonzero-level root, add vertices
                   and edges at all levels [1, 2, .., rl]
                   piercing through the graph. If mode=="in"
                   they pierce the other way */
                if (mode != IGRAPH_IN) {
                    VECTOR(newedges)[edgeptr++] = no_of_nodes;
                    VECTOR(newedges)[edgeptr++] = rn;
                    for (j = 0; j < rl - 1; j++) {
                        VECTOR(newedges)[edgeptr++] = no_of_nodes + 1;
                        VECTOR(newedges)[edgeptr++] = no_of_nodes;
                        no_of_nodes++;
                    }
                } else {
                    VECTOR(newedges)[edgeptr++] = rn;
                    VECTOR(newedges)[edgeptr++] = no_of_nodes;
                    for (j = 0; j < rl - 1; j++) {
                        VECTOR(newedges)[edgeptr++] = no_of_nodes;
                        VECTOR(newedges)[edgeptr++] = no_of_nodes + 1;
                        no_of_nodes++;
                    }
                }

                /* move on to the next root */
                VECTOR(*roots)[i] = no_of_nodes++;
            }

            /* actually add the edges to the graph */
            IGRAPH_CHECK(igraph_add_edges(&extended, &newedges, 0));
        }
    }

    /* We have root vertices now. If one or more nonzero-level roots were
       chosen by the user, we have copied the graph and added a few vertices
       and (directed) edges to connect those floating roots to nonfloating,
       zero-level equivalent roots.

       Below, the function

       igraph_i_layout_reingold_tilford(pextended, res, mode, real_root)

       calculates the actual rt coordinates of the graph. However, for
       simplicity that function requires a connected graph and a single root.
       For directed graphs, it needs not be strongly connected, however all
       nodes must be reachable from the root following the stream (i.e. the
       root must be a "mother vertex").

       So before we call that function we have to make sure the (copied) graph
       satisfies that condition. That requires:
         1. if there is more than one root, defining a single real_root
         2. if a real_root is defined, adding edges to connect all roots to it
         3. ensure real_root is mother of the whole graph. If it is not,
            add shortcut edges from real_root to any disconnected node for now.

      NOTE: 3. could be done better, e.g. by topological sorting of some kind.
      But for now it's ok like this.
    */
    /* if there is only one root, no need for real_root */
    if (igraph_vector_int_size(proots) == 1) {
        real_root = VECTOR(*proots)[0];
        if (real_root < 0 || real_root >= no_of_nodes) {
            IGRAPH_ERROR("Invalid vertex ID.", IGRAPH_EINVVID);
        }

        /* else, we need to make real_root */
    } else {
        igraph_integer_t no_of_newedges;

        /* Make copy of the graph unless it exists already */
        if (pextended == graph) {
            pextended = &extended;
            IGRAPH_CHECK(igraph_copy(&extended, graph));
            IGRAPH_FINALLY(igraph_destroy, &extended);
        }

        /* add real_root to the vertices */
        real_root = no_of_nodes;
        IGRAPH_CHECK(igraph_add_vertices(&extended, 1, 0));
        no_of_nodes++;

        /* add edges from the roots to real_root */
        no_of_newedges = igraph_vector_int_size(proots);
        igraph_vector_int_resize(&newedges, no_of_newedges * 2);
        for (igraph_integer_t i = 0; i < no_of_newedges; i++) {
            VECTOR(newedges)[2 * i] = no_of_nodes - 1;
            VECTOR(newedges)[2 * i + 1] = VECTOR(*proots)[i];
        }

        IGRAPH_CHECK(igraph_add_edges(&extended, &newedges, 0));
    }

    /* prepare edges to unreachable parts of the graph */
    IGRAPH_CHECK(igraph_i_layout_reingold_tilford_unreachable(pextended, mode, real_root, no_of_nodes, &newedges));

    if (igraph_vector_int_size(&newedges) != 0) {
        /* Make copy of the graph unless it exists already */
        if (pextended == graph) {
            pextended = &extended;
            IGRAPH_CHECK(igraph_copy(&extended, graph));
            IGRAPH_FINALLY(igraph_destroy, &extended);
        }

        IGRAPH_CHECK(igraph_add_edges(&extended, &newedges, 0));
    }
    igraph_vector_int_destroy(&newedges);
    IGRAPH_FINALLY_CLEAN(1);

    /* ----------------------------------------------------------------------- */
    /* Layout */
    IGRAPH_CHECK(igraph_i_layout_reingold_tilford(pextended, res, mode, real_root));

    /* Remove the new vertices from the layout */
    if (no_of_nodes != no_of_nodes_orig) {
        if (no_of_nodes - 1 == no_of_nodes_orig) {
            IGRAPH_CHECK(igraph_matrix_remove_row(res, no_of_nodes_orig));
        } else {
            igraph_matrix_t tmp;
            igraph_integer_t i;
            IGRAPH_MATRIX_INIT_FINALLY(&tmp, no_of_nodes_orig, 2);
            for (i = 0; i < no_of_nodes_orig; i++) {
                MATRIX(tmp, i, 0) = MATRIX(*res, i, 0);
                MATRIX(tmp, i, 1) = MATRIX(*res, i, 1);
            }
            IGRAPH_CHECK(igraph_matrix_update(res, &tmp));
            igraph_matrix_destroy(&tmp);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    if (pextended != graph) {
        igraph_destroy(&extended);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* Remove the roots vector if it was created by us */
    if (proots != roots) {
        igraph_vector_int_destroy(&myroots);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_layout_reingold_tilford_circular
 * \brief Circular Reingold-Tilford layout for trees.
 *
 * This layout is almost the same as \ref igraph_layout_reingold_tilford(), but
 * the tree is drawn in a circular way, with the root vertex in the center.
 *
 * \param graph The graph object.
 * \param res The result, the coordinates in a matrix. The parameter
 *   should point to an initialized matrix object and will be resized.
 * \param mode Specifies which edges to consider when building the tree.
 *   If it is \c IGRAPH_OUT then only the outgoing, if it is \c IGRAPH_IN
 *   then only the incoming edges of a parent are considered. If it is
 *   \c IGRAPH_ALL then all edges are used (this was the behavior in
 *   igraph 0.5 and before). This parameter also influences how the root
 *   vertices are calculated, if they are not given. See the \p roots parameter.
 * \param roots The index of the root vertex or root vertices. The set of roots
 *   should be specified so that all vertices of the graph are reachable from them.
 *   Simply put, in the undirected case, one root should be given from each
 *   connected component. If \p roots is \c NULL or a pointer to an empty vector,
 *   then the roots will be selected automatically. Currently, automatic root
 *   selection prefers low eccentricity vertices in graphs with fewer than
 *   500 vertices, and high degree vertices (according to \p mode) in larger graphs.
 *   The root selection heuristic may change without notice. To ensure a consistent
 *   output, please specify the roots manually.
 * \param rootlevel This argument can be useful when drawing forests which are
 *   not trees (i.e. they are unconnected and have tree components). It specifies
 *   the level of the root vertices for every tree in the forest. It is only
 *   considered if not a null pointer and the \p roots argument is also given
 *   (and it is not a null pointer or an empty vector).
 * \return Error code.
 *
 * \sa \ref igraph_layout_reingold_tilford().
 */
igraph_error_t igraph_layout_reingold_tilford_circular(const igraph_t *graph,
        igraph_matrix_t *res,
        igraph_neimode_t mode,
        const igraph_vector_int_t *roots,
        const igraph_vector_int_t *rootlevel) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_real_t ratio;
    igraph_real_t minx, maxx;

    IGRAPH_CHECK(igraph_layout_reingold_tilford(graph, res, mode, roots, rootlevel));

    if (no_of_nodes == 0) {
        return IGRAPH_SUCCESS;
    }

    ratio = 2 * M_PI * (no_of_nodes - 1.0) / no_of_nodes;

    minx = maxx = MATRIX(*res, 0, 0);
    for (igraph_integer_t i = 1; i < no_of_nodes; i++) {
        if (MATRIX(*res, i, 0) > maxx) {
            maxx = MATRIX(*res, i, 0);
        }
        if (MATRIX(*res, i, 0) < minx) {
            minx = MATRIX(*res, i, 0);
        }
    }
    if (maxx > minx) {
        ratio /= (maxx - minx);
    }
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_real_t phi = (MATRIX(*res, i, 0) - minx) * ratio;
        igraph_real_t r = MATRIX(*res, i, 1);
        MATRIX(*res, i, 0) = r * cos(phi);
        MATRIX(*res, i, 1) = r * sin(phi);
    }

    return IGRAPH_SUCCESS;
}
