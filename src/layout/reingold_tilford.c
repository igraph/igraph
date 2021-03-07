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
#include "igraph_progress.h"
#include "igraph_structural.h"
#include "igraph_topology.h"

#include "core/math.h"

static int igraph_i_layout_reingold_tilford_unreachable(
    const igraph_t *graph,
    igraph_neimode_t mode,
    long int real_root,
    long int no_of_nodes,
    igraph_vector_t *pnewedges) {

    long int no_of_newedges;
    igraph_vector_t visited;
    long int i, j, n;
    igraph_dqueue_t q = IGRAPH_DQUEUE_NULL;
    igraph_adjlist_t allneis;
    igraph_vector_int_t *neis;

    igraph_vector_resize(pnewedges, 0);

    /* traverse from real_root and see what nodes you cannot reach */
    no_of_newedges = 0;
    IGRAPH_VECTOR_INIT_FINALLY(&visited, no_of_nodes);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, mode, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

    /* start from real_root and go BFS */
    IGRAPH_CHECK(igraph_dqueue_push(&q, real_root));
    while (!igraph_dqueue_empty(&q)) {
        long int actnode = (long int) igraph_dqueue_pop(&q);
        neis = igraph_adjlist_get(&allneis, actnode);
        n = igraph_vector_int_size(neis);
        VECTOR(visited)[actnode] = 1;
        for (j = 0; j < n; j++) {
            long int neighbor = (long int) VECTOR(*neis)[j];
            if (!(long int)VECTOR(visited)[neighbor]) {
                IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
            }
        }
    }

    for (j = 0; j < no_of_nodes; j++) {
        no_of_newedges += 1 - VECTOR(visited)[j];
    }

    /* if any nodes are unreachable, add edges between them and real_root */
    if (no_of_newedges != 0) {

        igraph_vector_resize(pnewedges, no_of_newedges * 2);
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

    igraph_dqueue_destroy(&q);
    igraph_adjlist_destroy(&allneis);
    igraph_vector_destroy(&visited);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}


/* Internal structure for Reingold-Tilford layout */
struct igraph_i_reingold_tilford_vertex {
    long int parent;        /* Parent node index */
    long int level;         /* Level of the node */
    igraph_real_t offset;     /* X offset from parent node */
    long int left_contour;  /* Next left node of the contour
              of the subtree rooted at this node */
    long int right_contour; /* Next right node of the contour
              of the subtree rooted at this node */
    igraph_real_t offset_to_left_contour;  /* X offset when following the left contour */
    igraph_real_t offset_to_right_contour;  /* X offset when following the right contour */
    long int left_extreme;  /* Leftmost node on the deepest layer of the subtree rooted at this node */
    long int right_extreme; /* Rightmost node on the deepest layer of the subtree rooted at this node */
    igraph_real_t offset_to_left_extreme;  /* X offset when jumping to the left extreme node */
    igraph_real_t offset_to_right_extreme;  /* X offset when jumping to the right extreme node */
};

static int igraph_i_layout_reingold_tilford_postorder(struct igraph_i_reingold_tilford_vertex *vdata,
                                                      long int node, long int vcount);
static int igraph_i_layout_reingold_tilford_calc_coords(struct igraph_i_reingold_tilford_vertex *vdata,
                                                        igraph_matrix_t *res, long int node,
                                                        long int vcount, igraph_real_t xpos);

/* uncomment the next line for debugging the Reingold-Tilford layout */
/* #define LAYOUT_RT_DEBUG 1 */

static int igraph_i_layout_reingold_tilford(const igraph_t *graph,
                                            igraph_matrix_t *res,
                                            igraph_neimode_t mode,
                                            long int root) {
    long int no_of_nodes = igraph_vcount(graph);
    long int i, n, j;
    igraph_dqueue_t q = IGRAPH_DQUEUE_NULL;
    igraph_adjlist_t allneis;
    igraph_vector_int_t *neis;
    struct igraph_i_reingold_tilford_vertex *vdata;

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, mode, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

    vdata = IGRAPH_CALLOC(no_of_nodes, struct igraph_i_reingold_tilford_vertex);
    if (vdata == 0) {
        IGRAPH_ERROR("igraph_layout_reingold_tilford failed", IGRAPH_ENOMEM);
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
    IGRAPH_CHECK(igraph_dqueue_push(&q, root));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
    while (!igraph_dqueue_empty(&q)) {
        long int actnode = (long int) igraph_dqueue_pop(&q);
        long int actdist = (long int) igraph_dqueue_pop(&q);
        neis = igraph_adjlist_get(&allneis, actnode);
        n = igraph_vector_int_size(neis);

        for (j = 0; j < n; j++) {
            long int neighbor = (long int) VECTOR(*neis)[j];
            if (vdata[neighbor].parent >= 0) {
                continue;
            }
            MATRIX(*res, neighbor, 1) = actdist + 1;
            IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
            IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
            vdata[neighbor].parent = actnode;
            vdata[neighbor].level = actdist + 1;
        }
    }

    /* Step 2: postorder tree traversal, determines the appropriate X
     * offsets for every node */
    igraph_i_layout_reingold_tilford_postorder(vdata, root, no_of_nodes);

    /* Step 3: calculate real coordinates based on X offsets */
    igraph_i_layout_reingold_tilford_calc_coords(vdata, res, root, no_of_nodes, vdata[root].offset);

    igraph_dqueue_destroy(&q);
    igraph_adjlist_destroy(&allneis);
    igraph_free(vdata);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_PROGRESS("Reingold-Tilford tree layout", 100.0, NULL);

#ifdef LAYOUT_RT_DEBUG
    for (i = 0; i < no_of_nodes; i++) {
        printf(
            "%3ld: offset = %.2f, contours = [%ld, %ld], contour offsets = [%.2f, %.2f]\n",
            i, vdata[i].offset,
            vdata[i].left_contour, vdata[i].right_contour,
            vdata[i].offset_to_left_contour, vdata[i].offset_to_right_contour
        );
        if (vdata[i].left_extreme != i || vdata[i].right_extreme != i) {
            printf(
                "     extrema = [%ld, %ld], offsets to extrema = [%.2f, %.2f]\n",
                vdata[i].left_extreme, vdata[i].right_extreme,
                vdata[i].offset_to_left_extreme, vdata[i].offset_to_right_extreme
            );
        }
    }
#endif

    return 0;
}

static int igraph_i_layout_reingold_tilford_calc_coords(
        struct igraph_i_reingold_tilford_vertex *vdata,
        igraph_matrix_t *res, long int node,
        long int vcount, igraph_real_t xpos) {
    long int i;
    MATRIX(*res, node, 0) = xpos;
    for (i = 0; i < vcount; i++) {
        if (i == node) {
            continue;
        }
        if (vdata[i].parent == node) {
            igraph_i_layout_reingold_tilford_calc_coords(vdata, res, i, vcount,
                    xpos + vdata[i].offset);
        }
    }
    return 0;
}

static int igraph_i_layout_reingold_tilford_postorder(
        struct igraph_i_reingold_tilford_vertex *vdata,
        long int node, long int vcount) {
    long int i, j, childcount, leftroot, leftrootidx;
    const igraph_real_t minsep = 1;
    igraph_real_t avg;

#ifdef LAYOUT_RT_DEBUG
    printf("Starting visiting node %ld\n", node);
#endif

    /* Check whether this node is a leaf node */
    childcount = 0;
    for (i = 0; i < vcount; i++) {
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
        return 0;
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
    printf("Visited node %ld and arranged its subtrees\n", node);
#endif
    for (i = 0, j = 0; i < vcount; i++) {
        if (i == node) {
            continue;
        }
        if (vdata[i].parent == node) {
            if (leftroot >= 0) {
                /* Now we will follow the right contour of leftroot and the
                 * left contour of the subtree rooted at i */
                long lnode, rnode, auxnode;
                igraph_real_t loffset, roffset, rootsep, newoffset;

#ifdef LAYOUT_RT_DEBUG
                printf("  Placing child %ld on level %ld, to the right of %ld\n", i, vdata[i].level, leftroot);
#endif
                lnode = leftroot; rnode = i;
                rootsep = vdata[leftroot].offset + minsep;
                loffset = vdata[leftroot].offset; roffset = loffset + minsep;

                /* Keep on updating the right contour now that we have attached
                 * a new node to the subtree being built */
                vdata[node].right_contour = i;
                vdata[node].offset_to_right_contour = rootsep;

#ifdef LAYOUT_RT_DEBUG
                printf("    Contour: [%ld, %ld], offsets: [%lf, %lf], rootsep: %lf\n",
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
                            printf("      Left subtree ended earlier, continuing left subtree's left and right contour on right subtree (node %ld gets connected to node %ld)\n", auxnode, vdata[rnode].left_contour);
                            printf("      New contour following offset for node %ld is %lf\n", auxnode, vdata[auxnode].offset_to_left_contour);
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
                            printf("      Right subtree ended earlier, continuing right subtree's left and right contour on left subtree (node %ld gets connected to node %ld)\n", auxnode, lnode);
                            printf("      New contour following offset for node %ld is %lf\n", auxnode, vdata[auxnode].offset_to_left_contour);
#endif
                        }
                        rnode = -1;
                    }
#ifdef LAYOUT_RT_DEBUG
                    printf("    Contour: [%ld, %ld], offsets: [%lf, %lf], rootsep: %lf\n",
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
                printf("  Offset of subtree with root node %ld will be %lf\n", i, rootsep);
#endif
                vdata[i].offset = rootsep;
                vdata[node].offset_to_right_contour = rootsep;
                avg = (avg * j) / (j + 1) + rootsep / (j + 1);
                leftrootidx = j;
                leftroot = i;
            } else {
                /* This is the first child of the node being considered so we
                 * can simply place the subtree on our virtual canvas */
#ifdef LAYOUT_RT_DEBUG
                printf("  Placing child %ld on level %ld as first child\n", i, vdata[i].level);
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
    printf("Shifting node %ld to be centered above children. Shift amount: %lf\n", node, avg);
#endif
    vdata[node].offset_to_left_contour -= avg;
    vdata[node].offset_to_right_contour -= avg;
    vdata[node].offset_to_left_extreme -= avg;
    vdata[node].offset_to_right_extreme -= avg;
    for (i = 0, j = 0; i < vcount; i++) {
        if (i == node) {
            continue;
        }
        if (vdata[i].parent == node) {
            vdata[i].offset -= avg;
        }
    }

    return 0;
}

/**
 * \function igraph_layout_reingold_tilford
 * \brief Reingold-Tilford layout for tree graphs
 *
 * </para><para>
 * Arranges the nodes in a tree where the given node is used as the root.
 * The tree is directed downwards and the parents are centered above its
 * children. For the exact algorithm, see:
 *
 * </para><para>
 * Reingold, E and Tilford, J: Tidier drawing of trees.
 * IEEE Trans. Softw. Eng., SE-7(2):223--228, 1981
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
 * \param roots The index of the root vertex or root vertices.
 *   If this is a non-empty vector then the supplied vertex ids are used
 *   as the roots of the trees (or a single tree if the graph is connected).
 *   If it is a null pointer of a pointer to an empty vector, then the root
 *   vertices are automatically calculated based on topological sorting,
 *   performed with the opposite mode than the \p mode argument.
 *   After the vertices have been sorted, one is selected from each component.
 * \param rootlevel This argument can be useful when drawing forests which are
 *   not trees (i.e. they are unconnected and have tree components). It specifies
 *   the level of the root vertices for every tree in the forest. It is only
 *   considered if not a null pointer and the \p roots argument is also given
 *   (and it is not a null pointer of an empty vector).
 * \return Error code.
 *
 * Added in version 0.2.
 *
 * \sa \ref igraph_layout_reingold_tilford_circular().
 *
 * \example examples/simple/igraph_layout_reingold_tilford.c
 */
int igraph_layout_reingold_tilford(const igraph_t *graph,
                                   igraph_matrix_t *res,
                                   igraph_neimode_t mode,
                                   const igraph_vector_t *roots,
                                   const igraph_vector_t *rootlevel) {

    long int no_of_nodes_orig = igraph_vcount(graph);
    long int no_of_nodes = no_of_nodes_orig;
    long int real_root;
    igraph_t extended;
    const igraph_t *pextended = graph;
    igraph_vector_t myroots;
    const igraph_vector_t *proots = roots;
    igraph_neimode_t mode2;
    long int i;
    igraph_vector_t newedges;

    /* TODO: possible speedup could be achieved if we use a table for storing
     * the children of each node in the tree. (Now the implementation uses a
     * single array containing the parent of each node and a node's children
     * are determined by looking for other nodes that have this node as parent)
     */

    /* at various steps it might be necessary to add edges to the graph */
    IGRAPH_VECTOR_INIT_FINALLY(&newedges, 0);

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    if ( (!roots || igraph_vector_size(roots) == 0) &&
         rootlevel && igraph_vector_size(rootlevel) != 0 ) {
        IGRAPH_WARNING("Reingold-Tilford layout: 'rootlevel' ignored");
    }

    /* ----------------------------------------------------------------------- */
    /* If root vertices are not given, then do a topological sort and take
       the last element from every component for directed graphs and mode == out,
       or the first element from every component for directed graphs and mode ==
       in,or select the vertex with the maximum degree from each component for
       undirected graphs */

    if (!roots || igraph_vector_size(roots) == 0) {

        igraph_vector_t order, membership;
        igraph_integer_t no_comps;
        long int i, noseen = 0;

        IGRAPH_VECTOR_INIT_FINALLY(&myroots, 0);
        IGRAPH_VECTOR_INIT_FINALLY(&order, no_of_nodes);
        IGRAPH_VECTOR_INIT_FINALLY(&membership, no_of_nodes);

        if (mode != IGRAPH_ALL) {
            /* look for roots by swimming against the stream */
            mode2 = (mode == IGRAPH_IN) ? IGRAPH_OUT : IGRAPH_IN;

            IGRAPH_CHECK(igraph_topological_sorting(graph, &order, mode2));
            IGRAPH_CHECK(igraph_clusters(graph, &membership, /*csize=*/ 0,
                                         &no_comps, IGRAPH_WEAK));
        } else {
            IGRAPH_CHECK(igraph_sort_vertex_ids_by_degree(graph, &order,
                         igraph_vss_all(), IGRAPH_ALL, 0, IGRAPH_ASCENDING, 0));
            IGRAPH_CHECK(igraph_clusters(graph, &membership, /*csize=*/ 0,
                                         &no_comps, IGRAPH_WEAK));
        }

        IGRAPH_CHECK(igraph_vector_resize(&myroots, no_comps));

        /* go backwards and fill the roots vector with indices [1, no_of_nodes]
           The index 0 is used to signal this root has not been found yet:
           all indices are then decreased by one to [0, no_of_nodes - 1] */
        igraph_vector_null(&myroots);
        proots = &myroots;
        for (i = no_of_nodes - 1; noseen < no_comps && i >= 0; i--) {
            long int v = (long int) VECTOR(order)[i];
            long int mem = (long int) VECTOR(membership)[v];
            if (VECTOR(myroots)[mem] == 0) {
                noseen += 1;
                VECTOR(myroots)[mem] = v + 1;
            }
        }
        for (i = 0; i < no_comps; i++) {
            VECTOR(myroots)[i] -= 1;
        }

        igraph_vector_destroy(&membership);
        igraph_vector_destroy(&order);
        IGRAPH_FINALLY_CLEAN(2);

    } else if (rootlevel && igraph_vector_size(rootlevel) > 0 &&
               igraph_vector_size(roots) > 1) {

        /* ----------------------------------------------------------------------- */
        /* Many roots were given to us, check 'rootlevel' */

        long int plus_levels = 0;
        long int i;

        if (igraph_vector_size(roots) != igraph_vector_size(rootlevel)) {
            IGRAPH_ERROR("Reingold-Tilford: 'roots' and 'rootlevel' lengths differ",
                         IGRAPH_EINVAL);
        }

        /* count the rootlevels that are not zero */
        for (i = 0; i < igraph_vector_size(roots); i++) {
            plus_levels += VECTOR(*rootlevel)[i];
        }

        /* make copy of graph, add vertices/edges */
        if (plus_levels != 0) {
            long int edgeptr = 0;

            pextended = &extended;
            IGRAPH_CHECK(igraph_copy(&extended, graph));
            IGRAPH_FINALLY(igraph_destroy, &extended);
            IGRAPH_CHECK(igraph_add_vertices(&extended,
                                             (igraph_integer_t) plus_levels, 0));

            igraph_vector_resize(&newedges, plus_levels * 2);

            for (i = 0; i < igraph_vector_size(roots); i++) {
                long int rl = (long int) VECTOR(*rootlevel)[i];
                long int rn = (long int) VECTOR(*roots)[i];
                long int j;

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
    if (igraph_vector_size(proots) == 1) {
        real_root = (long int) VECTOR(*proots)[0];
        if (real_root < 0 || real_root >= no_of_nodes) {
            IGRAPH_ERROR("Invalid vertex id.", IGRAPH_EINVVID);
        }

        /* else, we need to make real_root */
    } else {
        long int no_of_newedges;

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
        no_of_newedges = igraph_vector_size(proots);
        igraph_vector_resize(&newedges, no_of_newedges * 2);
        for (i = 0; i < no_of_newedges; i++) {
            VECTOR(newedges)[2 * i] = no_of_nodes - 1;
            VECTOR(newedges)[2 * i + 1] = VECTOR(*proots)[i];
        }

        IGRAPH_CHECK(igraph_add_edges(&extended, &newedges, 0));
    }

    /* prepare edges to unreachable parts of the graph */
    IGRAPH_CHECK(igraph_i_layout_reingold_tilford_unreachable(pextended, mode, real_root, no_of_nodes, &newedges));

    if (igraph_vector_size(&newedges) != 0) {
        /* Make copy of the graph unless it exists already */
        if (pextended == graph) {
            pextended = &extended;
            IGRAPH_CHECK(igraph_copy(&extended, graph));
            IGRAPH_FINALLY(igraph_destroy, &extended);
        }

        IGRAPH_CHECK(igraph_add_edges(&extended, &newedges, 0));
    }
    igraph_vector_destroy(&newedges);
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
            long int i;
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
        igraph_vector_destroy(&myroots);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/**
 * \function igraph_layout_reingold_tilford_circular
 * \brief Circular Reingold-Tilford layout for trees
 *
 * </para><para>
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
 * \param roots The index of the root vertex or root vertices.
 *   If this is a non-empty vector then the supplied vertex ids are used
 *   as the roots of the trees (or a single tree if the graph is connected).
 *   If it is a null pointer of a pointer to an empty vector, then the root
 *   vertices are automatically calculated based on topological sorting,
 *   performed with the opposite mode than the \p mode argument.
 *   After the vertices have been sorted, one is selected from each component.
 * \param rootlevel This argument can be useful when drawing forests which are
 *   not trees (i.e. they are unconnected and have tree components). It specifies
 *   the level of the root vertices for every tree in the forest. It is only
 *   considered if not a null pointer and the \p roots argument is also given
 *   (and it is not a null pointer or an empty vector).
 * \return Error code.
 *
 * \sa \ref igraph_layout_reingold_tilford().
 */
int igraph_layout_reingold_tilford_circular(const igraph_t *graph,
        igraph_matrix_t *res,
        igraph_neimode_t mode,
        const igraph_vector_t *roots,
        const igraph_vector_t *rootlevel) {

    long int no_of_nodes = igraph_vcount(graph);
    long int i;
    igraph_real_t ratio;
    igraph_real_t minx, maxx;

    IGRAPH_CHECK(igraph_layout_reingold_tilford(graph, res, mode, roots, rootlevel));

    if (no_of_nodes == 0) {
        return IGRAPH_SUCCESS;
    }

    ratio = 2 * M_PI * (no_of_nodes - 1.0) / no_of_nodes;

    minx = maxx = MATRIX(*res, 0, 0);
    for (i = 1; i < no_of_nodes; i++) {
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
    for (i = 0; i < no_of_nodes; i++) {
        igraph_real_t phi = (MATRIX(*res, i, 0) - minx) * ratio;
        igraph_real_t r = MATRIX(*res, i, 1);
        MATRIX(*res, i, 0) = r * cos(phi);
        MATRIX(*res, i, 1) = r * sin(phi);
    }

    return IGRAPH_SUCCESS;
}
