/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_cohesive_blocks.h"

#include "igraph_constructors.h"
#include "igraph_dqueue.h"
#include "igraph_flow.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_operators.h"
#include "igraph_separators.h"
#include "igraph_statusbar.h"
#include "igraph_structural.h"

#include "core/interruption.h"

static void igraph_i_cohesive_blocks_free_graphs(igraph_vector_ptr_t *ptr) {
    igraph_integer_t i, n = igraph_vector_ptr_size(ptr);

    for (i = 0; i < n; i++) {
        igraph_t *g = VECTOR(*ptr)[i];
        if (g) {
            igraph_destroy(g);
            IGRAPH_FREE(VECTOR(*ptr)[i]); /* also sets it to NULL */
        }
    }
}

/* This is kind of a BFS to find the components of the graph, after
 * deleting the vertices marked in 'excluded'.
 * These vertices are not put in the BFS queue, but they are added to
 * all neighboring components.
 */

static igraph_error_t igraph_i_cb_components(igraph_t *graph,
                                  const igraph_vector_bool_t *excluded,
                                  igraph_vector_int_t *components,
                                  igraph_integer_t *no,
                                  /* working area follows */
                                  igraph_vector_int_t *compid,
                                  igraph_dqueue_int_t *Q,
                                  igraph_vector_int_t *neis) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t i;
    igraph_integer_t cno = 0;

    igraph_vector_int_clear(components);
    igraph_dqueue_int_clear(Q);
    IGRAPH_CHECK(igraph_vector_int_resize(compid, no_of_nodes));
    igraph_vector_int_null(compid);

    for (i = 0; i < no_of_nodes; i++) {

        if (VECTOR(*compid)[i])   {
            continue;
        }
        if (VECTOR(*excluded)[i]) {
            continue;
        }

        IGRAPH_CHECK(igraph_dqueue_int_push(Q, i));
        IGRAPH_CHECK(igraph_vector_int_push_back(components, i));
        VECTOR(*compid)[i] = ++cno;

        while (!igraph_dqueue_int_empty(Q)) {
            igraph_integer_t node = igraph_dqueue_int_pop(Q);
            igraph_integer_t j, n;
            IGRAPH_CHECK(igraph_neighbors(graph, neis, node, IGRAPH_ALL));
            n = igraph_vector_int_size(neis);
            for (j = 0; j < n; j++) {
                igraph_integer_t v = VECTOR(*neis)[j];
                if (VECTOR(*excluded)[v]) {
                    if (VECTOR(*compid)[v] != cno) {
                        VECTOR(*compid)[v] = cno;
                        IGRAPH_CHECK(igraph_vector_int_push_back(components, v));
                    }
                } else {
                    if (VECTOR(*compid)[v] == 0) {
                        VECTOR(*compid)[v] = cno; /* could be anything positive */
                        IGRAPH_CHECK(igraph_vector_int_push_back(components, v));
                        IGRAPH_CHECK(igraph_dqueue_int_push(Q, v));
                    }
                }
            }
        } /* while !igraph_dqueue_int_empty */

        IGRAPH_CHECK(igraph_vector_int_push_back(components, -1));

    } /* for i<no_of_nodes */

    *no = cno;

    return IGRAPH_SUCCESS;
}

static igraph_bool_t igraph_i_cb_isin(const igraph_vector_int_t *needle,
                                      const igraph_vector_int_t *haystack) {
    igraph_integer_t nlen = igraph_vector_int_size(needle);
    igraph_integer_t hlen = igraph_vector_int_size(haystack);
    igraph_integer_t np = 0, hp = 0;

    if (hlen < nlen) {
        return false;
    }

    while (np < nlen && hp < hlen) {
        if (VECTOR(*needle)[np] == VECTOR(*haystack)[hp]) {
            np++; hp++;
        } else if (VECTOR(*needle)[np] < VECTOR(*haystack)[hp]) {
            return false;
        } else {
            hp++;
        }
    }

    return np == nlen;
}

/**
 * \function igraph_cohesive_blocks
 * \brief Identifies the hierarchical cohesive block structure of a graph.
 *
 * Cohesive blocking is a method of determining hierarchical subsets of
 * graph vertices based on their structural cohesion (or vertex
 * connectivity). For a given graph G, a subset of its vertices
 * S is said to be maximally k-cohesive if there is
 * no superset of S with vertex connectivity greater than or equal to k.
 * Cohesive blocking is a process through which, given a
 * k-cohesive set of vertices, maximally l-cohesive subsets are
 * recursively identified with l>k. Thus a hiearchy of vertex subsets
 * is found, with the entire graph G at its root.
 *
 * </para><para>
 * This function implements cohesive blocking and
 * calculates the complete cohesive block hierarchy of a graph.
 *
 * </para><para>
 * See the following reference for details:
 *
 * </para><para>
 * J. Moody and D. R. White. Structural
 * cohesion and embeddedness: A hierarchical concept of social
 * groups. American Sociological Review, 68(1):103--127, Feb 2003.
 * https://doi.org/10.2307/3088904
 *
 * \param graph The input graph. It must be undirected and simple. See
 *    \ref igraph_is_simple().
 * \param blocks If not a null pointer, then it must be an initialized
 *    list of integers vectors; the cohesive blocks will be stored here.
 *    Each block is encoded with a vector of type \ref igraph_vector_int_t that
 *    contains the vertex IDs of the block.
 * \param cohesion If not a null pointer, then it must be an initialized
 *    vector and the cohesion of the blocks is stored here, in the same
 *    order as the blocks in the \p blocks vector list.
 * \param parent If not a null pointer, then it must be an initialized
 *    vector and the block hierarchy is stored here. For each block, the
 *    ID (i.e. the position in the \p blocks vector list) of its
 *    parent block is stored. For the top block in the hierarchy,
 *    <code>-1</code> is stored.
 * \param block_tree If not a null pointer, then it must be a pointer
 *    to an uninitialized graph, and the block hierarchy is stored
 *    here as an igraph graph. The vertex IDs correspond to the order
 *    of the blocks in the \p blocks vector.
 * \return Error code.
 *
 * Time complexity: TODO.
 *
 * \example examples/simple/cohesive_blocks.c
 */

igraph_error_t igraph_cohesive_blocks(const igraph_t *graph,
                           igraph_vector_int_list_t *blocks,
                           igraph_vector_int_t *cohesion,
                           igraph_vector_int_t *parent,
                           igraph_t *block_tree) {

    /* Some implementation comments. Everything is relatively
       straightforward, except, that we need to follow the vertex IDs
       of the various subgraphs, without having to store two-way
       mappings at each level. The subgraphs can overlap, this
       complicates things a bit.

       The 'Q' vector is used as a double ended queue and it contains
       the subgraphs to work on in the future. Some other vectors are
       associated with it. 'Qparent' gives the parent graph of a graph
       in Q. Qmapping gives the mapping of the vertices from the graph
       to the parent graph. Qcohesion is the vertex connectivity of the
       graph.

       Qptr is an integer and points to the next graph to work on.
    */

    /* In theory, Q could be an igraph_graph_list_t; however, in that case
     * we would not be able to pop off graphs from the front of the list as
     * all elements of an igraph_graph_list_t are expected to be initialized,
     * valid graphs. That's why we use an igraph_vector_ptr_t instead. */

    igraph_vector_ptr_t Q;
    igraph_vector_int_list_t Qmapping;
    igraph_vector_int_t Qparent;
    igraph_vector_int_t Qcohesion;
    igraph_vector_bool_t Qcheck;
    igraph_integer_t Qptr = 0;
    igraph_integer_t conn;
    igraph_bool_t is_simple;

    igraph_t *graph_copy;

    igraph_vector_int_list_t separators;
    igraph_vector_int_t compvertices;
    igraph_vector_int_t components;
    igraph_vector_int_t newmapping;
    igraph_vector_bool_t marked;

    igraph_vector_int_t compid;
    igraph_dqueue_int_t bfsQ;
    igraph_vector_int_t neis;

    if (igraph_is_directed(graph)) {
        IGRAPH_ERROR("Cohesive blocking only works on undirected graphs.",
                     IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_is_simple(graph, &is_simple));
    if (!is_simple) {
        IGRAPH_ERROR("Cohesive blocking only works on simple graphs.",
                     IGRAPH_EINVAL);
    }

    if (blocks)   {
        igraph_vector_int_list_clear(blocks);
    }
    if (cohesion) {
        igraph_vector_int_clear(cohesion);
    }
    if (parent)   {
        igraph_vector_int_clear(parent);
    }

    IGRAPH_CHECK(igraph_vector_ptr_init(&Q, 1));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &Q);
    IGRAPH_FINALLY(igraph_i_cohesive_blocks_free_graphs, &Q);

    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&Qmapping, 1);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&Qparent, 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&Qcohesion, 1);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&Qcheck, 1);

    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&separators, 0);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&compvertices, 0);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&marked, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_dqueue_int_init(&bfsQ, 100));
    IGRAPH_FINALLY(igraph_dqueue_int_destroy, &bfsQ);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&compid, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&components, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&newmapping, 0);

    /* Put the input graph in the queue */
    graph_copy = IGRAPH_CALLOC(1, igraph_t);
    IGRAPH_CHECK_OOM(graph_copy, "Insufficient memory for cohesive blocking.");

    IGRAPH_CHECK(igraph_copy(graph_copy, graph));
    VECTOR(Q)[0] = graph_copy;
    VECTOR(Qparent)[0] = -1;  /* Has no parent */
    IGRAPH_CHECK(igraph_vertex_connectivity(graph, &conn, /*checks=*/ true));
    VECTOR(Qcohesion)[0] = conn;
    VECTOR(Qcheck)[0] = false;

    /* Then work until the queue is empty */
    while (Qptr < igraph_vector_ptr_size(&Q)) {
        igraph_t *mygraph = VECTOR(Q)[Qptr];
        igraph_bool_t mycheck = VECTOR(Qcheck)[Qptr];
        igraph_integer_t mynodes = igraph_vcount(mygraph);
        igraph_integer_t i, nsep;
        igraph_integer_t no, kept = 0;
        igraph_integer_t cptr = 0;
        igraph_integer_t nsepv = 0;
        igraph_bool_t addedsep = false;

        IGRAPH_ALLOW_INTERRUPTION();

        /* Get the separators */
        IGRAPH_CHECK(igraph_minimum_size_separators(mygraph, &separators));
        nsep = igraph_vector_int_list_size(&separators);

        /* Remove them from the graph, also mark them */
        IGRAPH_CHECK(igraph_vector_bool_resize(&marked, mynodes));
        igraph_vector_bool_null(&marked);
        for (i = 0; i < nsep; i++) {
            igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(&separators, i);
            igraph_integer_t j, n = igraph_vector_int_size(v);
            for (j = 0; j < n; j++) {
                igraph_integer_t vv = VECTOR(*v)[j];
                if (!VECTOR(marked)[vv]) {
                    nsepv++;
                    VECTOR(marked)[vv] = true;
                }
            }
        }

        /* Find the connected components, omitting the separator vertices,
           but including the neighboring separator vertices
         */
        IGRAPH_CHECK(igraph_i_cb_components(mygraph, &marked,
                                            &components, &no,
                                            &compid, &bfsQ, &neis));

        /* Add the separator vertices themselves, as another component,
           but only if there is at least one vertex not included in any
           separator. */
        if (nsepv != mynodes) {
            addedsep = true;
            for (i = 0; i < mynodes; i++) {
                if (VECTOR(marked)[i]) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(&components, i));
                }
            }
            IGRAPH_CHECK(igraph_vector_int_push_back(&components, -1));
            no++;
        }

        for (i = 0; i < no; i++) {
            igraph_t *newgraph;
            igraph_integer_t maxdeg;

            igraph_vector_int_clear(&compvertices);

            while (true) {
                igraph_integer_t v = VECTOR(components)[cptr++];
                if (v < 0) {
                    break;
                }
                IGRAPH_CHECK(igraph_vector_int_push_back(&compvertices, v));
            }

            newgraph = IGRAPH_CALLOC(1, igraph_t);
            IGRAPH_CHECK_OOM(newgraph, "Insufficient memory for cohesive blocking.");
            IGRAPH_FINALLY(igraph_free, newgraph);

            IGRAPH_CHECK(igraph_induced_subgraph_map(mygraph, newgraph,
                                                     igraph_vss_vector(&compvertices),
                                                     IGRAPH_SUBGRAPH_AUTO,
                                                     /*map=*/ NULL,
                                                     /*invmap=*/ &newmapping));
            IGRAPH_FINALLY(igraph_destroy, newgraph);

            IGRAPH_CHECK(igraph_maxdegree(newgraph, &maxdeg, igraph_vss_all(),
                                          IGRAPH_ALL, IGRAPH_LOOPS));
            if (maxdeg > VECTOR(Qcohesion)[Qptr]) {
                igraph_integer_t newconn;
                kept++;
                IGRAPH_CHECK(igraph_vector_ptr_push_back(&Q, newgraph));
                IGRAPH_FINALLY_CLEAN(2);
                IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(&Qmapping, &newmapping));
                IGRAPH_CHECK(igraph_vertex_connectivity(newgraph, &newconn,
                                                        /*checks=*/ 1));
                IGRAPH_CHECK(igraph_vector_int_push_back(&Qcohesion, newconn));
                IGRAPH_CHECK(igraph_vector_int_push_back(&Qparent, Qptr));
                IGRAPH_CHECK(igraph_vector_bool_push_back(&Qcheck,
                             mycheck || addedsep));
            } else {
                igraph_destroy(newgraph);
                igraph_free(newgraph);
                IGRAPH_FINALLY_CLEAN(2);
            }
        }

        igraph_destroy(mygraph);
        igraph_free(mygraph);
        VECTOR(Q)[Qptr] = NULL;

        Qptr++;
    }

    igraph_vector_int_destroy(&newmapping);
    igraph_vector_int_destroy(&components);
    igraph_vector_int_destroy(&compid);
    igraph_dqueue_int_destroy(&bfsQ);
    igraph_vector_int_destroy(&neis);
    igraph_vector_bool_destroy(&marked);
    igraph_vector_int_destroy(&compvertices);
    igraph_vector_int_list_destroy(&separators);
    IGRAPH_FINALLY_CLEAN(8);

    if (blocks || cohesion || parent || block_tree) {
        igraph_integer_t noblocks = Qptr, badblocks = 0;
        igraph_vector_bool_t removed;
        igraph_integer_t i, resptr = 0;
        igraph_vector_int_t rewritemap;

        IGRAPH_CHECK(igraph_vector_bool_init(&removed, noblocks));
        IGRAPH_FINALLY(igraph_vector_bool_destroy, &removed);
        IGRAPH_CHECK(igraph_vector_int_init(&rewritemap, noblocks));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &rewritemap);

        for (i = 1; i < noblocks; i++) {
            igraph_integer_t p = VECTOR(Qparent)[i];
            while (VECTOR(removed)[p]) {
                p = VECTOR(Qparent)[p];
            }
            if (VECTOR(Qcohesion)[p] >= VECTOR(Qcohesion)[i]) {
                VECTOR(removed)[i] = true;
                badblocks++;
            }
        }

        /* Rewrite the mappings */
        for (i = 1; i < Qptr; i++) {
            igraph_integer_t j, n, p = VECTOR(Qparent)[i];
            igraph_vector_int_t *mapping, *pmapping;

            if (p == 0) {
                continue;
            }

            mapping = igraph_vector_int_list_get_ptr(&Qmapping, i);
            pmapping = igraph_vector_int_list_get_ptr(&Qmapping, p);

            n = igraph_vector_int_size(mapping);
            for (j = 0; j < n; j++) {
                igraph_integer_t v = VECTOR(*mapping)[j];
                VECTOR(*mapping)[j] = VECTOR(*pmapping)[v];
            }
        }

        /* Because we also put the separator vertices in the queue, it is
           not ensured that the found blocks are not subsets of each other.
           We check this now. */
        for (i = 1; i < noblocks; i++) {
            igraph_integer_t j, ic;
            igraph_vector_int_t *ivec;
            if (!VECTOR(Qcheck)[i] || VECTOR(removed)[i]) {
                continue;
            }
            ivec = igraph_vector_int_list_get_ptr(&Qmapping, i);
            ic = VECTOR(Qcohesion)[i];
            for (j = 1; j < noblocks; j++) {
                igraph_vector_int_t *jvec;
                igraph_integer_t jc;
                if (j == i || !VECTOR(Qcheck)[j] || VECTOR(removed)[j]) {
                    continue;
                }
                jvec = igraph_vector_int_list_get_ptr(&Qmapping, j);
                jc = VECTOR(Qcohesion)[j];
                if (igraph_i_cb_isin(ivec, jvec) && jc >= ic) {
                    badblocks++;
                    VECTOR(removed)[i] = true;
                    break;
                }
            }
        }

        noblocks -= badblocks;

        if (blocks) {
            IGRAPH_CHECK(igraph_vector_int_list_resize(blocks, noblocks));
        }
        if (cohesion) {
            IGRAPH_CHECK(igraph_vector_int_resize(cohesion, noblocks));
        }
        if (parent) {
            IGRAPH_CHECK(igraph_vector_int_resize(parent, noblocks));
        }

        for (i = 0; i < Qptr; i++) {
            if (VECTOR(removed)[i]) {
                continue;
            }
            VECTOR(rewritemap)[i] = resptr;
            if (cohesion) {
                VECTOR(*cohesion)[resptr] = VECTOR(Qcohesion)[i];
            }
            if (parent || block_tree) {
                igraph_integer_t p = VECTOR(Qparent)[i];
                while (p >= 0 && VECTOR(removed)[p]) {
                    p = VECTOR(Qparent)[p];
                }
                if (p >= 0) {
                    p = VECTOR(rewritemap)[p];
                }
                VECTOR(Qparent)[i] = p;
                if (parent) {
                    VECTOR(*parent)[resptr] = p;
                }
            }
            if (blocks) {
                IGRAPH_CHECK(
                    igraph_vector_int_update(
                        igraph_vector_int_list_get_ptr(blocks, resptr),
                        igraph_vector_int_list_get_ptr(&Qmapping, i)
                    )
                );
                igraph_vector_int_clear(igraph_vector_int_list_get_ptr(&Qmapping, i));
            }
            resptr++;
        }

        /* Plus the original graph */
        if (blocks) {
            igraph_integer_t num_vertices = igraph_vcount(graph);
            igraph_vector_int_t *orig = igraph_vector_int_list_get_ptr(blocks, 0);
            IGRAPH_CHECK(igraph_vector_int_resize(orig, num_vertices));
            for (i = 0; i < num_vertices; i++) {
                VECTOR(*orig)[i] = i;
            }
        }

        if (block_tree) {
            igraph_vector_int_t edges;
            igraph_integer_t eptr = 0;
            IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, noblocks * 2 - 2);
            for (i = 1; i < Qptr; i++) {
                if (VECTOR(removed)[i]) {
                    continue;
                }
                VECTOR(edges)[eptr++] = VECTOR(Qparent)[i];
                VECTOR(edges)[eptr++] = VECTOR(rewritemap)[i];
            }

            IGRAPH_CHECK(igraph_create(block_tree, &edges, noblocks,
                                       IGRAPH_DIRECTED));
            igraph_vector_int_destroy(&edges);
            IGRAPH_FINALLY_CLEAN(1);
        }

        igraph_vector_int_destroy(&rewritemap);
        igraph_vector_bool_destroy(&removed);
        IGRAPH_FINALLY_CLEAN(2);

    }

    igraph_vector_bool_destroy(&Qcheck);
    igraph_vector_int_destroy(&Qcohesion);
    igraph_vector_int_destroy(&Qparent);
    igraph_vector_int_list_destroy(&Qmapping);
    IGRAPH_FINALLY_CLEAN(4);

    igraph_vector_ptr_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(2); /* + the elements of Q, they were already destroyed */

    return IGRAPH_SUCCESS;
}
