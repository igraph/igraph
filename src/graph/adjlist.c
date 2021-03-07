/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_adjlist.h"
#include "igraph_memory.h"
#include "igraph_interface.h"

#include "core/interruption.h"

#include <string.h>   /* memset */
#include <stdio.h>

/**
 * Helper function that simplifies a sorted adjacency vector by removing
 * duplicate elements and optionally self-loops.
 */
static int igraph_i_simplify_sorted_int_adjacency_vector_in_place(
    igraph_vector_int_t *v, igraph_integer_t index, igraph_neimode_t mode,
    igraph_loops_t loops, igraph_multiple_t multiple
);

/**
 * Helper function that removes loops from an incidence vector (either both
 * occurrences or only one of them).
 */
static int igraph_i_remove_loops_from_incidence_vector_in_place(
    igraph_vector_int_t *v, const igraph_t *graph, igraph_loops_t loops
);

/**
 * \section about_adjlists
 * <para>Sometimes it is easier to work with a graph which is in
 * adjacency list format: a list of vectors; each vector contains the
 * neighbor vertices or incident edges of a given vertex. Typically,
 * this representation is good if we need to iterate over the neighbors
 * of all vertices many times. E.g. when finding the shortest paths
 * between all pairs of vertices or calculating closeness centrality
 * for all the vertices.</para>
 *
 * <para>The <type>igraph_adjlist_t</type> stores the adjacency lists
 * of a graph. After creation it is independent of the original graph,
 * it can be modified freely with the usual vector operations, the
 * graph is not affected. E.g. the adjacency list can be used to
 * rewire the edges of a graph efficiently. If one used the
 * straightforward \ref igraph_delete_edges() and \ref
 * igraph_add_edges() combination for this that needs O(|V|+|E|) time
 * for every single deletion and insertion operation, it is thus very
 * slow if many edges are rewired. Extracting the graph into an
 * adjacency list, do all the rewiring operations on the vectors of
 * the adjacency list and then creating a new graph needs (depending
 * on how exactly the rewiring is done) typically O(|V|+|E|) time for
 * the whole rewiring process.</para>
 *
 * <para>Lazy adjacency lists are a bit different. When creating a
 * lazy adjacency list, the neighbors of the vertices are not queried,
 * only some memory is allocated for the vectors. When \ref
 * igraph_lazy_adjlist_get() is called for vertex v the first time,
 * the neighbors of v are queried and stored in a vector of the
 * adjacency list, so they don't need to be queried again. Lazy
 * adjacency lists are handy if you have an at least linear operation
 * (because initialization is generally linear in terms of the number of
 * vertices), but you don't know how many vertices you will visit
 * during the computation.
 * </para>
 *
 * <para>
 * \example examples/simple/adjlist.c
 * </para>
 */

/**
 * \function igraph_adjlist_init
 * \brief Constructs an adjacency list of vertices from a given graph.
 *
 * Creates a list of vectors containing the neighbors of all vertices
 * in a graph. The adjacency list is independent of the graph after
 * creation, e.g. the graph can be destroyed and modified, the
 * adjacency list contains the state of the graph at the time of its
 * initialization.
 * \param graph The input graph.
 * \param al Pointer to an uninitialized <type>igraph_adjlist_t</type> object.
 * \param mode Constant specifying whether outgoing
 *   (<code>IGRAPH_OUT</code>), incoming (<code>IGRAPH_IN</code>),
 *   or both (<code>IGRAPH_ALL</code>) types of neighbors to include
 *   in the adjacency list. It is ignored for undirected networks.
 * \param loops Specifies how to treat loop edges. <code>IGRAPH_NO_LOOPS</code>
 *   removes loop edges from the adjacency list. <code>IGRAPH_LOOPS_ONCE</code>
 *   makes each loop edge appear only once in the adjacency list of the
 *   corresponding vertex. <code>IGRAPH_LOOPS_TWICE</code> makes loop edges
 *   appear \em twice in the adjacency list of the corresponding vertex,
 *   but only if the graph is undirected or <code>mode</code> is set to
 *   <code>IGRAPH_ALL</code>.
 * \param multiple Specifies how to treat multiple (parallel) edges.
 *   <code>IGRAPH_NO_MULTIPLE</code> collapses parallel edges into a single one;
 *   <code>IGRAPH_MULTIPLE</code> keeps the multiplicities of parallel edges
 *   so the same vertex will appear as many times in the adjacency list of
 *   another vertex as the number of parallel edges going between the two
 *   vertices.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 */

int igraph_adjlist_init(const igraph_t *graph, igraph_adjlist_t *al,
                        igraph_neimode_t mode, igraph_loops_t loops,
                        igraph_multiple_t multiple) {
    igraph_integer_t i;
    igraph_vector_t tmp;
    int j, n;

    if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Cannot create adjacency list view", IGRAPH_EINVMODE);
    }

    igraph_vector_init(&tmp, 0);
    IGRAPH_FINALLY(igraph_vector_destroy, &tmp);

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    al->length = igraph_vcount(graph);
    al->adjs = IGRAPH_CALLOC(al->length, igraph_vector_int_t);
    if (al->adjs == 0) {
        IGRAPH_ERROR("Cannot create adjacency list view", IGRAPH_ENOMEM);
    }

    IGRAPH_FINALLY(igraph_adjlist_destroy, al);

    for (i = 0; i < al->length; i++) {
        IGRAPH_ALLOW_INTERRUPTION();
        
        IGRAPH_CHECK(igraph_neighbors(graph, &tmp, i, mode));

        n = igraph_vector_size(&tmp);
        IGRAPH_CHECK(igraph_vector_int_init(&al->adjs[i], n));
        for (j = 0; j < n; j++) {
            VECTOR(al->adjs[i])[j] = VECTOR(tmp)[j];
        }

        IGRAPH_CHECK(igraph_i_simplify_sorted_int_adjacency_vector_in_place(
            &al->adjs[i], i, mode, loops, multiple
        ));
    }

    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

/**
 * \function igraph_adjlist_init_empty
 * \brief Initializes an empty adjacency list.
 *
 * Creates a list of vectors, one for each vertex. This is useful when you
 * are \em constructing a graph using an adjacency list representation as
 * it does not require your graph to exist yet.
 * \param no_of_nodes The number of vertices
 * \param al Pointer to an uninitialized <type>igraph_adjlist_t</type> object.
 * \return Error code.
 *
 * Time complexity: O(|V|), linear in the number of vertices.
 */
int igraph_adjlist_init_empty(igraph_adjlist_t *al, igraph_integer_t no_of_nodes) {
    long int i;

    al->length = no_of_nodes;
    al->adjs = IGRAPH_CALLOC(al->length, igraph_vector_int_t);
    if (al->adjs == 0) {
        IGRAPH_ERROR("Cannot create adjlist view", IGRAPH_ENOMEM);
    }

    IGRAPH_FINALLY(igraph_adjlist_destroy, al);
    for (i = 0; i < al->length; i++) {
        IGRAPH_CHECK(igraph_vector_int_init(&al->adjs[i], 0));
    }
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_adjlist_init_complementer
 * \brief Adjacency lists for the complementer graph.
 *
 * This function creates adjacency lists for the complementer
 * of the input graph. In the complementer graph all edges are present
 * which are not present in the original graph. Multiple edges in the
 * input graph are ignored.
 * \param graph The input graph.
 * \param al Pointer to a not yet initialized adjacency list.
 * \param mode Constant specifying whether outgoing
 *   (<code>IGRAPH_OUT</code>), incoming (<code>IGRAPH_IN</code>),
 *   or both (<code>IGRAPH_ALL</code>) types of neighbors (in the
 *   complementer graph) to include in the adjacency list. It is
 *   ignored for undirected networks.
 * \param loops Whether to consider loop edges.
 * \return Error code.
 *
 * Time complexity: O(|V|^2+|E|), quadratic in the number of vertices.
 */
int igraph_adjlist_init_complementer(const igraph_t *graph,
                                     igraph_adjlist_t *al,
                                     igraph_neimode_t mode,
                                     igraph_bool_t loops) {
    igraph_integer_t i, j, k, n;
    igraph_bool_t* seen;
    igraph_vector_t vec;

    if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Cannot create complementer adjlist view", IGRAPH_EINVMODE);
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    al->length = igraph_vcount(graph);
    al->adjs = IGRAPH_CALLOC(al->length, igraph_vector_int_t);
    if (al->adjs == 0) {
        IGRAPH_ERROR("Cannot create complementer adjlist view", IGRAPH_ENOMEM);
    }

    IGRAPH_FINALLY(igraph_adjlist_destroy, al);

    n = al->length;
    seen = IGRAPH_CALLOC(n, igraph_bool_t);
    if (seen == 0) {
        IGRAPH_ERROR("Cannot create complementer adjlist view", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, seen);

    IGRAPH_VECTOR_INIT_FINALLY(&vec, 0);

    for (i = 0; i < al->length; i++) {
        IGRAPH_ALLOW_INTERRUPTION();
        igraph_neighbors(graph, &vec, i, mode);
        memset(seen, 0, sizeof(igraph_bool_t) * (unsigned) al->length);
        n = al->length;
        if (!loops) {
            seen[i] = 1;
            n--;
        }
        for (j = 0; j < igraph_vector_size(&vec); j++) {
            if (! seen [ (long int) VECTOR(vec)[j] ] ) {
                n--;
                seen[ (long int) VECTOR(vec)[j] ] = 1;
            }
        }
        IGRAPH_CHECK(igraph_vector_int_init(&al->adjs[i], n));
        for (j = 0, k = 0; k < n; j++) {
            if (!seen[j]) {
                VECTOR(al->adjs[i])[k++] = j;
            }
        }
    }

    IGRAPH_FREE(seen);
    igraph_vector_destroy(&vec);
    IGRAPH_FINALLY_CLEAN(3);
    return 0;
}

/**
 * \function igraph_adjlist_destroy
 * \brief Deallocates an adjacency list.
 *
 * Free all memory allocated for an adjacency list.
 * \param al The adjacency list to destroy.
 *
 * Time complexity: depends on memory management.
 */
void igraph_adjlist_destroy(igraph_adjlist_t *al) {
    long int i;
    for (i = 0; i < al->length; i++) {
        if (&al->adjs[i]) {
            igraph_vector_int_destroy(&al->adjs[i]);
        }
    }
    IGRAPH_FREE(al->adjs);
}

/**
 * \function igraph_adjlist_clear
 * Removes all edges from an adjacency list.
 *
 * \param al The adjacency list.
 * Time complexity: depends on memory management, typically O(n), where n is
 * the total number of elements in the adjacency list.
 */
void igraph_adjlist_clear(igraph_adjlist_t *al) {
    long int i;
    for (i = 0; i < al->length; i++) {
        igraph_vector_int_clear(&al->adjs[i]);
    }
}

/**
 * \function igraph_adjlist_size
 * \brief Returns the number of vertices in an adjacency list.
 *
 * \param al The adjacency list.
 * \return The number of vertices in the adjacency list.
 *
 * Time complexity: O(1).
 */
igraph_integer_t igraph_adjlist_size(const igraph_adjlist_t *al) {
    return al->length;
}

/**
 * \function igraph_adjlist_sort
 * \brief Sorts each vector in an adjacency list.
 *
 * Sorts every vector of the adjacency list.
 * \param al The adjacency list.
 *
 * Time complexity: O(n log n), n is the total number of elements in
 * the adjacency list.
 */
void igraph_adjlist_sort(igraph_adjlist_t *al) {
    long int i;
    for (i = 0; i < al->length; i++) {
        igraph_vector_int_sort(&al->adjs[i]);
    }
}

/**
 * \function igraph_adjlist_simplify
 * \brief Simplifies an adjacency list.
 *
 * Simplifies an adjacency list, i.e. removes loop and multiple edges.
 *
 * \param al The adjacency list.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of edges and
 * vertices.
 */
int igraph_adjlist_simplify(igraph_adjlist_t *al) {
    long int i;
    long int n = al->length;
    igraph_vector_int_t mark;
    igraph_vector_int_init(&mark, n);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &mark);
    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = &al->adjs[i];
        long int j, l = igraph_vector_int_size(v);
        VECTOR(mark)[i] = i + 1;
        for (j = 0; j < l; /* nothing */) {
            long int e = (long int) VECTOR(*v)[j];
            if (VECTOR(mark)[e] != i + 1) {
                VECTOR(mark)[e] = i + 1;
                j++;
            } else {
                VECTOR(*v)[j] = igraph_vector_int_tail(v);
                igraph_vector_int_pop_back(v);
                l--;
            }
        }
    }

    igraph_vector_int_destroy(&mark);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

int igraph_adjlist_remove_duplicate(const igraph_t *graph,
                                    igraph_adjlist_t *al) {
    long int i, j, l, n, p;
    igraph_vector_int_t *v;
    
    IGRAPH_WARNING(
        "igraph_adjlist_remove_duplicate() is deprecated; use the constructor "
        "arguments of igraph_adjlist_init() to specify whether you want loop "
        "edges to appear once or twice in the adjacency list."
    );

    IGRAPH_UNUSED(graph);

    n = al->length;
    for (i = 0; i < n; i++) {
        v = &al->adjs[i];
        l = igraph_vector_int_size(v);
        if (l > 0) {
            p = 1;
            for (j = 1; j < l; j++) {
                long int e = (long int) VECTOR(*v)[j];
                /* Non-loop edges, and one end of loop edges are fine. */
                /* We assume that the vector is sorted and we also keep it sorted */
                if (e != i || VECTOR(*v)[j - 1] != e) {
                    VECTOR(*v)[p++] = e;
                }
            }
            igraph_vector_int_resize(v, p);
        }
    }

    return 0;
}

#ifndef USING_R
int igraph_adjlist_print(const igraph_adjlist_t *al) {
    long int i;
    long int n = al->length;
    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = &al->adjs[i];
        igraph_vector_int_print(v);
    }
    return 0;
}
#endif

int igraph_adjlist_fprint(const igraph_adjlist_t *al, FILE *outfile) {
    long int i;
    long int n = al->length;
    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = &al->adjs[i];
        igraph_vector_int_fprint(v, outfile);
    }
    return 0;
}

#define ADJLIST_CANON_EDGE(from, to, directed) \
    do {                     \
        igraph_integer_t temp;         \
        if((!directed) && from < to) {     \
            temp = to;               \
            to = from;               \
            from = temp;             \
        }                      \
    } while(0);

igraph_bool_t igraph_adjlist_has_edge(igraph_adjlist_t* al, igraph_integer_t from, igraph_integer_t to, igraph_bool_t directed) {
    igraph_vector_int_t* fromvec;
    ADJLIST_CANON_EDGE(from, to, directed);
    fromvec = igraph_adjlist_get(al, from);
    return igraph_vector_int_binsearch2(fromvec, to);

}

int igraph_adjlist_replace_edge(igraph_adjlist_t* al, igraph_integer_t from, igraph_integer_t oldto, igraph_integer_t newto, igraph_bool_t directed) {
    igraph_vector_int_t *oldfromvec, *newfromvec;
    int err1, err2;
    long int oldpos, newpos;
    igraph_integer_t oldfrom = from, newfrom = from;
    ADJLIST_CANON_EDGE(oldfrom, oldto, directed);
    ADJLIST_CANON_EDGE(newfrom, newto, directed);

    oldfromvec = igraph_adjlist_get(al, oldfrom);
    newfromvec = igraph_adjlist_get(al, newfrom);


    err1 = igraph_vector_int_binsearch(oldfromvec, oldto, &oldpos);
    err2 = igraph_vector_int_binsearch(newfromvec, newto, &newpos);

    /* oldfrom -> oldto should exist; newfrom -> newto should not. */
    if ((!err1) || err2) {
        return 1;
    }

    igraph_vector_int_remove(oldfromvec, oldpos);
    if (oldfromvec == newfromvec && oldpos < newpos) {
        --newpos;
    }
    IGRAPH_CHECK(igraph_vector_int_insert(newfromvec, newpos, newto));

    return 0;

}

static int igraph_i_remove_loops_from_incidence_vector_in_place(
    igraph_vector_int_t *v, const igraph_t *graph, igraph_loops_t loops
) {
    long int i, length, eid, write_ptr;
    igraph_vector_int_t *seen_loops = 0;

    /* In this function we make use of the fact that we are dealing with
     * _incidence_ lists, and the only way for an edge ID to appear twice
     * within an incidence list is if the edge is a loop edge; otherwise each
     * element will be unique.
     * 
     * Note that incidence vectors are not sorted by edge ID, so we need to
     * look up the edges in the graph to decide whether they are loops or not.
     * 
     * Also, it may be tempting to introduce a boolean in case of IGRAPH_LOOPS_ONCE,
     * and flip it every time we see a loop to get rid of half of the occurrences,
     * but the problem is that even if the same loop edge ID appears twice in
     * the input list, they are not guaranteed to be next to each other; it
     * may be the case that there are multiple loop edges, each edge appears
     * twice, and we want to keep exactly one of them for each ID. That's why
     * we have a "seen_loops" vector.
     */

    if (loops == IGRAPH_LOOPS_TWICE) {
        /* Loop edges appear twice by default, nothing to do. */
        return IGRAPH_SUCCESS;
    }

    length = igraph_vector_int_size(v);
    if (length == 0) {
        return IGRAPH_SUCCESS;
    }

    if (loops == IGRAPH_LOOPS_ONCE) {
        /* We need a helper vector */
        seen_loops = IGRAPH_CALLOC(1, igraph_vector_int_t);
        IGRAPH_FINALLY(igraph_free, seen_loops);
        IGRAPH_CHECK(igraph_vector_int_init(seen_loops, 0));
        IGRAPH_FINALLY(igraph_vector_int_destroy, seen_loops);
    } else if (loops != IGRAPH_NO_LOOPS) {
        IGRAPH_ERROR("Invalid value for 'loops' argument", IGRAPH_EINVAL);
    }

    for (i = 0, write_ptr = 0; i < length; i++) {
        eid = VECTOR(*v)[i];
        if (IGRAPH_FROM(graph, eid) == IGRAPH_TO(graph, eid)) {
            /* Loop edge */
            if (seen_loops && !igraph_vector_int_contains(seen_loops, eid)) {
                VECTOR(*v)[write_ptr++] = eid;
                IGRAPH_CHECK(igraph_vector_int_push_back(seen_loops, eid));
            }
        } else {
            /* Not a loop edge */
            VECTOR(*v)[write_ptr++] = eid;
        }
    }

    /* Always succeeds since we never grow the vector */
    igraph_vector_int_resize(v, write_ptr);

    /* Destroy the helper vector */
    if (seen_loops) {
        igraph_vector_int_destroy(seen_loops);
        IGRAPH_FREE(seen_loops);
        IGRAPH_FINALLY_CLEAN(2);
    }
    
    return IGRAPH_SUCCESS;
}

int igraph_inclist_remove_duplicate(const igraph_t *graph, igraph_inclist_t *il) {
    long int i, n;

    IGRAPH_WARNING(
        "igraph_inclist_remove_duplicate() is deprecated; use the constructor "
        "arguments of igraph_inclist_init() to specify whether you want loop "
        "edges to appear once or twice in the incidence list."
    );

    IGRAPH_UNUSED(graph);

    n = il->length;
    for (i = 0; i < n; i++) {
        IGRAPH_CHECK(
            igraph_i_remove_loops_from_incidence_vector_in_place(
                &il->incs[i], graph, IGRAPH_LOOPS_ONCE
            )
        );
    }

    return 0;
}

#ifndef USING_R
int igraph_inclist_print(const igraph_inclist_t *al) {
    long int i;
    long int n = al->length;
    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = &al->incs[i];
        igraph_vector_int_print(v);
    }
    return 0;
}
#endif

int igraph_inclist_fprint(const igraph_inclist_t *al, FILE *outfile) {
    long int i;
    long int n = al->length;
    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = &al->incs[i];
        igraph_vector_int_fprint(v, outfile);
    }
    return 0;
}

/**
 * \function igraph_inclist_init
 * \brief Initializes an incidence list.
 *
 * Creates a list of vectors containing the incident edges for all
 * vertices. The incidence list is independent of the graph after
 * creation, subsequent changes of the graph object do not update the
 * incidence list, and changes to the incidence list do not update the
 * graph.
 *
 * </para><para>
 * When \p mode is \c IGRAPH_IN or \c IGRAPH_OUT, each edge ID will appear
 * in the incidence list \em once. When \p mode is \c IGRAPH_ALL, each edge ID
 * will appear in the incidence list \em twice, once for the source vertex
 * and once for the target edge. It also means that the edge IDs of loop edges
 * will appear \em twice for the \em same vertex.
 *
 * \param graph The input graph.
 * \param il Pointer to an uninitialized incidence list.
 * \param mode Constant specifying whether incoming edges
 *   (<code>IGRAPH_IN</code>), outgoing edges (<code>IGRAPH_OUT</code>) or
 *   both (<code>IGRAPH_ALL</code>) to include in the incidence lists
 *   of directed graphs. It is ignored for undirected graphs.
 * \param loops Specifies how to treat loop edges. <code>IGRAPH_NO_LOOPS</code>
 *   removes loop edges from the incidence list. <code>IGRAPH_LOOPS_ONCE</code>
 *   makes each loop edge appear only once in the incidence list of the
 *   corresponding vertex. <code>IGRAPH_LOOPS_TWICE</code> makes loop edges
 *   appear \em twice in the incidence list of the corresponding vertex,
 *   but only if the graph is undirected or <code>mode</code> is set to
 *   <code>IGRAPH_ALL</code>.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 */
int igraph_inclist_init(const igraph_t *graph,
                        igraph_inclist_t *il,
                        igraph_neimode_t mode,
                        igraph_loops_t loops) {
    igraph_integer_t i;
    igraph_vector_t tmp;

    if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Cannot create incidence list view", IGRAPH_EINVMODE);
    }

    igraph_vector_init(&tmp, 0);
    IGRAPH_FINALLY(igraph_vector_destroy, &tmp);

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    il->length = igraph_vcount(graph);
    il->incs = IGRAPH_CALLOC(il->length, igraph_vector_int_t);
    if (il->incs == 0) {
        IGRAPH_ERROR("Cannot create incidence list view", IGRAPH_ENOMEM);
    }

    IGRAPH_FINALLY(igraph_inclist_destroy, il);
    for (i = 0; i < il->length; i++) {
        int j, n;

        IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_CHECK(igraph_incident(graph, &tmp, i, mode));

        n = igraph_vector_size(&tmp);
        IGRAPH_CHECK(igraph_vector_int_init(&il->incs[i], n));

        for (j = 0; j < n; j++) {
            VECTOR(il->incs[i])[j] = VECTOR(tmp)[j];
        }

        if (loops != IGRAPH_LOOPS_TWICE) {
            IGRAPH_CHECK(
                igraph_i_remove_loops_from_incidence_vector_in_place(&il->incs[i], graph, loops)
            );
        }
    }

    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

/**
 * \function igraph_inclist_init_empty
 * \brief Initializes an incidence list corresponding to an empty graph.
 *
 * This function essentially creates a list of empty vectors that may
 * be treated as an incidence list for a graph with a given number of
 * vertices.
 *
 * \param il Pointer to an uninitialized incidence list.
 * \param n  The number of vertices in the incidence list.
 * \return Error code.
 *
 * Time complexity: O(|V|), linear in the number of vertices.
 */

int igraph_inclist_init_empty(igraph_inclist_t *il, igraph_integer_t n) {
    long int i;

    il->length = n;
    il->incs = IGRAPH_CALLOC(il->length, igraph_vector_int_t);
    if (il->incs == 0) {
        IGRAPH_ERROR("Cannot create incidence list view", IGRAPH_ENOMEM);
    }

    IGRAPH_FINALLY(igraph_inclist_destroy, il);
    for (i = 0; i < n; i++) {
        IGRAPH_CHECK(igraph_vector_int_init(&il->incs[i], 0));
    }

    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_inclist_destroy
 * \brief Frees all memory allocated for an incidence list.
 *
 * \param eal The incidence list to destroy.
 *
 * Time complexity: depends on memory management.
 */

void igraph_inclist_destroy(igraph_inclist_t *il) {
    long int i;
    for (i = 0; i < il->length; i++) {
        /* This works if some igraph_vector_int_t's are 0,
           because igraph_vector_destroy can handle this. */
        igraph_vector_int_destroy(&il->incs[i]);
    }
    IGRAPH_FREE(il->incs);
}

/**
 * \function igraph_inclist_clear
 * \brief Removes all edges from an incidence list.
 *
 * \param il The incidence list.
 * 
 * Time complexity: depends on memory management, typically O(n), where n is
 * the total number of elements in the incidence list.
 */
void igraph_inclist_clear(igraph_inclist_t *il) {
    long int i;
    for (i = 0; i < il->length; i++) {
        igraph_vector_int_clear(&il->incs[i]);
    }
}

/**
 * \function igraph_inclist_size
 * \brief Returns the number of vertices in an incidence list.
 *
 * \param il The incidence list.
 * \return The number of vertices in the incidence list.
 *
 * Time complexity: O(1).
 */
igraph_integer_t igraph_inclist_size(const igraph_inclist_t *il) {
    return il->length;
}

static int igraph_i_simplify_sorted_int_adjacency_vector_in_place(
    igraph_vector_int_t *v, igraph_integer_t index, igraph_neimode_t mode,
    igraph_loops_t loops, igraph_multiple_t multiple
) {
    long int i, p = 0;
    long int n = igraph_vector_int_size(v);

    if (
        multiple == IGRAPH_MULTIPLE &&
        (
            loops == IGRAPH_LOOPS_TWICE ||
            (loops == IGRAPH_LOOPS_ONCE && (mode == IGRAPH_IN || mode == IGRAPH_OUT))
        )
    ) {
        /* nothing to simplify */
        return IGRAPH_SUCCESS;        
    }

    if (loops == IGRAPH_NO_LOOPS) {
        if (multiple == IGRAPH_NO_MULTIPLE) {
            /* We need to get rid of loops and multiple edges completely */
            for (i = 0; i < n; i++) {
                if (VECTOR(*v)[i] != index &&
                    (i == n - 1 || VECTOR(*v)[i + 1] != VECTOR(*v)[i])) {
                    VECTOR(*v)[p] = VECTOR(*v)[i];
                    p++;
                }
            }
        } else {
            /* We need to get rid of loops but keep multiple edges */
            for (i = 0; i < n; i++) {
                if (VECTOR(*v)[i] != index) {
                    VECTOR(*v)[p] = VECTOR(*v)[i];
                    p++;
                }
            }
        }
    } else if (loops == IGRAPH_LOOPS_ONCE) {
        if (multiple == IGRAPH_NO_MULTIPLE) {
            /* We need to get rid of multiple edges completely (including
             * multiple loop edges), but keep one edge from each loop edge */
            /* TODO(ntamas): think this through! */
            for (i = 0; i < n; i++) {
                if (i == n - 1 || VECTOR(*v)[i + 1] != VECTOR(*v)[i]) {
                    VECTOR(*v)[p] = VECTOR(*v)[i];
                    p++;
                }
            }
        } else {
            /* We need to keep one edge from each loop edge and we don't need to
             * touch multiple edges. Note that we can get here only if
             * mode == IGRAPH_ALL; if mode was IGRAPH_IN or IGRAPH_OUT, we would
             * have bailed out earlier */
            for (i = 0; i < n; i++) {
                VECTOR(*v)[p] = VECTOR(*v)[i];
                if (VECTOR(*v)[i] == index) {
                    /* this was a loop edge so if the next element is the same, we
                    * need to skip that */
                    if (i < n-1 && VECTOR(*v)[i + 1] == index) {
                        i++;
                    }
                }
                p++;
            }
        }
    } else if (loops == IGRAPH_LOOPS_TWICE && multiple == IGRAPH_NO_MULTIPLE) {
        /* We need to get rid of multiple edges completely (including
         * multiple loop edges), but keep both edge from each loop edge */
        /* TODO(ntamas): think this through! */
        for (i = 0; i < n; i++) {
            if (i == n - 1 || VECTOR(*v)[i + 1] != VECTOR(*v)[i]) {
                VECTOR(*v)[p] = VECTOR(*v)[i];
                p++;
            } else {
                /* Current item is the same as the next one, but if it is a
                 * loop edge, then the first one or two items are okay. We need
                 * to keep one if mode == IGRAPH_IN or mode == IGRAPH_OUT,
                 * otherwise we need to keep two */
                if (VECTOR(*v)[i] == index) {
                    VECTOR(*v)[p] = VECTOR(*v)[i];
                    p++;
                    if (mode == IGRAPH_ALL) {
                        VECTOR(*v)[p] = VECTOR(*v)[i];
                        p++;
                    }
                    /* skip over all the items corresponding to the loop edges */
                    while (i < n && VECTOR(*v)[i] == index) {
                        i++;
                    }
                    i--; /* because the for loop also increases i by 1 */
                }
            }
        }
    } else {
        /* TODO; we don't use this combination yet */
        return IGRAPH_UNIMPLEMENTED;
    }

    /* always succeeds since we are never growing the vector */
    igraph_vector_int_resize(v, p);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_lazy_adjlist_init
 * \brief Initialized a lazy adjacency list.
 *
 * Create a lazy adjacency list for vertices. This function only
 * allocates some memory for storing the vectors of an adjacency list,
 * but the neighbor vertices are not queried, only at the \ref
 * igraph_lazy_adjlist_get() calls.
 * \param graph The input graph.
 * \param al Pointer to an uninitialized adjacency list object.
 * \param mode Constant, it gives whether incoming edges
 *   (<code>IGRAPH_IN</code>), outgoing edges
 *   (<code>IGRPAH_OUT</code>) or both types of edges
 *   (<code>IGRAPH_ALL</code>) are considered. It is ignored for
 *   undirected graphs.
 * \param simplify Constant, it gives whether to simplify the vectors
 *   in the adjacency list (<code>IGRAPH_SIMPLIFY</code>) or not
 *   (<code>IGRAPH_DONT_SIMPLIFY</code>).
 * \return Error code.
 *
 * Time complexity: O(|V|), the number of vertices, possibly, but
 * depends on the underlying memory management too.
 */

int igraph_lazy_adjlist_init(const igraph_t *graph,
                             igraph_lazy_adjlist_t *al,
                             igraph_neimode_t mode,
                             igraph_loops_t loops,
                             igraph_multiple_t multiple) {
    if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Cannor create lazy adjacency list view", IGRAPH_EINVMODE);
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    al->mode = mode;
    al->loops = loops;
    al->multiple = multiple;
    al->graph = graph;

    al->length = igraph_vcount(graph);
    al->adjs = IGRAPH_CALLOC(al->length, igraph_vector_int_t*);

    if (al->adjs == 0) {
        IGRAPH_ERROR("Cannot create lazy adjacency list view", IGRAPH_ENOMEM);
    }

    IGRAPH_FINALLY(igraph_free, al->adjs);

    IGRAPH_CHECK(igraph_vector_init(&al->dummy, 0));
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_lazy_adjlist_destroy
 * \brief Deallocate a lazt adjacency list.
 *
 * Free all allocated memory for a lazy adjacency list.
 * \param al The adjacency list to deallocate.
 *
 * Time complexity: depends on the memory management.
 */

void igraph_lazy_adjlist_destroy(igraph_lazy_adjlist_t *al) {
    igraph_lazy_adjlist_clear(al);
    igraph_vector_destroy(&al->dummy);
    IGRAPH_FREE(al->adjs);
}

/**
 * \function igraph_lazy_adjlist_clear
 * \brief Removes all edges from a lazy adjacency list.
 *
 * \param al The lazy adjacency list.
 * Time complexity: depends on memory management, typically O(n), where n is
 * the total number of elements in the adjacency list.
 */
void igraph_lazy_adjlist_clear(igraph_lazy_adjlist_t *al) {
    long int i, n = al->length;
    for (i = 0; i < n; i++) {
        if (al->adjs[i] != 0) {
            igraph_vector_int_destroy(al->adjs[i]);
            IGRAPH_FREE(al->adjs[i]);
        }
    }
}

/**
 * \function igraph_lazy_adjlist_size
 * \brief Returns the number of vertices in a lazy adjacency list.
 *
 * \param al The lazy adjacency list.
 * \return The number of vertices in the lazy adjacency list.
 *
 * Time complexity: O(1).
 */
igraph_integer_t igraph_lazy_adjlist_size(const igraph_lazy_adjlist_t *al) {
    return al->length;
}

igraph_vector_int_t *igraph_i_lazy_adjlist_get_real(igraph_lazy_adjlist_t *al,
        igraph_integer_t pno) {
    igraph_integer_t no = pno;
    long int i, n;
    int ret;

    if (al->adjs[no] == 0) {
        ret = igraph_neighbors(al->graph, &al->dummy, no, al->mode);
        if (ret != 0) {
            igraph_error("", IGRAPH_FILE_BASENAME, __LINE__, ret);
            return 0;
        }

        al->adjs[no] = IGRAPH_CALLOC(1, igraph_vector_int_t);
        if (al->adjs[no] == 0) {
            igraph_error("Lazy adjlist failed", IGRAPH_FILE_BASENAME, __LINE__,
                         IGRAPH_ENOMEM);
            return 0;
        }

        n = igraph_vector_size(&al->dummy);
        ret = igraph_vector_int_init(al->adjs[no], n);
        if (ret != 0) {
            IGRAPH_FREE(al->adjs[no]);
            igraph_error("", IGRAPH_FILE_BASENAME, __LINE__, ret);
            return 0;
        }

        for (i = 0; i < n; i++) {
            VECTOR(*al->adjs[no])[i] = VECTOR(al->dummy)[i];
        }

        ret = igraph_i_simplify_sorted_int_adjacency_vector_in_place(
            al->adjs[no], no, al->mode, al->loops, al->multiple
        );
        if (ret != 0) {
            igraph_vector_int_destroy(al->adjs[no]);
            IGRAPH_FREE(al->adjs[no]);
            igraph_error("", IGRAPH_FILE_BASENAME, __LINE__, ret);
            return 0;
        }
    }

    return al->adjs[no];
}

/**
 * \function igraph_lazy_inclist_init
 * \brief Initializes a lazy incidence list of edges.
 *
 * Create a lazy incidence list for edges. This function only
 * allocates some memory for storing the vectors of an incidence list,
 * but the incident edges are not queried, only when \ref
 * igraph_lazy_inclist_get() is called.
 *
 * </para><para>
 * When \p mode is \c IGRAPH_IN or \c IGRAPH_OUT, each edge ID will appear
 * in the incidence list \em once. When \p mode is \c IGRAPH_ALL, each edge ID
 * will appear in the incidence list \em twice, once for the source vertex
 * and once for the target edge. It also means that the edge IDs of loop edges
 * will appear \em twice for the \em same vertex.
 *
 * \param graph The input graph.
 * \param al Pointer to an uninitialized incidence list.
 * \param mode Constant, it gives whether incoming edges
 *   (<code>IGRAPH_IN</code>), outgoing edges
 *   (<code>IGRAPH_OUT</code>) or both types of edges
 *   (<code>IGRAPH_ALL</code>) are considered. It is ignored for
 *   undirected graphs.
 * \param loops Specifies how to treat loop edges. <code>IGRAPH_NO_LOOPS</code>
 *   removes loop edges from the incidence list. <code>IGRAPH_LOOPS_ONCE</code>
 *   makes each loop edge appear only once in the incidence list of the
 *   corresponding vertex. <code>IGRAPH_LOOPS_TWICE</code> makes loop edges
 *   appear \em twice in the incidence list of the corresponding vertex,
 *   but only if the graph is undirected or <code>mode</code> is set to
 *   <code>IGRAPH_ALL</code>.
 * \return Error code.
 *
 * Time complexity: O(|V|), the number of vertices, possibly. But it
 * also depends on the underlying memory management.
 */

int igraph_lazy_inclist_init(const igraph_t *graph,
                             igraph_lazy_inclist_t *il,
                             igraph_neimode_t mode,
                             igraph_loops_t loops) {

    if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Cannot create lazy incidence list view", IGRAPH_EINVMODE);
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    il->graph = graph;
    il->loops = loops;
    il->mode = mode;

    il->length = igraph_vcount(graph);
    il->incs = IGRAPH_CALLOC(il->length, igraph_vector_int_t*);
    if (il->incs == 0) {
     
        IGRAPH_ERROR("Cannot create lazy incidence list view", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, il->incs);

    IGRAPH_CHECK(igraph_vector_init(&il->dummy, 0));
    IGRAPH_FINALLY_CLEAN(1);

    return 0;

}

/**
 * \function igraph_lazy_inclist_destroy
 * \brief Deallocates a lazy incidence list.
 *
 * Frees all allocated memory for a lazy incidence list.
 * \param al The incidence list to deallocate.
 *
 * Time complexity: depends on memory management.
 */

void igraph_lazy_inclist_destroy(igraph_lazy_inclist_t *il) {
    igraph_lazy_inclist_clear(il);
    igraph_vector_destroy(&il->dummy);
    IGRAPH_FREE(il->incs);
}

/**
 * \function igraph_lazy_inclist_clear
 * \brief Removes all edges from a lazy incidence list.
 *
 * \param il The lazy incidence list.
 *
 * Time complexity: depends on memory management, typically O(n), where n is
 * the total number of elements in the incidence list.
 */
void igraph_lazy_inclist_clear(igraph_lazy_inclist_t *il) {
    long int i, n = il->length;
    for (i = 0; i < n; i++) {
        if (il->incs[i] != 0) {
            igraph_vector_int_destroy(il->incs[i]);
            IGRAPH_FREE(il->incs[i]);
        }
    }
}

/**
 * \function igraph_lazy_inclist_size
 * \brief Returns the number of vertices in a lazy incidence list.
 *
 * \param il The lazy incidence list.
 * \return The number of vertices in the lazy incidence list.
 *
 * Time complexity: O(1).
 */
igraph_integer_t igraph_lazy_inclist_size(const igraph_lazy_inclist_t *il) {
    return il->length;
}

igraph_vector_int_t *igraph_i_lazy_inclist_get_real(igraph_lazy_inclist_t *il,
        igraph_integer_t pno) {
    igraph_integer_t no = pno;
    int ret;
    long int i, n;

    if (il->incs[no] == 0) {
        ret = igraph_incident(il->graph, &il->dummy, no, il->mode);
        if (ret != 0) {
            igraph_error("", IGRAPH_FILE_BASENAME, __LINE__, ret);
            return 0;
        }

        il->incs[no] = IGRAPH_CALLOC(1, igraph_vector_int_t);
        if (il->incs[no] == 0) {
            igraph_error("Lazy incidence list query failed", IGRAPH_FILE_BASENAME, __LINE__,
                         IGRAPH_ENOMEM);
            return 0;
        }

        n = igraph_vector_size(&il->dummy);
        ret = igraph_vector_int_init(il->incs[no], n);
        if (ret != 0) {
            IGRAPH_FREE(il->incs[no]);
            igraph_error("", IGRAPH_FILE_BASENAME, __LINE__, ret);
            return 0;
        }

        for (i = 0; i < n; i++) {
            VECTOR(*il->incs[no])[i] = VECTOR(il->dummy)[i];
        }

        if (il->loops != IGRAPH_LOOPS_TWICE) {
            ret = igraph_i_remove_loops_from_incidence_vector_in_place(il->incs[no], il->graph, il->loops);
            if (ret != 0) {
                igraph_vector_int_destroy(il->incs[no]);
                IGRAPH_FREE(il->incs[no]);
                return 0;
            }
        }
    }

    return il->incs[no];
}
