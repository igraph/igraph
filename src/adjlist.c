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
#include "igraph_interrupt_internal.h"
#include "config.h"

#include <string.h>   /* memset */
#include <stdio.h>

/**
 * \section about_adjlists
 * <para>Sometimes it is easier to work with a graph which is in
 * adjacency list format: a list of vectors; each vector contains the
 * neighbor vertices or incident edges of a given vertex. Typically,
 * this representation is good if we need to iterate over the neighbors
 * of all vertices many times. E.g. when finding the shortest paths
 * between every pairs of vertices or calculating closeness centrality
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
 * (because initialization is generally linear in terms of number of
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
 * Initialize an adjacency list of vertices from a given graph
 *
 * Create a list of vectors containing the neighbors of all vertices
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
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 */

int igraph_adjlist_init(const igraph_t *graph, igraph_adjlist_t *al,
                        igraph_neimode_t mode) {
    igraph_integer_t i;
    igraph_vector_t tmp;

    if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Cannot create adjlist view", IGRAPH_EINVMODE);
    }

    igraph_vector_init(&tmp, 0);
    IGRAPH_FINALLY(igraph_vector_destroy, &tmp);

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    al->length = igraph_vcount(graph);
    al->adjs = igraph_Calloc(al->length, igraph_vector_int_t);
    if (al->adjs == 0) {
        IGRAPH_ERROR("Cannot create adjlist view", IGRAPH_ENOMEM);
    }

    IGRAPH_FINALLY(igraph_adjlist_destroy, al);
    for (i = 0; i < al->length; i++) {
        int j, n;
        IGRAPH_ALLOW_INTERRUPTION();
        IGRAPH_CHECK(igraph_neighbors(graph, &tmp, i, mode));
        n = igraph_vector_size(&tmp);
        IGRAPH_CHECK(igraph_vector_int_init(&al->adjs[i], n));
        for (j = 0; j < n; j++) {
            VECTOR(al->adjs[i])[j] = VECTOR(tmp)[j];
        }
    }

    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

/**
 * \function igraph_adjlist_init_empty
 * Initialize an empty adjacency list
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
    al->adjs = igraph_Calloc(al->length, igraph_vector_int_t);
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
 * Adjacency lists for the complementer graph
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
    al->adjs = igraph_Calloc(al->length, igraph_vector_int_t);
    if (al->adjs == 0) {
        IGRAPH_ERROR("Cannot create complementer adjlist view", IGRAPH_ENOMEM);
    }

    IGRAPH_FINALLY(igraph_adjlist_destroy, al);

    n = al->length;
    seen = igraph_Calloc(n, igraph_bool_t);
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

    igraph_Free(seen);
    igraph_vector_destroy(&vec);
    IGRAPH_FINALLY_CLEAN(3);
    return 0;
}

/**
 * \function igraph_adjlist_destroy
 * Deallocate memory
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
    igraph_Free(al->adjs);
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
 * Number of vertices in an adjacency list.
 *
 * \param al The adjacency list.
 * \return The number of elements.
 *
 * Time complexity: O(1).
 */

igraph_integer_t igraph_adjlist_size(const igraph_adjlist_t *al) {
    return al->length;
}

/* igraph_vector_int_t *igraph_adjlist_get(igraph_adjlist_t *al, igraph_integer_t no) { */
/*   return &al->adjs[(long int)no]; */
/* } */

/**
 * \function igraph_adjlist_sort
 * Sort each vector in an adjacency list.
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
 * Simplify
 *
 * Simplify an adjacency list, ie. remove loop and multiple edges.
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
    long int i;
    long int n = al->length;
    IGRAPH_UNUSED(graph);
    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = &al->adjs[i];
        long int j, p = 1, l = igraph_vector_int_size(v);
        for (j = 1; j < l; j++) {
            long int e = (long int) VECTOR(*v)[j];
            /* Non-loop edges, and one end of loop edges are fine. */
            /* We use here, that the vector is sorted and we also keep it sorted */
            if (e != i || VECTOR(*v)[j - 1] != e) {
                VECTOR(*v)[p++] = e;
            }
        }
        igraph_vector_int_resize(v, p);
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

int igraph_adjedgelist_remove_duplicate(const igraph_t *graph,
                                        igraph_inclist_t *al) {
    IGRAPH_WARNING("igraph_adjedgelist_remove_duplicate() is deprecated, use "
                   "igraph_inclist_remove_duplicate() instead");
    return igraph_inclist_remove_duplicate(graph, al);
}

#ifndef USING_R
int igraph_adjedgelist_print(const igraph_inclist_t *al, FILE *outfile) {
    IGRAPH_WARNING("igraph_adjedgelist_print() is deprecated, use "
                   "igraph_inclist_print() instead");
    return igraph_inclist_fprint(al, outfile);
}
#endif

/**
 * \function igraph_adjedgelist_init
 * Initialize an incidence list of edges
 *
 * This function was superseded by \ref igraph_inclist_init() in igraph 0.6.
 * Please use \ref igraph_inclist_init() instead of this function.
 *
 * </para><para>
 * Deprecated in version 0.6.
 */
int igraph_adjedgelist_init(const igraph_t *graph,
                            igraph_inclist_t *il,
                            igraph_neimode_t mode) {
    IGRAPH_WARNING("igraph_adjedgelist_init() is deprecated, use "
                   "igraph_inclist_init() instead");
    return igraph_inclist_init(graph, il, mode);
}

/**
 * \function igraph_adjedgelist_destroy
 * Frees all memory allocated for an incidence list.
 *
 * This function was superseded by \ref igraph_inclist_destroy() in igraph 0.6.
 * Please use \ref igraph_inclist_destroy() instead of this function.
 *
 * </para><para>
 * Deprecated in version 0.6.
 */
void igraph_adjedgelist_destroy(igraph_inclist_t *il) {
    IGRAPH_WARNING("igraph_adjedgelist_destroy() is deprecated, use "
                   "igraph_inclist_destroy() instead");
    igraph_inclist_destroy(il);
}

int igraph_inclist_remove_duplicate(const igraph_t *graph,
                                    igraph_inclist_t *al) {
    long int i;
    long int n = al->length;
    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = &al->incs[i];
        long int j, p = 1, l = igraph_vector_int_size(v);
        for (j = 1; j < l; j++) {
            long int e = (long int) VECTOR(*v)[j];
            /* Non-loop edges and one end of loop edges are fine. */
            /* We use here, that the vector is sorted and we also keep it sorted */
            if (IGRAPH_FROM(graph, e) != IGRAPH_TO(graph, e) ||
                VECTOR(*v)[j - 1] != e) {
                VECTOR(*v)[p++] = e;
            }
        }
        igraph_vector_int_resize(v, p);
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
 * Initialize an incidence list of edges
 *
 * Create a list of vectors containing the incident edges for all
 * vertices. The incidence list is independent of the graph after
 * creation, subsequent changes of the graph object do not update the
 * incidence list, and changes to the incidence list do not update the
 * graph.
 * \param graph The input graph.
 * \param il Pointer to an uninitialized incidence list.
 * \param mode Constant specifying whether incoming edges
 *   (<code>IGRAPH_IN</code>), outgoing edges (<code>IGRAPH_OUT</code>) or
 *   both (<code>IGRAPH_ALL</code>) to include in the incidence lists
 *   of directed graphs. It is ignored for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 */

int igraph_inclist_init(const igraph_t *graph,
                        igraph_inclist_t *il,
                        igraph_neimode_t mode) {
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
    il->incs = igraph_Calloc(il->length, igraph_vector_int_t);
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
    }

    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

/**
 * \function igraph_inclist_init_empty
 * \brief Initialize an incidence list corresponding to an empty graph.
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
    il->incs = igraph_Calloc(il->length, igraph_vector_int_t);
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
 * Frees all memory allocated for an incidence list.
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
    igraph_Free(il->incs);
}

/**
 * \function igraph_inclist_clear
 * Removes all edges from an incidence list.
 *
 * \param il The incidence list.
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
 * \function igraph_lazy_adjlist_init
 * Constructor
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
                             igraph_lazy_adlist_simplify_t simplify) {
    if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Cannor create adjlist view", IGRAPH_EINVMODE);
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }
    al->mode = mode;
    al->simplify = simplify;
    al->graph = graph;

    al->length = igraph_vcount(graph);
    al->adjs = igraph_Calloc(al->length, igraph_vector_t*);
    if (al->adjs == 0) {
        IGRAPH_ERROR("Cannot create lazy adjlist view", IGRAPH_ENOMEM);
    }

    return 0;
}

/**
 * \function igraph_lazy_adjlist_destroy
 * Deallocate memory
 *
 * Free all allocated memory for a lazy adjacency list.
 * \param al The adjacency list to deallocate.
 *
 * Time complexity: depends on the memory management.
 */

void igraph_lazy_adjlist_destroy(igraph_lazy_adjlist_t *al) {
    igraph_lazy_adjlist_clear(al);
    igraph_Free(al->adjs);
}

/**
 * \function igraph_lazy_adjlist_clear
 * Removes all edges from a lazy adjacency list.
 *
 * \param al The lazy adjacency list.
 * Time complexity: depends on memory management, typically O(n), where n is
 * the total number of elements in the adjacency list.
 */
void igraph_lazy_adjlist_clear(igraph_lazy_adjlist_t *al) {
    long int i, n = al->length;
    for (i = 0; i < n; i++) {
        if (al->adjs[i] != 0) {
            igraph_vector_destroy(al->adjs[i]);
            igraph_Free(al->adjs[i]);
        }
    }
}

igraph_vector_t *igraph_lazy_adjlist_get_real(igraph_lazy_adjlist_t *al,
        igraph_integer_t pno) {
    igraph_integer_t no = pno;
    int ret;
    if (al->adjs[no] == 0) {
        al->adjs[no] = igraph_Calloc(1, igraph_vector_t);
        if (al->adjs[no] == 0) {
            igraph_error("Lazy adjlist failed", __FILE__, __LINE__,
                         IGRAPH_ENOMEM);
        }
        ret = igraph_vector_init(al->adjs[no], 0);
        if (ret != 0) {
            igraph_error("", __FILE__, __LINE__, ret);
        }
        ret = igraph_neighbors(al->graph, al->adjs[no], no, al->mode);
        if (ret != 0) {
            igraph_error("", __FILE__, __LINE__, ret);
        }

        if (al->simplify == IGRAPH_SIMPLIFY) {
            igraph_vector_t *v = al->adjs[no];
            long int i, p = 0, n = igraph_vector_size(v);
            for (i = 0; i < n; i++) {
                if (VECTOR(*v)[i] != no &&
                    (i == n - 1 || VECTOR(*v)[i + 1] != VECTOR(*v)[i])) {
                    VECTOR(*v)[p] = VECTOR(*v)[i];
                    p++;
                }
            }
            igraph_vector_resize(v, p);
        }
    }

    return al->adjs[no];
}

/**
 * \function igraph_lazy_adjedgelist_init
 * Initializes a lazy incidence list of edges
 *
 * This function was superseded by \ref igraph_lazy_inclist_init() in igraph 0.6.
 * Please use \ref igraph_lazy_inclist_init() instead of this function.
 *
 * </para><para>
 * Deprecated in version 0.6.
 */
int igraph_lazy_adjedgelist_init(const igraph_t *graph,
                                 igraph_lazy_inclist_t *il,
                                 igraph_neimode_t mode) {
    IGRAPH_WARNING("igraph_lazy_adjedgelist_init() is deprecated, use "
                   "igraph_lazy_inclist_init() instead");
    return igraph_lazy_inclist_init(graph, il, mode);
}

/**
 * \function igraph_lazy_adjedgelist_destroy
 * Frees all memory allocated for an incidence list.
 *
 * This function was superseded by \ref igraph_lazy_inclist_destroy() in igraph 0.6.
 * Please use \ref igraph_lazy_inclist_destroy() instead of this function.
 *
 * </para><para>
 * Deprecated in version 0.6.
 */
void igraph_lazy_adjedgelist_destroy(igraph_lazy_inclist_t *il) {
    IGRAPH_WARNING("igraph_lazy_adjedgelist_destroy() is deprecated, use "
                   "igraph_lazy_inclist_destroy() instead");
    igraph_lazy_inclist_destroy(il);
}

igraph_vector_t *igraph_lazy_adjedgelist_get_real(igraph_lazy_adjedgelist_t *il,
        igraph_integer_t pno) {
    IGRAPH_WARNING("igraph_lazy_adjedgelist_get_real() is deprecated, use "
                   "igraph_lazy_inclist_get_real() instead");
    return igraph_lazy_inclist_get_real(il, pno);
}

/**
 * \function igraph_lazy_inclist_init
 * Initializes a lazy incidence list of edges
 *
 * Create a lazy incidence list for edges. This function only
 * allocates some memory for storing the vectors of an incidence list,
 * but the incident edges are not queried, only when \ref
 * igraph_lazy_inclist_get() is called.
 * \param graph The input graph.
 * \param al Pointer to an uninitialized incidence list.
 * \param mode Constant, it gives whether incoming edges
 *   (<code>IGRAPH_IN</code>), outgoing edges
 *   (<code>IGRPAH_OUT</code>) or both types of edges
 *   (<code>IGRAPH_ALL</code>) are considered. It is ignored for
 *   undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|V|), the number of vertices, possibly. But it
 * also depends on the underlying memory management.
 */

int igraph_lazy_inclist_init(const igraph_t *graph,
                             igraph_lazy_inclist_t *al,
                             igraph_neimode_t mode) {

    if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Cannot create lazy incidence list view", IGRAPH_EINVMODE);
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    al->mode = mode;
    al->graph = graph;

    al->length = igraph_vcount(graph);
    al->incs = igraph_Calloc(al->length, igraph_vector_t*);
    if (al->incs == 0) {
        IGRAPH_ERROR("Cannot create lazy incidence list view", IGRAPH_ENOMEM);
    }

    return 0;

}

/**
 * \function igraph_lazy_inclist_destroy
 * Deallocates memory
 *
 * Frees all allocated memory for a lazy incidence list.
 * \param al The incidence list to deallocate.
 *
 * Time complexity: depends on memory management.
 */

void igraph_lazy_inclist_destroy(igraph_lazy_inclist_t *il) {
    igraph_lazy_inclist_clear(il);
    igraph_Free(il->incs);
}

/**
 * \function igraph_lazy_inclist_clear
 * Removes all edges from a lazy incidence list.
 *
 * \param il The lazy incidence list.
 * Time complexity: depends on memory management, typically O(n), where n is
 * the total number of elements in the incidence list.
 */
void igraph_lazy_inclist_clear(igraph_lazy_inclist_t *il) {
    long int i, n = il->length;
    for (i = 0; i < n; i++) {
        if (il->incs[i] != 0) {
            igraph_vector_destroy(il->incs[i]);
            igraph_Free(il->incs[i]);
        }
    }
}

igraph_vector_t *igraph_lazy_inclist_get_real(igraph_lazy_inclist_t *il,
        igraph_integer_t pno) {
    igraph_integer_t no = pno;
    int ret;
    if (il->incs[no] == 0) {
        il->incs[no] = igraph_Calloc(1, igraph_vector_t);
        if (il->incs[no] == 0) {
            igraph_error("Lazy incidence list query failed", __FILE__, __LINE__,
                         IGRAPH_ENOMEM);
        }
        ret = igraph_vector_init(il->incs[no], 0);
        if (ret != 0) {
            igraph_error("", __FILE__, __LINE__, ret);
        }
        ret = igraph_incident(il->graph, il->incs[no], no, il->mode);
        if (ret != 0) {
            igraph_error("", __FILE__, __LINE__, ret);
        }
    }
    return il->incs[no];
}
