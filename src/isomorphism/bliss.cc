/*
 Copyright (C) 2003-2006 Tommi Junttila

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License version 2
 as published by the Free Software Foundation.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* FSF address fixed in the above notice on 1 Oct 2009 by Tamas Nepusz */

#include "bliss/graph.hh"

#include "igraph_topology.h"
#include "igraph_conversion.h"
#include "igraph_interface.h"
#include "igraph_interrupt.h"
#include "igraph_memory.h"
#include "igraph_vector.h"
#include "igraph_vector_ptr.h"

#include "core/exceptions.h"

using namespace bliss;
using namespace std;

/**
 * \section about_bliss
 *
 * <para>
 * Bliss is a successor of the famous NAUTY algorithm and
 * implementation. While using the same ideas in general, with better
 * heuristics and data structures Bliss outperforms NAUTY on most
 * graphs.
 * </para>
 *
 * <para>
 * Bliss was developed and implemented by Tommi Junttila and Petteri Kaski at
 * Helsinki University of Technology, Finland. For more information,
 * see the Bliss homepage at https://users.aalto.fi/~tjunttil/bliss/ and the following
 * publication:
 * </para>
 *
 * <para>
 * Tommi Junttila and Petteri Kaski: "Engineering an Efficient Canonical Labeling
 * Tool for Large and Sparse Graphs" In ALENEX 2007, pages 135–149, 2007
 * https://doi.org/10.1137/1.9781611972870.13
 * </para>
 *
 * <para>
 * Tommi Junttila and Petteri Kaski: "Conflict Propagation and Component Recursion
 * for Canonical Labeling" in TAPAS 2011, pages 151–162, 2011.
 * https://doi.org/10.1007/978-3-642-19754-3_16
 * </para>
 *
 * <para>
 * Bliss works with both directed graphs and undirected graphs. It supports graphs with
 * self-loops, but not graphs with multi-edges.
 * </para>
 *
 * <para>
 * Bliss version 0.75 is included in igraph.
 * </para>
 */

namespace { // unnamed namespace

inline AbstractGraph *bliss_from_igraph(const igraph_t *graph) {
    unsigned int nof_vertices = (unsigned int)igraph_vcount(graph);
    unsigned int nof_edges = (unsigned int)igraph_ecount(graph);

    AbstractGraph *g;

    if (igraph_is_directed(graph)) {
        g = new Digraph(nof_vertices);
    } else {
        g = new Graph(nof_vertices);
    }

    /* g->set_verbose_level(0); */

    for (unsigned int i = 0; i < nof_edges; i++) {
        g->add_edge((unsigned int)IGRAPH_FROM(graph, i), (unsigned int)IGRAPH_TO(graph, i));
    }

    return g;
}


void bliss_free_graph(AbstractGraph *g) {
    delete g;
}


inline int bliss_set_sh(AbstractGraph *g, igraph_bliss_sh_t sh, bool directed) {
    if (directed) {
        Digraph::SplittingHeuristic gsh = Digraph::shs_fsm;
        switch (sh) {
        case IGRAPH_BLISS_F:    gsh = Digraph::shs_f;   break;
        case IGRAPH_BLISS_FL:   gsh = Digraph::shs_fl;  break;
        case IGRAPH_BLISS_FS:   gsh = Digraph::shs_fs;  break;
        case IGRAPH_BLISS_FM:   gsh = Digraph::shs_fm;  break;
        case IGRAPH_BLISS_FLM:  gsh = Digraph::shs_flm; break;
        case IGRAPH_BLISS_FSM:  gsh = Digraph::shs_fsm; break;
        default: IGRAPH_ERROR("Invalid splitting heuristic.", IGRAPH_EINVAL);
        }
        static_cast<Digraph *>(g)->set_splitting_heuristic(gsh);
    } else {
        Graph::SplittingHeuristic gsh = Graph::shs_fsm;
        switch (sh) {
        case IGRAPH_BLISS_F:    gsh = Graph::shs_f;   break;
        case IGRAPH_BLISS_FL:   gsh = Graph::shs_fl;  break;
        case IGRAPH_BLISS_FS:   gsh = Graph::shs_fs;  break;
        case IGRAPH_BLISS_FM:   gsh = Graph::shs_fm;  break;
        case IGRAPH_BLISS_FLM:  gsh = Graph::shs_flm; break;
        case IGRAPH_BLISS_FSM:  gsh = Graph::shs_fsm; break;
        default: IGRAPH_ERROR("Invalid splitting heuristic.", IGRAPH_EINVAL);
        }
        static_cast<Graph *>(g)->set_splitting_heuristic(gsh);
    }
    return IGRAPH_SUCCESS;
}


inline int bliss_set_colors(AbstractGraph *g, const igraph_vector_int_t *colors) {
    if (colors == NULL) {
        return IGRAPH_SUCCESS;
    }
    const int n = g->get_nof_vertices();
    if (n != igraph_vector_int_size(colors)) {
        IGRAPH_ERROR("Invalid vertex color vector length.", IGRAPH_EINVAL);
    }
    for (int i = 0; i < n; ++i) {
        g->change_color(i, VECTOR(*colors)[i]);
    }
    return IGRAPH_SUCCESS;
}


inline int bliss_info_to_igraph(igraph_bliss_info_t *info, const Stats &stats) {
    if (info) {
        size_t group_size_strlen;

        info->max_level      = stats.get_max_level();
        info->nof_nodes      = stats.get_nof_nodes();
        info->nof_leaf_nodes = stats.get_nof_leaf_nodes();
        info->nof_bad_nodes  = stats.get_nof_bad_nodes();
        info->nof_canupdates = stats.get_nof_canupdates();
        info->nof_generators = stats.get_nof_generators();

        mpz_t group_size;
        mpz_init(group_size);
        stats.get_group_size().get(group_size);
        group_size_strlen = mpz_sizeinbase(group_size, /* base */ 10) + 2;
        info->group_size = IGRAPH_CALLOC(group_size_strlen, char);
        if (! info->group_size) {
            IGRAPH_ERROR("Insufficient memory to retrieve automotphism group size.", IGRAPH_ENOMEM);
        }
        mpz_get_str(info->group_size, /* base */ 10, group_size);
        mpz_clear(group_size);
    }

    return IGRAPH_SUCCESS;
}


// This is the callback function that can tell Bliss to terminate early.
struct AbortChecker {
    bool aborted;

    AbortChecker() : aborted(false) { }
    bool operator()() {
        if (igraph_allow_interruption(NULL) != IGRAPH_SUCCESS) {
            aborted = true;
            return true;
        }
        return false;
    }
};


// This is the callback function used with AbstractGraph::find_automorphisms().
// It collects the automorphism group generators into a pointer vector.
class AutCollector {
    igraph_vector_ptr_t *generators;

public:
    AutCollector(igraph_vector_ptr_t *generators_) : generators(generators_) { }

    void operator ()(unsigned int n, const unsigned int *aut) {
        int err;
        igraph_vector_t *newvector = IGRAPH_CALLOC(1, igraph_vector_t);
        if (! newvector) {
            throw bad_alloc();
        }
        err = igraph_vector_init(newvector, n);
        if (err) {
            throw bad_alloc();
        }
        copy(aut, aut + n, newvector->stor_begin); // takes care of unsigned int -> double conversion
        err = igraph_vector_ptr_push_back(generators, newvector);
        if (err) {
            throw bad_alloc();
        }
    }
};

} // end unnamed namespace


/**
 * \function igraph_canonical_permutation
 * Canonical permutation using Bliss
 *
 * This function computes the canonical permutation which transforms
 * the graph into a canonical form by using the Bliss algorithm.
 *
 * \param graph The input graph. Multiple edges between the same nodes
 *   are not supported and will cause an incorrect result to be returned.
 * \param colors An optional vertex color vector for the graph. Supply a
 *   null pointer is the graph is not colored.
 * \param labeling Pointer to a vector, the result is stored here. The
 *    permutation takes vertex 0 to the first element of the vector,
 *    vertex 1 to the second, etc. The vector will be resized as
 *    needed.
 * \param sh The splitting heuristics to be used in Bliss. See \ref
 *    igraph_bliss_sh_t.
 * \param info If not \c NULL then information on Bliss internals is
 *    stored here. See \ref igraph_bliss_info_t.
 * \return Error code.
 *
 * Time complexity: exponential, in practice it is fast for many graphs.
 */
int igraph_canonical_permutation(const igraph_t *graph, const igraph_vector_int_t *colors,
                                 igraph_vector_t *labeling, igraph_bliss_sh_t sh, igraph_bliss_info_t *info) {
    IGRAPH_HANDLE_EXCEPTIONS(
        AbstractGraph *g = bliss_from_igraph(graph);
        IGRAPH_FINALLY(bliss_free_graph, g);
        const unsigned int N = g->get_nof_vertices();

        IGRAPH_CHECK(bliss_set_sh(g, sh, igraph_is_directed(graph)));
        IGRAPH_CHECK(bliss_set_colors(g, colors));

        Stats stats;
        AbortChecker checker;
        const unsigned int *cl = g->canonical_form(stats, /* report */ nullptr, /* terminate */ checker);
        if (checker.aborted) {
            return IGRAPH_INTERRUPTED;
        }

        IGRAPH_CHECK(igraph_vector_resize(labeling, N));
        for (unsigned int i = 0; i < N; i++) {
            VECTOR(*labeling)[i] = cl[i];
        }

        IGRAPH_CHECK(bliss_info_to_igraph(info, stats));

        delete g;
        IGRAPH_FINALLY_CLEAN(1);
    );

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_automorphisms
 * Number of automorphisms using Bliss
 *
 * The number of automorphisms of a graph is computed using Bliss. The
 * result is returned as part of the \p info structure, in tag \c
 * group_size. It is returned as a string, as it can be very high even
 * for relatively small graphs. If the GNU MP library is used then
 * this number is exact, otherwise a <type>long double</type> is used
 * and it is only approximate. See also \ref igraph_bliss_info_t.
 *
 * \param graph The input graph. Multiple edges between the same nodes
 *   are not supported and will cause an incorrect result to be returned.
 * \param colors An optional vertex color vector for the graph. Supply a
 *   null pointer is the graph is not colored.
 * \param sh The splitting heuristics to be used in Bliss. See \ref
 *    igraph_bliss_sh_t.
 * \param info The result is stored here, in particular in the \c
 *    group_size tag of \p info.
 * \return Error code.
 *
 * Time complexity: exponential, in practice it is fast for many graphs.
 */
int igraph_automorphisms(const igraph_t *graph, const igraph_vector_int_t *colors,
                         igraph_bliss_sh_t sh, igraph_bliss_info_t *info) {
    IGRAPH_HANDLE_EXCEPTIONS(
        AbstractGraph *g = bliss_from_igraph(graph);
        IGRAPH_FINALLY(bliss_free_graph, g);

        IGRAPH_CHECK(bliss_set_sh(g, sh, igraph_is_directed(graph)));
        IGRAPH_CHECK(bliss_set_colors(g, colors));

        Stats stats;
        AbortChecker checker;
        g->find_automorphisms(stats, /* report */ nullptr, /* terminate */ checker);
        if (checker.aborted) {
            return IGRAPH_INTERRUPTED;
        }

        IGRAPH_CHECK(bliss_info_to_igraph(info, stats));

        delete g;
        IGRAPH_FINALLY_CLEAN(1);
    );

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_automorphism_group
 * Automorphism group generators using Bliss
 *
 * The generators of the automorphism group of a graph are computed
 * using Bliss. The generator set may not be minimal and may depend on
 * the splitting heuristics.
 *
 * \param graph The input graph. Multiple edges between the same nodes
 *   are not supported and will cause an incorrect result to be returned.
 * \param colors An optional vertex color vector for the graph. Supply a
 *   null pointer is the graph is not colored.
 * \param generators Must be an initialized pointer vector. It will
 *    contain pointers to \ref igraph_vector_t objects
 *    representing generators of the automorphism group.
 * \param sh The splitting heuristics to be used in Bliss. See \ref
 *    igraph_bliss_sh_t.
 * \param info If not \c NULL then information on Bliss internals is
 *    stored here. See \ref igraph_bliss_info_t.
 * \return Error code.
 *
 * Time complexity: exponential, in practice it is fast for many graphs.
 */
int igraph_automorphism_group(
    const igraph_t *graph, const igraph_vector_int_t *colors, igraph_vector_ptr_t *generators,
    igraph_bliss_sh_t sh, igraph_bliss_info_t *info) {
    IGRAPH_HANDLE_EXCEPTIONS(
        AbstractGraph *g = bliss_from_igraph(graph);
        IGRAPH_FINALLY(bliss_free_graph, g);

        IGRAPH_CHECK(bliss_set_sh(g, sh, igraph_is_directed(graph)));
        IGRAPH_CHECK(bliss_set_colors(g, colors));

        Stats stats;
        igraph_vector_ptr_resize(generators, 0);
        AutCollector collector(generators);
        AbortChecker checker;
        g->find_automorphisms(stats, collector, checker);
        if (checker.aborted) {
            return IGRAPH_INTERRUPTED;
        }
        IGRAPH_CHECK(bliss_info_to_igraph(info, stats));

        delete g;
        IGRAPH_FINALLY_CLEAN(1);
    );

    return IGRAPH_SUCCESS;
}


/* The following license notice applies to the rest of this file */

/*
   IGraph library.
   Copyright (C) 2006-2021  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

/**
 * \function igraph_isomorphic_bliss
 * Graph isomorphism via Bliss
 *
 * This function uses the Bliss graph isomorphism algorithm, a
 * successor of the famous NAUTY algorithm and implementation. Bliss
 * is open source and licensed according to the GNU LGPL. See
 * https://users.aalto.fi/~tjunttil/bliss/ for
 * details. Currently the 0.75 version of Bliss is included in igraph.
 *
 * </para><para>
 *
 * \param graph1 The first input graph. Multiple edges between the same nodes
 *   are not supported and will cause an incorrect result to be returned.
 * \param graph2 The second input graph. Multiple edges between the same nodes
 *   are not supported and will cause an incorrect result to be returned.
 * \param colors1 An optional vertex color vector for the first graph. Supply a
 *   null pointer if your graph is not colored.
 * \param colors2 An optional vertex color vector for the second graph. Supply a
 *   null pointer if your graph is not colored.
 * \param iso Pointer to a boolean, the result is stored here.
 * \param map12 A vector or \c NULL pointer. If not \c NULL then an
 *   isomorphic mapping from \p graph1 to \p graph2 is stored here.
 *   If the input graphs are not isomorphic then this vector is
 *   cleared, i.e. it will have length zero.
 * \param map21 Similar to \p map12, but for the mapping from \p
 *   graph2 to \p graph1.
 * \param sh Splitting heuristics to be used for the graphs. See
 *   \ref igraph_bliss_sh_t.
 * \param info1 If not \c NULL, information about the canonization of
 *    the first input graph is stored here. See \ref igraph_bliss_info_t
 *    for details. Note that if the two graphs have different number
 *    of vertices or edges, then this is not filled.
 * \param info2 Same as \p info1, but for the second graph.
 * \return Error code.
 *
 * Time complexity: exponential, but in practice it is quite fast.
 */
int igraph_isomorphic_bliss(const igraph_t *graph1, const igraph_t *graph2,
                            const igraph_vector_int_t *colors1, const igraph_vector_int_t *colors2,
                            igraph_bool_t *iso, igraph_vector_t *map12,
                            igraph_vector_t *map21, igraph_bliss_sh_t sh,
                            igraph_bliss_info_t *info1, igraph_bliss_info_t *info2) {

    long int no_of_nodes = igraph_vcount(graph1);
    long int no_of_edges = igraph_ecount(graph1);
    igraph_vector_t perm1, perm2;
    igraph_vector_t vmap12, *mymap12 = &vmap12;
    igraph_vector_t from, to, index;
    igraph_vector_t from2, to2, index2;
    igraph_bool_t directed;
    long int i, j;

    *iso = 0;
    if (info1) {
        info1->nof_nodes = info1->nof_leaf_nodes = info1->nof_bad_nodes =
                               info1->nof_canupdates = info1->max_level = info1->nof_generators = 0;
        info1->group_size = 0;
    }
    if (info2) {
        info2->nof_nodes = info2->nof_leaf_nodes = info2->nof_bad_nodes =
                               info2->nof_canupdates = info2->max_level = info2->nof_generators = 0;
        info2->group_size = 0;
    }

    directed = igraph_is_directed(graph1);
    if (igraph_is_directed(graph2) != directed) {
        IGRAPH_ERROR("Cannot compare directed and undirected graphs.",
                     IGRAPH_EINVAL);
    }
    if ((colors1 == NULL || colors2 == NULL) && colors1 != colors2) {
        IGRAPH_WARNING("Only one of the graphs is vertex colored, colors will be ignored.");
        colors1 = NULL; colors2 = NULL;
    }

    if (no_of_nodes != igraph_vcount(graph2) ||
        no_of_edges != igraph_ecount(graph2)) {
        if (map12) {
            igraph_vector_clear(map12);
        }
        if (map21) {
            igraph_vector_clear(map21);
        }
        return 0;
    }

    if (map12) {
        mymap12 = map12;
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(mymap12, 0);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&perm1, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&perm2, no_of_nodes);

    IGRAPH_CHECK(igraph_canonical_permutation(graph1, colors1, &perm1, sh, info1));
    IGRAPH_CHECK(igraph_canonical_permutation(graph2, colors2, &perm2, sh, info2));

    IGRAPH_CHECK(igraph_vector_resize(mymap12, no_of_nodes));

    /* The inverse of perm2 is produced in mymap12 */
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(*mymap12)[ (long int)VECTOR(perm2)[i] ] = i;
    }
    /* Now we produce perm2^{-1} o perm1 in perm2 */
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(perm2)[i] = VECTOR(*mymap12)[ (long int) VECTOR(perm1)[i] ];
    }
    /* Copy it to mymap12 */
    igraph_vector_update(mymap12, &perm2);

    igraph_vector_destroy(&perm1);
    igraph_vector_destroy(&perm2);
    IGRAPH_FINALLY_CLEAN(2);

    /* Check isomorphism, we apply the permutation in mymap12 to graph1
       and should get graph2 */

    IGRAPH_VECTOR_INIT_FINALLY(&from, no_of_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&to, no_of_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&index, no_of_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&from2, no_of_edges * 2);
    IGRAPH_VECTOR_INIT_FINALLY(&to2, no_of_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&index2, no_of_edges);

    for (i = 0; i < no_of_edges; i++) {
        VECTOR(from)[i] = VECTOR(*mymap12)[ (long int) IGRAPH_FROM(graph1, i) ];
        VECTOR(to)[i]   = VECTOR(*mymap12)[ (long int) IGRAPH_TO  (graph1, i) ];
        if (! directed && VECTOR(from)[i] < VECTOR(to)[i]) {
            igraph_real_t tmp = VECTOR(from)[i];
            VECTOR(from)[i] = VECTOR(to)[i];
            VECTOR(to)[i] = tmp;
        }
    }
    igraph_vector_order(&from, &to, &index, no_of_nodes);

    igraph_get_edgelist(graph2, &from2, /*bycol=*/ 1);
    for (i = 0, j = no_of_edges; i < no_of_edges; i++, j++) {
        VECTOR(to2)[i] = VECTOR(from2)[j];
        if (! directed && VECTOR(from2)[i] < VECTOR(to2)[i]) {
            igraph_real_t tmp = VECTOR(from2)[i];
            VECTOR(from2)[i] = VECTOR(to2)[i];
            VECTOR(to2)[i] = tmp;
        }
    }
    igraph_vector_resize(&from2, no_of_edges);
    igraph_vector_order(&from2, &to2, &index2, no_of_nodes);

    *iso = 1;
    for (i = 0; i < no_of_edges; i++) {
        long int i1 = (long int) VECTOR(index)[i];
        long int i2 = (long int) VECTOR(index2)[i];
        if (VECTOR(from)[i1] != VECTOR(from2)[i2] ||
            VECTOR(to)[i1] != VECTOR(to2)[i2]) {
            *iso = 0;
            break;
        }
    }

    /* If the graphs are coloured, we also need to check that applying the
       permutation mymap12 to colors1 gives colors2. */

    if (*iso && colors1 != NULL) {
        for (i = 0; i < no_of_nodes; i++) {
            if (VECTOR(*colors1)[i] != VECTOR(*colors2)[(long int) VECTOR(*mymap12)[i] ]) {
                *iso = 0;
                break;
            }
        }
    }

    igraph_vector_destroy(&index2);
    igraph_vector_destroy(&to2);
    igraph_vector_destroy(&from2);
    igraph_vector_destroy(&index);
    igraph_vector_destroy(&to);
    igraph_vector_destroy(&from);
    IGRAPH_FINALLY_CLEAN(6);

    if (*iso) {
        /* The inverse of mymap12 */
        if (map21) {
            IGRAPH_CHECK(igraph_vector_resize(map21, no_of_nodes));
            for (i = 0; i < no_of_nodes; i++) {
                VECTOR(*map21)[ (long int) VECTOR(*mymap12)[i] ] = i;
            }
        }
    } else {
        if (map12) {
            igraph_vector_clear(map12);
        }
        if (map21) {
            igraph_vector_clear(map21);
        }
    }

    if (!map12) {
        igraph_vector_destroy(mymap12);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}
