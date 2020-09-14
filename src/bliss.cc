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

#include "igraph_types.h"
#include "igraph_topology.h"

#include "igraph_datatype.h"
#include "igraph_interface.h"

#include "igraph_handle_exceptions.h"

using namespace bliss;
using namespace std;


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

    g->set_verbose_level(0);

    for (unsigned int i = 0; i < nof_edges; i++) {
        g->add_edge((unsigned int)IGRAPH_FROM(graph, i), (unsigned int)IGRAPH_TO(graph, i));
    }
    return g;
}


static void bliss_free_graph(AbstractGraph *g) {
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
        default: IGRAPH_ERROR("Invalid splitting heuristic", IGRAPH_EINVAL);
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
        default: IGRAPH_ERROR("Invalid splitting heuristic", IGRAPH_EINVAL);
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
        IGRAPH_ERROR("Invalid vertex color vector length", IGRAPH_EINVAL);
    }
    for (int i = 0; i < n; ++i) {
        g->change_color(i, VECTOR(*colors)[i]);
    }
    return IGRAPH_SUCCESS;
}


inline void bliss_info_to_igraph(igraph_bliss_info_t *info, const Stats &stats) {
    if (info) {
        info->max_level      = stats.get_max_level();
        info->nof_nodes      = stats.get_nof_nodes();
        info->nof_leaf_nodes = stats.get_nof_leaf_nodes();
        info->nof_bad_nodes  = stats.get_nof_bad_nodes();
        info->nof_canupdates = stats.get_nof_canupdates();
        info->nof_generators = stats.get_nof_generators();
        stats.group_size.tostring(&info->group_size);
    }
}


// this is the callback function used with AbstractGraph::find_automorphisms()
// it collects the group generators into a pointer vector
static void collect_generators(void *generators, unsigned int n, const unsigned int *aut) {
    igraph_vector_ptr_t *gen = static_cast<igraph_vector_ptr_t *>(generators);
    igraph_vector_t *newvector = igraph_Calloc(1, igraph_vector_t);
    igraph_vector_init(newvector, n);
    copy(aut, aut + n, newvector->stor_begin); // takes care of unsigned int -> double conversion
    igraph_vector_ptr_push_back(gen, newvector);
}

} // end unnamed namespace

/**
 * \function igraph_canonical_permutation
 * Canonical permutation using BLISS
 *
 * This function computes the canonical permutation which transforms
 * the graph into a canonical form by using the BLISS algorithm.
 *
 * \param graph The input graph. Multiple edges between the same nodes
 *   are not supported and will cause an incorrect result to be returned.
 * \param colors An optional vertex color vector for the graph. Supply a
 *   null pointer is the graph is not colored.
 * \param labeling Pointer to a vector, the result is stored here. The
 *    permutation takes vertex 0 to the first element of the vector,
 *    vertex 1 to the second, etc. The vector will be resized as
 *    needed.
 * \param sh The splitting heuristics to be used in BLISS. See \ref
 *    igraph_bliss_sh_t.
 * \param info If not \c NULL then information on BLISS internals is
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
        const unsigned int *cl = g->canonical_form(stats, NULL, NULL);
        IGRAPH_CHECK(igraph_vector_resize(labeling, N));
        for (unsigned int i = 0; i < N; i++) {
            VECTOR(*labeling)[i] = cl[i];
        }

        bliss_info_to_igraph(info, stats);

        delete g;
        IGRAPH_FINALLY_CLEAN(1);
    );

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_automorphisms
 * Number of automorphisms using BLISS
 *
 * The number of automorphisms of a graph is computed using BLISS. The
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
 * \param sh The splitting heuristics to be used in BLISS. See \ref
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
        g->find_automorphisms(stats, NULL, NULL);

        bliss_info_to_igraph(info, stats);

        delete g;
        IGRAPH_FINALLY_CLEAN(1);
    );
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_automorphism_group
 * Automorphism group generators using BLISS
 *
 * The generators of the automorphism group of a graph are computed
 * using BLISS. The generator set may not be minimal and may depend on
 * the splitting heuristics.
 *
 * \param graph The input graph. Multiple edges between the same nodes
 *   are not supported and will cause an incorrect result to be returned.
 * \param colors An optional vertex color vector for the graph. Supply a
 *   null pointer is the graph is not colored.
 * \param generators Must be an initialized pointer vector. It will
 *    contain pointers to \ref igraph_vector_t objects
 *    representing generators of the automorphism group.
 * \param sh The splitting heuristics to be used in BLISS. See \ref
 *    igraph_bliss_sh_t.
 * \param info If not \c NULL then information on BLISS internals is
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
        g->find_automorphisms(stats, collect_generators, generators);

        bliss_info_to_igraph(info, stats);

        delete g;
        IGRAPH_FINALLY_CLEAN(1);
    );

    return IGRAPH_SUCCESS;
}



