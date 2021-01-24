/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_cliques.h"
#include "igraph_constants.h"
#include "igraph_interface.h"
#include "igraph_community.h"
#include "igraph_adjlist.h"
#include "core/interruption.h"
#include "igraph_memory.h"
#include "igraph_progress.h"
#include "core/math.h"

#define CONCAT2x(a,b) a ## b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define FUNCTION(name,sfx) CONCAT2(name,sfx)

static igraph_long_t igraph_i_maximal_cliques_reorder_adjlists(
        const igraph_vector_long_t *PX,
        igraph_long_t PS, igraph_long_t PE, igraph_long_t XS, igraph_long_t XE,
        const igraph_vector_long_t *pos,
        igraph_adjlist_t *adjlist);

static igraph_long_t igraph_i_maximal_cliques_select_pivot(
        const igraph_vector_long_t *PX,
        igraph_long_t PS, igraph_long_t PE, igraph_long_t XS, igraph_long_t XE,
        const igraph_vector_long_t *pos,
        const igraph_adjlist_t *adjlist,
        igraph_long_t *pivot,
        igraph_vector_long_t *nextv,
        igraph_long_t oldPS, igraph_long_t oldXE);

static igraph_long_t igraph_i_maximal_cliques_down(
        igraph_vector_long_t *PX,
        igraph_long_t PS, igraph_long_t PE, igraph_long_t XS, igraph_long_t XE,
        igraph_vector_long_t *pos,
        igraph_adjlist_t *adjlist, igraph_long_t mynextv,
        igraph_vector_long_t *R,
        igraph_long_t *newPS, igraph_long_t *newXE);

static igraph_long_t igraph_i_maximal_cliques_PX(
        igraph_vector_long_t *PX, igraph_long_t PS, igraph_long_t *PE,
        igraph_long_t *XS, igraph_long_t XE, igraph_vector_long_t *pos,
        igraph_adjlist_t *adjlist, igraph_long_t v,
        igraph_vector_long_t *H);

static igraph_long_t igraph_i_maximal_cliques_up(
        igraph_vector_long_t *PX, igraph_long_t PS, igraph_long_t PE,
        igraph_long_t XS, igraph_long_t XE, igraph_vector_long_t *pos,
        igraph_adjlist_t *adjlist,
        igraph_vector_long_t *R,
        igraph_vector_long_t *H);

#define PRINT_PX do {                              \
        igraph_long_t j;                                 \
        printf("PX=");                             \
        for (j=0; j<PS; j++) {                         \
            printf("%i ", VECTOR(*PX)[j]);                       \
        }                                      \
        printf("( ");                              \
        for (; j<=PE; j++) {                           \
            printf("%i ", VECTOR(*PX)[j]);                       \
        }                                      \
        printf("| ");                              \
        for (; j<=XE; j++) {                           \
            printf("%i ", VECTOR(*PX)[j]);                       \
        }                                      \
        printf(") ");                              \
        for (; j<igraph_vector_long_size(PX); j++) {                \
            printf("%i ", VECTOR(*PX)[j]);                       \
        }                                      \
        printf("\n");                              \
    } while (0);

#define PRINT_PX1 do {                             \
        igraph_long_t j;                                 \
        printf("PX=");                             \
        for (j=0; j<PS; j++) {                         \
            printf("%i ", VECTOR(*PX)[j]);                       \
        }                                      \
        printf("( ");                              \
        for (; j<=*PE; j++) {                          \
            printf("%i ", VECTOR(*PX)[j]);                       \
        }                                      \
        printf("| ");                              \
        for (; j<=XE; j++) {                           \
            printf("%i ", VECTOR(*PX)[j]);                       \
        }                                      \
        printf(") ");                              \
        for (; j<igraph_vector_long_size(PX); j++) {                \
            printf("%i ", VECTOR(*PX)[j]);                       \
        }                                      \
        printf("\n");                              \
    } while (0)

static igraph_long_t igraph_i_maximal_cliques_reorder_adjlists(
        const igraph_vector_long_t *PX,
        igraph_long_t PS, igraph_long_t PE, igraph_long_t XS, igraph_long_t XE,
        const igraph_vector_long_t *pos,
        igraph_adjlist_t *adjlist) {
    igraph_long_t j;
    igraph_long_t sPS = PS + 1, sPE = PE + 1;

    IGRAPH_UNUSED(XS);

    for (j = PS; j <= XE; j++) {
        igraph_long_t av = VECTOR(*PX)[j];
        igraph_vector_long_t *avneis = igraph_adjlist_get(adjlist, av);
        igraph_long_t *avp = VECTOR(*avneis);
        igraph_long_t avlen = igraph_vector_long_size(avneis);
        igraph_long_t *ave = avp + avlen;
        igraph_long_t *avnei = avp, *pp = avp;

        for (; avnei < ave; avnei++) {
            igraph_long_t avneipos = VECTOR(*pos)[(igraph_long_t)(*avnei)];
            if (avneipos >= sPS && avneipos <= sPE) {
                if (pp != avnei) {
                    igraph_long_t tmp = *avnei;
                    *avnei = *pp;
                    *pp = tmp;
                }
                pp++;
            }
        }
    }
    return 0;
}

static igraph_long_t igraph_i_maximal_cliques_select_pivot(
        const igraph_vector_long_t *PX,
        igraph_long_t PS, igraph_long_t PE, igraph_long_t XS, igraph_long_t XE,
        const igraph_vector_long_t *pos,
        const igraph_adjlist_t *adjlist,
        igraph_long_t *pivot,
        igraph_vector_long_t *nextv,
        igraph_long_t oldPS, igraph_long_t oldXE) {
    igraph_vector_long_t *pivotvectneis;
    igraph_long_t i, pivotvectlen, j, usize = -1;
    igraph_long_t soldPS = oldPS + 1, soldXE = oldXE + 1, sPS = PS + 1, sPE = PE + 1;

    IGRAPH_UNUSED(XS);

    /* Choose a pivotvect, and bring up P vertices at the same time */
    for (i = PS; i <= XE; i++) {
        igraph_long_t av = VECTOR(*PX)[i];
        igraph_vector_long_t *avneis = igraph_adjlist_get(adjlist, av);
        igraph_long_t *avp = VECTOR(*avneis);
        igraph_long_t avlen = igraph_vector_long_size(avneis);
        igraph_long_t *ave = avp + avlen;
        igraph_long_t *avnei = avp, *pp = avp;

        for (; avnei < ave; avnei++) {
            igraph_long_t avneipos = VECTOR(*pos)[(igraph_long_t)(*avnei)];
            if (avneipos < soldPS || avneipos > soldXE) {
                break;
            }
            if (avneipos >= sPS && avneipos <= sPE) {
                if (pp != avnei) {
                    igraph_long_t tmp = *avnei;
                    *avnei = *pp;
                    *pp = tmp;
                }
                pp++;
            }
        }
        if ((j = pp - avp) > usize) {
            *pivot = av;
            usize = j;
        }
    }

    igraph_vector_long_push_back(nextv, -1);
    pivotvectneis = igraph_adjlist_get(adjlist, *pivot);
    pivotvectlen = igraph_vector_long_size(pivotvectneis);

    for (j = PS; j <= PE; j++) {
        igraph_long_t vcand = VECTOR(*PX)[j];
        igraph_bool_t nei = 0;
        igraph_long_t k = 0;
        for (k = 0; k < pivotvectlen; k++) {
            igraph_long_t unv = VECTOR(*pivotvectneis)[k];
            igraph_long_t unvpos = VECTOR(*pos)[unv];
            if (unvpos < sPS || unvpos > sPE) {
                break;
            }
            if (unv == vcand) {
                nei = 1;
                break;
            }
        }
        if (!nei) {
            igraph_vector_long_push_back(nextv, vcand);
        }
    }

    return 0;
}

#define SWAP(p1,p2) do {            \
        igraph_long_t v1=VECTOR(*PX)[p1];         \
        igraph_long_t v2=VECTOR(*PX)[p2];         \
        VECTOR(*PX)[p1] = v2;           \
        VECTOR(*PX)[p2] = v1;           \
        VECTOR(*pos)[v1] = (p2)+1;          \
        VECTOR(*pos)[v2] = (p1)+1;          \
    } while (0)

static igraph_long_t igraph_i_maximal_cliques_down(igraph_vector_long_t *PX,
                                         igraph_long_t PS, igraph_long_t PE, igraph_long_t XS, igraph_long_t XE,
                                         igraph_vector_long_t *pos,
                                         igraph_adjlist_t *adjlist, igraph_long_t mynextv,
                                         igraph_vector_long_t *R,
                                         igraph_long_t *newPS, igraph_long_t *newXE) {

    igraph_vector_long_t *vneis = igraph_adjlist_get(adjlist, mynextv);
    igraph_long_t j, vneislen = igraph_vector_long_size(vneis);
    igraph_long_t sPS = PS + 1, sPE = PE + 1, sXS = XS + 1, sXE = XE + 1;

    *newPS = PE + 1; *newXE = XS - 1;
    for (j = 0; j < vneislen; j++) {
        igraph_long_t vnei = VECTOR(*vneis)[j];
        igraph_long_t vneipos = VECTOR(*pos)[vnei];
        if (vneipos >= sPS && vneipos <= sPE) {
            (*newPS)--;
            SWAP(vneipos - 1, *newPS);
        } else if (vneipos >= sXS && vneipos <= sXE) {
            (*newXE)++;
            SWAP(vneipos - 1, *newXE);
        }
    }

    igraph_vector_long_push_back(R, mynextv);

    return 0;
}

#undef SWAP

static igraph_long_t igraph_i_maximal_cliques_PX(igraph_vector_long_t *PX, igraph_long_t PS, igraph_long_t *PE,
                                       igraph_long_t *XS, igraph_long_t XE, igraph_vector_long_t *pos,
                                       igraph_adjlist_t *adjlist, igraph_long_t v,
                                       igraph_vector_long_t *H) {

    igraph_long_t vpos = VECTOR(*pos)[v] - 1;
    igraph_long_t tmp = VECTOR(*PX)[*PE];

    IGRAPH_UNUSED(PS);
    IGRAPH_UNUSED(XE);
    IGRAPH_UNUSED(adjlist);

    VECTOR(*PX)[vpos] = tmp;
    VECTOR(*PX)[*PE] = v;
    VECTOR(*pos)[v] = (*PE) + 1;
    VECTOR(*pos)[tmp] = vpos + 1;
    (*PE)--; (*XS)--;
    igraph_vector_long_push_back(H, v);

    return 0;
}

static igraph_long_t igraph_i_maximal_cliques_up(igraph_vector_long_t *PX, igraph_long_t PS, igraph_long_t PE,
                                       igraph_long_t XS, igraph_long_t XE, igraph_vector_long_t *pos,
                                       igraph_adjlist_t *adjlist,
                                       igraph_vector_long_t *R,
                                       igraph_vector_long_t *H) {
    igraph_long_t vv;

    IGRAPH_UNUSED(PS);
    IGRAPH_UNUSED(XE);
    IGRAPH_UNUSED(adjlist);

    igraph_vector_long_pop_back(R);

    while ((vv = igraph_vector_long_pop_back(H)) != -1) {
        igraph_long_t vvpos = VECTOR(*pos)[vv];
        igraph_long_t tmp = VECTOR(*PX)[XS];
        VECTOR(*PX)[XS] = vv;
        VECTOR(*PX)[vvpos - 1] = tmp;
        VECTOR(*pos)[vv] = XS + 1;
        VECTOR(*pos)[tmp] = vvpos;
        PE++; XS++;
    }

    return 0;
}

/**
 * \function igraph_maximal_cliques
 * \brief Finds all maximal cliques in a graph.
 *
 * </para><para>
 * A maximal clique is a clique which can't be extended any more by
 * adding a new vertex to it.
 *
 * </para><para>
 * If you are only interested in the size of the largest clique in the
 * graph, use \ref igraph_clique_number() instead.
 *
 * </para><para>
 * The current implementation uses a modified Bron-Kerbosch
 * algorithm to find the maximal cliques, see: David Eppstein,
 * Maarten Löffler, Darren Strash: Listing All Maximal Cliques in
 * Sparse Graphs in Near-Optimal Time. Algorithms and Computation,
 * Lecture Notes in Computer Science Volume 6506, 2010, pp 403-414.
 *
 * </para><para>The implementation of this function changed between
 * igraph 0.5 and 0.6 and also between 0.6 and 0.7, so the order of
 * the cliques and the order of vertices within the cliques will
 * almost surely be different between these three versions.
 *
 * \param graph The input graph.
 * \param res Pointer to a pointer vector, the result will be stored
 *   here, i.e. \p res will contain pointers to \ref igraph_vector_t
 *   objects which contain the indices of vertices involved in a clique.
 *   The pointer vector will be resized if needed but note that the
 *   objects in the pointer vector will not be freed. Note that vertices
 *   of a clique may be returned in arbitrary order.
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_independent_vertex_sets(), \ref
 * igraph_clique_number()
 *
 * Time complexity: O(d(n-d)3^(d/3)) worst case, d is the degeneracy
 * of the graph, this is typically small for sparse graphs.
 *
 * \example examples/simple/igraph_maximal_cliques.c
 */

igraph_long_t igraph_maximal_cliques(const igraph_t *graph,
                           igraph_vector_ptr_t *res,
                           igraph_long_t min_size,
                           igraph_long_t max_size);

#define IGRAPH_MC_ORIG
#include "maximal_cliques_template.h"
#undef IGRAPH_MC_ORIG

/**
 * \function igraph_maximal_cliques_count
 * Count the number of maximal cliques in a graph
 *
 * </para><para>
 * The current implementation uses a modified Bron-Kerbosch
 * algorithm to find the maximal cliques, see: David Eppstein,
 * Maarten Löffler, Darren Strash: Listing All Maximal Cliques in
 * Sparse Graphs in Near-Optimal Time. Algorithms and Computation,
 * Lecture Notes in Computer Science Volume 6506, 2010, pp 403-414.
 *
 * \param graph The input graph.
 * \param res Pointer to an \c igraph_long_t; the number of maximal
 *   cliques will be stored here.
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_cliques().
 *
 * Time complexity: O(d(n-d)3^(d/3)) worst case, d is the degeneracy
 * of the graph, this is typically small for sparse graphs.
 *
 * \example examples/simple/igraph_maximal_cliques.c
 */

igraph_long_t igraph_maximal_cliques_count(const igraph_t *graph,
                                 igraph_long_t *res,
                                 igraph_long_t min_size,
                                 igraph_long_t max_size);

#define IGRAPH_MC_COUNT
#include "maximal_cliques_template.h"
#undef IGRAPH_MC_COUNT

/**
 * \function igraph_maximal_cliques_file
 * Find maximal cliques and write them to a file
 *
 * This function enumerates all maximal cliques and writes them to file.
 *
 * </para><para>
 *
 * Edge directions are ignored.
 *
 * </para><para>
 *
 * \param graph The input graph.
 * \param outfile Pointer to the output file, it should be writable.
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_cliques().
 *
 * Time complexity: O(d(n-d)3^(d/3)) worst case, d is the degeneracy
 * of the graph, this is typically small for sparse graphs.*
 *
 */

igraph_long_t igraph_maximal_cliques_file(const igraph_t *graph,
                                FILE *outfile,
                                igraph_long_t min_size,
                                igraph_long_t max_size);

#define IGRAPH_MC_FILE
#include "maximal_cliques_template.h"
#undef IGRAPH_MC_FILE

/**
 * \function igraph_maximal_cliques_subset
 * Maximal cliques for a subset of initial vertices
 *
 * This function enumerates all maximal cliques for a subset of initial
 * vertices and writes them to file.
 *
 * </para><para>
 *
 * Edge directions are ignored.
 *
 * </para><para>
 *
 * \param graph The input graph.
 * \param subset Pointer to an \c  igraph_vector_long_t containing the
 *   subset of initial vertices
 * \param res Pointer to an \c igraph_ptr_t; the cliques will be
 *   stored here
 * \param no Pointer to an \c igraph_long_t; the number of maximal
 *   cliques will be stored here.
 * \param outfile Pointer to the output file, it should be writable.
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_cliques().
 *
 * Time complexity: O(d(n-d)3^(d/3)) worst case, d is the degeneracy
 * of the graph, this is typically small for sparse graphs.
 *
 */

igraph_long_t igraph_maximal_cliques_subset(const igraph_t *graph,
                                  igraph_vector_long_t *subset,
                                  igraph_vector_ptr_t *res,
                                  igraph_long_t *no,
                                  FILE *outfile,
                                  igraph_long_t min_size,
                                  igraph_long_t max_size);

#define IGRAPH_MC_FULL
#include "maximal_cliques_template.h"
#undef IGRAPH_MC_FULL


/**
 * \function igraph_maximal_cliques_callback
 * \brief Finds maximal cliques in a graph and calls a function for each one.
 *
 * This function enumerates all maximal cliques within the given size range
 * and calls \p cliquehandler_fn for each of them. The cliques are passed to the
 * callback function as a pointer to an \ref igraph_vector_t.  Destroying and
 * freeing this vector is left up to the user.  Use \ref igraph_vector_destroy()
 * to destroy it first, then free it using \ref igraph_free().
 *
 * </para><para>
 *
 * Edge directions are ignored.
 *
 * </para><para>
 *
 * \param graph The input graph.
 * \param cliquehandler_fn Callback function to be called for each clique.
 * See also \ref igraph_clique_handler_t.
 * \param arg Extra argument to supply to \p cliquehandler_fn.
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_cliques().
 *
 * Time complexity: O(d(n-d)3^(d/3)) worst case, d is the degeneracy
 * of the graph, this is typically small for sparse graphs.
 *
 */

igraph_long_t igraph_maximal_cliques_callback(const igraph_t *graph,
                                    igraph_clique_handler_t *cliquehandler_fn, void *arg,
                                    igraph_long_t min_size, igraph_long_t max_size);

#define IGRAPH_MC_CALLBACK
#include "maximal_cliques_template.h"
#undef IGRAPH_MC_CALLBACK


/**
 * \function igraph_maximal_cliques_hist
 * \brief Counts the number of maximal cliques of each size in a graph.
 *
 * This function counts how many maximal cliques of each size are present in
 * the graph. Size-1 maximal cliques are simply isolated vertices.
 *
 * </para><para>
 *
 * Edge directions are ignored.
 *
 * </para><para>
 *
 * \param graph The input graph.
 * \param hist Pointer to an initialized vector. The result will be stored
 * here. The first element will store the number of size-1 maximal cliques,
 * the second element the number of size-2 maximal cliques, etc.
 * For cliques smaller than \p min_size, zero counts will be returned.
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_cliques().
 *
 * Time complexity: O(d(n-d)3^(d/3)) worst case, d is the degeneracy
 * of the graph, this is typically small for sparse graphs.
 *
 */

igraph_long_t igraph_maximal_cliques_hist(const igraph_t *graph,
                                igraph_vector_t *hist,
                                igraph_long_t min_size,
                                igraph_long_t max_size);

#define IGRAPH_MC_HIST
#include "maximal_cliques_template.h"
#undef IGRAPH_MC_HIST
