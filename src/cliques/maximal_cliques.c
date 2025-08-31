/*
   igraph library.
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

#include "igraph_adjlist.h"
#include "igraph_constants.h"
#include "igraph_community.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_progress.h"

#include "core/interruption.h"

#define CONCAT2x(a,b) a ## b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define FUNCTION(name,sfx) CONCAT2(name,sfx)

static igraph_error_t igraph_i_maximal_cliques_reorder_adjlists(
        const igraph_vector_int_t *PX,
        igraph_int_t PS, igraph_int_t PE, igraph_int_t XS, igraph_int_t XE,
        const igraph_vector_int_t *pos,
        igraph_adjlist_t *adjlist);

static igraph_error_t igraph_i_maximal_cliques_select_pivot(
        const igraph_vector_int_t *PX,
        igraph_int_t PS, igraph_int_t PE, igraph_int_t XS, igraph_int_t XE,
        const igraph_vector_int_t *pos,
        const igraph_adjlist_t *adjlist,
        igraph_int_t *pivot,
        igraph_vector_int_t *nextv,
        igraph_int_t oldPS, igraph_int_t oldXE);

static igraph_error_t igraph_i_maximal_cliques_down(
        igraph_vector_int_t *PX,
        igraph_int_t PS, igraph_int_t PE, igraph_int_t XS, igraph_int_t XE,
        igraph_vector_int_t *pos,
        igraph_adjlist_t *adjlist, igraph_int_t mynextv,
        igraph_vector_int_t *R,
        igraph_int_t *newPS, igraph_int_t *newXE);

static igraph_error_t igraph_i_maximal_cliques_PX(
        igraph_vector_int_t *PX, igraph_int_t PS, igraph_int_t *PE,
        igraph_int_t *XS, igraph_int_t XE, igraph_vector_int_t *pos,
        igraph_adjlist_t *adjlist, igraph_int_t v,
        igraph_vector_int_t *H);

static igraph_error_t igraph_i_maximal_cliques_up(
        igraph_vector_int_t *PX, igraph_int_t PS, igraph_int_t PE,
        igraph_int_t XS, igraph_int_t XE, igraph_vector_int_t *pos,
        igraph_adjlist_t *adjlist,
        igraph_vector_int_t *R,
        igraph_vector_int_t *H);

#define PRINT_PX do { \
        igraph_int_t j; \
        printf("PX="); \
        for (j=0; j<PS; j++) { \
            printf("%" IGRAPH_PRId " ", VECTOR(*PX)[j]); \
        } \
        printf("( "); \
        for (; j<=PE; j++) { \
            printf("%" IGRAPH_PRId " ", VECTOR(*PX)[j]); \
        } \
        printf("| "); \
        for (; j<=XE; j++) { \
            printf("%" IGRAPH_PRId " ", VECTOR(*PX)[j]); \
        } \
        printf(") "); \
        for (; j<igraph_vector_int_size(PX); j++) { \
            printf("%" IGRAPH_PRId " ", VECTOR(*PX)[j]); \
        } \
        printf("\n"); \
    } while (0);

#define PRINT_PX1 do { \
        igraph_int_t j; \
        printf("PX="); \
        for (j=0; j<PS; j++) { \
            printf("%" IGRAPH_PRId " ", VECTOR(*PX)[j]); \
        } \
        printf("( "); \
        for (; j<=*PE; j++) { \
            printf("%" IGRAPH_PRId " ", VECTOR(*PX)[j]); \
        } \
        printf("| "); \
        for (; j<=XE; j++) { \
            printf("%" IGRAPH_PRId " ", VECTOR(*PX)[j]); \
        } \
        printf(") "); \
        for (; j<igraph_vector_int_size(PX); j++) { \
            printf("%" IGRAPH_PRId " ", VECTOR(*PX)[j]); \
        } \
        printf("\n"); \
    } while (0)

static igraph_error_t igraph_i_maximal_cliques_reorder_adjlists(
        const igraph_vector_int_t *PX,
        igraph_int_t PS, igraph_int_t PE,
        igraph_int_t XS, igraph_int_t XE,
        const igraph_vector_int_t *pos,
        igraph_adjlist_t *adjlist) {
    igraph_int_t j;
    igraph_int_t sPS = PS + 1, sPE = PE + 1;

    IGRAPH_UNUSED(XS);

    for (j = PS; j <= XE; j++) {
        igraph_int_t av = VECTOR(*PX)[j];
        igraph_vector_int_t *avneis = igraph_adjlist_get(adjlist, av);
        igraph_int_t *avp = VECTOR(*avneis);
        igraph_int_t avlen = igraph_vector_int_size(avneis);
        igraph_int_t *ave = avp + avlen;
        igraph_int_t *avnei = avp, *pp = avp;

        for (; avnei < ave; avnei++) {
            igraph_int_t avneipos = VECTOR(*pos)[(*avnei)];
            if (avneipos >= sPS && avneipos <= sPE) {
                if (pp != avnei) {
                    igraph_int_t tmp = *avnei;
                    *avnei = *pp;
                    *pp = tmp;
                }
                pp++;
            }
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_maximal_cliques_select_pivot(
        const igraph_vector_int_t *PX,
        igraph_int_t PS, igraph_int_t PE,
        igraph_int_t XS, igraph_int_t XE,
        const igraph_vector_int_t *pos,
        const igraph_adjlist_t *adjlist,
        igraph_int_t *pivot,
        igraph_vector_int_t *nextv,
        igraph_int_t oldPS, igraph_int_t oldXE) {
    igraph_vector_int_t *pivotvectneis;
    igraph_int_t j, pivotvectlen;
    igraph_int_t i, usize = -1;
    igraph_int_t soldPS = oldPS + 1, soldXE = oldXE + 1, sPS = PS + 1, sPE = PE + 1;

    IGRAPH_UNUSED(XS);

    /* Choose a pivotvect, and bring up P vertices at the same time */
    for (i = PS; i <= XE; i++) {
        igraph_int_t av = VECTOR(*PX)[i];
        igraph_vector_int_t *avneis = igraph_adjlist_get(adjlist, av);
        igraph_int_t *avp = VECTOR(*avneis);
        igraph_int_t avlen = igraph_vector_int_size(avneis);
        igraph_int_t *ave = avp + avlen;
        igraph_int_t *avnei = avp, *pp = avp;

        for (; avnei < ave; avnei++) {
            igraph_int_t avneipos = VECTOR(*pos)[(*avnei)];
            if (avneipos < soldPS || avneipos > soldXE) {
                break;
            }
            if (avneipos >= sPS && avneipos <= sPE) {
                if (pp != avnei) {
                    igraph_int_t tmp = *avnei;
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

    IGRAPH_CHECK(igraph_vector_int_push_back(nextv, -1));
    pivotvectneis = igraph_adjlist_get(adjlist, *pivot);
    pivotvectlen = igraph_vector_int_size(pivotvectneis);

    for (j = PS; j <= PE; j++) {
        igraph_int_t vcand = VECTOR(*PX)[j];
        igraph_bool_t nei = false;
        igraph_int_t k = 0;
        for (k = 0; k < pivotvectlen; k++) {
            igraph_int_t unv = VECTOR(*pivotvectneis)[k];
            igraph_int_t unvpos = VECTOR(*pos)[unv];
            if (unvpos < sPS || unvpos > sPE) {
                break;
            }
            if (unv == vcand) {
                nei = true;
                break;
            }
        }
        if (!nei) {
            IGRAPH_CHECK(igraph_vector_int_push_back(nextv, vcand));
        }
    }

    return IGRAPH_SUCCESS;
}

#define SWAP(p1,p2) do { \
        igraph_int_t v1=VECTOR(*PX)[p1]; \
        igraph_int_t v2=VECTOR(*PX)[p2]; \
        VECTOR(*PX)[p1] = v2; \
        VECTOR(*PX)[p2] = v1; \
        VECTOR(*pos)[v1] = (p2)+1; \
        VECTOR(*pos)[v2] = (p1)+1; \
    } while (0)

static igraph_error_t igraph_i_maximal_cliques_down(igraph_vector_int_t *PX,
                                         igraph_int_t PS, igraph_int_t PE,
                                         igraph_int_t XS, igraph_int_t XE,
                                         igraph_vector_int_t *pos,
                                         igraph_adjlist_t *adjlist, igraph_int_t mynextv,
                                         igraph_vector_int_t *R,
                                         igraph_int_t *newPS, igraph_int_t *newXE) {

    igraph_vector_int_t *vneis = igraph_adjlist_get(adjlist, mynextv);
    igraph_int_t j, vneislen = igraph_vector_int_size(vneis);
    igraph_int_t sPS = PS + 1, sPE = PE + 1, sXS = XS + 1, sXE = XE + 1;

    *newPS = PE + 1; *newXE = XS - 1;
    for (j = 0; j < vneislen; j++) {
        igraph_int_t vnei = VECTOR(*vneis)[j];
        igraph_int_t vneipos = VECTOR(*pos)[vnei];
        if (vneipos >= sPS && vneipos <= sPE) {
            (*newPS)--;
            SWAP(vneipos - 1, *newPS);
        } else if (vneipos >= sXS && vneipos <= sXE) {
            (*newXE)++;
            SWAP(vneipos - 1, *newXE);
        }
    }

    IGRAPH_CHECK(igraph_vector_int_push_back(R, mynextv));

    return IGRAPH_SUCCESS;
}

#undef SWAP

static igraph_error_t igraph_i_maximal_cliques_PX(igraph_vector_int_t *PX,
    igraph_int_t PS, igraph_int_t *PE, igraph_int_t *XS, igraph_int_t XE,
    igraph_vector_int_t *pos, igraph_adjlist_t *adjlist, igraph_int_t v,
    igraph_vector_int_t *H
) {

    igraph_int_t vpos = VECTOR(*pos)[v] - 1;
    igraph_int_t tmp = VECTOR(*PX)[*PE];

    IGRAPH_UNUSED(PS);
    IGRAPH_UNUSED(XE);
    IGRAPH_UNUSED(adjlist);

    VECTOR(*PX)[vpos] = tmp;
    VECTOR(*PX)[*PE] = v;
    VECTOR(*pos)[v] = (*PE) + 1;
    VECTOR(*pos)[tmp] = vpos + 1;
    (*PE)--; (*XS)--;
    IGRAPH_CHECK(igraph_vector_int_push_back(H, v));

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_maximal_cliques_up(
    igraph_vector_int_t *PX, igraph_int_t PS, igraph_int_t PE,
    igraph_int_t XS, igraph_int_t XE, igraph_vector_int_t *pos,
    igraph_adjlist_t *adjlist,
    igraph_vector_int_t *R,
    igraph_vector_int_t *H
) {
    igraph_int_t vv;

    IGRAPH_UNUSED(PS);
    IGRAPH_UNUSED(PE);
    IGRAPH_UNUSED(XE);
    IGRAPH_UNUSED(adjlist);

    igraph_vector_int_pop_back(R);

    while ((vv = igraph_vector_int_pop_back(H)) != -1) {
        igraph_int_t vvpos = VECTOR(*pos)[vv];
        igraph_int_t tmp = VECTOR(*PX)[XS];
        VECTOR(*PX)[XS] = vv;
        VECTOR(*PX)[vvpos - 1] = tmp;
        VECTOR(*pos)[vv] = XS + 1;
        VECTOR(*pos)[tmp] = vvpos;
        PE++; XS++;
    }

    return IGRAPH_SUCCESS;
}


/* igraph_maximal_cliques */

igraph_error_t igraph_i_maximal_cliques(
        const igraph_t *graph,
        igraph_vector_int_list_t *res,
        igraph_int_t min_size, igraph_int_t max_size,
        igraph_int_t max_results);

#define IGRAPH_MC_ORIG
#include "maximal_cliques_template.h"
#undef IGRAPH_MC_ORIG

/**
 * \function igraph_maximal_cliques
 * \brief Finds all maximal cliques in a graph.
 *
 * This function lists maximal cliques within a size range, ignoring edge
 * directions. A clique is a subset of vertices in which all vertex pairs are
 * connected. A \em maximal clique is a clique which is not a strict subset
 * of any larger clique.
 *
 * </para><para>
 * No guarantees are given about the order in which cliques are returned.
 *
 * </para><para>
 * The current implementation uses a modified Bron-Kerbosch algorithm due to
 * Eppstein, Löffler and Strash.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * David Eppstein, Maarten Löffler, Darren Strash:
 * Listing All Maximal Cliques in Sparse Graphs in Near-Optimal Time.
 * Algorithms and Computation, Lecture Notes in Computer Science,
 * volume 6506, pp 403-414 (2010).
 * https://doi.org/10.1007/978-3-642-17517-6_36
 * https://arxiv.org/abs/1006.5440
 *
 * \param graph The input graph. Edge directions are ignored.
 * \param res Pointer to list of integer vectors. The maximal cliques
 *   will be returned here as vectors of vertex IDs. Note that vertices
 *   of a clique may be returned in arbitrary order.
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \param max_results At most this many cliques will be recorded. If
 *   negative, or \ref IGRAPH_UNLIMITED, no limit is applied.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_independent_vertex_sets() to find maximal
 * independent sets, which are cliques of the complement graph;
 * \ref igraph_clique_number() to find the size of the largest clique;
 * \ref igraph_cliques() to find all cliques.
 *
 * Time complexity: O(d(n-d)3^(d/3)) worst case, d is the degeneracy
 * of the graph, this is typically small for sparse graphs.
 *
 * \example examples/simple/igraph_maximal_cliques.c
 */

igraph_error_t igraph_maximal_cliques(
        const igraph_t *graph,
        igraph_vector_int_list_t *res,
        igraph_int_t min_size, igraph_int_t max_size,
        igraph_int_t max_results) {
    return igraph_i_maximal_cliques(graph, res, min_size, max_size, max_results);
}


/* igraph_maximal_cliques_count */

igraph_error_t igraph_i_maximal_cliques_count(
        const igraph_t *graph,
        igraph_int_t *res,
        igraph_int_t min_size, igraph_int_t max_size,
        igraph_int_t max_results);

#define IGRAPH_MC_COUNT
#include "maximal_cliques_template.h"
#undef IGRAPH_MC_COUNT

/**
 * \function igraph_maximal_cliques_count
 * \brief Count the number of maximal cliques in a graph.
 *
 * See \ref igraph_maximal_cliques() for details.
 *
 * \param graph The input graph. Edge directions are ignored.
 * \param res Pointer to an \c igraph_int_t; the number of maximal
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

igraph_error_t igraph_maximal_cliques_count(
        const igraph_t *graph,
        igraph_int_t *res,
        igraph_int_t min_size, igraph_int_t max_size) {
    return igraph_i_maximal_cliques_count(graph, res, min_size, max_size, IGRAPH_UNLIMITED);
}

/* igraph_maximal_cliques_file */

igraph_error_t igraph_i_maximal_cliques_file(
        const igraph_t *graph,
        FILE *outfile,
        igraph_int_t min_size, igraph_int_t max_size,
        igraph_int_t max_results);

#define IGRAPH_MC_FILE
#include "maximal_cliques_template.h"
#undef IGRAPH_MC_FILE

/**
 * \function igraph_maximal_cliques_file
 * \brief Find maximal cliques and write them to a file.
 *
 * This function enumerates all maximal cliques within a size range
 * and writes them to file. See \ref igraph_maximal_cliques() for
 * details
 *
 * \param graph The input graph. Edge directions are ignored.
 * \param outfile Pointer to the output file, it should be writable.
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \param max_results At most this many cliques will be output. If
 *   negative, or \ref IGRAPH_UNLIMITED, no limit is applied.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_cliques().
 *
 * Time complexity: O(d(n-d)3^(d/3)) worst case, d is the degeneracy
 * of the graph, this is typically small for sparse graphs.*
 *
 */

igraph_error_t igraph_maximal_cliques_file(
        const igraph_t *graph,
        FILE *outfile,
        igraph_int_t min_size, igraph_int_t max_size,
        igraph_int_t max_results) {
    return igraph_i_maximal_cliques_file(graph, outfile, min_size, max_size, max_results);
}


/* igraph_maximal_cliques_subset */

igraph_error_t igraph_i_maximal_cliques_subset(
        const igraph_t *graph, const igraph_vector_int_t *subset,
        igraph_vector_int_list_t *res, igraph_int_t *no,
        FILE *outfile, igraph_int_t min_size, igraph_int_t max_size,
        igraph_int_t max_results);

#define IGRAPH_MC_FULL
#include "maximal_cliques_template.h"
#undef IGRAPH_MC_FULL

/**
 * \function igraph_maximal_cliques_subset
 * \brief Maximal cliques for a subset of initial vertices.
 *
 * This function enumerates all maximal cliques for a subset of initial
 * vertices and writes them to file. See \ref igraph_maximal_cliques()
 * for details.
 *
 * \param graph The input graph. Edge directions are ignored.
 * \param subset Pointer to an \c  igraph_vector_int_t containing the
 *   subset of initial vertices.
 * \param res Pointer to a list of integer vectors; the cliques will be
 *   stored here.
 * \param no Pointer to an \c igraph_int_t; the number of maximal
 *   cliques will be stored here.
 * \param outfile Pointer to an output file or \c NULL.
 *   When not \c NULL, the file should be writable.
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \param max_results At most this many cliques will be recorded. If
 *   negative, or \ref IGRAPH_UNLIMITED, no limit is applied.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_cliques().
 *
 * Time complexity: O(d(n-d)3^(d/3)) worst case, d is the degeneracy
 * of the graph, this is typically small for sparse graphs.
 *
 */

igraph_error_t igraph_maximal_cliques_subset(
        const igraph_t *graph, const igraph_vector_int_t *subset,
        igraph_vector_int_list_t *res, igraph_int_t *no, FILE *outfile,
        igraph_int_t min_size, igraph_int_t max_size,
        igraph_int_t max_results) {
    return igraph_i_maximal_cliques_subset(graph, subset, res, no, outfile, min_size, max_size, max_results);
}


/* igraph_maximal_cliques_callback */

igraph_error_t igraph_i_maximal_cliques_callback(
        const igraph_t *graph,
        igraph_clique_handler_t *cliquehandler_fn, void *arg,
        igraph_int_t min_size, igraph_int_t max_size,
        igraph_int_t max_results);

#define IGRAPH_MC_CALLBACK
#include "maximal_cliques_template.h"
#undef IGRAPH_MC_CALLBACK

/**
 * \function igraph_maximal_cliques_callback
 * \brief Finds maximal cliques in a graph and calls a function for each one.
 *
 * This function enumerates all maximal cliques within the given size range
 * and calls \p cliquehandler_fn for each of them. The cliques are passed to the
 * callback function as a pointer to an \ref igraph_vector_int_t. The vector is
 * owned by the maximal clique search routine so users are expected to make a
 * copy of the vector using \ref igraph_vector_int_init_copy() if they want to
 * hold on to it.
 *
 * \param graph The input graph. Edge directions are ignored.
 * \param cliquehandler_fn Callback function to be called for each clique.
 * See also \ref igraph_clique_handler_t.
 * \param arg Extra argument to supply to \p cliquehandler_fn.
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_cliques(), \ref igraph_cliques_callback().
 *
 * Time complexity: O(d(n-d)3^(d/3)) worst case, d is the degeneracy
 * of the graph, this is typically small for sparse graphs.
 *
 */

igraph_error_t igraph_maximal_cliques_callback(
        const igraph_t *graph,
        igraph_int_t min_size, igraph_int_t max_size,
        igraph_clique_handler_t *cliquehandler_fn, void *arg) {
    return igraph_i_maximal_cliques_callback(graph, cliquehandler_fn, arg, min_size, max_size, IGRAPH_UNLIMITED);
}


/* igraph_maximal_cliques_hist */

igraph_error_t igraph_i_maximal_cliques_hist(
        const igraph_t *graph,
        igraph_vector_t *hist,
        igraph_int_t min_size, igraph_int_t max_size,
        igraph_int_t max_results);

#define IGRAPH_MC_HIST
#include "maximal_cliques_template.h"
#undef IGRAPH_MC_HIST

/**
 * \function igraph_maximal_cliques_hist
 * \brief Counts the number of maximal cliques of each size in a graph.
 *
 * This function counts how many maximal cliques of each size are present in
 * the graph. Maximal cliques of size one are simply isolated vertices.
 *
 * \param graph The input graph. Edge directions are ignored.
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
 * \sa \ref igraph_maximal_cliques(), \ref igraph_clique_size_hist().
 *
 * Time complexity: O(d(n-d)3^(d/3)) worst case, d is the degeneracy
 * of the graph, this is typically small for sparse graphs.
 *
 */

igraph_error_t igraph_maximal_cliques_hist(
        const igraph_t *graph,
        igraph_vector_t *hist,
        igraph_int_t min_size, igraph_int_t max_size) {
    return igraph_i_maximal_cliques_hist(graph, hist, min_size, max_size, IGRAPH_UNLIMITED);
}
