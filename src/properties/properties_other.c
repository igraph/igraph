/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_structural.h"
#include "igraph_transitivity.h"
#include "igraph_paths.h"
#include "core/math.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_centrality.h"
#include "igraph_components.h"
#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_dqueue.h"
#include "igraph_attributes.h"
#include "igraph_neighborhood.h"
#include "igraph_operators.h"
#include "igraph_stack.h"
#include "igraph_topology.h"
#include "igraph_qsort.h"
#include "igraph_error.h"
#include "config.h"

#include <string.h>
#include <limits.h>

#include "core/fixed_vectorlist.h"
#include "core/indheap.h"
#include "core/interruption.h"
#include "properties/properties_internal.h"

/**
 * \section about_structural
 *
 * <para>These functions usually calculate some structural property
 * of a graph, like its diameter, the degree of the nodes, etc.</para>
 */

/**
 * \ingroup structural
 * \function igraph_subcomponent
 * \brief The vertices in the same component as a given vertex.
 *
 * \param graph The graph object.
 * \param res The result, vector with the ids of the vertices in the
 *        same component.
 * \param vertex The id of the vertex of which the component is
 *        searched.
 * \param mode Type of the component for directed graphs, possible
 *        values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the set of vertices reachable \em from the
 *          \p vertex,
 *        \cli IGRAPH_IN
 *          the set of vertices from which the
 *          \p vertex is reachable.
 *        \cli IGRAPH_ALL
 *          the graph is considered as an
 *          undirected graph. Note that this is \em not the same
 *          as the union of the previous two.
 *        \endclist
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *          not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           \p vertex is an invalid vertex id
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument passed.
 *        \endclist
 *
 * Time complexity: O(|V|+|E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph.
 *
 * \sa \ref igraph_induced_subgraph() if you want a graph object consisting only
 * a given set of vertices and the edges between them.
 */

int igraph_subcomponent(const igraph_t *graph, igraph_vector_t *res, igraph_real_t vertex,
                        igraph_neimode_t mode) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_t q = IGRAPH_DQUEUE_NULL;
    char *already_added;
    long int i, vsize;
    igraph_vector_t tmp = IGRAPH_VECTOR_NULL;

    if (!IGRAPH_FINITE(vertex) || vertex < 0 || vertex >= no_of_nodes) {
        IGRAPH_ERROR("subcomponent failed", IGRAPH_EINVVID);
    }
    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("invalid mode argument", IGRAPH_EINVMODE);
    }

    already_added = igraph_Calloc(no_of_nodes, char);
    if (already_added == 0) {
        IGRAPH_ERROR("subcomponent failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, already_added);

    igraph_vector_clear(res);

    IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_dqueue_push(&q, vertex));
    IGRAPH_CHECK(igraph_vector_push_back(res, vertex));
    already_added[(long int)vertex] = 1;

    while (!igraph_dqueue_empty(&q)) {
        long int actnode = (long int) igraph_dqueue_pop(&q);

        IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_CHECK(igraph_neighbors(graph, &tmp, (igraph_integer_t) actnode,
                                      mode));
        vsize = igraph_vector_size(&tmp);
        for (i = 0; i < vsize; i++) {
            long int neighbor = (long int) VECTOR(tmp)[i];

            if (already_added[neighbor]) {
                continue;
            }
            already_added[neighbor] = 1;
            IGRAPH_CHECK(igraph_vector_push_back(res, neighbor));
            IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
        }
    }

    igraph_dqueue_destroy(&q);
    igraph_vector_destroy(&tmp);
    igraph_Free(already_added);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

/**
 * \ingroup structural
 * \function igraph_reciprocity
 * \brief Calculates the reciprocity of a directed graph.
 *
 * </para><para>
 * The measure of reciprocity defines the proportion of mutual
 * connections, in a directed graph. It is most commonly defined as
 * the probability that the opposite counterpart of a directed edge is
 * also included in the graph. In adjacency matrix notation:
 * <code>sum(i, j, (A.*A')ij) / sum(i, j, Aij)</code>, where
 * <code>A.*A'</code> is the element-wise product of matrix
 * <code>A</code> and its transpose. This measure is
 * calculated if the \p mode argument is \c
 * IGRAPH_RECIPROCITY_DEFAULT.
 *
 * </para><para>
 * Prior to igraph version 0.6, another measure was implemented,
 * defined as the probability of mutual connection between a vertex
 * pair if we know that there is a (possibly non-mutual) connection
 * between them. In other words, (unordered) vertex pairs are
 * classified into three groups: (1) disconnected, (2)
 * non-reciprocally connected, (3) reciprocally connected.
 * The result is the size of group (3), divided by the sum of group
 * sizes (2)+(3). This measure is calculated if \p mode is \c
 * IGRAPH_RECIPROCITY_RATIO.
 *
 * \param graph The graph object.
 * \param res Pointer to an \c igraph_real_t which will contain the result.
 * \param ignore_loops Whether to ignore loop edges.
 * \param mode Type of reciprocity to calculate, possible values are
 *    \c IGRAPH_RECIPROCITY_DEFAULT and \c IGRAPH_RECIPROCITY_RATIO,
 *    please see their description above.
 * \return Error code:
 *         \c IGRAPH_EINVAL: graph has no edges
 *         \c IGRAPH_ENOMEM: not enough memory for
 *         temporary data.
 *
 * Time complexity: O(|V|+|E|), |V| is the number of vertices,
 * |E| is the number of edges.
 *
 * \example examples/simple/igraph_reciprocity.c
 */

int igraph_reciprocity(const igraph_t *graph, igraph_real_t *res,
                       igraph_bool_t ignore_loops,
                       igraph_reciprocity_t mode) {

    igraph_integer_t nonrec = 0, rec = 0, loops = 0;
    igraph_vector_t inneis, outneis;
    long int i;
    long int no_of_nodes = igraph_vcount(graph);

    if (mode != IGRAPH_RECIPROCITY_DEFAULT &&
        mode != IGRAPH_RECIPROCITY_RATIO) {
        IGRAPH_ERROR("Invalid reciprocity type", IGRAPH_EINVAL);
    }

    /* THIS IS AN EXIT HERE !!!!!!!!!!!!!! */
    if (!igraph_is_directed(graph)) {
        *res = 1.0;
        return 0;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&inneis, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&outneis, 0);

    for (i = 0; i < no_of_nodes; i++) {
        long int ip, op;
        igraph_neighbors(graph, &inneis, (igraph_integer_t) i, IGRAPH_IN);
        igraph_neighbors(graph, &outneis, (igraph_integer_t) i, IGRAPH_OUT);

        ip = op = 0;
        while (ip < igraph_vector_size(&inneis) &&
               op < igraph_vector_size(&outneis)) {
            if (VECTOR(inneis)[ip] < VECTOR(outneis)[op]) {
                nonrec += 1;
                ip++;
            } else if (VECTOR(inneis)[ip] > VECTOR(outneis)[op]) {
                nonrec += 1;
                op++;
            } else {

                /* loop edge? */
                if (VECTOR(inneis)[ip] == i) {
                    loops += 1;
                    if (!ignore_loops) {
                        rec += 1;
                    }
                } else {
                    rec += 1;
                }

                ip++;
                op++;
            }
        }
        nonrec += (igraph_vector_size(&inneis) - ip) +
                  (igraph_vector_size(&outneis) - op);
    }

    if (mode == IGRAPH_RECIPROCITY_DEFAULT) {
        if (ignore_loops) {
            *res = (igraph_real_t) rec / (igraph_ecount(graph) - loops);
        } else {
            *res = (igraph_real_t) rec / (igraph_ecount(graph));
        }
    } else if (mode == IGRAPH_RECIPROCITY_RATIO) {
        *res = (igraph_real_t) rec / (rec + nonrec);
    }

    igraph_vector_destroy(&inneis);
    igraph_vector_destroy(&outneis);
    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

/**
 * \function igraph_constraint
 * \brief Burt's constraint scores.
 *
 * </para><para>
 * This function calculates Burt's constraint scores for the given
 * vertices, also known as structural holes.
 *
 * </para><para>
 * Burt's constraint is higher if ego has less, or mutually stronger
 * related (i.e. more redundant) contacts. Burt's measure of
 * constraint, C[i], of vertex i's ego network V[i], is defined for
 * directed and valued graphs,
 * <blockquote><para>
 * C[i] = sum( sum( (p[i,q] p[q,j])^2, q in V[i], q != i,j ), j in
 * V[], j != i)
 * </para></blockquote>
 * for a graph of order (i.e. number of vertices) N, where proportional
 * tie strengths are defined as
 * <blockquote><para>
 * p[i,j]=(a[i,j]+a[j,i]) / sum(a[i,k]+a[k,i], k in V[i], k != i),
 * </para></blockquote>
 * a[i,j] are elements of A and
 * the latter being the graph adjacency matrix. For isolated vertices,
 * constraint is undefined.
 *
 * </para><para>
 * Burt, R.S. (2004). Structural holes and good ideas. American
 * Journal of Sociology 110, 349-399.
 *
 * </para><para>
 * The first R version of this function was contributed by Jeroen
 * Bruggeman.
 * \param graph A graph object.
 * \param res Pointer to an initialized vector, the result will be
 *        stored here. The vector will be resized to have the
 *        appropriate size for holding the result.
 * \param vids Vertex selector containing the vertices for which the
 *        constraint should be calculated.
 * \param weights Vector giving the weights of the edges. If it is
 *        \c NULL then each edge is supposed to have the same weight.
 * \return Error code.
 *
 * Time complexity: O(|V|+E|+n*d^2), n is the number of vertices for
 * which the constraint is calculated and d is the average degree, |V|
 * is the number of vertices, |E| the number of edges in the
 * graph. If the weights argument is \c NULL then the time complexity
 * is O(|V|+n*d^2).
 */

int igraph_constraint(const igraph_t *graph, igraph_vector_t *res,
                      igraph_vs_t vids, const igraph_vector_t *weights) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_vit_t vit;
    long int nodes_to_calc;
    long int a, b, c, i, j, q, vsize, vsize2;
    igraph_integer_t edge, from, to, edge2;

    igraph_vector_t contrib;
    igraph_vector_t degree;
    igraph_vector_t ineis_in, ineis_out, jneis_in, jneis_out;

    if (weights != 0 && igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid length of weight vector", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&contrib, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&ineis_in, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&ineis_out, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&jneis_in, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&jneis_out, 0);

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    nodes_to_calc = IGRAPH_VIT_SIZE(vit);

    if (weights == 0) {
        IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
                                   IGRAPH_ALL, IGRAPH_NO_LOOPS));
    } else {
        for (a = 0; a < no_of_edges; a++) {
            igraph_edge(graph, (igraph_integer_t) a, &from, &to);
            if (from != to) {
                VECTOR(degree)[(long int) from] += VECTOR(*weights)[a];
                VECTOR(degree)[(long int) to  ] += VECTOR(*weights)[a];
            }
        }
    }

    IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
    igraph_vector_null(res);

    for (a = 0; a < nodes_to_calc; a++, IGRAPH_VIT_NEXT(vit)) {
        i = IGRAPH_VIT_GET(vit);

        /* get neighbors of i */
        IGRAPH_CHECK(igraph_incident(graph, &ineis_in, (igraph_integer_t) i,
                                     IGRAPH_IN));
        IGRAPH_CHECK(igraph_incident(graph, &ineis_out, (igraph_integer_t) i,
                                     IGRAPH_OUT));

        /* NaN for isolates */
        if (igraph_vector_size(&ineis_in) == 0 &&
            igraph_vector_size(&ineis_out) == 0) {
            VECTOR(*res)[a] = IGRAPH_NAN;
        }

        /* zero their contribution */
        vsize = igraph_vector_size(&ineis_in);
        for (b = 0; b < vsize; b++) {
            edge = (igraph_integer_t) VECTOR(ineis_in)[b];
            j = (long int) IGRAPH_OTHER(graph, edge, i);
            VECTOR(contrib)[j] = 0.0;
        }
        vsize = igraph_vector_size(&ineis_out);
        for (b = 0; b < vsize; b++) {
            edge = (igraph_integer_t) VECTOR(ineis_out)[b];
            j = (long int) IGRAPH_OTHER(graph, edge, i);
            VECTOR(contrib)[j] = 0.0;
        }

        /* add the direct contributions, in-neighbors and out-neighbors */
        vsize = igraph_vector_size(&ineis_in);
        for (b = 0; b < vsize; b++) {
            edge = (igraph_integer_t) VECTOR(ineis_in)[b];
            j = (long int) IGRAPH_OTHER(graph, edge, i);
            if (i != j) {     /* excluding loops */
                if (weights) {
                    VECTOR(contrib)[j] +=
                        VECTOR(*weights)[(long int)edge] / VECTOR(degree)[i];
                } else {
                    VECTOR(contrib)[j] += 1.0 / VECTOR(degree)[i];
                }
            }
        }
        if (igraph_is_directed(graph)) {
            vsize = igraph_vector_size(&ineis_out);
            for (b = 0; b < vsize; b++) {
                edge = (igraph_integer_t) VECTOR(ineis_out)[b];
                j = (long int) IGRAPH_OTHER(graph, edge, i);
                if (i != j) {
                    if (weights) {
                        VECTOR(contrib)[j] +=
                            VECTOR(*weights)[(long int)edge] / VECTOR(degree)[i];
                    } else {
                        VECTOR(contrib)[j] += 1.0 / VECTOR(degree)[i];
                    }
                }
            }
        }

        /* add the indirect contributions, in-in, in-out, out-in, out-out */
        vsize = igraph_vector_size(&ineis_in);
        for (b = 0; b < vsize; b++) {
            edge = (igraph_integer_t) VECTOR(ineis_in)[b];
            j = (long int) IGRAPH_OTHER(graph, edge, i);
            if (i == j) {
                continue;
            }
            IGRAPH_CHECK(igraph_incident(graph, &jneis_in, (igraph_integer_t) j,
                                         IGRAPH_IN));
            IGRAPH_CHECK(igraph_incident(graph, &jneis_out, (igraph_integer_t) j,
                                         IGRAPH_OUT));
            vsize2 = igraph_vector_size(&jneis_in);
            for (c = 0; c < vsize2; c++) {
                edge2 = (igraph_integer_t) VECTOR(jneis_in)[c];
                q = (long int) IGRAPH_OTHER(graph, edge2, j);
                if (j != q) {
                    if (weights) {
                        VECTOR(contrib)[q] +=
                            VECTOR(*weights)[(long int)edge] *
                            VECTOR(*weights)[(long int)edge2] /
                            VECTOR(degree)[i] / VECTOR(degree)[j];
                    } else {
                        VECTOR(contrib)[q] += 1 / VECTOR(degree)[i] / VECTOR(degree)[j];
                    }
                }
            }
            if (igraph_is_directed(graph)) {
                vsize2 = igraph_vector_size(&jneis_out);
                for (c = 0; c < vsize2; c++) {
                    edge2 = (igraph_integer_t) VECTOR(jneis_out)[c];
                    q = (long int) IGRAPH_OTHER(graph, edge2, j);
                    if (j != q) {
                        if (weights) {
                            VECTOR(contrib)[q] +=
                                VECTOR(*weights)[(long int)edge] *
                                VECTOR(*weights)[(long int)edge2] /
                                VECTOR(degree)[i] / VECTOR(degree)[j];
                        } else {
                            VECTOR(contrib)[q] += 1 / VECTOR(degree)[i] / VECTOR(degree)[j];
                        }
                    }
                }
            }
        }
        if (igraph_is_directed(graph)) {
            vsize = igraph_vector_size(&ineis_out);
            for (b = 0; b < vsize; b++) {
                edge = (igraph_integer_t) VECTOR(ineis_out)[b];
                j = (long int) IGRAPH_OTHER(graph, edge, i);
                if (i == j) {
                    continue;
                }
                IGRAPH_CHECK(igraph_incident(graph, &jneis_in, (igraph_integer_t) j,
                                             IGRAPH_IN));
                IGRAPH_CHECK(igraph_incident(graph, &jneis_out, (igraph_integer_t) j,
                                             IGRAPH_OUT));
                vsize2 = igraph_vector_size(&jneis_in);
                for (c = 0; c < vsize2; c++) {
                    edge2 = (igraph_integer_t) VECTOR(jneis_in)[c];
                    q = (long int) IGRAPH_OTHER(graph, edge2, j);
                    if (j != q) {
                        if (weights) {
                            VECTOR(contrib)[q] +=
                                VECTOR(*weights)[(long int)edge] *
                                VECTOR(*weights)[(long int)edge2] /
                                VECTOR(degree)[i] / VECTOR(degree)[j];
                        } else {
                            VECTOR(contrib)[q] += 1 / VECTOR(degree)[i] / VECTOR(degree)[j];
                        }
                    }
                }
                vsize2 = igraph_vector_size(&jneis_out);
                for (c = 0; c < vsize2; c++) {
                    edge2 = (igraph_integer_t) VECTOR(jneis_out)[c];
                    q = (long int) IGRAPH_OTHER(graph, edge2, j);
                    if (j != q) {
                        if (weights) {
                            VECTOR(contrib)[q] +=
                                VECTOR(*weights)[(long int)edge] *
                                VECTOR(*weights)[(long int)edge2] /
                                VECTOR(degree)[i] / VECTOR(degree)[j];
                        } else {
                            VECTOR(contrib)[q] += 1 / VECTOR(degree)[i] / VECTOR(degree)[j];
                        }
                    }
                }
            }
        }

        /* squared sum of the contributions */
        vsize = igraph_vector_size(&ineis_in);
        for (b = 0; b < vsize; b++) {
            edge = (igraph_integer_t) VECTOR(ineis_in)[b];
            j = (long int) IGRAPH_OTHER(graph, edge, i);
            if (i == j) {
                continue;
            }
            VECTOR(*res)[a] += VECTOR(contrib)[j] * VECTOR(contrib)[j];
            VECTOR(contrib)[j] = 0.0;
        }
        if (igraph_is_directed(graph)) {
            vsize =  igraph_vector_size(&ineis_out);
            for (b = 0; b < vsize; b++) {
                edge = (igraph_integer_t) VECTOR(ineis_out)[b];
                j = (long int) IGRAPH_OTHER(graph, edge, i);
                if (i == j) {
                    continue;
                }
                VECTOR(*res)[a] += VECTOR(contrib)[j] * VECTOR(contrib)[j];
                VECTOR(contrib)[j] = 0.0;
            }
        }
    }

    igraph_vit_destroy(&vit);
    igraph_vector_destroy(&jneis_out);
    igraph_vector_destroy(&jneis_in);
    igraph_vector_destroy(&ineis_out);
    igraph_vector_destroy(&ineis_in);
    igraph_vector_destroy(&degree);
    igraph_vector_destroy(&contrib);
    IGRAPH_FINALLY_CLEAN(7);

    return 0;
}

/**
 * \function igraph_density
 * Calculate the density of a graph.
 *
 * </para><para>The density of a graph is simply the ratio number of
 * edges and the number of possible edges. Note that density is
 * ill-defined for graphs with multiple and/or loop edges, so consider
 * calling \ref igraph_simplify() on the graph if you know that it
 * contains multiple or loop edges.
 * \param graph The input graph object.
 * \param res Pointer to a real number, the result will be stored
 *   here.
 * \param loops Logical constant, whether to include loops in the
 *   calculation. If this constant is TRUE then
 *   loop edges are thought to be possible in the graph (this does not
 *   necessarily mean that the graph really contains any loops). If
 *   this is FALSE then the result is only correct if the graph does not
 *   contain loops.
 * \return Error code.
 *
 * Time complexity: O(1).
 */

int igraph_density(const igraph_t *graph, igraph_real_t *res,
                   igraph_bool_t loops) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_real_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);

    if (no_of_nodes == 0) {
        *res = IGRAPH_NAN;
        return 0;
    }

    if (!loops) {
        if (no_of_nodes == 1) {
            *res = IGRAPH_NAN;
        } else if (directed) {
            *res = no_of_edges / no_of_nodes / (no_of_nodes - 1);
        } else {
            *res = no_of_edges / no_of_nodes * 2.0 / (no_of_nodes - 1);
        }
    } else {
        if (directed) {
            *res = no_of_edges / no_of_nodes / no_of_nodes;
        } else {
            *res = no_of_edges / no_of_nodes * 2.0 / (no_of_nodes + 1);
        }
    }

    return 0;
}

/**
 * \function igraph_neighborhood_size
 * \brief Calculates the size of the neighborhood of a given vertex.
 *
 * The neighborhood of a given order of a vertex includes all vertices
 * which are closer to the vertex than the order. I.e., order 0 is
 * always the vertex itself, order 1 is the vertex plus its immediate
 * neighbors, order 2 is order 1 plus the immediate neighbors of the
 * vertices in order 1, etc.
 *
 * </para><para> This function calculates the size of the neighborhood
 * of the given order for the given vertices.
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the result will be
 *    stored here. It will be resized as needed.
 * \param vids The vertices for which the calculation is performed.
 * \param order Integer giving the order of the neighborhood.
 * \param mode Specifies how to use the direction of the edges if a
 *   directed graph is analyzed. For \c IGRAPH_OUT only the outgoing
 *   edges are followed, so all vertices reachable from the source
 *   vertex in at most \c order steps are counted. For \c IGRAPH_IN
 *   all vertices from which the source vertex is reachable in at most
 *   \c order steps are counted. \c IGRAPH_ALL ignores the direction
 *   of the edges. This argument is ignored for undirected graphs.
 * \param mindist The minimum distance to include a vertex in the counting.
 *   If this is one, then the starting vertex is not counted. If this is
 *   two, then its neighbors are not counted, either, etc.
 * \return Error code.
 *
 * \sa \ref igraph_neighborhood() for calculating the actual neighborhood,
 * \ref igraph_neighborhood_graphs() for creating separate graphs from
 * the neighborhoods.
 *
 * Time complexity: O(n*d*o), where n is the number vertices for which
 * the calculation is performed, d is the average degree, o is the order.
 */

int igraph_neighborhood_size(const igraph_t *graph, igraph_vector_t *res,
                             igraph_vs_t vids, igraph_integer_t order,
                             igraph_neimode_t mode,
                             igraph_integer_t mindist) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_t q;
    igraph_vit_t vit;
    long int i, j;
    long int *added;
    igraph_vector_t neis;

    if (order < 0) {
        IGRAPH_ERROR("Negative order in neighborhood size", IGRAPH_EINVAL);
    }

    if (mindist < 0 || mindist > order) {
        IGRAPH_ERROR("Minimum distance should be between zero and order",
                     IGRAPH_EINVAL);
    }

    added = igraph_Calloc(no_of_nodes, long int);
    if (added == 0) {
        IGRAPH_ERROR("Cannot calculate neighborhood size", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));

    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        long int node = IGRAPH_VIT_GET(vit);
        long int size = mindist == 0 ? 1 : 0;
        added[node] = i + 1;
        igraph_dqueue_clear(&q);
        if (order > 0) {
            igraph_dqueue_push(&q, node);
            igraph_dqueue_push(&q, 0);
        }

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);
            long int actdist = (long int) igraph_dqueue_pop(&q);
            long int n;
            igraph_neighbors(graph, &neis, (igraph_integer_t) actnode, mode);
            n = igraph_vector_size(&neis);

            if (actdist < order - 1) {
                /* we add them to the q */
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
                        if (actdist + 1 >= mindist) {
                            size++;
                        }
                    }
                }
            } else {
                /* we just count them, but don't add them */
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (actdist + 1 >= mindist) {
                            size++;
                        }
                    }
                }
            }

        } /* while q not empty */

        VECTOR(*res)[i] = size;
    } /* for VIT, i */

    igraph_vector_destroy(&neis);
    igraph_vit_destroy(&vit);
    igraph_dqueue_destroy(&q);
    igraph_Free(added);
    IGRAPH_FINALLY_CLEAN(4);

    return 0;
}

/**
 * \function igraph_neighborhood
 * Calculate the neighborhood of vertices.
 *
 * The neighborhood of a given order of a vertex includes all vertices
 * which are closer to the vertex than the order. I.e., order 0 is
 * always the vertex itself, order 1 is the vertex plus its immediate
 * neighbors, order 2 is order 1 plus the immediate neighbors of the
 * vertices in order 1, etc.
 *
 * </para><para> This function calculates the vertices within the
 * neighborhood of the specified vertices.
 * \param graph The input graph.
 * \param res An initialized pointer vector. Note that the objects
 *    (pointers) in the vector will \em not be freed, but the pointer
 *    vector will be resized as needed. The result of the calculation
 *    will be stored here in \ref igraph_vector_t objects.
 * \param vids The vertices for which the calculation is performed.
 * \param order Integer giving the order of the neighborhood.
 * \param mode Specifies how to use the direction of the edges if a
 *   directed graph is analyzed. For \c IGRAPH_OUT only the outgoing
 *   edges are followed, so all vertices reachable from the source
 *   vertex in at most \p order steps are included. For \c IGRAPH_IN
 *   all vertices from which the source vertex is reachable in at most
 *   \p order steps are included. \c IGRAPH_ALL ignores the direction
 *   of the edges. This argument is ignored for undirected graphs.
 * \param mindist The minimum distance to include a vertex in the counting.
 *   If this is one, then the starting vertex is not counted. If this is
 *   two, then its neighbors are not counted, either, etc.
 * \return Error code.
 *
 * \sa \ref igraph_neighborhood_size() to calculate the size of the
 * neighborhood, \ref igraph_neighborhood_graphs() for creating
 * graphs from the neighborhoods.
 *
 * Time complexity: O(n*d*o), n is the number of vertices for which
 * the calculation is performed, d is the average degree, o is the
 * order.
 */

int igraph_neighborhood(const igraph_t *graph, igraph_vector_ptr_t *res,
                        igraph_vs_t vids, igraph_integer_t order,
                        igraph_neimode_t mode, igraph_integer_t mindist) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_t q;
    igraph_vit_t vit;
    long int i, j;
    long int *added;
    igraph_vector_t neis;
    igraph_vector_t tmp;
    igraph_vector_t *newv;

    if (order < 0) {
        IGRAPH_ERROR("Negative order in neighborhood size", IGRAPH_EINVAL);
    }

    if (mindist < 0 || mindist > order) {
        IGRAPH_ERROR("Minimum distance should be between zero and order",
                     IGRAPH_EINVAL);
    }

    added = igraph_Calloc(no_of_nodes, long int);
    if (added == 0) {
        IGRAPH_ERROR("Cannot calculate neighborhood size", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
    IGRAPH_CHECK(igraph_vector_ptr_resize(res, IGRAPH_VIT_SIZE(vit)));

    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        long int node = IGRAPH_VIT_GET(vit);
        added[node] = i + 1;
        igraph_vector_clear(&tmp);
        if (mindist == 0) {
            IGRAPH_CHECK(igraph_vector_push_back(&tmp, node));
        }
        if (order > 0) {
            igraph_dqueue_push(&q, node);
            igraph_dqueue_push(&q, 0);
        }

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);
            long int actdist = (long int) igraph_dqueue_pop(&q);
            long int n;
            igraph_neighbors(graph, &neis, (igraph_integer_t) actnode, mode);
            n = igraph_vector_size(&neis);

            if (actdist < order - 1) {
                /* we add them to the q */
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_push_back(&tmp, nei));
                        }
                    }
                }
            } else {
                /* we just count them but don't add them to q */
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_push_back(&tmp, nei));
                        }
                    }
                }
            }

        } /* while q not empty */

        newv = igraph_Calloc(1, igraph_vector_t);
        if (newv == 0) {
            IGRAPH_ERROR("Cannot calculate neighborhood", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, newv);
        IGRAPH_CHECK(igraph_vector_copy(newv, &tmp));
        VECTOR(*res)[i] = newv;
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_destroy(&tmp);
    igraph_vector_destroy(&neis);
    igraph_vit_destroy(&vit);
    igraph_dqueue_destroy(&q);
    igraph_Free(added);
    IGRAPH_FINALLY_CLEAN(5);

    return 0;
}

/**
 * \function igraph_neighborhood_graphs
 * Create graphs from the neighborhood(s) of some vertex/vertices.
 *
 * The neighborhood of a given order of a vertex includes all vertices
 * which are closer to the vertex than the order. Ie. order 0 is
 * always the vertex itself, order 1 is the vertex plus its immediate
 * neighbors, order 2 is order 1 plus the immediate neighbors of the
 * vertices in order 1, etc.
 *
 * </para><para> This function finds every vertex in the neighborhood
 * of a given parameter vertex and creates a graph from these
 * vertices.
 *
 * </para><para> The first version of this function was written by
 * Vincent Matossian, thanks Vincent.
 * \param graph The input graph.
 * \param res Pointer to a pointer vector, the result will be stored
 *   here, ie. \p res will contain pointers to \c igraph_t
 *   objects. It will be resized if needed but note that the
 *   objects in the pointer vector will not be freed.
 * \param vids The vertices for which the calculation is performed.
 * \param order Integer giving the order of the neighborhood.
 * \param mode Specifies how to use the direction of the edges if a
 *   directed graph is analyzed. For \c IGRAPH_OUT only the outgoing
 *   edges are followed, so all vertices reachable from the source
 *   vertex in at most \p order steps are counted. For \c IGRAPH_IN
 *   all vertices from which the source vertex is reachable in at most
 *   \p order steps are counted. \c IGRAPH_ALL ignores the direction
 *   of the edges. This argument is ignored for undirected graphs.
 * \param mindist The minimum distance to include a vertex in the counting.
 *   If this is one, then the starting vertex is not counted. If this is
 *   two, then its neighbors are not counted, either, etc.
 * \return Error code.
 *
 * \sa \ref igraph_neighborhood_size() for calculating the neighborhood
 * sizes only, \ref igraph_neighborhood() for calculating the
 * neighborhoods (but not creating graphs).
 *
 * Time complexity: O(n*(|V|+|E|)), where n is the number vertices for
 * which the calculation is performed, |V| and |E| are the number of
 * vertices and edges in the original input graph.
 */

int igraph_neighborhood_graphs(const igraph_t *graph, igraph_vector_ptr_t *res,
                               igraph_vs_t vids, igraph_integer_t order,
                               igraph_neimode_t mode,
                               igraph_integer_t mindist) {
    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_t q;
    igraph_vit_t vit;
    long int i, j;
    long int *added;
    igraph_vector_t neis;
    igraph_vector_t tmp;
    igraph_t *newg;

    if (order < 0) {
        IGRAPH_ERROR("Negative order in neighborhood size", IGRAPH_EINVAL);
    }

    if (mindist < 0 || mindist > order) {
        IGRAPH_ERROR("Minimum distance should be between zero and order",
                     IGRAPH_EINVAL);
    }

    added = igraph_Calloc(no_of_nodes, long int);
    if (added == 0) {
        IGRAPH_ERROR("Cannot calculate neighborhood size", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
    IGRAPH_CHECK(igraph_vector_ptr_resize(res, IGRAPH_VIT_SIZE(vit)));

    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        long int node = IGRAPH_VIT_GET(vit);
        added[node] = i + 1;
        igraph_vector_clear(&tmp);
        if (mindist == 0) {
            IGRAPH_CHECK(igraph_vector_push_back(&tmp, node));
        }
        if (order > 0) {
            igraph_dqueue_push(&q, node);
            igraph_dqueue_push(&q, 0);
        }

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);
            long int actdist = (long int) igraph_dqueue_pop(&q);
            long int n;
            igraph_neighbors(graph, &neis, (igraph_integer_t) actnode, mode);
            n = igraph_vector_size(&neis);

            if (actdist < order - 1) {
                /* we add them to the q */
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_push_back(&tmp, nei));
                        }
                    }
                }
            } else {
                /* we just count them but don't add them to q */
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_push_back(&tmp, nei));
                        }
                    }
                }
            }

        } /* while q not empty */

        newg = igraph_Calloc(1, igraph_t);
        if (newg == 0) {
            IGRAPH_ERROR("Cannot create neighborhood graph", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, newg);
        if (igraph_vector_size(&tmp) < no_of_nodes) {
            IGRAPH_CHECK(igraph_induced_subgraph(graph, newg,
                                                 igraph_vss_vector(&tmp),
                                                 IGRAPH_SUBGRAPH_AUTO));
        } else {
            IGRAPH_CHECK(igraph_copy(newg, graph));
        }
        VECTOR(*res)[i] = newg;
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_destroy(&tmp);
    igraph_vector_destroy(&neis);
    igraph_vit_destroy(&vit);
    igraph_dqueue_destroy(&q);
    igraph_Free(added);
    IGRAPH_FINALLY_CLEAN(5);

    return 0;
}

/**
 * \function igraph_girth
 * \brief The girth of a graph is the length of the shortest cycle in it.
 *
 * </para><para>
 * The current implementation works for undirected graphs only,
 * directed graphs are treated as undirected graphs. Self-loops and
 * multiple edges are ignored.
 *
 * </para><para>
 * For graphs that contain no cycles, and only for such graphs, 
 * zero is returned. Note that in some applications, it is customary 
 * to define the girth of acyclic graphs to be infinity. However, infinity
 * is not representable as an \c igraph_integer_t, therefore zero is used
 * for this case.
 *
 * </para><para>
 * This implementation is based on Alon Itai and Michael Rodeh:
 * Finding a minimum circuit in a graph
 * \emb Proceedings of the ninth annual ACM symposium on Theory of
 * computing \eme, 1-10, 1977. The first implementation of this
 * function was done by Keith Briggs, thanks Keith.
 * \param graph The input graph.
 * \param girth Pointer to an integer, if not \c NULL then the result
 *     will be stored here.
 * \param circle Pointer to an initialized vector, the vertex ids in
 *     the shortest circle will be stored here. If \c NULL then it is
 *     ignored.
 * \return Error code.
 *
 * Time complexity: O((|V|+|E|)^2), |V| is the number of vertices, |E|
 * is the number of edges in the general case. If the graph has no
 * cycles at all then the function needs O(|V|+|E|) time to realize
 * this and then it stops.
 *
 * \example examples/simple/igraph_girth.c
 */

int igraph_girth(const igraph_t *graph, igraph_integer_t *girth,
                 igraph_vector_t *circle) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_t q;
    igraph_lazy_adjlist_t adjlist;
    long int mincirc = LONG_MAX, minvertex = 0;
    long int node;
    igraph_bool_t triangle = 0;
    igraph_vector_t *neis;
    igraph_vector_long_t level;
    long int stoplevel = no_of_nodes + 1;
    igraph_bool_t anycircle = 0;
    long int t1 = 0, t2 = 0;

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, IGRAPH_ALL,
                                          IGRAPH_SIMPLIFY));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vector_long_init(&level, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &level);

    for (node = 0; !triangle && node < no_of_nodes; node++) {

        /* Are there circles in this graph at all? */
        if (node == 1 && anycircle == 0) {
            igraph_bool_t conn;
            IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_WEAK));
            if (conn) {
                /* No, there are none */
                break;
            }
        }

        anycircle = 0;
        igraph_dqueue_clear(&q);
        igraph_vector_long_null(&level);
        IGRAPH_CHECK(igraph_dqueue_push(&q, node));
        VECTOR(level)[node] = 1;

        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);
            long int actlevel = VECTOR(level)[actnode];
            long int i, n;

            if (actlevel >= stoplevel) {
                break;
            }

            neis = igraph_lazy_adjlist_get(&adjlist, (igraph_integer_t) actnode);
            n = igraph_vector_size(neis);
            for (i = 0; i < n; i++) {
                long int nei = (long int) VECTOR(*neis)[i];
                long int neilevel = VECTOR(level)[nei];
                if (neilevel != 0) {
                    if (neilevel == actlevel - 1) {
                        continue;
                    } else {
                        /* found circle */
                        stoplevel = neilevel;
                        anycircle = 1;
                        if (actlevel < mincirc) {
                            /* Is it a minimum circle? */
                            mincirc = actlevel + neilevel - 1;
                            minvertex = node;
                            t1 = actnode; t2 = nei;
                            if (neilevel == 2) {
                                /* Is it a triangle? */
                                triangle = 1;
                            }
                        }
                        if (neilevel == actlevel) {
                            break;
                        }
                    }
                } else {
                    igraph_dqueue_push(&q, nei);
                    VECTOR(level)[nei] = actlevel + 1;
                }
            }

        } /* while q !empty */
    } /* node */

    if (girth) {
        if (mincirc == LONG_MAX) {
            *girth = mincirc = 0;
        } else {
            *girth = (igraph_integer_t) mincirc;
        }
    }

    /* Store the actual circle, if needed */
    if (circle) {
        IGRAPH_CHECK(igraph_vector_resize(circle, mincirc));
        if (mincirc != 0) {
            long int i, n, idx = 0;
            igraph_dqueue_clear(&q);
            igraph_vector_long_null(&level); /* used for father pointers */
#define FATHER(x) (VECTOR(level)[(x)])
            IGRAPH_CHECK(igraph_dqueue_push(&q, minvertex));
            FATHER(minvertex) = minvertex;
            while (FATHER(t1) == 0 || FATHER(t2) == 0) {
                long int actnode = (long int) igraph_dqueue_pop(&q);
                neis = igraph_lazy_adjlist_get(&adjlist, (igraph_integer_t) actnode);
                n = igraph_vector_size(neis);
                for (i = 0; i < n; i++) {
                    long int nei = (long int) VECTOR(*neis)[i];
                    if (FATHER(nei) == 0) {
                        FATHER(nei) = actnode + 1;
                        igraph_dqueue_push(&q, nei);
                    }
                }
            }  /* while q !empty */
            /* Ok, now use FATHER to create the path */
            while (t1 != minvertex) {
                VECTOR(*circle)[idx++] = t1;
                t1 = FATHER(t1) - 1;
            }
            VECTOR(*circle)[idx] = minvertex;
            idx = mincirc - 1;
            while (t2 != minvertex) {
                VECTOR(*circle)[idx--] = t2;
                t2 = FATHER(t2) - 1;
            }
        } /* anycircle */
    } /* circle */
#undef FATHER

    igraph_vector_long_destroy(&level);
    igraph_dqueue_destroy(&q);
    igraph_lazy_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

/**
 * \function igraph_add_edge
 * \brief Adds a single edge to a graph.
 *
 * </para><para>
 * For directed graphs the edge points from \p from to \p to.
 *
 * </para><para>
 * Note that if you want to add many edges to a big graph, then it is
 * inefficient to add them one by one, it is better to collect them into
 * a vector and add all of them via a single \ref igraph_add_edges() call.
 * \param igraph The graph.
 * \param from The id of the first vertex of the edge.
 * \param to The id of the second vertex of the edge.
 * \return Error code.
 *
 * \sa \ref igraph_add_edges() to add many edges, \ref
 * igraph_delete_edges() to remove edges and \ref
 * igraph_add_vertices() to add vertices.
 *
 * Time complexity: O(|V|+|E|), the number of edges plus the number of
 * vertices.
 */

int igraph_add_edge(igraph_t *graph, igraph_integer_t from, igraph_integer_t to) {

    igraph_vector_t edges;
    int ret;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 2);

    VECTOR(edges)[0] = from;
    VECTOR(edges)[1] = to;
    IGRAPH_CHECK(ret = igraph_add_edges(graph, &edges, 0));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return ret;
}

/*
 * \example examples/simple/graph_convergence_degree.c
 */

int igraph_convergence_degree(const igraph_t *graph, igraph_vector_t *result,
                              igraph_vector_t *ins, igraph_vector_t *outs) {
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    long int i, j, k, n;
    long int *geodist;
    igraph_vector_int_t *eids;
    igraph_vector_t *ins_p, *outs_p, ins_v, outs_v;
    igraph_dqueue_t q;
    igraph_inclist_t inclist;
    igraph_bool_t directed = igraph_is_directed(graph);

    if (result != 0) {
        IGRAPH_CHECK(igraph_vector_resize(result, no_of_edges));
    }
    IGRAPH_CHECK(igraph_dqueue_init(&q, 100));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &q);

    if (ins == 0) {
        ins_p = &ins_v;
        IGRAPH_VECTOR_INIT_FINALLY(ins_p, no_of_edges);
    } else {
        ins_p = ins;
        IGRAPH_CHECK(igraph_vector_resize(ins_p, no_of_edges));
        igraph_vector_null(ins_p);
    }

    if (outs == 0) {
        outs_p = &outs_v;
        IGRAPH_VECTOR_INIT_FINALLY(outs_p, no_of_edges);
    } else {
        outs_p = outs;
        IGRAPH_CHECK(igraph_vector_resize(outs_p, no_of_edges));
        igraph_vector_null(outs_p);
    }

    geodist = igraph_Calloc(no_of_nodes, long int);
    if (geodist == 0) {
        IGRAPH_ERROR("Cannot calculate convergence degrees", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, geodist);

    /* Collect shortest paths originating from/to every node to correctly
     * determine input field sizes */
    for (k = 0; k < (directed ? 2 : 1); k++) {
        igraph_neimode_t neimode = (k == 0) ? IGRAPH_OUT : IGRAPH_IN;
        igraph_real_t *vec;
        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, neimode));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
        vec = (k == 0) ? VECTOR(*ins_p) : VECTOR(*outs_p);
        for (i = 0; i < no_of_nodes; i++) {
            igraph_dqueue_clear(&q);
            memset(geodist, 0, sizeof(long int) * (size_t) no_of_nodes);
            geodist[i] = 1;
            IGRAPH_CHECK(igraph_dqueue_push(&q, i));
            IGRAPH_CHECK(igraph_dqueue_push(&q, 0.0));
            while (!igraph_dqueue_empty(&q)) {
                long int actnode = (long int) igraph_dqueue_pop(&q);
                long int actdist = (long int) igraph_dqueue_pop(&q);
                IGRAPH_ALLOW_INTERRUPTION();
                eids = igraph_inclist_get(&inclist, actnode);
                n = igraph_vector_int_size(eids);
                for (j = 0; j < n; j++) {
                    long int neighbor = IGRAPH_OTHER(graph, VECTOR(*eids)[j], actnode);
                    if (geodist[neighbor] != 0) {
                        /* we've already seen this node, another shortest path? */
                        if (geodist[neighbor] - 1 == actdist + 1) {
                            /* Since this edge is in the BFS tree rooted at i, we must
                             * increase either the size of the infield or the outfield */
                            if (!directed) {
                                if (actnode < neighbor) {
                                    VECTOR(*ins_p)[(long int)VECTOR(*eids)[j]] += 1;
                                } else {
                                    VECTOR(*outs_p)[(long int)VECTOR(*eids)[j]] += 1;
                                }
                            } else {
                                vec[(long int)VECTOR(*eids)[j]] += 1;
                            }
                        } else if (geodist[neighbor] - 1 < actdist + 1) {
                            continue;
                        }
                    } else {
                        /* we haven't seen this node yet */
                        IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
                        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
                        /* Since this edge is in the BFS tree rooted at i, we must
                         * increase either the size of the infield or the outfield */
                        if (!directed) {
                            if (actnode < neighbor) {
                                VECTOR(*ins_p)[(long int)VECTOR(*eids)[j]] += 1;
                            } else {
                                VECTOR(*outs_p)[(long int)VECTOR(*eids)[j]] += 1;
                            }
                        } else {
                            vec[(long int)VECTOR(*eids)[j]] += 1;
                        }
                        geodist[neighbor] = actdist + 2;
                    }
                }
            }
        }

        igraph_inclist_destroy(&inclist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (result != 0) {
        for (i = 0; i < no_of_edges; i++)
            VECTOR(*result)[i] = (VECTOR(*ins_p)[i] - VECTOR(*outs_p)[i]) /
                                 (VECTOR(*ins_p)[i] + VECTOR(*outs_p)[i]);
        if (!directed) {
            for (i = 0; i < no_of_edges; i++)
                if (VECTOR(*result)[i] < 0) {
                    VECTOR(*result)[i] = -VECTOR(*result)[i];
                }
        }
    }

    if (ins == 0) {
        igraph_vector_destroy(ins_p);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (outs == 0) {
        igraph_vector_destroy(outs_p);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_Free(geodist);
    igraph_dqueue_destroy(&q);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

/**
 * \function igraph_unfold_tree
 * Unfolding a graph into a tree, by possibly multiplicating its vertices.
 *
 * A graph is converted into a tree (or forest, if it is unconnected),
 * by performing a breadth-first search on it, and replicating
 * vertices that were found a second, third, etc. time.
 * \param graph The input graph, it can be either directed or
 *   undirected.
 * \param tree Pointer to an uninitialized graph object, the result is
 *   stored here.
 * \param mode For directed graphs; whether to follow paths along edge
 *    directions (\c IGRAPH_OUT), or the opposite (\c IGRAPH_IN), or
 *    ignore edge directions completely (\c IGRAPH_ALL). It is ignored
 *    for undirected graphs.
 * \param roots A numeric vector giving the root vertex, or vertices
 *   (if the graph is not connected), to start from.
 * \param vertex_index Pointer to an initialized vector, or a null
 *   pointer. If not a null pointer, then a mapping from the vertices
 *   in the new graph to the ones in the original is created here.
 * \return Error code.
 *
 * Time complexity: O(n+m), linear in the number vertices and edges.
 *
 */

int igraph_unfold_tree(const igraph_t *graph, igraph_t *tree,
                       igraph_neimode_t mode, const igraph_vector_t *roots,
                       igraph_vector_t *vertex_index) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    long int no_of_roots = igraph_vector_size(roots);
    long int tree_vertex_count = no_of_nodes;

    igraph_vector_t edges;
    igraph_vector_bool_t seen_vertices;
    igraph_vector_bool_t seen_edges;

    igraph_dqueue_t Q;
    igraph_vector_t neis;

    long int i, n, r, v_ptr = no_of_nodes;

    /* TODO: handle not-connected graphs, multiple root vertices */

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    igraph_vector_reserve(&edges, no_of_edges * 2);
    IGRAPH_DQUEUE_INIT_FINALLY(&Q, 100);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&seen_vertices, no_of_nodes);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&seen_edges, no_of_edges);

    if (vertex_index) {
        IGRAPH_CHECK(igraph_vector_resize(vertex_index, no_of_nodes));
        for (i = 0; i < no_of_nodes; i++) {
            VECTOR(*vertex_index)[i] = i;
        }
    }

    for (r = 0; r < no_of_roots; r++) {

        long int root = (long int) VECTOR(*roots)[r];
        VECTOR(seen_vertices)[root] = 1;
        igraph_dqueue_push(&Q, root);

        while (!igraph_dqueue_empty(&Q)) {
            long int actnode = (long int) igraph_dqueue_pop(&Q);

            IGRAPH_CHECK(igraph_incident(graph, &neis, (igraph_integer_t) actnode, mode));
            n = igraph_vector_size(&neis);
            for (i = 0; i < n; i++) {

                long int edge = (long int) VECTOR(neis)[i];
                long int from = IGRAPH_FROM(graph, edge);
                long int to = IGRAPH_TO(graph, edge);
                long int nei = IGRAPH_OTHER(graph, edge, actnode);

                if (! VECTOR(seen_edges)[edge]) {

                    VECTOR(seen_edges)[edge] = 1;

                    if (! VECTOR(seen_vertices)[nei]) {

                        igraph_vector_push_back(&edges, from);
                        igraph_vector_push_back(&edges, to);

                        VECTOR(seen_vertices)[nei] = 1;
                        IGRAPH_CHECK(igraph_dqueue_push(&Q, nei));

                    } else {

                        tree_vertex_count++;
                        if (vertex_index) {
                            IGRAPH_CHECK(igraph_vector_push_back(vertex_index, nei));
                        }

                        if (from == nei) {
                            igraph_vector_push_back(&edges, v_ptr++);
                            igraph_vector_push_back(&edges, to);
                        } else {
                            igraph_vector_push_back(&edges, from);
                            igraph_vector_push_back(&edges, v_ptr++);
                        }
                    }
                }

            } /* for i<n */

        } /* ! igraph_dqueue_empty(&Q) */

    } /* r < igraph_vector_size(roots) */

    igraph_vector_bool_destroy(&seen_edges);
    igraph_vector_bool_destroy(&seen_vertices);
    igraph_vector_destroy(&neis);
    igraph_dqueue_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(4);

    IGRAPH_CHECK(igraph_create(tree, &edges, tree_vertex_count, igraph_is_directed(graph)));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_diversity
 * Structural diversity index of the vertices
 *
 * This measure was defined in Nathan Eagle, Michael Macy and Rob
 * Claxton: Network Diversity and Economic Development, Science 328,
 * 1029--1031, 2010.
 *
 * </para><para>
 * It is simply the (normalized) Shannon entropy of the
 * incident edges' weights. D(i)=H(i)/log(k[i]), and
 * H(i) = -sum(p[i,j] log(p[i,j]), j=1..k[i]),
 * where p[i,j]=w[i,j]/sum(w[i,l], l=1..k[i]),  k[i] is the (total)
 * degree of vertex i, and w[i,j] is the weight of the edge(s) between
 * vertex i and j.
 * \param graph The input graph, edge directions are ignored.
 * \param weights The edge weights, in the order of the edge ids, must
 *    have appropriate length.
 * \param res An initialized vector, the results are stored here.
 * \param vids Vector with the vertex ids for which to calculate the
 *    measure.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear.
 *
 */

int igraph_diversity(igraph_t *graph, const igraph_vector_t *weights,
                     igraph_vector_t *res, const igraph_vs_t vids) {

    int no_of_nodes = igraph_vcount(graph);
    int no_of_edges = igraph_ecount(graph);
    igraph_vector_t incident;
    igraph_vit_t vit;
    igraph_real_t s, ent, w;
    int i, j, k;

    if (!weights) {
        IGRAPH_ERROR("Edge weights must be given", IGRAPH_EINVAL);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid edge weight vector length", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&incident, 10);

    if (igraph_vs_is_all(&vids)) {
        IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
        for (i = 0; i < no_of_nodes; i++) {
            s = ent = 0.0;
            IGRAPH_CHECK(igraph_incident(graph, &incident, i, /*mode=*/ IGRAPH_ALL));
            for (j = 0, k = (int) igraph_vector_size(&incident); j < k; j++) {
                w = VECTOR(*weights)[(long int)VECTOR(incident)[j]];
                s += w;
                ent += (w * log(w));
            }
            VECTOR(*res)[i] = (log(s) - ent / s) / log(k);
        }
    } else {
        IGRAPH_CHECK(igraph_vector_resize(res, 0));
        IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
        IGRAPH_FINALLY(igraph_vit_destroy, &vit);

        for (IGRAPH_VIT_RESET(vit), i = 0;
             !IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), i++) {
            long int v = IGRAPH_VIT_GET(vit);
            s = ent = 0.0;
            IGRAPH_CHECK(igraph_incident(graph, &incident, (igraph_integer_t) v,
                                         /*mode=*/ IGRAPH_ALL));
            for (j = 0, k = (int) igraph_vector_size(&incident); j < k; j++) {
                w = VECTOR(*weights)[(long int)VECTOR(incident)[j]];
                s += w;
                ent += (w * log(w));
            }
            IGRAPH_CHECK(igraph_vector_push_back(res, (log(s) - ent / s) / log(k)));
        }

        igraph_vit_destroy(&vit);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_destroy(&incident);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}


/* igraph_is_tree -- check if a graph is a tree */

/* count the number of vertices reachable from the root */
static int igraph_i_is_tree_visitor(igraph_integer_t root, const igraph_adjlist_t *al, igraph_integer_t *visited_count) {
    igraph_stack_int_t stack;
    igraph_vector_bool_t visited;
    long i;

    IGRAPH_CHECK(igraph_vector_bool_init(&visited, igraph_adjlist_size(al)));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &visited);

    IGRAPH_CHECK(igraph_stack_int_init(&stack, 0));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &stack);

    *visited_count = 0;

    /* push the root into the stack */
    IGRAPH_CHECK(igraph_stack_int_push(&stack, root));

    while (! igraph_stack_int_empty(&stack)) {
        igraph_integer_t u;
        igraph_vector_int_t *neighbors;
        long ncount;

        /* take a vertex from the stack, mark it as visited */
        u = igraph_stack_int_pop(&stack);
        if (IGRAPH_LIKELY(! VECTOR(visited)[u])) {
            VECTOR(visited)[u] = 1;
            *visited_count += 1;
        }

        /* register all its yet-unvisited neighbours for future processing */
        neighbors = igraph_adjlist_get(al, u);
        ncount = igraph_vector_int_size(neighbors);
        for (i = 0; i < ncount; ++i) {
            igraph_integer_t v = VECTOR(*neighbors)[i];
            if (! VECTOR(visited)[v]) {
                IGRAPH_CHECK(igraph_stack_int_push(&stack, v));
            }
        }
    }

    igraph_stack_int_destroy(&stack);
    igraph_vector_bool_destroy(&visited);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup structural
 * \function igraph_is_tree
 * \brief Decides whether the graph is a tree.
 *
 * An undirected graph is a tree if it is connected and has no cycles.
 * </para><para>
 *
 * In the directed case, a possible additional requirement is that all
 * edges are oriented away from a root (out-tree or arborescence) or all edges
 * are oriented towards a root (in-tree or anti-arborescence).
 * This test can be controlled using the \p mode parameter.
 * </para><para>
 *
 * By convention, the null graph (i.e. the graph with no vertices) is considered not to be a tree.
 *
 * \param graph The graph object to analyze.
 * \param res Pointer to a logical variable, the result will be stored
 *        here.
 * \param root If not \c NULL, the root node will be stored here. When \p mode
 *        is \c IGRAPH_ALL or the graph is undirected, any vertex can be the root
 *        and \p root is set to 0 (the first vertex). When \p mode is \c IGRAPH_OUT
 *        or \c IGRAPH_IN, the root is set to the vertex with zero in- or out-degree,
 *        respectively.
 * \param mode For a directed graph this specifies whether to test for an
 *        out-tree, an in-tree or ignore edge directions. The respective
 *        possible values are:
 *        \c IGRAPH_OUT, \c IGRAPH_IN, \c IGRAPH_ALL. This argument is
 *        ignored for undirected graphs.
 * \return Error code:
 *        \c IGRAPH_EINVAL: invalid mode argument.
 *
 * Time complexity: At most O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa igraph_is_weakly_connected()
 *
 * \example examples/simple/igraph_tree.c
 */

int igraph_is_tree(const igraph_t *graph, igraph_bool_t *res, igraph_integer_t *root, igraph_neimode_t mode) {
    igraph_adjlist_t al;
    igraph_integer_t iroot = 0;
    igraph_integer_t visited_count;
    igraph_integer_t vcount, ecount;

    vcount = igraph_vcount(graph);
    ecount = igraph_ecount(graph);

    /* A tree must have precisely vcount-1 edges. */
    /* By convention, the zero-vertex graph will not be considered a tree. */
    if (ecount != vcount - 1) {
        *res = 0;
        return IGRAPH_SUCCESS;
    }

    /* The single-vertex graph is a tree, provided it has no edges (checked in the previous if (..)) */
    if (vcount == 1) {
        *res = 1;
        if (root) {
            *root = 0;
        }
        return IGRAPH_SUCCESS;
    }

    /* For higher vertex counts we cannot short-circuit due to the possibility
     * of loops or multi-edges even when the edge count is correct. */

    /* Ignore mode for undirected graphs. */
    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    IGRAPH_CHECK(igraph_adjlist_init(graph, &al, mode));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);

    /* The main algorithm:
     * We find a root and check that all other vertices are reachable from it.
     * We have already checked the number of edges, so with the additional
     * reachability condition we can verify if the graph is a tree.
     *
     * For directed graphs, the root is the node with no incoming/outgoing
     * connections, depending on 'mode'. For undirected, it is arbitrary, so
     * we choose 0.
     */

    *res = 1; /* assume success */

    switch (mode) {
    case IGRAPH_ALL:
        iroot = 0;
        break;

    case IGRAPH_IN:
    case IGRAPH_OUT: {
        igraph_vector_t degree;
        igraph_integer_t i;

        IGRAPH_CHECK(igraph_vector_init(&degree, 0));
        IGRAPH_FINALLY(igraph_vector_destroy, &degree);

        IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), mode == IGRAPH_IN ? IGRAPH_OUT : IGRAPH_IN, /* loops = */ 1));

        for (i = 0; i < vcount; ++i)
            if (VECTOR(degree)[i] == 0) {
                break;
            }

        /* if no suitable root is found, the graph is not a tree */
        if (i == vcount) {
            *res = 0;
        } else {
            iroot = i;
        }

        igraph_vector_destroy(&degree);
        IGRAPH_FINALLY_CLEAN(1);
    }

    break;
    default:
        IGRAPH_ERROR("Invalid mode", IGRAPH_EINVMODE);
    }

    /* if no suitable root was found, skip visiting vertices */
    if (*res) {
        IGRAPH_CHECK(igraph_i_is_tree_visitor(iroot, &al, &visited_count));
        *res = visited_count == vcount;
    }

    if (root) {
        *root = iroot;
    }

    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
