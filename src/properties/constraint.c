/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2021 The igraph development team

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

#include "igraph_centrality.h"

#include "igraph_interface.h"

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
