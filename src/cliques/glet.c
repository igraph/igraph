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

#include "igraph_graphlets.h"

#include "igraph_conversion.h"
#include "igraph_constructors.h"
#include "igraph_cliques.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_operators.h"
#include "igraph_qsort.h"
#include "igraph_structural.h"

/**
 * \section graphlets_intro Introduction
 *
 * <para>
 * Graphlet decomposition models a weighted undirected graph
 * via the union of potentially overlapping dense social groups.
 * This is done by a two-step algorithm. In the first step, a candidate
 * set of groups (a candidate basis) is created by finding cliques
 * in the thresholded input graph. In the second step,
 * the graph is projected onto the candidate basis, resulting in a
 * weight coefficient for each clique in the candidate basis.
 * </para>
 *
 * <para>
 * For more information on graphlet decomposition, see
 * Hossein Azari Soufiani and Edoardo M Airoldi: "Graphlet decomposition of a weighted network",
 * https://arxiv.org/abs/1203.2821 and http://proceedings.mlr.press/v22/azari12/azari12.pdf
 * </para>
 *
 * <para>
 * igraph contains three functions for performing the graphlet
 * decomponsition of a graph. The first is \ref igraph_graphlets(), which
 * performs both steps of the method and returns a list of subgraphs
 * with their corresponding weights. The other two functions
 * correspond to the first and second steps of the algorithm, and they are
 * useful if the user wishes to perform them individually:
 * \ref igraph_graphlets_candidate_basis() and
 * \ref igraph_graphlets_project().
 * </para>
 *
 * <para>
 * <remark>
 * Note: The term "graphlet" is used for several unrelated concepts
 * in the literature. If you are looking to count induced subgraphs, see
 * \ref igraph_motifs_randesu() and \ref igraph_subisomorphic_lad().
 * </remark>
 * </para>
 */

typedef struct {
    igraph_vector_int_t *resultids;
    igraph_t *result;
    igraph_vector_t *resultweights;
    igraph_integer_t nc;
} igraph_i_subclique_next_free_t;

static void igraph_i_subclique_next_free(void *ptr) {
    igraph_i_subclique_next_free_t *data = ptr;
    igraph_integer_t i;
    if (data->resultids) {
        for (i = 0; i < data->nc; i++) {
            igraph_vector_int_destroy(&data->resultids[i]);
        }
        IGRAPH_FREE(data->resultids);
    }
    if (data->result) {
        for (i = 0; i < data->nc; i++) {
            igraph_destroy(&data->result[i]);
        }
        IGRAPH_FREE(data->result);
    }
    if (data->resultweights) {
        for (i = 0; i < data->nc; i++) {
            igraph_vector_destroy(&data->resultweights[i]);
        }
        IGRAPH_FREE(data->resultweights);
    }
}

/**
 * \function igraph_i_subclique_next
 * Calculate subcliques of the cliques found at the previous level
 *
 * \param graph Input graph.
 * \param weight Edge weights.
 * \param ids The IDs of the vertices in the input graph.
 * \param cliques A list of \ref igraph_vector_int_t, vertex IDs for cliques.
 * \param result The result is stored here, a list of graphs is stored
 *        here.
 * \param resultids The IDs of the vertices in the result graphs is
 *        stored here.
 * \param clique_thr The thresholds for the cliques are stored here,
 *        if not a null pointer.
 * \param next_thr The next thresholds for the cliques are stored
 *        here, if not a null pointer.
 *
 */

static igraph_error_t igraph_i_subclique_next(const igraph_t *graph,
                                   const igraph_vector_t *weights,
                                   const igraph_vector_int_t *ids,
                                   const igraph_vector_int_list_t *cliques,
                                   igraph_t **result,
                                   igraph_vector_t **resultweights,
                                   igraph_vector_int_t **resultids,
                                   igraph_vector_t *clique_thr,
                                   igraph_vector_t *next_thr) {

    /* The input is a set of cliques, that were found at a previous level.
       For each clique, we calculate the next threshold, drop the isolate
       vertices, and create a new graph from them. */

    igraph_vector_int_t mark, map;
    igraph_vector_int_t edges;
    igraph_vector_int_t neis, newedges;
    igraph_integer_t c, nc = igraph_vector_int_list_size(cliques);
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_i_subclique_next_free_t freedata = { NULL, NULL, NULL, nc };

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid length of weight vector", IGRAPH_EINVAL);
    }

    if (igraph_vector_int_size(ids) != no_of_nodes) {
        IGRAPH_ERROR("Invalid length of ID vector", IGRAPH_EINVAL);
    }

    IGRAPH_FINALLY(igraph_i_subclique_next_free, &freedata);

    *resultids = IGRAPH_CALLOC(nc, igraph_vector_int_t);
    IGRAPH_CHECK_OOM(*resultids, "Cannot calculate next cliques.");
    freedata.resultids = *resultids;

    *resultweights = IGRAPH_CALLOC(nc, igraph_vector_t);
    IGRAPH_CHECK_OOM(*resultweights, "Cannot calculate next cliques.");
    freedata.resultweights = *resultweights;

    *result = IGRAPH_CALLOC(nc, igraph_t);
    IGRAPH_CHECK_OOM(*result, "Cannot calculate next cliques.");
    freedata.result = *result;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&newedges, 100);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&mark, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&map, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 100);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 10);

    if (clique_thr) {
        IGRAPH_CHECK(igraph_vector_resize(clique_thr, nc));
    }
    if (next_thr) {
        IGRAPH_CHECK(igraph_vector_resize(next_thr, nc));
    }

    /* Iterate over all cliques. We will create graphs for all
       subgraphs defined by the cliques. */

    for (c = 0; c < nc; c++) {
        igraph_vector_int_t *clique = igraph_vector_int_list_get_ptr(cliques, c);
        igraph_real_t minweight = IGRAPH_INFINITY, nextweight = IGRAPH_INFINITY;
        igraph_integer_t e, v, clsize = igraph_vector_int_size(clique);
        igraph_integer_t noe, nov = 0;
        igraph_vector_int_t *newids = (*resultids) + c;
        igraph_vector_t *neww = (*resultweights) + c;
        igraph_t *newgraph = (*result) + c;

        igraph_vector_int_clear(&edges);
        igraph_vector_int_clear(&newedges);

        /* --------------------------------------------------- */

        /* Iterate over the vertices of a clique and find the
           edges within the clique, put them in a list.
           At the same time, search for the minimum edge weight within
           the clique and the next edge weight if any. */

        for (v = 0; v < clsize; v++) {
            igraph_integer_t i, neilen, node = VECTOR(*clique)[v];
            IGRAPH_CHECK(igraph_incident(graph, &neis, node, IGRAPH_ALL));
            neilen = igraph_vector_int_size(&neis);
            VECTOR(mark)[node] = c + 1;
            for (i = 0; i < neilen; i++) {
                igraph_integer_t edge = VECTOR(neis)[i];
                igraph_integer_t nei = IGRAPH_OTHER(graph, edge, node);
                if (VECTOR(mark)[nei] == c + 1) {
                    igraph_real_t w = VECTOR(*weights)[edge];
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, edge));
                    if (w < minweight) {
                        nextweight = minweight;
                        minweight = w;
                    } else if (w > minweight && w < nextweight) {
                        nextweight = w;
                    }
                }
            }
        } /* v < clsize */

        /* --------------------------------------------------- */

        /* OK, we have stored the edges and found the weight of
           the clique and the next weight to consider */

        if (clique_thr) {
            VECTOR(*clique_thr)[c] = minweight;
        }
        if (next_thr)   {
            VECTOR(*next_thr  )[c] = nextweight;
        }

        /* --------------------------------------------------- */

        /* Now we create the subgraph from the edges above the next
           threshold, and their incident vertices. */

        IGRAPH_CHECK(igraph_vector_int_init(newids, 0));
        IGRAPH_CHECK(igraph_vector_init(neww, 0));

        /* We use mark[] to denote the vertices already mapped to
           the new graph. If this is -(c+1), then the vertex was
           mapped, otherwise it was not. The mapping itself is in
           map[]. */

        noe = igraph_vector_int_size(&edges);
        for (e = 0; e < noe; e++) {
            igraph_integer_t edge = VECTOR(edges)[e];
            igraph_integer_t from, to;
            igraph_real_t w = VECTOR(*weights)[edge];
            IGRAPH_CHECK(igraph_edge(graph, edge, &from, &to));
            if (w >= nextweight) {
                if (VECTOR(mark)[from] == c + 1) {
                    VECTOR(map)[from] = nov++;
                    VECTOR(mark)[from] = -(c + 1);
                    IGRAPH_CHECK(igraph_vector_int_push_back(newids, VECTOR(*ids)[from]));
                }
                if (VECTOR(mark)[to] == c + 1) {
                    VECTOR(map)[to] = nov++;
                    VECTOR(mark)[to] = -(c + 1);
                    IGRAPH_CHECK(igraph_vector_int_push_back(newids, VECTOR(*ids)[to]));
                }
                IGRAPH_CHECK(igraph_vector_push_back(neww, w));
                IGRAPH_CHECK(igraph_vector_int_push_back(&newedges, VECTOR(map)[from]));
                IGRAPH_CHECK(igraph_vector_int_push_back(&newedges, VECTOR(map)[to]));
            }
        }

        IGRAPH_CHECK(igraph_create(newgraph, &newedges, nov, IGRAPH_UNDIRECTED));

        /* --------------------------------------------------- */

    } /* c < nc */

    igraph_vector_int_destroy(&neis);
    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&mark);
    igraph_vector_int_destroy(&map);
    igraph_vector_int_destroy(&newedges);
    IGRAPH_FINALLY_CLEAN(6);  /* + freedata */

    return IGRAPH_SUCCESS;
}

static void igraph_i_graphlets_destroy_clique_list(igraph_vector_ptr_t *vl) {
    igraph_integer_t i, n = igraph_vector_ptr_size(vl);
    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = (igraph_vector_int_t*) VECTOR(*vl)[i];
        if (v) {
            igraph_vector_int_destroy(v);
            IGRAPH_FREE(v);
        }
    }
    igraph_vector_ptr_destroy(vl);
}

static igraph_error_t igraph_i_graphlets(const igraph_t *graph,
                              const igraph_vector_t *weights,
                              igraph_vector_ptr_t *cliques,
                              igraph_vector_t *thresholds,
                              const igraph_vector_int_t *ids,
                              igraph_real_t startthr) {

    /* This version is different from the main function, and is
       appropriate to use in recursive calls, because it _adds_ the
       results to 'cliques' and 'thresholds' and uses the supplied
       'startthr' */

    igraph_vector_int_list_t mycliques;
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vector_int_t subv;
    igraph_t subg;
    igraph_t *newgraphs = NULL;
    igraph_vector_t *newweights = NULL;
    igraph_vector_int_t *newids = NULL;
    igraph_vector_t clique_thr, next_thr;
    igraph_i_subclique_next_free_t freedata = { NULL, NULL, NULL, 0 };

    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&mycliques, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&subv, 0);

    /* We start by finding cliques at the lowest threshold */
    for (igraph_integer_t i = 0; i < no_of_edges; i++) {
        if (VECTOR(*weights)[i] >= startthr) {
            IGRAPH_CHECK(igraph_vector_int_push_back(&subv, i));
        }
    }
    IGRAPH_CHECK(igraph_subgraph_from_edges(graph, &subg, igraph_ess_vector(&subv), /*delete_vertices=*/ 0));
    IGRAPH_FINALLY(igraph_destroy, &subg);
    IGRAPH_CHECK(igraph_maximal_cliques(&subg, &mycliques, /*min_size=*/ 0, /*max_size=*/ 0));
    igraph_destroy(&subg);
    igraph_vector_int_destroy(&subv);
    IGRAPH_FINALLY_CLEAN(2);

    const igraph_integer_t nocliques = igraph_vector_int_list_size(&mycliques);

    /* Get the next cliques and thresholds */
    IGRAPH_VECTOR_INIT_FINALLY(&next_thr, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&clique_thr, 0);

    IGRAPH_CHECK(igraph_i_subclique_next(
        graph, weights, ids, &mycliques, &newgraphs, &newweights, &newids,
        &clique_thr, &next_thr
    ));

    freedata.result = newgraphs;
    freedata.resultids = newids;
    freedata.resultweights = newweights;
    freedata.nc = nocliques;
    IGRAPH_FINALLY(igraph_i_subclique_next_free, &freedata);

    /* Store cliques at the current level */

    IGRAPH_CHECK(igraph_vector_append(thresholds, &clique_thr));

    igraph_vector_destroy(&clique_thr);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_vector_ptr_resize(cliques, igraph_vector_ptr_size(cliques) + nocliques));
    for (igraph_integer_t i = 0, j = igraph_vector_ptr_size(cliques) - 1; i < nocliques; i++, j--) {
        igraph_vector_int_t *cl = IGRAPH_CALLOC(1, igraph_vector_int_t);
        IGRAPH_CHECK_OOM(cl, "Cannot find graphlets.");
        IGRAPH_FINALLY(igraph_free, cl);

        *cl = igraph_vector_int_list_pop_back(&mycliques);

        /* From this point onwards, _we_ own the clique and not `mycliques'.
         * We pass on the ownership to `cliques' */
        VECTOR(*cliques)[j] = cl;
        IGRAPH_FINALLY_CLEAN(1);

        const igraph_integer_t n = igraph_vector_int_size(cl);
        for (igraph_integer_t k = 0; k < n; k++) {
            igraph_integer_t node = VECTOR(*cl)[k];
            VECTOR(*cl)[k] = VECTOR(*ids)[node];
        }
        igraph_vector_int_sort(cl);
    }

    igraph_vector_int_list_destroy(&mycliques); /* contents was copied over to `cliques' */
    IGRAPH_FINALLY_CLEAN(1);

    /* Recursive calls for cliques found */
    for (igraph_integer_t i = 0; i < nocliques; i++) {
        igraph_t *g = newgraphs + i;
        if (igraph_vcount(g) > 1) {
            igraph_vector_t *w_sub = newweights + i;
            igraph_vector_int_t *ids_sub = newids + i;
            IGRAPH_CHECK(igraph_i_graphlets(g, w_sub, cliques, thresholds, ids_sub, VECTOR(next_thr)[i]));
        }
    }

    igraph_vector_destroy(&next_thr);
    igraph_i_subclique_next_free(&freedata);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

typedef struct {
    const igraph_vector_ptr_t *cliques;
    const igraph_vector_t *thresholds;
} igraph_i_graphlets_filter_t;

static int igraph_i_graphlets_filter_cmp(void *data, const void *a, const void *b) {
    igraph_i_graphlets_filter_t *ddata = (igraph_i_graphlets_filter_t *) data;
    igraph_integer_t *aa = (igraph_integer_t*) a;
    igraph_integer_t *bb = (igraph_integer_t*) b;
    igraph_real_t t_a = VECTOR(*ddata->thresholds)[*aa];
    igraph_real_t t_b = VECTOR(*ddata->thresholds)[*bb];
    igraph_vector_int_t *v_a, *v_b;
    igraph_integer_t s_a, s_b;

    if (t_a < t_b) {
        return -1;
    } else if (t_a > t_b) {
        return 1;
    }

    v_a = (igraph_vector_int_t*) VECTOR(*ddata->cliques)[*aa];
    v_b = (igraph_vector_int_t*) VECTOR(*ddata->cliques)[*bb];
    s_a = igraph_vector_int_size(v_a);
    s_b = igraph_vector_int_size(v_b);

    if (s_a < s_b) {
        return -1;
    } else if (s_a > s_b) {
        return 1;
    } else {
        return 0;
    }
}

static igraph_error_t igraph_i_graphlets_filter(igraph_vector_ptr_t *cliques,
                                     igraph_vector_t *thresholds) {

    /* Filter out non-maximal cliques. Every non-maximal clique is
       part of a maximal clique, at the same threshold.

       First we order the cliques, according to their threshold, and
       then according to their size. So when we look for a candidate
       superset, we only need to check the cliques next in the list,
       until their threshold is different. */

    igraph_integer_t i, iptr, nocliques = igraph_vector_ptr_size(cliques);
    igraph_vector_int_t order;
    igraph_i_graphlets_filter_t sortdata = { cliques, thresholds };

    IGRAPH_CHECK(igraph_vector_int_init_range(&order, 0, nocliques));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &order);

    igraph_qsort_r(VECTOR(order), nocliques, sizeof(VECTOR(order)[0]), &sortdata,
                   igraph_i_graphlets_filter_cmp);

    for (i = 0; i < nocliques - 1; i++) {
        igraph_integer_t ri = VECTOR(order)[i];
        igraph_vector_int_t *needle = VECTOR(*cliques)[ri];
        igraph_real_t thr_i = VECTOR(*thresholds)[ri];
        igraph_integer_t n_i = igraph_vector_int_size(needle);

        for (igraph_integer_t j = i + 1; j < nocliques; j++) {
            igraph_integer_t rj = VECTOR(order)[j];
            igraph_real_t thr_j = VECTOR(*thresholds)[rj];
            igraph_vector_int_t *hay;
            igraph_integer_t n_j, pi = 0, pj = 0;

            /* Done, not found */
            if (thr_j != thr_i) {
                break;
            }

            /* Check size of hay */
            hay = VECTOR(*cliques)[rj];
            n_j = igraph_vector_int_size(hay);
            if (n_i > n_j) {
                continue;
            }

            /* Check if hay is a superset */
            while (pi < n_i && pj < n_j && n_i - pi <= n_j - pj) {
                igraph_integer_t ei = VECTOR(*needle)[pi];
                igraph_integer_t ej = VECTOR(*hay)[pj];
                if (ei < ej) {
                    break;
                } else if (ei > ej) {
                    pj++;
                } else {
                    pi++; pj++;
                }
            }
            if (pi == n_i) {
                /* Found, delete immediately */
                igraph_vector_int_destroy(needle);
                igraph_free(needle);
                VECTOR(*cliques)[ri] = 0;
                break;
            }
        }
    }

    /* Remove null pointers from the list of cliques */
    for (i = 0, iptr = 0; i < nocliques; i++) {
        igraph_vector_int_t *v = VECTOR(*cliques)[i];
        if (v) {
            VECTOR(*cliques)[iptr] = v;
            VECTOR(*thresholds)[iptr] = VECTOR(*thresholds)[i];
            iptr++;
        }
    }
    IGRAPH_CHECK(igraph_vector_ptr_resize(cliques, iptr));
    IGRAPH_CHECK(igraph_vector_resize(thresholds, iptr));

    igraph_vector_int_destroy(&order);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_graphlets_candidate_basis
 * Calculate a candidate graphlets basis
 *
 * \param graph The input graph, it must be a simple graph, edge directions are
 *        ignored.
 * \param weights Weights of the edges, a vector.
 * \param cliques An initialized list of integer vectors. The graphlet basis is
 *        stored here. Each element of the list is an integer vector of
 *        vertex IDs, encoding a single basis subgraph.
 * \param thresholds An initialized vector, the (highest possible)
 *        weight thresholds for finding the basis subgraphs are stored
 *        here.
 * \return Error code.
 *
 * See also: \ref igraph_graphlets() and \ref igraph_graphlets_project().
 */

igraph_error_t igraph_graphlets_candidate_basis(const igraph_t *graph,
                                     const igraph_vector_t *weights,
                                     igraph_vector_int_list_t *cliques,
                                     igraph_vector_t *thresholds) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_real_t minthr;
    igraph_vector_int_t ids;
    igraph_bool_t simple;
    igraph_integer_t i, no_of_cliques;
    igraph_vector_ptr_t mycliques;

    /* Some checks */
    if (weights == NULL) {
        IGRAPH_ERROR("Graphlet functions require weighted graphs", IGRAPH_EINVAL);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_is_simple(graph, &simple));
    if (!simple) {
        IGRAPH_ERROR("Graphlets work on simple graphs only", IGRAPH_EINVAL);
    }
    if (igraph_is_directed(graph)) {
        /* When the graph is directed, mutual edges are effectively multi-edges as we
         * are ignoring edge directions. */
        igraph_bool_t has_mutual;
        IGRAPH_CHECK(igraph_has_mutual(graph, &has_mutual, false));
        if (has_mutual) {
            IGRAPH_ERROR("Graphlets work on simple graphs only", IGRAPH_EINVAL);
        }
    }

    /* Internally, we will still use igraph_vector_ptr_t instead of
     * igraph_vector_int_list_t to manage the list of cliques; this is because
     * we are going to append & filter the list and it's more complicated to
     * do with an igraph_vector_int_list_t */
    IGRAPH_CHECK(igraph_vector_ptr_init(&mycliques, 0));
    IGRAPH_FINALLY(igraph_i_graphlets_destroy_clique_list, &mycliques);

    igraph_vector_int_list_clear(cliques);
    igraph_vector_clear(thresholds);

    minthr = igraph_vector_min(weights);

    IGRAPH_CHECK(igraph_vector_int_init_range(&ids, 0, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &ids);

    IGRAPH_CHECK(igraph_i_graphlets(graph, weights, &mycliques, thresholds, &ids, minthr));

    igraph_vector_int_destroy(&ids);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_i_graphlets_filter(&mycliques, thresholds));

    /* Pass ownership of cliques in `mycliques' to `cliques' so the user does
     * not have to work with igraph_vector_ptr_t */
    no_of_cliques = igraph_vector_ptr_size(&mycliques);
    for (i = 0; i < no_of_cliques; i++) {
        IGRAPH_CHECK(igraph_vector_int_list_push_back(
            cliques, VECTOR(mycliques)[i]
        ));
        IGRAPH_FREE(VECTOR(mycliques)[i]);
    }

    /* `mycliques' is now empty so we can clear and destroy */
    igraph_vector_ptr_clear(&mycliques);
    igraph_vector_ptr_destroy(&mycliques);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* TODO: not made static because it is used by the R interface */
igraph_error_t igraph_i_graphlets_project(
    const igraph_t *graph, const igraph_vector_t *weights,
    const igraph_vector_int_list_t *cliques, igraph_vector_t *Mu, igraph_bool_t startMu,
    igraph_integer_t niter, igraph_integer_t vid1
) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_cliques = igraph_vector_int_list_size(cliques);
    igraph_vector_int_t vcl, vclidx, ecl, eclidx, cel, celidx;
    igraph_vector_int_t edgelist;
    igraph_vector_t newweights, normfact;
    igraph_integer_t i, total_vertices, e, ptr, total_edges;
    igraph_bool_t simple;

    /* Check arguments */
    if (weights == NULL) {
        IGRAPH_ERROR("Graphlet functions require weighted graphs", IGRAPH_EINVAL);
    }
    if (no_of_edges != igraph_vector_size(weights)) {
        IGRAPH_ERROR("Invalid weight vector size", IGRAPH_EINVAL);
    }
    if (startMu && igraph_vector_size(Mu) != no_cliques) {
        IGRAPH_ERROR("Invalid start coefficient vector size", IGRAPH_EINVAL);
    }
    if (niter < 0) {
        IGRAPH_ERROR("Number of iterations must be non-negative", IGRAPH_EINVAL);
    }
    IGRAPH_CHECK(igraph_is_simple(graph, &simple));
    if (!simple) {
        IGRAPH_ERROR("Graphlets work on simple graphs only", IGRAPH_EINVAL);
    }
    if (igraph_is_directed(graph)) {
        /* When the graph is directed, mutual edges are effectively multi-edges as we
         * are ignoring edge directions. */
        igraph_bool_t has_mutual;
        IGRAPH_CHECK(igraph_has_mutual(graph, &has_mutual, false));
        if (has_mutual) {
            IGRAPH_ERROR("Graphlets work on simple graphs only", IGRAPH_EINVAL);
        }
    }

    if (!startMu) {
        IGRAPH_CHECK(igraph_vector_resize(Mu, no_cliques));
        igraph_vector_fill(Mu, 1);
    }

    /* Count # cliques per vertex. Also, create an index
       for the edges per clique. */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vclidx, no_of_nodes + 2);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&celidx, no_cliques + 3);
    for (i = 0, total_vertices = 0, total_edges = 0; i < no_cliques; i++) {
        igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(cliques, i);
        igraph_integer_t j, n = igraph_vector_int_size(v);
        total_vertices += n;
        total_edges += n * (n - 1) / 2;
        VECTOR(celidx)[i + 2] = total_edges;
        for (j = 0; j < n; j++) {
            igraph_integer_t vv = VECTOR(*v)[j] - vid1;
            VECTOR(vclidx)[vv + 2] += 1;
        }
    }
    VECTOR(celidx)[i + 2] = total_edges;

    /* Finalize index vector */
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(vclidx)[i + 2] += VECTOR(vclidx)[i + 1];
    }

    /* Create vertex-clique list, the cliques for each vertex. */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vcl, total_vertices);
    for (i = 0; i < no_cliques; i++) {
        igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(cliques, i);
        igraph_integer_t j, n = igraph_vector_int_size(v);
        for (j = 0; j < n; j++) {
            igraph_integer_t vv = VECTOR(*v)[j] - vid1;
            igraph_integer_t p = VECTOR(vclidx)[vv + 1];
            VECTOR(vcl)[p] = i;
            VECTOR(vclidx)[vv + 1] += 1;
        }
    }

    /* Create an edge-clique list, the cliques of each edge */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&ecl, total_edges);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eclidx, no_of_edges + 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edgelist, no_of_edges * 2);
    IGRAPH_CHECK(igraph_get_edgelist(graph, &edgelist, /*by_col=*/ 0));
    for (i = 0, e = 0, ptr = 0; e < no_of_edges; e++) {
        igraph_integer_t from = VECTOR(edgelist)[i++];
        igraph_integer_t to = VECTOR(edgelist)[i++];
        igraph_integer_t from_s = VECTOR(vclidx)[from];
        igraph_integer_t from_e = VECTOR(vclidx)[from + 1];
        igraph_integer_t to_s = VECTOR(vclidx)[to];
        igraph_integer_t to_e = VECTOR(vclidx)[to + 1];
        VECTOR(eclidx)[e] = ptr;
        while (from_s < from_e && to_s < to_e) {
            igraph_integer_t from_v = VECTOR(vcl)[from_s];
            igraph_integer_t to_v = VECTOR(vcl)[to_s];
            if (from_v == to_v) {
                VECTOR(ecl)[ptr++] = from_v;
                from_s++; to_s++;
            } else if (from_v < to_v) {
                from_s++;
            } else {
                to_s++;
            }
        }
    }
    VECTOR(eclidx)[e] = ptr;

    igraph_vector_int_destroy(&edgelist);
    IGRAPH_FINALLY_CLEAN(1);

    /* Convert the edge-clique list to a clique-edge list */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&cel, total_edges);
    for (i = 0; i < no_of_edges; i++) {
        igraph_integer_t ecl_s = VECTOR(eclidx)[i], ecl_e = VECTOR(eclidx)[i + 1], j;
        for (j = ecl_s; j < ecl_e; j++) {
            igraph_integer_t cl = VECTOR(ecl)[j];
            igraph_integer_t epos = VECTOR(celidx)[cl + 1];
            VECTOR(cel)[epos] = i;
            VECTOR(celidx)[cl + 1] += 1;
        }
    }

    /* Normalizing factors for the iteration */
    IGRAPH_VECTOR_INIT_FINALLY(&normfact, no_cliques);
    for (i = 0; i < no_cliques; i++) {
        igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(cliques, i);
        igraph_integer_t n = igraph_vector_int_size(v);
        VECTOR(normfact)[i] = n * (n + 1) / 2;
    }

    /* We have the clique-edge list, so do the projection now */
    IGRAPH_VECTOR_INIT_FINALLY(&newweights, no_of_edges);
    for (i = 0; i < niter; i++) {
        for (e = 0; e < no_of_edges; e++) {
            igraph_integer_t start = VECTOR(eclidx)[e];
            igraph_integer_t end = VECTOR(eclidx)[e + 1];
            VECTOR(newweights)[e] = 0.0001;
            while (start < end) {
                igraph_integer_t clique = VECTOR(ecl)[start++];
                VECTOR(newweights)[e] += VECTOR(*Mu)[clique];
            }
        }
        for (e = 0; e < no_cliques; e++) {
            igraph_real_t sumratio = 0;
            igraph_integer_t start = VECTOR(celidx)[e];
            igraph_integer_t end = VECTOR(celidx)[e + 1];
            while (start < end) {
                igraph_integer_t edge = VECTOR(cel)[start++];
                sumratio += VECTOR(*weights)[edge] / VECTOR(newweights)[edge];
            }
            VECTOR(*Mu)[e] *= sumratio / VECTOR(normfact)[e];
        }
    }

    igraph_vector_destroy(&newweights);
    igraph_vector_destroy(&normfact);
    igraph_vector_int_destroy(&cel);
    igraph_vector_int_destroy(&eclidx);
    igraph_vector_int_destroy(&ecl);
    igraph_vector_int_destroy(&vcl);
    igraph_vector_int_destroy(&celidx);
    igraph_vector_int_destroy(&vclidx);
    IGRAPH_FINALLY_CLEAN(8);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_graphlets_project
 * Project a graph on a graphlets basis
 *
 * Note that the graph projected does not have to be the same that
 * was used to calculate the graphlet basis, but it is assumed that
 * it has the same number of vertices, and the vertex IDs of the two
 * graphs match.
 * \param graph The input graph, it must be a simple graph, edge directions are
 *        ignored.
 * \param weights Weights of the edges in the input graph, a vector.
 * \param cliques An initialized list of integer vectors. The graphlet basis is
 *        stored here. Each element of the list is an integer vector of
 *        vertex IDs, encoding a single basis subgraph.
 * \param Mu An initialized vector, the weights of the graphlets will
 *        be stored here. This vector is also used to initialize the
 *        the weight vector for the iterative algorithm, if the
 *        \c startMu argument is true.
 * \param startMu If true, then the supplied Mu vector is
 *        used as the starting point of the iteration. Otherwise a
 *        constant 1 vector is used.
 * \param niter Integer scalar, the number of iterations to perform.
 * \return Error code.
 *
 * See also: \ref igraph_graphlets() and
 * \ref igraph_graphlets_candidate_basis().
 */

igraph_error_t igraph_graphlets_project(const igraph_t *graph,
                             const igraph_vector_t *weights,
                             const igraph_vector_int_list_t *cliques,
                             igraph_vector_t *Mu, igraph_bool_t startMu,
                             igraph_integer_t niter) {

    return igraph_i_graphlets_project(graph, weights, cliques, Mu, startMu,
                                      niter, /*vid1=*/ 0);
}

typedef struct igraph_i_graphlets_order_t {
    const igraph_vector_int_list_t *cliques;
    const igraph_vector_t *Mu;
} igraph_i_graphlets_order_t;

static int igraph_i_graphlets_order_cmp(void *data, const void *a, const void *b) {
    igraph_i_graphlets_order_t *ddata = (igraph_i_graphlets_order_t*) data;
    igraph_integer_t *aa = (igraph_integer_t*) a;
    igraph_integer_t *bb = (igraph_integer_t*) b;
    igraph_real_t Mu_a = VECTOR(*ddata->Mu)[*aa];
    igraph_real_t Mu_b = VECTOR(*ddata->Mu)[*bb];

    if (Mu_a < Mu_b) {
        return 1;
    } else if (Mu_a > Mu_b) {
        return -1;
    } else {
        return 0;
    }
}

/**
 * \function igraph_graphlets
 * Calculate graphlets basis and project the graph on it
 *
 * This function simply calls \ref igraph_graphlets_candidate_basis()
 * and \ref igraph_graphlets_project(), and then orders the graphlets
 * according to decreasing weights.
 * \param graph The input graph, it must be a simple graph, edge directions are
 *        ignored.
 * \param weights Weights of the edges, a vector.
 * \param cliques An initialized list of integer vectors. The graphlet basis is
 *        stored here. Each element of the list is an integer vector of
 *        vertex IDs, encoding a single basis subgraph.
 * \param Mu An initialized vector, the weights of the graphlets will
 *        be stored here.
 * \param niter Integer scalar, the number of iterations to perform
 *        for the projection step.
 * \return Error code.
 *
 * See also: \ref igraph_graphlets_candidate_basis() and
 * \ref igraph_graphlets_project().
 */

igraph_error_t igraph_graphlets(const igraph_t *graph,
                     const igraph_vector_t *weights,
                     igraph_vector_int_list_t *cliques,
                     igraph_vector_t *Mu, igraph_integer_t niter) {

    igraph_integer_t nocliques;
    igraph_vector_t thresholds;
    igraph_vector_int_t order;
    igraph_i_graphlets_order_t sortdata = { cliques, Mu };

    IGRAPH_VECTOR_INIT_FINALLY(&thresholds, 0);
    IGRAPH_CHECK(igraph_graphlets_candidate_basis(graph, weights, cliques, &thresholds));
    igraph_vector_destroy(&thresholds);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_graphlets_project(graph, weights, cliques, Mu, /*startMu=*/ false, niter));

    nocliques = igraph_vector_int_list_size(cliques);
    IGRAPH_CHECK(igraph_vector_int_init_range(&order, 0, nocliques));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &order);

    igraph_qsort_r(VECTOR(order), nocliques, sizeof(VECTOR(order)[0]), &sortdata,
                   igraph_i_graphlets_order_cmp);

    IGRAPH_CHECK(igraph_vector_int_list_permute(cliques, &order));
    IGRAPH_CHECK(igraph_vector_index_int(Mu, &order));

    igraph_vector_int_destroy(&order);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
