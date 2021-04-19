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
    int nc;
} igraph_i_subclique_next_free_t;

static void igraph_i_subclique_next_free(void *ptr) {
    igraph_i_subclique_next_free_t *data = ptr;
    int i;
    if (data->resultids) {
        for (i = 0; i < data->nc; i++) {
            if (data->resultids + i) {
                igraph_vector_int_destroy(data->resultids + i);
            }
        }
        IGRAPH_FREE(data->resultids);
    }
    if (data->result) {
        for (i = 0; i < data->nc; i++) {
            if (data->result + i) {
                igraph_destroy(data->result + i);
            }
        }
        IGRAPH_FREE(data->result);
    }
    if (data->resultweights) {
        for (i = 0; i < data->nc; i++) {
            if (data->resultweights + i) {
                igraph_vector_destroy(data->resultweights + i);
            }
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
 * \param ids The ids of the vertices in the input graph.
 * \param cliques A list of vectors, vertex ids for cliques.
 * \param result The result is stored here, a list of graphs is stored
 *        here.
 * \param resultids The ids of the vertices in the result graphs is
 *        stored here.
 * \param clique_thr The thresholds for the cliques are stored here,
 *        if not a null pointer.
 * \param next_thr The next thresholds for the cliques are stored
 *        here, if not a null pointer.
 *
 */

static int igraph_i_subclique_next(const igraph_t *graph,
                                   const igraph_vector_t *weights,
                                   const igraph_vector_int_t *ids,
                                   const igraph_vector_ptr_t *cliques,
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
    igraph_vector_t neis, newedges;
    igraph_integer_t c, nc = igraph_vector_ptr_size(cliques);
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_i_subclique_next_free_t freedata = { 0, 0, 0, nc };

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid length of weight vector", IGRAPH_EINVAL);
    }

    if (igraph_vector_int_size(ids) != no_of_nodes) {
        IGRAPH_ERROR("Invalid length of ID vector", IGRAPH_EINVAL);
    }

    IGRAPH_FINALLY(igraph_i_subclique_next_free, &freedata);
    *resultids = IGRAPH_CALLOC(nc, igraph_vector_int_t);
    if (!*resultids) {
        IGRAPH_ERROR("Cannot calculate next cliques", IGRAPH_ENOMEM);
    }
    freedata.resultids = *resultids;
    *resultweights = IGRAPH_CALLOC(nc, igraph_vector_t);
    if (!*resultweights) {
        IGRAPH_ERROR("Cannot calculate next cliques", IGRAPH_ENOMEM);
    }
    freedata.resultweights = *resultweights;
    *result = IGRAPH_CALLOC(nc, igraph_t);
    if (!*result) {
        IGRAPH_ERROR("Cannot calculate next cliques", IGRAPH_ENOMEM);
    }
    freedata.result = *result;

    igraph_vector_init(&newedges, 100);
    IGRAPH_FINALLY(igraph_vector_destroy, &newedges);
    igraph_vector_int_init(&mark, no_of_nodes);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &mark);
    igraph_vector_int_init(&map, no_of_nodes);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &map);
    igraph_vector_int_init(&edges, 100);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &edges);
    igraph_vector_init(&neis, 10);
    IGRAPH_FINALLY(igraph_vector_destroy, &neis);

    if (clique_thr) {
        igraph_vector_resize(clique_thr, nc);
    }
    if (next_thr)   {
        igraph_vector_resize(next_thr,   nc);
    }

    /* Iterate over all cliques. We will create graphs for all
       subgraphs defined by the cliques. */

    for (c = 0; c < nc; c++) {
        igraph_vector_t *clique = VECTOR(*cliques)[c];
        igraph_real_t minweight = IGRAPH_INFINITY, nextweight = IGRAPH_INFINITY;
        igraph_integer_t e, v, clsize = igraph_vector_size(clique);
        igraph_integer_t noe, nov = 0;
        igraph_vector_int_t *newids = (*resultids) + c;
        igraph_vector_t *neww = (*resultweights) + c;
        igraph_t *newgraph = (*result) + c;
        igraph_vector_int_clear(&edges);
        igraph_vector_clear(&newedges);

        /* --------------------------------------------------- */

        /* Iterate over the vertices of a clique and find the
           edges within the clique, put them in a list.
           At the same time, search for the minimum edge weight within
           the clique and the next edge weight if any. */

        for (v = 0; v < clsize; v++) {
            igraph_integer_t i, neilen, node = VECTOR(*clique)[v];
            igraph_incident(graph, &neis, node, IGRAPH_ALL);
            neilen = igraph_vector_size(&neis);
            VECTOR(mark)[node] = c + 1;
            for (i = 0; i < neilen; i++) {
                igraph_integer_t edge = VECTOR(neis)[i];
                igraph_integer_t nei = IGRAPH_OTHER(graph, edge, node);
                if (VECTOR(mark)[nei] == c + 1) {
                    igraph_real_t w = VECTOR(*weights)[edge];
                    igraph_vector_int_push_back(&edges, edge);
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

        igraph_vector_int_init(newids, 0);
        igraph_vector_init(neww, 0);

        /* We use mark[] to denote the vertices already mapped to
           the new graph. If this is -(c+1), then the vertex was
           mapped, otherwise it was not. The mapping itself is in
           map[]. */

        noe = igraph_vector_int_size(&edges);
        for (e = 0; e < noe; e++) {
            igraph_integer_t edge = VECTOR(edges)[e];
            igraph_integer_t from, to;
            igraph_real_t w = VECTOR(*weights)[edge];
            igraph_edge(graph, edge, &from, &to);
            if (w >= nextweight) {
                if (VECTOR(mark)[from] == c + 1) {
                    VECTOR(map)[from] = nov++;
                    VECTOR(mark)[from] = -(c + 1);
                    igraph_vector_int_push_back(newids, VECTOR(*ids)[from]);
                }
                if (VECTOR(mark)[to] == c + 1) {
                    VECTOR(map)[to] = nov++;
                    VECTOR(mark)[to] = -(c + 1);
                    igraph_vector_int_push_back(newids, VECTOR(*ids)[to]);
                }
                igraph_vector_push_back(neww, w);
                igraph_vector_push_back(&newedges, VECTOR(map)[from]);
                igraph_vector_push_back(&newedges, VECTOR(map)[to]);
            }
        }

        igraph_create(newgraph, &newedges, nov, IGRAPH_UNDIRECTED);

        /* --------------------------------------------------- */

    } /* c < nc */

    igraph_vector_destroy(&neis);
    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&mark);
    igraph_vector_int_destroy(&map);
    igraph_vector_destroy(&newedges);
    IGRAPH_FINALLY_CLEAN(6);  /* + freedata */

    return 0;
}

static void igraph_i_graphlets_destroy_vectorlist(igraph_vector_ptr_t *vl) {
    int i, n = igraph_vector_ptr_size(vl);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = (igraph_vector_t*) VECTOR(*vl)[i];
        if (v) {
            igraph_vector_destroy(v);
        }
    }
    igraph_vector_ptr_destroy(vl);
}

static int igraph_i_graphlets(const igraph_t *graph,
                              const igraph_vector_t *weights,
                              igraph_vector_ptr_t *cliques,
                              igraph_vector_t *thresholds,
                              const igraph_vector_int_t *ids,
                              igraph_real_t startthr) {

    /* This version is different from the main function, and is
       appropriate to use in recursive calls, because it _adds_ the
       results to 'cliques' and 'thresholds' and uses the supplied
       'startthr' */

    igraph_vector_ptr_t mycliques;
    int no_of_edges = igraph_ecount(graph);
    igraph_vector_t subv;
    igraph_t subg;
    int i, nographs, nocliques;
    igraph_t *newgraphs = 0;
    igraph_vector_t *newweights = 0;
    igraph_vector_int_t *newids = 0;
    igraph_vector_t clique_thr, next_thr;
    igraph_i_subclique_next_free_t freedata = { 0, 0, 0, 0 };

    IGRAPH_CHECK(igraph_vector_ptr_init(&mycliques, 0));
    IGRAPH_FINALLY(igraph_i_graphlets_destroy_vectorlist, &mycliques);
    IGRAPH_VECTOR_INIT_FINALLY(&subv, 0);

    /* We start by finding cliques at the lowest threshold */
    for (i = 0; i < no_of_edges; i++) {
        if (VECTOR(*weights)[i] >= startthr) {
            IGRAPH_CHECK(igraph_vector_push_back(&subv, i));
        }
    }
    igraph_subgraph_edges(graph, &subg, igraph_ess_vector(&subv),
                          /*delete_vertices=*/ 0);
    IGRAPH_FINALLY(igraph_destroy, &subg);
    igraph_maximal_cliques(&subg, &mycliques, /*min_size=*/ 0, /*max_size=*/ 0);
    igraph_destroy(&subg);
    IGRAPH_FINALLY_CLEAN(1);
    nocliques = igraph_vector_ptr_size(&mycliques);

    igraph_vector_destroy(&subv);
    IGRAPH_FINALLY_CLEAN(1);

    /* Get the next cliques and thresholds */
    IGRAPH_VECTOR_INIT_FINALLY(&next_thr, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&clique_thr, 0);

    igraph_i_subclique_next(graph, weights, ids, &mycliques,
                            &newgraphs, &newweights, &newids,
                            &clique_thr, &next_thr);

    freedata.result = newgraphs;
    freedata.resultids = newids;
    freedata.resultweights = newweights;
    freedata.nc = nocliques;
    IGRAPH_FINALLY(igraph_i_subclique_next_free, &freedata);

    /* Store cliques at the current level */
    igraph_vector_append(thresholds, &clique_thr);
    for (i = 0; i < nocliques; i++) {
        igraph_vector_t *cl = (igraph_vector_t*) VECTOR(mycliques)[i];
        int j, n = igraph_vector_size(cl);
        for (j = 0; j < n; j++) {
            int node = VECTOR(*cl)[j];
            VECTOR(*cl)[j] = VECTOR(*ids)[node];
        }
        igraph_vector_sort(cl);
    }
    igraph_vector_ptr_append(cliques, &mycliques);

    /* Recursive calls for cliques found */
    nographs = igraph_vector_ptr_size(&mycliques);
    for (i = 0; i < nographs; i++) {
        igraph_t *g = newgraphs + i;
        if (igraph_vcount(g) > 1) {
            igraph_vector_t *w_sub = newweights + i;
            igraph_vector_int_t *ids_sub = newids + i;
            igraph_i_graphlets(g, w_sub, cliques, thresholds, ids_sub, VECTOR(next_thr)[i]);
        }
    }

    igraph_vector_destroy(&clique_thr);
    igraph_vector_destroy(&next_thr);
    igraph_i_subclique_next_free(&freedata);
    igraph_vector_ptr_destroy(&mycliques); /* contents was copied over */
    IGRAPH_FINALLY_CLEAN(4);

    return 0;
}

typedef struct {
    const igraph_vector_ptr_t *cliques;
    const igraph_vector_t *thresholds;
} igraph_i_graphlets_filter_t;

static int igraph_i_graphlets_filter_cmp(void *data, const void *a, const void *b) {
    igraph_i_graphlets_filter_t *ddata = (igraph_i_graphlets_filter_t *) data;
    int *aa = (int*) a;
    int *bb = (int*) b;
    igraph_real_t t_a = VECTOR(*ddata->thresholds)[*aa];
    igraph_real_t t_b = VECTOR(*ddata->thresholds)[*bb];
    igraph_vector_t *v_a, *v_b;
    int s_a, s_b;

    if (t_a < t_b) {
        return -1;
    } else if (t_a > t_b) {
        return 1;
    }

    v_a = (igraph_vector_t*) VECTOR(*ddata->cliques)[*aa];
    v_b = (igraph_vector_t*) VECTOR(*ddata->cliques)[*bb];
    s_a = igraph_vector_size(v_a);
    s_b = igraph_vector_size(v_b);

    if (s_a < s_b) {
        return -1;
    } else if (s_a > s_b) {
        return 1;
    } else {
        return 0;
    }
}

static int igraph_i_graphlets_filter(igraph_vector_ptr_t *cliques,
                                     igraph_vector_t *thresholds) {

    /* Filter out non-maximal cliques. Every non-maximal clique is
       part of a maximal clique, at the same threshold.

       First we order the cliques, according to their threshold, and
       then according to their size. So when we look for a candidate
       superset, we only need to check the cliques next in the list,
       until their threshold is different. */

    int i, iptr, nocliques = igraph_vector_ptr_size(cliques);
    igraph_vector_int_t order;
    igraph_i_graphlets_filter_t sortdata = { cliques, thresholds };

    igraph_vector_int_init(&order, nocliques);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &order);
    for (i = 0; i < nocliques; i++) {
        VECTOR(order)[i] = i;
    }

    igraph_qsort_r(VECTOR(order), nocliques, sizeof(int), &sortdata,
                   igraph_i_graphlets_filter_cmp);

    for (i = 0; i < nocliques - 1; i++) {
        int ri = VECTOR(order)[i];
        igraph_vector_t *needle = VECTOR(*cliques)[ri];
        igraph_real_t thr_i = VECTOR(*thresholds)[ri];
        int n_i = igraph_vector_size(needle);
        int j = i + 1;

        for (j = i + 1; j < nocliques; j++) {
            int rj = VECTOR(order)[j];
            igraph_real_t thr_j = VECTOR(*thresholds)[rj];
            igraph_vector_t *hay;
            int n_j, pi = 0, pj = 0;

            /* Done, not found */
            if (thr_j != thr_i) {
                break;
            }

            /* Check size of hay */
            hay = VECTOR(*cliques)[rj];
            n_j = igraph_vector_size(hay);
            if (n_i > n_j) {
                continue;
            }

            /* Check if hay is a superset */
            while (pi < n_i && pj < n_j && n_i - pi <= n_j - pj) {
                int ei = VECTOR(*needle)[pi];
                int ej = VECTOR(*hay)[pj];
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
                igraph_vector_destroy(needle);
                igraph_free(needle);
                VECTOR(*cliques)[ri] = 0;
                break;
            }
        }
    }

    /* Remove null pointers from the list of cliques */
    for (i = 0, iptr = 0; i < nocliques; i++) {
        igraph_vector_t *v = VECTOR(*cliques)[i];
        if (v) {
            VECTOR(*cliques)[iptr] = v;
            VECTOR(*thresholds)[iptr] = VECTOR(*thresholds)[i];
            iptr++;
        }
    }
    igraph_vector_ptr_resize(cliques, iptr);
    igraph_vector_resize(thresholds, iptr);

    igraph_vector_int_destroy(&order);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_graphlets_candidate_basis
 * Calculate a candidate graphlets basis
 *
 * \param graph The input graph, it must be a simple graph, edge directions are
 *        ignored.
 * \param weights Weights of the edges, a vector.
 * \param cliques An initialized vector of pointers.
 *        The graphlet basis is stored here. Each element of the pointer
 *        vector will be a vector of vertex ids. Each elements must be
 *        destroyed using \ref igraph_vector_destroy() and \ref igraph_free().
 * \param thresholds An initialized vector, the (highest possible)
 *        weight thresholds for finding the basis subgraphs are stored
 *        here.
 * \return Error code.
 *
 * See also: \ref igraph_graphlets() and \ref igraph_graphlets_project().
 */

int igraph_graphlets_candidate_basis(const igraph_t *graph,
                                     const igraph_vector_t *weights,
                                     igraph_vector_ptr_t *cliques,
                                     igraph_vector_t *thresholds) {

    int no_of_nodes = igraph_vcount(graph);
    int no_of_edges = igraph_ecount(graph);
    igraph_real_t minthr;
    igraph_vector_int_t ids;
    igraph_bool_t simple;
    int i;

    /* Some checks */
    if (weights == NULL) {
        IGRAPH_ERROR("Graphlet functions require weighted graphs", IGRAPH_EINVAL);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }

    igraph_is_simple(graph, &simple);
    if (!simple) {
        IGRAPH_ERROR("Graphlets work on simple graphs only", IGRAPH_EINVAL);
    }

    minthr = igraph_vector_min(weights);
    igraph_vector_ptr_clear(cliques);
    igraph_vector_clear(thresholds);
    igraph_vector_int_init(&ids, no_of_nodes);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &ids);
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(ids)[i] = i;
    }

    igraph_i_graphlets(graph, weights, cliques, thresholds, &ids, minthr);

    igraph_vector_int_destroy(&ids);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_i_graphlets_filter(cliques, thresholds);

    return 0;
}

/* TODO: not made static because it is used by the R interface */
int igraph_i_graphlets_project(const igraph_t *graph,
                               const igraph_vector_t *weights,
                               const igraph_vector_ptr_t *cliques,
                               igraph_vector_t *Mu, igraph_bool_t startMu,
                               int niter, int vid1) {

    int no_of_nodes = igraph_vcount(graph);
    int no_of_edges = igraph_ecount(graph);
    int no_cliques = igraph_vector_ptr_size(cliques);
    igraph_vector_int_t vcl, vclidx, ecl, eclidx, cel, celidx;
    igraph_vector_t edgelist, newweights, normfact;
    int i, total_vertices, e, ptr, total_edges;
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
    igraph_is_simple(graph, &simple);
    if (!simple) {
        IGRAPH_ERROR("Graphlets work on simple graphs only", IGRAPH_EINVAL);
    }

    if (!startMu) {
        igraph_vector_resize(Mu, no_cliques);
        igraph_vector_fill(Mu, 1);
    }

    /* Count # cliques per vertex. Also, create an index
       for the edges per clique. */
    IGRAPH_CHECK(igraph_vector_int_init(&vclidx, no_of_nodes + 2));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &vclidx);
    IGRAPH_CHECK(igraph_vector_int_init(&celidx, no_cliques + 3));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &celidx);
    for (i = 0, total_vertices = 0, total_edges = 0; i < no_cliques; i++) {
        igraph_vector_t *v = VECTOR(*cliques)[i];
        int j, n = igraph_vector_size(v);
        total_vertices += n;
        total_edges += n * (n - 1) / 2;
        VECTOR(celidx)[i + 2] = total_edges;
        for (j = 0; j < n; j++) {
            int vv = VECTOR(*v)[j] - vid1;
            VECTOR(vclidx)[vv + 2] += 1;
        }
    }
    VECTOR(celidx)[i + 2] = total_edges;

    /* Finalize index vector */
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(vclidx)[i + 2] += VECTOR(vclidx)[i + 1];
    }

    /* Create vertex-clique list, the cliques for each vertex. */
    IGRAPH_CHECK(igraph_vector_int_init(&vcl, total_vertices));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &vcl);
    for (i = 0; i < no_cliques; i++) {
        igraph_vector_t *v = VECTOR(*cliques)[i];
        int j, n = igraph_vector_size(v);
        for (j = 0; j < n; j++) {
            int vv = VECTOR(*v)[j] - vid1;
            int p = VECTOR(vclidx)[vv + 1];
            VECTOR(vcl)[p] = i;
            VECTOR(vclidx)[vv + 1] += 1;
        }
    }

    /* Create an edge-clique list, the cliques of each edge */
    IGRAPH_CHECK(igraph_vector_int_init(&ecl, total_edges));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &ecl);
    IGRAPH_CHECK(igraph_vector_int_init(&eclidx, no_of_edges + 1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &eclidx);
    IGRAPH_CHECK(igraph_vector_init(&edgelist, no_of_edges * 2));
    IGRAPH_FINALLY(igraph_vector_destroy, &edgelist);
    IGRAPH_CHECK(igraph_get_edgelist(graph, &edgelist, /*by_col=*/ 0));
    for (i = 0, e = 0, ptr = 0; e < no_of_edges; e++) {
        int from = VECTOR(edgelist)[i++];
        int to = VECTOR(edgelist)[i++];
        int from_s = VECTOR(vclidx)[from];
        int from_e = VECTOR(vclidx)[from + 1];
        int to_s = VECTOR(vclidx)[to];
        int to_e = VECTOR(vclidx)[to + 1];
        VECTOR(eclidx)[e] = ptr;
        while (from_s < from_e && to_s < to_e) {
            int from_v = VECTOR(vcl)[from_s];
            int to_v = VECTOR(vcl)[to_s];
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

    igraph_vector_destroy(&edgelist);
    IGRAPH_FINALLY_CLEAN(1);

    /* Convert the edge-clique list to a clique-edge list */
    IGRAPH_CHECK(igraph_vector_int_init(&cel, total_edges));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &cel);
    for (i = 0; i < no_of_edges; i++) {
        int ecl_s = VECTOR(eclidx)[i], ecl_e = VECTOR(eclidx)[i + 1], j;
        for (j = ecl_s; j < ecl_e; j++) {
            int cl = VECTOR(ecl)[j];
            int epos = VECTOR(celidx)[cl + 1];
            VECTOR(cel)[epos] = i;
            VECTOR(celidx)[cl + 1] += 1;
        }
    }

    /* Normalizing factors for the iteration */
    IGRAPH_CHECK(igraph_vector_init(&normfact, no_cliques));
    IGRAPH_FINALLY(igraph_vector_destroy, &normfact);
    for (i = 0; i < no_cliques; i++) {
        igraph_vector_t *v = VECTOR(*cliques)[i];
        int n = igraph_vector_size(v);
        VECTOR(normfact)[i] = n * (n + 1) / 2;
    }

    /* We have the clique-edge list, so do the projection now */
    IGRAPH_CHECK(igraph_vector_init(&newweights, no_of_edges));
    IGRAPH_FINALLY(igraph_vector_destroy, &newweights);
    for (i = 0; i < niter; i++) {
        for (e = 0; e < no_of_edges; e++) {
            int start = VECTOR(eclidx)[e];
            int end = VECTOR(eclidx)[e + 1];
            VECTOR(newweights)[e] = 0.0001;
            while (start < end) {
                int clique = VECTOR(ecl)[start++];
                VECTOR(newweights)[e] += VECTOR(*Mu)[clique];
            }
        }
        for (e = 0; e < no_cliques; e++) {
            igraph_real_t sumratio = 0;
            int start = VECTOR(celidx)[e];
            int end = VECTOR(celidx)[e + 1];
            while (start < end) {
                int edge = VECTOR(cel)[start++];
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

    return 0;
}

/**
 * \function igraph_graphlets_project
 * Project a graph on a graphlets basis
 *
 * Note that the graph projected does not have to be the same that
 * was used to calculate the graphlet basis, but it is assumed that
 * it has the same number of vertices, and the vertex ids of the two
 * graphs match.
 * \param graph The input graph, it must be a simple graph, edge directions are
 *        ignored.
 * \param weights Weights of the edges in the input graph, a vector.
 * \param cliques The graphlet basis, a pointer vector, in which each
 *        element is a vector of vertex ids.
 * \param Mu An initialized vector, the weights of the graphlets will
 *        be stored here. This vector is also used to initialize the
 *        the weight vector for the iterative algorithm, if the
 *        \c startMu argument is true (non-zero).
 * \param startMu If true (non-zero), then the supplied Mu vector is
 *        used as the starting point of the iteration. Otherwise a
 *        constant 1 vector is used.
 * \param niter Integer scalar, the number of iterations to perform.
 * \return Error code.
 *
 * See also: \ref igraph_graphlets() and
 * \ref igraph_graphlets_candidate_basis().
 */

int igraph_graphlets_project(const igraph_t *graph,
                             const igraph_vector_t *weights,
                             const igraph_vector_ptr_t *cliques,
                             igraph_vector_t *Mu, igraph_bool_t startMu,
                             int niter) {

    return igraph_i_graphlets_project(graph, weights, cliques, Mu, startMu,
                                      niter, /*vid1=*/ 0);
}

typedef struct igraph_i_graphlets_order_t {
    const igraph_vector_ptr_t *cliques;
    const igraph_vector_t *Mu;
} igraph_i_graphlets_order_t;

static int igraph_i_graphlets_order_cmp(void *data, const void *a, const void *b) {
    igraph_i_graphlets_order_t *ddata = (igraph_i_graphlets_order_t*) data;
    int *aa = (int*) a;
    int *bb = (int*) b;
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
 * \param cliques An initialized vector of pointers.
 *        The graphlet basis is stored here. Each element of the pointer
 *        vector will be a vector of vertex ids.
 * \param Mu An initialized vector, the weights of the graphlets will
 *        be stored here.
 * \param niter Integer scalar, the number of iterations to perform
 *        for the projection step.
 * \return Error code.
 *
 * See also: \ref igraph_graphlets_candidate_basis() and
 * \ref igraph_graphlets_project().
 */

int igraph_graphlets(const igraph_t *graph,
                     const igraph_vector_t *weights,
                     igraph_vector_ptr_t *cliques,
                     igraph_vector_t *Mu, int niter) {

    int i, nocliques;
    igraph_vector_t thresholds;
    igraph_vector_int_t order;
    igraph_i_graphlets_order_t sortdata = { cliques, Mu };

    igraph_vector_init(&thresholds, 0);
    IGRAPH_FINALLY(igraph_vector_destroy, &thresholds);
    igraph_graphlets_candidate_basis(graph, weights, cliques, &thresholds);
    igraph_vector_destroy(&thresholds);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_graphlets_project(graph, weights, cliques, Mu, /*startMu=*/ 0, niter);

    nocliques = igraph_vector_ptr_size(cliques);
    igraph_vector_int_init(&order, nocliques);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &order);
    for (i = 0; i < nocliques; i++) {
        VECTOR(order)[i] = i;
    }
    igraph_qsort_r(VECTOR(order), nocliques, sizeof(int), &sortdata,
                   igraph_i_graphlets_order_cmp);

    igraph_vector_ptr_index_int(cliques, &order);
    igraph_vector_index_int(Mu, &order);

    igraph_vector_int_destroy(&order);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}
