/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2008-2021  The igraph development team <igraph@igraph.org>

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

#include "igraph_layout.h"

#include "igraph_interface.h"
#include "igraph_lapack.h"
#include "igraph_matrix.h"
#include "igraph_nongraph.h"
#include "igraph_random.h"
#include "igraph_vector_list.h"

#include "layout/layout_internal.h"
#include "core/interruption.h"

#include <math.h>

/* This file contains the implementation of the UMAP algorithm.
 *
 * UMAP is typically used as a to reduce dimensionality of vectors, embedding them in
 * 2D (or, less commonly, in 3D). Despite this geometric flair, UMAP heavily relies on
 * graphs as intermediate data structures and is therefore a useful graph layout
 * algorithm in its own right. Conceptually, there are three steps:
 *
 * 1. Compute a sparse graph with edges connecting similar vectors, e.g. a k-nearest
 *    neighbor graph. A vector of distances is associated with the graph edges. This
 *    file does *not* perform this part of the computation since there are many
 *    libraries out there that can compute knn or other sparse graphs efficiently
 *    starting from vector spaces (e.g. faiss).
 * 2. Convert the distances into weights, which are weights between 0 and 1
 *    that are larger for short-distance edges. This step is exposed via
 *    igraph_layout_umap_compute_weights.
 * 3. Compute a layout for the graph, using its associated weights as edge
 *    weights. This step is exposed via igraph_layout_umap and its 3D counterpart.
 *    These two fuctions can also compute steps 2 and 3 in one go, since that's the
 *    most common use case: the argument "distances_are_weights" should be
 *    set to false.
 *
 * A few more details w/r/t steps 2 and 3, since they are computed in detail below.
 *
 * STEP 2
 * For each vertex, the distance to its closest neighbor, called rho, is "forfeited":
 * that edge begets weight 1 (in principle, at least). Farther neighbors beget
 * lower weights according to an exponential decay. The scale factor of this
 * decay is called sigma and is computed from the graph itself.
 *
 * STEP 3
 * The layout is computed via stochastic gradient descent, i.e. applying stochastic
 * forces along high-weight edges and, more rarely, low-weight edges.
 * To compute the stochastic forces, one needs a smooth function that approximates
 * weights but in the embedded space:
 *                        Q(d) = ( 1 + a*d^2b )^-1
 * where d is the 2D/3D distance between the vertices and a and b are constants that
 * are computed globally based on a user-chosen fudge parameter called min_dist.
 * Smaller min_dist will give rise to slightly more compact embeddings. We find a
 * and b via gradient descent, which is implemented de novo below.
 *
 * Repulsion is computed via negative sampling, typically a few nodes are picked
 * at random as repulsive sources each time an attractive force is computed.
 *
 * During the stochastic gradient descent, the learning rate - a multiplicative factor
 * on top of the stochastic forces themselves - is reduced linearly from 1 to 0. At
 * the end, the stochastic forces can be strong but their effect is reduced to almost
 * nothing by the small learning rate. Notice that UMAP does not formally converge:
 * instead, we reduce the forces' impact steadily to a trickle and finally quench it
 * altogether.
 *
 * FINAL COMMENTS
 * This implementation uses a few more tricks to improve the result:
 * - a few constants are defined to limit the force applied to vertices at each step
 *   and other geometric corrections
 * - the layout is centered at the end of the computation.
 * - a seed layout can be used. Notice that since UMAP runs for an essentially fixed
 *   time rather than until convergence, using a good/bad seed does not affect
 *   runtimes significantly.
 * */
#define UMAP_FORCE_LIMIT 4
#define UMAP_MIN_DISTANCE_ATTRACTION 0.0001
#define UMAP_CORRECT_DISTANCE_REPULSION 0.01

/* Find sigma for this vertex by binary search */
static igraph_error_t igraph_i_umap_find_sigma(const igraph_vector_t *distances,
        const igraph_vector_int_t *eids,
        igraph_real_t rho, igraph_real_t *sigma_p, igraph_real_t target) {

    igraph_real_t sigma = 1;
    igraph_real_t sum;
    igraph_real_t tol = 0.01;
    igraph_integer_t maxiter = 100;
    igraph_integer_t no_of_neis = igraph_vector_int_size(eids);
    igraph_integer_t eid;
    igraph_real_t step = sigma;
    igraph_integer_t seen_max = 0;

    /* Binary search */
    for (igraph_integer_t iter = 0; iter < maxiter; iter++) {
        sum = 0;
        for (igraph_integer_t j = 0; j < no_of_neis; j++) {
            eid = VECTOR(*eids)[j];
            sum += exp(-(VECTOR(*distances)[eid] - rho) / sigma);
        }

#ifdef UMAP_DEBUG
        printf("SIGMA function (no_of_neis = %" IGRAPH_PRId ")- sum: %g, "
               "target: %g, rho: %g, sigma: %g\n", no_of_neis, sum, target, rho, sigma);
#endif

        if (sum < target) {
            /* going back up after having seen an upper bound */
            if (seen_max == 1) {
                step /= 2;
            /* we need to go up but have not seen an upper bound yet
             * first iteration we want to increase by sigma, else we must come from
             * below, so we are sitting at 2 * step, we want to move to 4 * step */
            } else if (iter > 0) {
                step *= 2;
            }
            sigma += step;
        /* overshooting, we have definitely seen the max */
        } else {
            seen_max = 1;
            step /= 2;
            sigma -= step;
        }

        /* Check for convergence */
        if (fabs(sum - target) < tol) {
            break;
        }
    }

    *sigma_p = sigma;

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_layout_umap_compute_weights
 * \brief Compute weights for a UMAP layout starting from distances.
 *
 * \experimental
 *
 * UMAP is used to embed high-dimensional vectors in a low-dimensional space
 * (most commonly 2D). It uses a distance graph as an intermediate data structure,
 * making it also a useful graph layout algorithm. See \ref igraph_layout_umap()
 * for more information.
 *
 * </para><para>
 *
 * An early step in UMAP is to compute exponentially decaying "weights" from the
 * distance graph. Connectivities can also be viewed as edge weights that quantify
 * similarity between two vertices. This function computes weights from the
 * distance graph. To compute the layout from precomputed weights, call
 * \ref igraph_layout_umap() with the \p distances_are_weights argument set to \c true.
 *
 * </para><para>
 *
 * While the distance graph can be directed (e.g. in a k-nearest neighbors, it is
 * clear *who* you are a neighbor of), the weights are usually undirected. Whenever two
 * vertices are doubly connected in the distance graph, the resulting weight W is set as:
 *
 * W = W1 + W2 - W1 * W2
 *
 * Because UMAP weights are interpreted as probabilities, this is just the probability
 * that either edge is present, without double counting. It is called "fuzzy union" in
 * the original UMAP implementation and is the default. One could also require that both
 * edges are there, i.e. W = W1 * W2: this would represent the fuzzy intersection and is
 * not implemented in igraph. As a consequence of this symmetrization, information is lost,
 * i.e. one needs fewer weights than one had distances. To keep things efficient, here
 * we set the weight for one of the two edges as above and the weight for its opposite edge
 * as 0, so that it will be skipped in the UMAP gradient descent later on.
 *
 * </para><para>
 *
 * Technical note: For each vertex, this function computes its scale factor (sigma),
 * its connectivity correction (rho), and finally the weights themselves.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Leland McInnes, John Healy, and James Melville. https://arxiv.org/abs/1802.03426
 *
 * \param graph Pointer to the distance graph. This can be directed (e.g. connecting
 *   each vertex to its neighbors in a k-nearest neighbor) or undirected, but must
 *   have no loops nor parallel edges. The only exception is: if the graph is directed,
 *   having pairs of edges with opposite direction is accepted.
 * \param distances Pointer to the vector with the vertex-to-vertex distance associated with
 *   each edge. This argument can be NULL, in which case all edges are assumed to have the
 *   same distance.
 * \param weights Pointer to an initialized vector where the result will be stored. If the
 *   input graph is directed, the weights represent a symmetrized version which contains
 *   less information. Therefore, whenever two edges between the same vertices and opposite
 *   direction are present in the input graph, only one of the weights is set and the other
 *   is fixed to zero. That format is accepted by \ref igraph_layout_umap(), which skips
 *   all zero-weight edges from the layout optimization.
 *
 * \return Error code.
 */
igraph_error_t igraph_layout_umap_compute_weights(
        const igraph_t *graph,
        const igraph_vector_t *distances,
        igraph_vector_t *weights) {

    igraph_integer_t no_of_vertices = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_neis, eid, i, j, k, l;
    igraph_vector_int_t eids;
    igraph_vector_int_list_t neighbors_seen;
    igraph_vector_list_t weights_seen;
    igraph_vector_int_t* neighbors_seen_elt;
    igraph_vector_t* weights_seen_elt;
    igraph_real_t rho, dist_max, dist, sigma, weight, weight_inv, sigma_target, dist_min;

    /* reserve memory for the weights */
    IGRAPH_CHECK(igraph_vector_resize(weights, no_of_edges));

    /* UMAP is sometimes used on unweighted graphs, otherwise check distance vector. */
    if (distances != NULL) {
        if (igraph_vector_size(distances) != no_of_edges) {
            IGRAPH_ERROR("Distances must be the same number as the edges in the graph.", IGRAPH_EINVAL);
        }
        if (no_of_edges > 0) {
            dist_min = igraph_vector_min(distances);
            if (dist_min < 0) {
                IGRAPH_ERROR("Distance values must not be negative.", IGRAPH_EINVAL);
            } else if (isnan(dist_min)) {
                IGRAPH_ERROR("Distance values must not be NaN.", IGRAPH_EINVAL);
            }
        }
    }

    /* Initialize auxiliary vectors */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eids, 0);
    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&neighbors_seen, no_of_vertices);
    IGRAPH_VECTOR_LIST_INIT_FINALLY(&weights_seen, no_of_vertices);

    /* Iterate over vertices x, like in the paper */
    for (i = 0; i < no_of_vertices; i++) {
        /* Edges out of this vertex, e.g. to its k-nearest neighbors */
        IGRAPH_CHECK(igraph_incident(graph, &eids, i, IGRAPH_OUT));
        no_of_neis = igraph_vector_int_size(&eids);

        /* Vertex has no neighbors */
        if (no_of_neis == 0) {
            continue;
        }

        /* Find rho for this vertex, i.e. the minimal non-self distance */
        if (distances != NULL) {
            rho = VECTOR(*distances)[VECTOR(eids)[0]];
            dist_max = rho;
            for (j = 1; j < no_of_neis; j++) {
                eid = VECTOR(eids)[j];
                dist = VECTOR(*distances)[eid];
                rho = fmin(rho, dist);
                dist_max = fmax(dist_max, dist);
            }
        } else {
            rho = dist_max = 0;
        }

        /* If the maximal distance is rho, all neighbors are identical to
         * each other. This can happen e.g. if distances == NULL. */
        if (dist_max == rho) {
            /* This is a special flag for later on */
            sigma = -1;

        /* Else, find sigma for this vertex, from its rho plus binary search */
        } else {
            sigma_target = log2(no_of_neis);
            IGRAPH_CHECK(igraph_i_umap_find_sigma(distances,
                        &eids, rho, &sigma,
                        sigma_target));
        }

        /* Convert to weights */
        for (j = 0; j < no_of_neis; j++) {
            eid = VECTOR(eids)[j];

            /* Basically, nodes closer than rho have probability 1, the rest is
             * exponentially penalized keeping rough cardinality */
            weight = sigma < 0 ? 1 : exp(-(VECTOR(*distances)[eid] - rho) / sigma);

            #ifdef UMAP_DEBUG
            if (distances != NULL)
                printf("distance: %g\n", VECTOR(*distances)[eid]);
            printf("weight: %g\n", weight);
            #endif

            /* Store in vector lists for later symmetrization */
            k = IGRAPH_OTHER(graph, eid, i);
            if (k == i) {
                IGRAPH_ERROR("Input graph must contain no self-loops.", IGRAPH_EINVAL);
            }

            neighbors_seen_elt = igraph_vector_int_list_get_ptr(&neighbors_seen, i);
            IGRAPH_CHECK(igraph_vector_int_push_back(neighbors_seen_elt, k));

            weights_seen_elt = igraph_vector_list_get_ptr(&weights_seen, i);
            IGRAPH_CHECK(igraph_vector_push_back(weights_seen_elt, weight));
        }

    }

    /* Symmetrize the weights. UMAP weights are probabilities of that edge being a
     * "real" connection. Unlike the distances, which can represent a directed graph,
     * weights are usually symmetric. We symmetrize via fuzzy union. */
    for (eid=0; eid < no_of_edges; eid++) {
        i = IGRAPH_FROM(graph, eid);
        k = IGRAPH_TO(graph, eid);

        /* Direct weight, if found */
        /* NOTE: this and the subsequent loop could be faster if we sorted the vectors
         * beforehand. Probably not such a big deal. */
        weight = 0;
        neighbors_seen_elt = igraph_vector_int_list_get_ptr(&neighbors_seen, i);
        weights_seen_elt = igraph_vector_list_get_ptr(&weights_seen, i);
        no_of_neis = igraph_vector_int_size(neighbors_seen_elt);
        for (l=0; l < no_of_neis; l++) {
            if (VECTOR(*neighbors_seen_elt)[l] == k) {
                weight = VECTOR(*weights_seen_elt)[l];
                /* Tag this weight so we can ignore it later on if the opposite
                 * directed edge is found. It's ok to retag */
                VECTOR(*weights_seen_elt)[l] = -1;
                break;
            }
        }

        /* The opposite edge has already been union-ed, set this one to -1 */
        if (weight < 0) {
            VECTOR(*weights)[eid] = 0;
            continue;
        }

        /* Weight of the opposite edge, if found */
        weight_inv = 0;
        neighbors_seen_elt = igraph_vector_int_list_get_ptr(&neighbors_seen, k);
        weights_seen_elt = igraph_vector_list_get_ptr(&weights_seen, k);
        no_of_neis = igraph_vector_int_size(neighbors_seen_elt);
        for (l=0; l < no_of_neis; l++) {
            if (VECTOR(*neighbors_seen_elt)[l] == i) {
                weight_inv = VECTOR(*weights_seen_elt)[l];
                /* Tag this weight so we can ignore it later on if the opposite
                 * directed edge is found. It's ok to retag */
                VECTOR(*weights_seen_elt)[l] = -1;
                break;
            }
        }

        /* The opposite edge has already been union-ed, set this one to -1 */
        if (weight_inv < 0) {
            VECTOR(*weights)[eid] = 0;
            continue;
        }

        /* First time this edge or its opposite are seen, set the W */
        VECTOR(*weights)[eid] = weight + weight_inv - weight * weight_inv;
    }

    igraph_vector_list_destroy(&weights_seen);
    igraph_vector_int_list_destroy(&neighbors_seen);
    igraph_vector_int_destroy(&eids);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}


/* Helper function to compute a and b parameters (smoothing probability metric in embedding space) */
static igraph_error_t igraph_i_umap_get_ab_residuals(igraph_vector_t *residuals,
        igraph_real_t *squared_sum_res, igraph_integer_t nr_points, igraph_real_t a,
        igraph_real_t b, igraph_vector_t *powb, const igraph_vector_t *x, igraph_real_t min_dist)
{
    igraph_real_t tmp;

    *squared_sum_res = 0;
    for (igraph_integer_t i = 0; i < nr_points; i++) {
        /* The ideal probability is:
         *
         *     P(d) = d < min_dist ? 1 : e^{-(d - min_dist)}
         *
         * which is the same as the high-dimensional probability, except
         * min_dist plays the role of rho and sigma is fixed at 1. However,
         * this function has a kink at min_dist (first derivative is not
         * continuous). So we smoothen it with:
         *
         *     Q(d) = ( 1 + a*d^2b )^-1
         *
         * which is quite similar throughout for appropriate a and b. Notice
         * that we do not need to smoothen the high-dimensional probability
         * function because the vertices are not moved in the high-dimensional
         * space, so there is no need for differentiating that function.
         *
         * The residual is of course:
         *
         *    Q(d) - P(d) = ( 1 + a*d^2b )^-1 - [ d < min_dist ? 1 : e^{-(d - min_dist)} ]
         *
         * This function also sets the auxiliary vector powb.
         * */
        VECTOR(*powb)[i] = pow(VECTOR(*x)[i], 2 * b);
        tmp = 1 / (1 + a * VECTOR(*powb)[i]);
        tmp -= VECTOR(*x)[i] <= min_dist ? 1 : exp(-(VECTOR(*x)[i] - min_dist));
        VECTOR(*residuals)[i] = tmp;
        *squared_sum_res += tmp * tmp;
    }
    return IGRAPH_SUCCESS;
}

/* UMAP minimizes the cross-entropy between probability of being a true edge in
 * high and low dimensions. For the low-dimensional computation, it uses a smooth
 * function of the Euclidean distance between two vertices:
 *
 * P(d) = (1 + a*d^2b)^-1
 *
 * where d is the distance and a and b are hyperparameters that basically determine
 * the cutoff distance at which the probability starts to decrease.
 *
 * We fit these two parameters using nonlinear least squares (Gauss-Newton + line search)
 * on a grid of artificial distances. There is only one user-chosen input argument that
 * determines this fit, called min_dist, which is approximately the cutoff distance we
 * are trying to achieve.
 *
 * ADVANCED NOTE:
 * In a way, the whole UMAP layout is invariant upon scaling transformations, of course,
 * so min_dist is basically meaningless. Another way to see this is that for any pair
 * (a,b) that minimize the least squares for dist_min, we can easily find a solution for
 * a new dist_min2 := alpha * dist_min:
 *
 * P(d, a, b) = (1 + a*d^2b)^-1
 *
 * P(alpha * d, a', b') = (1 + a'*(alpha * d)^2b' )^-1
 *
 * that is:
 *
 * a*d^2b = a'*alpha^2b'*d^2b'   for each  d >= 0.
 *
 * So for d = 1        ->  a = a'*alpha^2b'
 * and for d = sqrt(2) ->  a*2^b = a'*alpha^2b'*2^b'
 *
 * which solves as:
 *
 * b' = b
 * a' = a / alpha^2b
 *
 * For instance, if b = 1, a -> 0.01*a moves the fit a decade towards larger min_dist,
 * and a -> 100*a moves the fit a decade towards smaller min_dist.
 * */
igraph_error_t igraph_i_umap_fit_ab(igraph_real_t min_dist, igraph_real_t *a_p, igraph_real_t *b_p)
{
    /* Grid points */
    igraph_vector_t x;
     /* Make a lattice from 0 to 3 * sigma with 300 points. This is what
      * umap.umap_.fit_ab_params does, but sigma is fixed to 1.0 here since
      * that's the default value used in scanpy and by virtually everyone */
    igraph_integer_t nr_points = 300;
    igraph_real_t end_point = 3.0;
    /* Initial values takes as reasonable assumptions from typical min_dist values */
    igraph_real_t b = 0.8;
    igraph_real_t a = 1.8;
    /* deltas */
    igraph_real_t da, db;
    /* Residuals */
    igraph_vector_t residuals;
    igraph_real_t squared_sum_res, squared_sum_res_old, squared_sum_res_tmp;
    /* Needed for the Gauss-Newton search */
    igraph_matrix_t jacobian, jTj, jTr;
    igraph_real_t tol = 0.001;
    igraph_real_t maxiter = 100;
    /* Auxiliary vars */
    igraph_real_t tmp;
    igraph_vector_t powb;
    int lapack_info;

    /* Distance lattice */
    IGRAPH_VECTOR_INIT_FINALLY(&x, nr_points);
    /* Residuals */
    IGRAPH_VECTOR_INIT_FINALLY(&residuals, nr_points);
    /* First derivatives, for the fitting (direction) */
    IGRAPH_MATRIX_INIT_FINALLY(&jacobian, nr_points, 2);
    /* Composite matrices/vectors for linear least squares at each iteration */
    IGRAPH_MATRIX_INIT_FINALLY(&jTj, 2, 2);
    IGRAPH_MATRIX_INIT_FINALLY(&jTr, 2, 1);
    /* Auxiliary vars for convenience */
    IGRAPH_VECTOR_INIT_FINALLY(&powb, nr_points);

    /* Distance |x-y| (this is a lattice, there are no actual x and y) */
    for (igraph_integer_t i = 0; i < nr_points; i++) {
        VECTOR(x)[i] = (end_point / nr_points) * i + 0.001; /* added a 0.001 to prevent NaNs */
    }

    /* Initialize squared_sum_res_old to a dummy value to prevent some compilers
     * from complaining about uninitialized values */
    squared_sum_res_old = IGRAPH_INFINITY;

#ifdef UMAP_DEBUG
    printf("start fit_ab\n");
#endif
    for (igraph_integer_t iter = 0; iter < maxiter; iter++) {
        IGRAPH_CHECK(igraph_i_umap_get_ab_residuals(&residuals, &squared_sum_res, nr_points, a, b,
                    &powb, &x, min_dist));

        /* break if good fit (conergence to truth) */
        if (squared_sum_res < tol * tol) {
#ifdef UMAP_DEBUG
            printf("convergence to zero (wow!)\n");
#endif
            break;
        }
        /* break if no change (convergence) */
        if ((iter > 0) && fabs(sqrt(squared_sum_res_old) - sqrt(squared_sum_res)) < tol) {
#ifdef UMAP_DEBUG
            printf("no-change absolute convergence\n");
#endif
            break;
        }

        /* Jacobian (first derivatives) of squared residuals at (a, b) */
        for (igraph_integer_t i = 0; i < nr_points; i++) {
            tmp = 1 + a * VECTOR(powb)[i];
            MATRIX(jacobian, i, 0) = - 2 * VECTOR(powb)[i] / tmp / tmp;
            MATRIX(jacobian, i, 1) = MATRIX(jacobian, i, 0) * a * log(VECTOR(x)[i]) * 2;
        }

        /* At each iteration, we want to minimize the linear approximation of the sum of squared
         * residuals:
         *
         * sum_i (Ji @ d(a,b) -r_i)^2
         *
         * Putting the first derivative to zero results in a linear system of 2 equations
         * (for a and b):
         *
         * sum_i J_i^T @ J_i @ d(a,b) = sum_i J_i^T r_i
         * *
         * or more compactly:
         *
         * J^T @ J @ d(a,b) = J^T @ r
         *
         * where J_T is the transpose of the Jacobian. Defining A := J^T @ J, B = J^T @ r:
         *
         * A @ d(a,b) = B
         *
         * This can be solved for d(a,b) using LAPACK within igraph
         * */
        /* Compute A and B, i.e. J^T @ J and J^T @ r */
        MATRIX(jTj, 0, 0) = MATRIX(jTj, 0, 1) = MATRIX(jTj, 1, 0) = MATRIX(jTj, 1, 1) = 0;
        MATRIX(jTr, 0, 0) = MATRIX(jTr, 1, 0) = 0;
        for (igraph_integer_t i = 0; i < nr_points; i++) {
            for (igraph_integer_t j1 = 0; j1 < 2; j1++) {
                for (igraph_integer_t j2 = 0; j2 < 2; j2++) {
                    MATRIX(jTj, j1, j2) += MATRIX(jacobian, i, j1) * MATRIX(jacobian, i, j2);
                }
                MATRIX(jTr, j1, 0) += MATRIX(jacobian, i, j1) * VECTOR(residuals)[i];
            }
        }
        /* LAPACK puts solution into jTr */
        IGRAPH_CHECK(igraph_lapack_dgesv(&jTj, 0, &jTr, &lapack_info));

        /* This might go wrong, in which case we should fail graciously */
        if (lapack_info != 0) {
            IGRAPH_ERROR("Singular matrix in the estimation of a and b for UMAP",  IGRAPH_EINVAL);
        }

        da = -MATRIX(jTr, 0, 0);
        db = -MATRIX(jTr, 1, 0);

        /* Improvement over GN: rough exponential line search for best delta
         * start from largest change, and keep shrinking as long as we are going down
         * */
        squared_sum_res_old = squared_sum_res;
        IGRAPH_CHECK(igraph_i_umap_get_ab_residuals(&residuals, &squared_sum_res, nr_points, a + da,
                    b + db, &powb, &x, min_dist));

#ifdef UMAP_DEBUG
        printf("start line search, SSR before delta: %g, current SSR:, %g\n", squared_sum_res_old,
                squared_sum_res);
#endif
        for (igraph_integer_t k = 0; k < 30; k++) {
            /* Try new parameters */
            da /= 2.0;
            db /= 2.0;
            squared_sum_res_tmp = squared_sum_res;
            IGRAPH_CHECK(igraph_i_umap_get_ab_residuals(&residuals, &squared_sum_res, nr_points,
                        a + da, b + db, &powb, &x, min_dist));

            /* Compare and if we are going back uphill, undo last step and break */
#ifdef UMAP_DEBUG
            printf("during line search, k = %d, old SSR:, %g, new SSR (half a,b):, %g\n", k,
                    squared_sum_res_tmp, squared_sum_res);
#endif
            if (squared_sum_res > squared_sum_res_tmp - tol) {
                da *= 2;
                db *= 2;
                break;
            }
        }
#ifdef UMAP_DEBUG
        printf("end of line search and iteration, squared_sum_res: %g \n\n", squared_sum_res_tmp);
#endif

        /* assign a, b*/
        a += da;
        b += db;

    }

    /* Free memory and tidy up stack */
    igraph_vector_destroy(&x);
    igraph_vector_destroy(&residuals);
    igraph_matrix_destroy(&jacobian);
    igraph_matrix_destroy(&jTj);
    igraph_matrix_destroy(&jTr);
    igraph_vector_destroy(&powb);
    IGRAPH_FINALLY_CLEAN(6);

#ifdef UMAP_DEBUG
    printf("a, b: %g %g\n", a, b);
#endif

    *a_p = a;
    *b_p = b;

    return IGRAPH_SUCCESS;

}

/* cross-entropy */
#ifdef UMAP_DEBUG
static igraph_error_t igraph_i_umap_compute_cross_entropy(const igraph_t *graph,
       const igraph_vector_t *umap_weights, const igraph_matrix_t *layout, igraph_real_t a, igraph_real_t b,
       igraph_real_t *cross_entropy) {

    igraph_real_t mu, nu, xd, yd, sqd;
    igraph_integer_t from, to;
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_vertices = igraph_vcount(graph);
    igraph_matrix_t edge_seen;

    IGRAPH_MATRIX_INIT_FINALLY(&edge_seen, no_of_vertices, no_of_vertices);

    /* Measure the (variable part of the) cross-entropy terms for debugging:
     * 1. - sum_edge_e mu(e) * log(nu(e))
     * 2. - sum_edge_e (1 - mu(e)) * log(1 - nu(e))
     * NOTE: the sum goes over the whole adjacency matrix, i.e. all potential edges,
     * not just the actual edges. That is because otherwise there's no benefit from
     * repelling unconnected edges.
     * */
    *cross_entropy = 0;
    for (igraph_integer_t eid = 0; eid < no_of_edges; eid++) {
        mu = VECTOR(*umap_weights)[eid];

        /* Find vertices */
        from = IGRAPH_FROM(graph, eid);
        to = IGRAPH_TO(graph, eid);
        /* Find distance in layout space */
        xd = (MATRIX(*layout, from, 0) - MATRIX(*layout, to, 0));
        yd = (MATRIX(*layout, from, 1) - MATRIX(*layout, to, 1));
        sqd = xd * xd + yd * yd;
        /* Find probability associated with distance using fitted Phi */
        nu = 1.0 / (1 + a * pow(sqd, b));

        /* Term 1: entropy from the edges */
        if (mu > 0)
            *cross_entropy -= mu * log(nu);
        /* Term 2: entropy from the missing edges */
        if (mu < 1)
            *cross_entropy -= (1 - mu) * log(1 - nu);

        MATRIX(edge_seen, from, to) = MATRIX(edge_seen, to, from) = 1;
    }
    /* Add the entropy from the missing edges */
    for (igraph_integer_t from = 0; from < no_of_vertices; from++) {
        for (igraph_integer_t to = 0; to < from; to++) {
            if (MATRIX(edge_seen, from, to) > 0) {
                continue;
            }

            /* Find distance in layout space */
            xd = (MATRIX(*layout, from, 0) - MATRIX(*layout, to, 0));
            yd = (MATRIX(*layout, from, 1) - MATRIX(*layout, to, 1));
            sqd = xd * xd + yd * yd;

            /* Find probability associated with distance using fitted Phi */
            nu = 1.0 / (1 + a * pow(sqd, b));

            /* Term 2*/
            *cross_entropy -= log(1 - nu);

            MATRIX(edge_seen, from, to) = MATRIX(edge_seen, to, from) = 1;
        }
    }

    igraph_matrix_destroy(&edge_seen);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
#endif /* UMAP_DEBUG */


/* clip forces to avoid too rapid shifts */
static igraph_real_t igraph_i_umap_clip_force(igraph_real_t force, igraph_real_t limit) {
    return force > limit ? limit : (force < -limit ? -limit : force);
}

static igraph_real_t igraph_i_umap_attract(
        igraph_real_t dsq,
        igraph_real_t a,
        igraph_real_t b)
{
    return - (2 * a * b * pow(dsq, b - 1.)) / (1. + a * pow(dsq, b));
}

static igraph_real_t igraph_i_umap_repel(
        igraph_real_t dsq,
        igraph_real_t a,
        igraph_real_t b)
{
    igraph_real_t dsq_min = UMAP_CORRECT_DISTANCE_REPULSION * UMAP_CORRECT_DISTANCE_REPULSION;

    return (2 * b) / (dsq_min + dsq) / (1. + a * pow(dsq, b));
}

static igraph_error_t igraph_i_umap_apply_forces(
        const igraph_t *graph,
        const igraph_vector_t *umap_weights,
        igraph_matrix_t *layout,
        igraph_real_t a,
        igraph_real_t b,
        igraph_real_t learning_rate,
        igraph_bool_t avoid_neighbor_repulsion,
        igraph_integer_t negative_sampling_rate,
        igraph_integer_t epoch,
        igraph_vector_t *next_epoch_sample_per_edge)
{
    igraph_integer_t no_of_vertices = igraph_matrix_nrow(layout);
    igraph_integer_t ndim = igraph_matrix_ncol(layout);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t from, to, nneis, eid;
    igraph_vector_t from_emb, to_emb, delta;
    igraph_real_t force = 0, dsq, force_d;
    /* The following is only used for small graphs, to avoid repelling your neighbors
     * For large sparse graphs, it's not necessary. For large dense graphs, you should
     * not be doing UMAP.
     * */
    igraph_vector_int_t neis, negative_vertices;
    igraph_integer_t n_negative_vertices = (no_of_vertices - 1 < negative_sampling_rate) ? (no_of_vertices - 1) : negative_sampling_rate;

    /* Initialize vectors */
    IGRAPH_VECTOR_INIT_FINALLY(&from_emb, ndim);
    IGRAPH_VECTOR_INIT_FINALLY(&to_emb, ndim);
    IGRAPH_VECTOR_INIT_FINALLY(&delta, ndim);

    if (avoid_neighbor_repulsion) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    }
    IGRAPH_VECTOR_INT_INIT_FINALLY(&negative_vertices, 0);

    /* Iterate over edges. Stronger edges are sampled more often */
    for (eid = 0; eid < no_of_edges; eid++) {
        /* Zero-weight edges do not affect vertex positions. They can
         * also emerge during the weight symmetrization. */
        if (VECTOR(*umap_weights)[eid] <= 0) {
            continue;
        }

        /* We sample all and only edges that are supposed to be moved at this time */
        if ((VECTOR(*next_epoch_sample_per_edge)[eid] - epoch) >= 1) {
            continue;
        }

        /* set next epoch at which this edge will be sampled */
        VECTOR(*next_epoch_sample_per_edge)[eid] += 1.0 / VECTOR(*umap_weights)[eid];

        /* we move all vertices on one end of the edges, then we come back for
         * the vertices on the other end. This way we don't move both ends at the
         * same time, which is almost a wasted move since they attract each other */
        int swapflag = (int)(RNG_UNIF01() > 0.5);
        int swapflag_end = swapflag + 2;
        for (; swapflag < swapflag_end; swapflag++) {

            /* half the time, swap the from/to, otherwise some vertices are never moved.
             * This has to do with the graph representation within igraph */
            if (swapflag % 2) {
                from = IGRAPH_FROM(graph, eid);
                to = IGRAPH_TO(graph, eid);
            } else {
                to = IGRAPH_FROM(graph, eid);
                from = IGRAPH_TO(graph, eid);
            }


            /* Current coordinates of both vertices */
            dsq = 0;
            for (igraph_integer_t d = 0; d != ndim; d++) {
                VECTOR(from_emb)[d] = MATRIX(*layout, from, d);
                VECTOR(to_emb)[d] = MATRIX(*layout, to, d);
                VECTOR(delta)[d] = MATRIX(*layout, from, d) - MATRIX(*layout, to, d);
                dsq += VECTOR(delta)[d] * VECTOR(delta)[d];
            }

            /* Apply attractive force since they are neighbors */
            /* NOTE: If they are already together, no force needed */
            if (dsq >= UMAP_MIN_DISTANCE_ATTRACTION * UMAP_MIN_DISTANCE_ATTRACTION) {
                force = igraph_i_umap_attract(dsq, a, b);
                for (igraph_integer_t d = 0; d != ndim; d++) {
                    force_d = force * VECTOR(delta)[d];
                    /* clip force to avoid too rapid change */
                    force_d = igraph_i_umap_clip_force(force_d, UMAP_FORCE_LIMIT);

            #ifdef UMAP_DEBUG
                    fprintf(stderr, "force attractive: delta[%ld] = %g, forces[%ld] = %g\n", d, VECTOR(delta)[d], d, force_d);
            #endif

                    MATRIX(*layout, from, d) += learning_rate * force_d;
                }
            }

            /* Random other nodes repel the focal vertex */
            IGRAPH_CHECK(igraph_random_sample(&negative_vertices,
                        0, no_of_vertices - 2, n_negative_vertices));
            for (igraph_integer_t j = 0; j < n_negative_vertices; j++) {

                IGRAPH_ALLOW_INTERRUPTION();

                /* Get random neighbor */
                to = VECTOR(negative_vertices)[j];
                /* obviously you cannot repel yourself */
                if (to >= from) {
                    to++;
                }

                /* do not repel neighbors for small graphs, for big graphs this
                 * does not matter as long as the k in knn << number of vertices */
                if (avoid_neighbor_repulsion) {
                    /* NOTE: the efficiency of this step could be improved but it
                     * should be only used for small graphs anyway, so it's fine */
                    igraph_bool_t skip = 0;
                    IGRAPH_CHECK(igraph_incident(graph, &neis, from, IGRAPH_ALL));
                    nneis = igraph_vector_int_size(&neis);
                    for (igraph_integer_t k = 0; k < nneis; k++) {
                        igraph_integer_t eid2 = VECTOR(neis)[k];
                        igraph_integer_t from2, to2;
                        from2 = IGRAPH_FROM(graph, eid2);
                        to2 = IGRAPH_TO(graph, eid2);
                        if (((from2 == from) && (to2 == to)) || ((from2 == to) && (from == to2))) {
                            skip = 1;
                            break;
                        }
                    }
                    if (skip == 1) {
                        continue;
                    }
                }

                /* Get layout of random neighbor and gradient in embedding */
                dsq = 0;
                for (igraph_integer_t d = 0; d != ndim; d++) {
                    VECTOR(to_emb)[d] = MATRIX(*layout, to, d);
                    VECTOR(delta)[d] = MATRIX(*layout, from, d) - MATRIX(*layout, to, d);
                    dsq += VECTOR(delta)[d] * VECTOR(delta)[d];
                }

                /* This repels the other vertex assuming it's a negative example
                 * that is no weight, no edge */
                force = igraph_i_umap_repel(dsq, a, b);
                /* The repulsive force is already *away* from the other (non-neighbor) vertex */
                for (igraph_integer_t d = 0; d != ndim; d++) {
                    force_d = force * VECTOR(delta)[d];

                    /* clip force to avoid too rapid change */
                    force_d = igraph_i_umap_clip_force(force_d, UMAP_FORCE_LIMIT);

                #ifdef UMAP_DEBUG
                    fprintf(stderr, "force repulsive: delta[%ld] = %g, forces[%ld] = %g\n", d, VECTOR(delta)[d], d, force_d);
                #endif

                    MATRIX(*layout, from, d) += learning_rate * force_d;
                }
            }
        }
    }

    /* Free vectors */
    igraph_vector_int_destroy(&negative_vertices);
    igraph_vector_destroy(&from_emb);
    igraph_vector_destroy(&to_emb);
    igraph_vector_destroy(&delta);
    IGRAPH_FINALLY_CLEAN(4);

    /* Free vector of neighbors if needed */
    if (avoid_neighbor_repulsion) {
        igraph_vector_int_destroy(&neis);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/* Edges with heavier weight/higher probability should be sampled more often. In
 * other words, vertices at each end of those edges should be moved more often. If the
 * edge weight is 1.0, which happens to each nearest neighbor due to the correction via
 * rho, that vertices at the end of that edge are moved each single epoch. Conversely,
 * vertices at the end of weak edges can be moved only once in a while. */
static igraph_error_t igraph_i_umap_optimize_layout_stochastic_gradient(
        const igraph_t *graph,
        const igraph_vector_t *umap_weights,
        igraph_real_t a,
        igraph_real_t b,
        igraph_matrix_t *layout,
        igraph_integer_t epochs,
        igraph_integer_t negative_sampling_rate) {

    igraph_real_t learning_rate = 1;
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vector_t next_epoch_sample_per_edge;

#ifdef UMAP_DEBUG
    igraph_real_t cross_entropy, cross_entropy_old;
#endif

    IGRAPH_VECTOR_INIT_FINALLY(&next_epoch_sample_per_edge, no_of_edges);

    /* Explicit avoidance of neighbor repulsion, only useful in small graphs
     * which are never very sparse. This is because negative sampling as implemented
     * relies on an approximation that only works if the graph is sparse, which is never
     * quite true for small graphs (i.e. |V| << |E| << |V|^2 is hard to judge if
     * |V| is small) */
    igraph_bool_t avoid_neighbor_repulsion = 0;
    if (igraph_vcount(graph) < 100) {
        avoid_neighbor_repulsion = 1;
    }

    /* Measure the (variable part of the) cross-entropy terms for debugging:
     * 1. - sum_edge_e mu(e) * log(nu(e))
     * 2. + sum_edge_e (1 - mu(e)) * log(1 - nu(e))
     * The latter is approximated by negative sampling as:
     * 2b. + sum_random_ij 1 * log(1 - nu_ij)
     * whereby the mu = 0 because we assume there's no edge between i and j, and nu_ij
     * is basically their distance in embedding space, lensed through the probability
     * function Phi.
     * */
#ifdef UMAP_DEBUG
    igraph_umap_compute_cross_entropy(
            graph, umap_weights, layout, a, b, &cross_entropy);
#endif

    for (igraph_integer_t e = 0; e < epochs; e++) {
        /* Apply (stochastic) forces */
        igraph_i_umap_apply_forces(
                graph,
                umap_weights,
                layout,
                a, b,
                learning_rate,
                avoid_neighbor_repulsion,
                negative_sampling_rate,
                e,
                &next_epoch_sample_per_edge);

#ifdef UMAP_DEBUG
        /* Recompute CE and check how it's going*/
        cross_entropy_old = cross_entropy;
        igraph_umap_compute_cross_entropy(
                graph, umap_weights, layout, a, b, &cross_entropy);

        printf("Cross-entropy before shift: %g, after shift: %g\n", cross_entropy_old, cross_entropy);
#endif

         /* Adjust learning rate */
        learning_rate = 1.0 - (igraph_real_t)(e + 1) / epochs;
    }

    igraph_vector_destroy(&next_epoch_sample_per_edge);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* Center layout around (0,0) at the end, just for convenience */
static igraph_error_t igraph_i_umap_center_layout(igraph_matrix_t *layout) {
    igraph_integer_t no_of_vertices = igraph_matrix_nrow(layout);
    igraph_real_t xm = 0, ym = 0;

    /* Compute center */
    xm = 0;
    ym = 0;
    for (igraph_integer_t i = 0; i < no_of_vertices; i++) {
        xm += MATRIX(*layout, i, 0);
        ym += MATRIX(*layout, i, 1);
    }
    xm /= no_of_vertices;
    ym /= no_of_vertices;

    /* Shift vertices */
    for (igraph_integer_t i = 0; i < no_of_vertices; i++) {
        MATRIX(*layout, i, 0) -= xm;
        MATRIX(*layout, i, 1) -= ym;
    }

    return IGRAPH_SUCCESS;
}


/* This is the main function that works for any dimensionality of the embedding
 * (currently hard-constrained to 2 or 3 ONLY in the initialization). */
static igraph_error_t igraph_i_layout_umap(
        const igraph_t *graph,
        igraph_matrix_t *res,
        igraph_bool_t use_seed,
        const igraph_vector_t *distances,
        igraph_real_t min_dist,
        igraph_integer_t epochs,
        igraph_integer_t ndim,
        igraph_bool_t distances_are_weights) {

    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_vertices = igraph_vcount(graph);
    /* probabilities of each edge being a real connection */
    igraph_vector_t weights;
    igraph_vector_t *weightsp;
    /* The smoothing parameters given min_dist */
    igraph_real_t a, b;
    /* How many repulsions for each attraction */
    igraph_integer_t negative_sampling_rate = 5;

    /* Check input arguments */
    if (min_dist < 0) {
        IGRAPH_ERRORF("Minimum distance must not be negative, got %g.",
                IGRAPH_EINVAL, min_dist);
    }

    if (epochs < 0) {
        IGRAPH_ERRORF("Number of epochs must be non-negative, got %" IGRAPH_PRId ".",
                IGRAPH_EINVAL, epochs);
    }

    if ((ndim != 2) && (ndim != 3)) {
        IGRAPH_ERRORF("Number of dimensions must be 2 or 3, got %" IGRAPH_PRId ".",
                IGRAPH_EINVAL, ndim);

    }

    /* Compute weights (exponential weights) from distances if required.
     * If the weights have already been computed, they are stored in
     * the "distances" vector and we can recycle the pointer. */
    if (distances_are_weights) {
        weightsp = (igraph_vector_t *) distances;
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(&weights, no_of_edges);
        IGRAPH_CHECK(igraph_layout_umap_compute_weights(
                    graph, distances, &weights));
        weightsp = &weights;
    }
    /* From now on everything lives in probability space, it does not matter whether
     * the original graph was weighted/distanced or unweighted */

    /* Compute initial layout if required. If a seed layout is used, then just
     * check that the dimensions of the layout make sense. */
    if (use_seed) {
        if ((igraph_matrix_nrow(res) != no_of_vertices) || (igraph_matrix_ncol(res) != ndim)) {
            if (!distances_are_weights) {
                igraph_vector_destroy(&weights);
                IGRAPH_FINALLY_CLEAN(1);
            }
            IGRAPH_ERRORF("Seed layout should have %" IGRAPH_PRId " points in %" IGRAPH_PRId " dimensions, got %" IGRAPH_PRId " points in %" IGRAPH_PRId " dimensions.",
                          IGRAPH_EINVAL, no_of_vertices, ndim,
                          igraph_matrix_nrow(res),
                          igraph_matrix_ncol(res));
        }

        /* Trivial graphs (0 or 1 nodes) with seed - do nothing */
        if (no_of_vertices <= 1) {
            if (!distances_are_weights) {
                igraph_vector_destroy(&weights);
                IGRAPH_FINALLY_CLEAN(1);
            }
            return IGRAPH_SUCCESS;
        }
    } else {
        /* Trivial graphs (0 or 1 nodes) beget trivial - but valid - layouts */
        if (no_of_vertices <= 1) {
            IGRAPH_CHECK(igraph_matrix_resize(res, no_of_vertices, ndim));
            igraph_matrix_null(res);
            if (!distances_are_weights) {
                igraph_vector_destroy(&weights);
                IGRAPH_FINALLY_CLEAN(1);
            }
            return IGRAPH_SUCCESS;
        }

        /* Skip spectral embedding for now (see #1971), initialize at random */
        if (ndim == 2) {
            igraph_layout_random(graph, res);
        } else {
            igraph_layout_random_3d(graph, res);
        }
    }

    RNG_BEGIN();

    /* Fit a and b parameter to find smooth approximation to
     * probability distribution in embedding space */
    IGRAPH_CHECK(igraph_i_umap_fit_ab(min_dist, &a, &b));

    /* Minimize cross-entropy between high-d and low-d probability
     * distributions */
    IGRAPH_CHECK(igraph_i_umap_optimize_layout_stochastic_gradient(
                graph,
                weightsp,
                a, b,
                res,
                epochs,
                negative_sampling_rate));

    if (!distances_are_weights) {
        igraph_vector_destroy(&weights);
        IGRAPH_FINALLY_CLEAN(1);
    }
    RNG_END();

    /* Center layout */
    IGRAPH_CHECK(igraph_i_umap_center_layout(res));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_layout_umap
 * \brief Layout using Uniform Manifold Approximation and Projection (UMAP).
 *
 * \experimental
 *
 * UMAP is mostly used to embed high-dimensional vectors in a low-dimensional space
 * (most commonly 2D). The algorithm is probabilistic and introduces nonlinearities,
 * unlike e.g. PCA and similar to T-distributed Stochastic Neighbor Embedding (t-SNE).
 * Nonlinearity helps "cluster" very similar vectors together without imposing a
 * global geometry on the embedded space (e.g. a rigid rotation + compression in PCA).
 * UMAP uses graphs as intermediate data structures, hence it can be used as a
 * graph layout algorithm as well.
 *
 * </para><para>
 *
 * The general UMAP workflow is to start from vectors, compute a sparse distance
 * graph that only contains edges between simiar points (e.g. a k-nearest neighbors
 * graph), and then convert these distances into exponentially decaying weights
 * between 0 and 1 that are larger for points that are closest neighbors in the
 * distance graph. If a graph without any distances associated to the edges is used,
 * all weights will be set to 1.
 *
 * </para><para>
 *
 * If you are trying to use this function to embed high-dimensional vectors, you should
 * first compute a k-nearest neighbors graph between your vectors and compute the
 * associated distances, and then call this function on that graph. If you already
 * have a distance graph, or you have a graph with no distances, you can call this
 * function directly. If you already have a graph with meaningful weights
 * associated to each edge, you can also call this function, but set the argument
 * distances_are_weights to true. To compute weights from distances
 * without computing the layout, see \ref igraph_layout_umap_compute_weights().
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Leland McInnes, John Healy, and James Melville. https://arxiv.org/abs/1802.03426

 * \param graph Pointer to the graph to find a layout for (i.e. to embed). This is
 *   typically a sparse graph with only edges for the shortest distances stored, e.g.
 *   a k-nearest neighbors graph.
 * \param res Pointer to the n by 2 matrix where the layout coordinates will be stored.
 * \param use_seed Logical, if true the supplied values in the \p res argument are used
 *   as an initial layout, if false a random initial layout is used.
 * \param distances Pointer to a vector of distances associated with the graph edges.
 *   If this argument is \c NULL, all weights will be set to 1.
 * \param min_dist A fudge parameter that decides how close two unconnected vertices
 *   can be in the embedding before feeling a repulsive force. It must not be
 *   negative. Typical values are between 0 and 1.
 * \param epochs Number of iterations of the main stochastic gradient descent loop on
 *   the cross-entropy. Typical values are between 30 and 500.
 * \param distances_are_weights Whether to use precomputed weights. If
 *   true, the "distances" vector contains precomputed weights. If false (the
 *   typical use case), this function will compute weights from distances and
 *   then use them to compute the layout.
 * \return Error code.
 */
igraph_error_t igraph_layout_umap(const igraph_t *graph,
                                  igraph_matrix_t *res,
                                  igraph_bool_t use_seed,
                                  const igraph_vector_t *distances,
                                  igraph_real_t min_dist,
                                  igraph_integer_t epochs,
                                  igraph_bool_t distances_are_weights) {
    return igraph_i_layout_umap(graph, res, use_seed,
            distances, min_dist, epochs, 2, distances_are_weights);
}


/**
 * \function igraph_layout_umap_3d
 * \brief 3D layout using UMAP.
 *
 * \experimental
 *
 * </para><para>
 *
 * This is the 3D version of the UMAP algorithm
 * (see \ref igraph_layout_umap() for the 2D version).
 *
 * \param graph Pointer to the graph to find a layout for (i.e. to embed). This is
 *   typically a directed, sparse graph with only edges for the shortest distances
 *   stored, e.g. a k-nearest neighbors graph with the edges going from each focal
 *   vertex to its neighbors. However, it can also be an undirected graph. If the
 *   \p distances_are_weights is \c true, this is treated as an undirected graph.
 * \param res Pointer to the n by 3 matrix where the layout coordinates will be stored.
 * \param use_seed Logical, if true the supplied values in the \p res argument are used
 *   as an initial layout, if false a random initial layout is used.
 * \param distances Pointer to a vector of distances associated with the graph edges.
 *   If this argument is \c NULL, all edges are assumed to have the same distance.
 * \param min_dist A fudge parameter that decides how close two unconnected vertices
 *   can be in the embedding before feeling a repulsive force. It must not be
 *   negative. Typical values are between 0 and 1.
 * \param epochs Number of iterations of the main stochastic gradient descent loop on
 *   the cross-entropy. Typical values are between 30 and 500.
 * \param distances_are_weights Whether to use precomputed weights. If \c false (the
 *   typical use case), this function will compute weights from distances and
 *   then use them to compute the layout.  If \c true, the "distances" vector contains
 *   precomputed weights, including possibly some weights equal to zero that are
 *   inconsequential for the layout optimization.
 * \return Error code.
 */
igraph_error_t igraph_layout_umap_3d(const igraph_t *graph,
                                     igraph_matrix_t *res,
                                     igraph_bool_t use_seed,
                                     const igraph_vector_t *distances,
                                     igraph_real_t min_dist,
                                     igraph_integer_t epochs,
                                     igraph_bool_t distances_are_weights) {
    return igraph_i_layout_umap(graph, res, use_seed,
            distances, min_dist, epochs, 3, distances_are_weights);
}
