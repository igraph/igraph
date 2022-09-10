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
#include "igraph_random.h"
#include "igraph_nongraph.h"

#include "layout/layout_internal.h"

#include <math.h>


/* rho is just the size of the distance from each vertex and its closest neighbor */
/* sigma is the the decay from each vertex, depends on its rho and the rest of its neighbor
 * distances */

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


/* Convert the graph with distances into a probability graph with exponential decay */
/* NOTE: this is funny because the distance was a correlation to start with... ?? */
static igraph_error_t igraph_i_umap_find_prob_graph(const igraph_t *graph,
        const igraph_vector_t *distances, igraph_vector_t *umap_weights) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_neis, eid;
    igraph_vector_int_t eids, weight_seen;
    igraph_real_t rho, dist_max, dist, sigma, weight, weight_inv, sigma_target;

    /* if the original graph is unweighted, probabilities are 1 throughout */
    if (distances == NULL) {
        for (igraph_integer_t j = 0; j < no_of_edges; j++) {
            VECTOR(*umap_weights)[j] = 1;
        }
        return IGRAPH_SUCCESS;
    }
    /* alright, the graph is weighted */

    /* Initialize vectors and matrices */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eids, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&weight_seen, no_of_edges);

    /* Iterate over vertices x, like in the paper */
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        /* Edges into this vertex */
        IGRAPH_CHECK(igraph_incident(graph, &eids, i, IGRAPH_ALL));
        no_of_neis = igraph_vector_int_size(&eids);

        /* Vertex has no neighbors */
        if (no_of_neis == 0) {
            continue;
        }

        /* Find rho for this vertex, i.e. the minimal non-self distance */
        rho = VECTOR(*distances)[VECTOR(eids)[0]];
        dist_max = rho;
        for (igraph_integer_t j = 1; j < no_of_neis; j++) {
            dist = VECTOR(*distances)[VECTOR(eids)[j]];
            rho = fmin(rho, dist);
            dist_max = fmax(dist_max, dist);
        }

        /* If the maximal distance is rho, all neighbors are identical to
         * each other. */
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

        /* Convert to umap_weight
         * Each edge is seen twice, from each of its two vertices. Because this weight
         * is a probability and the probability of the two vertices to be close are set
         * as the probability of either "edge direction" being legit, the final weight
         * is:
         *
         * P_{-> | <-} = P_{->} + P_{<-} - P_{->} * P_{<-}
         *
         * to avoid double counting.
         *
         * To implement that, we keep track of whether we have "seen" this edge before.
         * If no, set the P_{->}. If yes, it contains P_{->} and now we can substitute
         * it with the final result.
         *
         * */
        for (igraph_integer_t j = 1; j < no_of_neis; j++) {
            eid = VECTOR(eids)[j];
            /* Basically, nodes closer than rho have probability 1, but nothing disappears */
            weight = sigma < 0 ? 1: exp(-(VECTOR(*distances)[eid] - rho) / sigma);

            /* Compute the probability of either edge direction if you can */
            if (VECTOR(weight_seen)[eid] != 0) {
                weight_inv = VECTOR(*umap_weights)[eid];
                weight = weight + weight_inv - weight * weight_inv;
            }

#ifdef UMAP_DEBUG
            printf("distance: %g\n", VECTOR(*distances)[eid]);
            printf("weight: %g\n", weight);
#endif
            VECTOR(*umap_weights)[eid] = weight;
            VECTOR(weight_seen)[eid] += 1;
        }

    }

    igraph_vector_int_destroy(&weight_seen);
    igraph_vector_int_destroy(&eids);
    IGRAPH_FINALLY_CLEAN(2);

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
            igraph_vector_destroy(&x);
            igraph_vector_destroy(&residuals);
            igraph_matrix_destroy(&jacobian);
            igraph_matrix_destroy(&jTj);
            igraph_matrix_destroy(&jTr);
            igraph_vector_destroy(&powb);
            IGRAPH_FINALLY_CLEAN(6);
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
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_matrix_t edge_seen;

    IGRAPH_MATRIX_INIT_FINALLY(&edge_seen, no_of_nodes, no_of_nodes);

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
        /* NOT 2 * b since it's already squared */
        nu = 1.0 / (1 + a * pow(sqd, b));

        /* Term 1*/
        *cross_entropy -= mu * log(nu);
        /* Term 2*/
        *cross_entropy -= (1 - mu) * log(1 - nu);

        MATRIX(edge_seen, from, to) = MATRIX(edge_seen, to, from) = 1;
    }
    /* Add the entropy from the missing edges */
    for (igraph_integer_t from = 0; from < no_of_nodes; from++) {
        for (igraph_integer_t to = 0; to < from; to++) {
            if (MATRIX(edge_seen, from, to) > 0) {
                continue;
            }

            /* Find distance in layout space */
            xd = (MATRIX(*layout, from, 0) - MATRIX(*layout, to, 0));
            yd = (MATRIX(*layout, from, 1) - MATRIX(*layout, to, 1));
            sqd = xd * xd + yd * yd;

            /* Find probability associated with distance using fitted Phi */
            /* NOT 2 * b since it's already squared */
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
static igraph_error_t igraph_i_umap_clip_force(igraph_real_t *force, igraph_real_t limit) {

    *force  = fmax(fmin(*force, limit), -limit);

    return IGRAPH_SUCCESS;
}

/*xd is difference in x direction, mu is a weight */
/* NOTE: mu is the probability of a true edge in high dimensions, not affected
 * by the embedding (in particular, xd and yd), so it's a constant for the
 * derivative/force. Same applies for the repulsion */
static igraph_error_t igraph_i_umap_attract(
        igraph_vector_t *delta, igraph_real_t mu,
        igraph_real_t a, igraph_real_t b, igraph_vector_t *forces)
{
    igraph_real_t dsq, phi, force;
    igraph_integer_t ndim = igraph_vector_size(delta);

    dsq = 0;
    for (igraph_integer_t d = 0; d != ndim; d++) {
        dsq += VECTOR(*delta)[d] * VECTOR(*delta)[d];
    }

    phi = 1. / (1. + a * pow(dsq, b));
    force = - mu * (2 * a * b * pow(dsq, b - 1)) * phi;
    for (igraph_integer_t d = 0; d != ndim; d++) {
        VECTOR(*forces)[d] = force * VECTOR(*delta)[d];
        /* clip force to avoid too rapid change */
        igraph_i_umap_clip_force(&(VECTOR(*forces)[d]), 3);

#ifdef UMAP_DEBUG
        printf("force attractive: delta[%d] = %g, forces[%d] = %g\n", d, VECTOR(*delta)[d], d, VECTOR(*forces)[d]);
#endif
    }

    return IGRAPH_SUCCESS;
}

/*xd is difference in x direction, mu is a weight */
static igraph_error_t igraph_i_umap_repel(
        igraph_vector_t *delta, igraph_real_t mu,
        igraph_real_t a, igraph_real_t b, igraph_vector_t *forces)
{
    igraph_real_t dsq, force;
    igraph_real_t min_dist = 0.01;
    igraph_integer_t ndim = igraph_vector_size(delta);

    dsq = 0;
    for (igraph_integer_t d = 0; d != ndim; d++) {
        dsq += VECTOR(*delta)[d] * VECTOR(*delta)[d];
    }
    dsq = fmax(min_dist * min_dist, dsq);

    /* NOTE: in practice, in negative sampling mu is always zero because we
     * *assume* the sample to be negative i.e. never a true edge */
    force = (1 - mu) * (2 * b) / dsq / (1 + a * pow(dsq, b));
    for (igraph_integer_t d = 0; d != ndim; d++) {
        VECTOR(*forces)[d] = force * VECTOR(*delta)[d];

        /* clip force to avoid too rapid change */
        igraph_i_umap_clip_force(&(VECTOR(*forces)[d]), 3);

#ifdef UMAP_DEBUG
        printf("force repulsive: delta[%d] = %g, forces[%d] = %g\n", d, VECTOR(*delta)[d], d, VECTOR(*forces)[d]);
#endif
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_umap_apply_forces(
        const igraph_t *graph,  const igraph_vector_t *umap_weights,
        igraph_matrix_t *layout, igraph_real_t a, igraph_real_t b, igraph_real_t prob,
        igraph_real_t learning_rate, igraph_bool_t avoid_neighbor_repulsion)
{
    igraph_integer_t no_of_nodes = igraph_matrix_nrow(layout);
    igraph_integer_t ndim = igraph_matrix_ncol(layout);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t from, to, nneis;
    igraph_vector_t from_emb, to_emb, delta, forces;
    /* The following is only used for small graphs, to avoid repelling your neighbors
     * For large sparse graphs, it's not necessary. For large dense graphs, you should
     * not be doing UMAP.
     * */
    igraph_vector_int_t neis, negative_vertices;

    /* Initialize vectors */
    IGRAPH_VECTOR_INIT_FINALLY(&from_emb, ndim);
    IGRAPH_VECTOR_INIT_FINALLY(&to_emb, ndim);
    IGRAPH_VECTOR_INIT_FINALLY(&delta, ndim);
    IGRAPH_VECTOR_INIT_FINALLY(&forces, ndim);

    if (avoid_neighbor_repulsion) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    }
    IGRAPH_VECTOR_INT_INIT_FINALLY(&negative_vertices, 0);

    igraph_integer_t n_random_vertices = sqrt(no_of_nodes);

    /* iterate over a random subsample of edges. We have an igraph_random_sample() function to
     * do that, but that one would require the allocation of a separate vector and we don't
     * need that, especially if `prob` is high */
    /* TODO: check possible improvement with sampling from a geometric distribution if prob is
     * known to be small; it would require the generation of only ceil(no_of_edges * prob)
     * random numbers and not no_of_edges */
    for (igraph_integer_t eid = 0; eid < no_of_edges; eid++) {
        if (RNG_UNIF01() > prob) {
            continue;
        }

        /* half the time, swap the from/to, otherwise some vertices are never moved */
        if (RNG_UNIF01() > 0.5) {
            from = IGRAPH_FROM(graph, eid);
            to = IGRAPH_TO(graph, eid);
        } else {
            to = IGRAPH_FROM(graph, eid);
            from = IGRAPH_TO(graph, eid);
        }

        /* Current coordinates of both vertices */
        for (igraph_integer_t d = 0; d != ndim; d++) {
            VECTOR(from_emb)[d] = MATRIX(*layout, from, d);
            VECTOR(to_emb)[d] = MATRIX(*layout, to, d);
            VECTOR(delta)[d] = MATRIX(*layout, from, d) - MATRIX(*layout, to, d);
        }

        /* Apply attractive force since they are neighbors */
        IGRAPH_CHECK(igraph_i_umap_attract(&delta, VECTOR(*umap_weights)[eid], a, b, &forces));
        for (igraph_integer_t d = 0; d != ndim; d++) {
            MATRIX(*layout, from, d) += learning_rate * VECTOR(forces)[d];
        }

        /* Random other nodes are repelled from one (the first) vertex */
        IGRAPH_CHECK(igraph_random_sample(&negative_vertices, 0, no_of_nodes - 2, n_random_vertices));
        for (igraph_integer_t j = 0; j < n_random_vertices; j++) {
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
                igraph_bool_t skip = false;
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
            for (igraph_integer_t d = 0; d != ndim; d++) {
                VECTOR(to_emb)[d] = MATRIX(*layout, to, d);
                VECTOR(delta)[d] = MATRIX(*layout, from, d) - MATRIX(*layout, to, d);
            }

            /* This repels the other vertex assuming it's a negative example
             * that is no weight, no edge */
            IGRAPH_CHECK(igraph_i_umap_repel(&delta, 0, a, b, &forces));
            /* The repulsive force is already *away* from the other (non-neighbor) vertex */
            for (igraph_integer_t d = 0; d != ndim; d++) {
                MATRIX(*layout, from, d) += learning_rate * VECTOR(forces)[d];
            }
        }
    }

    /* Free vectors */
    igraph_vector_int_destroy(&negative_vertices);
    igraph_vector_destroy(&from_emb);
    igraph_vector_destroy(&to_emb);
    igraph_vector_destroy(&delta);
    igraph_vector_destroy(&forces);
    IGRAPH_FINALLY_CLEAN(5);

    /* Free vector of neighbors if needed */
    if (avoid_neighbor_repulsion) {
        igraph_vector_int_destroy(&neis);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_umap_optimize_layout_stochastic_gradient(const igraph_t *graph,
       const igraph_vector_t *umap_weights, igraph_real_t a, igraph_real_t b,
       igraph_matrix_t *layout, igraph_integer_t epochs, igraph_real_t sampling_prob) {

    igraph_real_t learning_rate = 1;

#ifdef UMAP_DEBUG
    igraph_real_t cross_entropy, cross_entropy_old;
#endif

    /* Explicit avoidance of neighbor repulsion, only useful in small graphs
     * which are never very sparse. This is because negative sampling as implemented
     * relies on an approximation that only works if the graph is sparse, which is never
     * quite true for small graphs (i.e. |V| << |E| << |V|^2 is hard to judge if
     * |V| is small) */
    igraph_bool_t avoid_neighbor_repulsion = false;
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
    igraph_i_umap_compute_cross_entropy(graph, umap_weights, layout, a, b, &cross_entropy);
#endif

    for (igraph_integer_t e = 0; e < epochs; e++) {
        /* Apply (stochastic) forces */
        igraph_i_umap_apply_forces(graph, umap_weights, layout, a, b, sampling_prob, learning_rate,
                avoid_neighbor_repulsion);

#ifdef UMAP_DEBUG
        /* Recompute CE and check how it's going*/
        cross_entropy_old = cross_entropy;
        igraph_i_umap_compute_cross_entropy(graph, umap_weights, layout, a, b, &cross_entropy);

        printf("Cross-entropy before shift: %g, after shift: %g\n", cross_entropy_old, cross_entropy);
#endif

         /* Adjust learning rate */
        learning_rate = 1.0 - (igraph_real_t)(e + 1) / epochs;
    }

    return IGRAPH_SUCCESS;
}

/* Center layout around (0,0) at the end, just for convenience */
static igraph_error_t igraph_i_umap_center_layout(igraph_matrix_t *layout) {
    igraph_integer_t no_of_nodes = igraph_matrix_nrow(layout);
    igraph_real_t xm = 0, ym = 0;

    /* Compute center */
    xm = 0;
    ym = 0;
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        xm += MATRIX(*layout, i, 0);
        ym += MATRIX(*layout, i, 1);
    }
    xm /= no_of_nodes;
    ym /= no_of_nodes;

    /* Shift vertices */
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        MATRIX(*layout, i, 0) -= xm;
        MATRIX(*layout, i, 1) -= ym;
    }

    return IGRAPH_SUCCESS;
}


/* Check the distances argument, which should be NULL or a vector of nonnegative numbers */
static igraph_error_t igraph_i_umap_check_distances(const igraph_vector_t *distances, igraph_integer_t no_of_edges) {

    if (distances == NULL) {
        return IGRAPH_SUCCESS;
    }

    if (igraph_vector_size(distances) != no_of_edges) {
        IGRAPH_ERROR("Distances must be the same number as the edges in the graph.", IGRAPH_EINVAL);
    }

    for (igraph_integer_t eid = 0; eid != no_of_edges; eid++) {
        if (VECTOR(*distances)[eid] < 0) {
            IGRAPH_ERROR("Distances cannot be negative.", IGRAPH_EINVAL);
        } else if (igraph_is_nan(VECTOR(*distances)[eid])) {
            IGRAPH_ERROR("Distances cannot contain NaN values.", IGRAPH_EINVAL);
        }
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
        igraph_real_t sampling_prob,
        igraph_integer_t ndim) {

    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    /* probabilities of each edge being a real connection */
    igraph_vector_t umap_weights;
    /* The smoothing parameters given min_dist */
    igraph_real_t a, b;

    /* Check input arguments */
    if (min_dist <= 0) {
        IGRAPH_ERRORF("Minimum distance must be positive, got %g.",
                IGRAPH_EINVAL, min_dist);
    }

    if (epochs < 0) {
        IGRAPH_ERRORF("Number of epochs must be non-negative, got %" IGRAPH_PRId ".",
                IGRAPH_EINVAL, epochs);
    }

    if ((sampling_prob <= 0) || (sampling_prob > 1)) {
        IGRAPH_ERRORF("Sampling probability must be in (0, 1], got %g.",
                IGRAPH_EINVAL, sampling_prob);
    }

    if ((ndim != 2) && (ndim != 3)) {
        IGRAPH_ERRORF("Number of dimensions must be 2 or 3, got %" IGRAPH_PRId ".",
                IGRAPH_EINVAL, ndim);

    }

    /* UMAP is sometimes used on unweighted graphs, that means distances are always zero */
    IGRAPH_CHECK(igraph_i_umap_check_distances(distances, no_of_edges));

    if (use_seed) {
        if ((igraph_matrix_nrow(res) != no_of_nodes) || (igraph_matrix_ncol(res) != ndim)) {
            IGRAPH_ERRORF("Seed layout should have %" IGRAPH_PRId " points in %" IGRAPH_PRId " dimensions, got %" IGRAPH_PRId " points in %" IGRAPH_PRId " dimensions.",
                          IGRAPH_EINVAL, no_of_nodes, ndim,
                          igraph_matrix_nrow(res),
                          igraph_matrix_ncol(res));
        }

        /* Trivial graphs (0 or 1 nodes) with seed - do nothing */
        if (no_of_nodes <= 1) {
            return IGRAPH_SUCCESS;
        }
    } else {
         /* Trivial graphs (0 or 1 nodes) beget trivial - but valid - layouts */
         if (no_of_nodes <= 1) {
             IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, ndim));
             igraph_matrix_null(res);
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
    IGRAPH_VECTOR_INIT_FINALLY(&umap_weights, no_of_edges);

    /* Make combined graph with smoothed probabilities */
    IGRAPH_CHECK(igraph_i_umap_find_prob_graph(graph, distances, &umap_weights));

    /* From now on everything lives in probability space, it does not matter whether
     * the original graph was weighted/distanced or unweighted */

    /* Fit a and b parameter to find smooth approximation to
     * probability distribution in embedding space */
    IGRAPH_CHECK(igraph_i_umap_fit_ab(min_dist, &a, &b));

    /* Minimize cross-entropy between high-d and low-d probability
     * distributions */
    IGRAPH_CHECK(igraph_i_umap_optimize_layout_stochastic_gradient(graph, &umap_weights, a, b,
                res, epochs, sampling_prob));

    igraph_vector_destroy(&umap_weights);
    IGRAPH_FINALLY_CLEAN(1);
    RNG_END();

    /* Center layout */
    IGRAPH_CHECK(igraph_i_umap_center_layout(res));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_layout_umap
 * \brief Layout using Uniform Manifold Approximation and Projection for Dimension Reduction.
 *
 * \experimental
 *
 * UMAP is mostly used to embed high-dimensional vectors in a low-dimensional space
 * (most commonly by far, 2D). The algorithm is probabilistic and introduces
 * nonlinearities, unlike e.g. PCA and similar to T-distributed Stochastic Neighbor
 * Embedding (t-SNE). Nonlinearity helps "cluster" very similar vectors together without
 * imposing a global geometry on the embedded space (e.g. a rigid rotation + compression
 * in PCA).
 *
 * </para><para>
 *
 * However, UMAP uses a graph with distances associated to the edges as a key
 * intermediate representation of the high-dimensional space, so it is also useful as
 * a general graph layouting algorithm, hence its inclusion in igraph.
 *
 * </para><para>
 *
 * Importantly, the edge-associated distances are derived from a similarity metric
 * between the high-dimensional vectors, often Pearson correlation:
 *
 * </para><para>
 *
 * <code>corr(v1, v2) = v1 . v2 / [ sqrt(v1 . v1) * sqrt(v2 . v2) ]</code>,
 *
 * </para><para>
 *
 * where <code>.</code> denotes the dot product.
 * In this case, the associated distance is usually defined as:
 *
 * </para><para>
 *
 * <code>d(v1, v2) = 1 - corr(v1, v2)</code>
 *
 * </para><para>
 *
 * This implementation can also work with unweighted similarity graphs, in which case
 * the distance parameter should be a null pointer and all edges beget a similarity
 * score of 1 (a distance of 0).
 *
 * </para><para>
 *
 * While all similarity graphs are theoretically embeddable, UMAP's stochastic gradient
 * descent approach really shines when the graph is sparse. In practice, most people
 * feed a k-nearest neighbor (either computed exactly or approximated) similarity graph
 * with some additional cutoff to exclude "quasi-neighbors" that lie beyond a certain
 * distance (e.g. correlation less than 0.2).
 *
 * </para><para>
 *
 * Therefore, if you are trying to use this function to embed high-dimensional vectors,
 * the steps are:
 *
 * </para><para>
 *
 * 1. Compute a sparse similarity graph (either exact or approximate) from your vectors,
 *    weighted or unweighted. If unsure, compute a k-nearest neighbors graph.
 *
 * </para><para>
 *
 * 2. If you keep the weights, convert them into distances or store them as a "weight"
 *    edge attribute and use a null pointer for the distances. If using similarity
 *    weights instead of distances, make sure they do not exceed 1.
 *
 * </para><para>
 *
 * 3. Feed the graph (and distances, if you have them) into this function.
 *
 * </para><para>
 *
 * Note: Step 1 above involves deciding if two high-dimensional vectors "look similar"
 *       which, because of the curse of dimensionality, is in many cases a highly
 *       subjective and potentially controversial operation: thread with care and at
 *       your own risk. Two high-dimensional vectors might look similar or extremely
 *       different depending on the point of view/angle, and there are a lot of
 *       viewpoints when the dimensionality ramps up.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Leland McInnes, John Healy, and James Melville. https://arxiv.org/abs/1802.03426
 *
 * \param graph Pointer to the similarity graph to find a layout for (i.e. to embed).
 * \param res Pointer to the n by 2 matrix where the layout coordinates will be stored.
 * \param use_seed Logical, if true the supplied values in the \p res argument are used
 *   as an initial layout, if false a random initial layout is used.
 * \param distances Pointer to a vector of edge lengths. Similarity graphs for
 *   UMAP are often originally meant in terms of similarity weights (e.g. correlation between
 *   high-dimensional vectors) and converted into distances by crude <code>dist = 1 - corr</code>.
 *   That is fine here too. If this argument \c NULL, all lengths are assumed to be the same.
 * \param min_dist A fudge parameter that decides how close two unconnected vertices can be in the
 *   embedding before feeling a repulsive force. It should be positive. Typically, 0.01 is a good
 *   number.
 * \param epochs Number of iterations of the main stochastic gradient descent loop on the
 *   cross-entropy. Usually, 500 epochs can be used if the graph is the graph is small
 *   (less than 50000 edges), 50 epochs are used for larger graphs.
 * \param sampling_prob The fraction of vertices moved at each iteration of the stochastic gradient
 *   descent (epoch). At fixed number of epochs, a higher fraction makes the algorithm slower.
 *   Vice versa, a too low number will converge very slowly, possibly too slowly.
 *
 * \return Error code.
 */
igraph_error_t igraph_layout_umap(const igraph_t *graph,
                                  igraph_matrix_t *res,
                                  igraph_bool_t use_seed,
                                  const igraph_vector_t *distances,
                                  igraph_real_t min_dist,
                                  igraph_integer_t epochs,
                                  igraph_real_t sampling_prob) {
    return igraph_i_layout_umap(graph, res, use_seed,
            distances, min_dist, epochs, sampling_prob, 2);
}


/**
 * \function igraph_layout_umap_3d
 * \brief 3D layout using UMAP.
 *
 * \experimental
 *
 * This is the 3D version of the UMAP algorithm
 * (see \ref igraph_layout_umap() for the 2D version).
 *
 * \param graph Pointer to the similarity graph to find a layout for (i.e. to embed).
 * \param res Pointer to the n by 3 matrix where the layout coordinates will be stored.
 * \param use_seed Logical, if true the supplied values in the \p res argument are used
 *   as an initial layout, if false a random initial layout is used.
 * \param distances Pointer to a vector of edge lengths. Similarity graphs for
 *   UMAP are often originally meant in terms of similarity weights (e.g. correlation between
 *   high-dimensional vectors) and converted into distances by crude <code>dist = 1 - corr</code>.
 *   That is fine here too. If this argument is \c NULL, all lengths are assumed to be the same.
 * \param min_dist A fudge parameter that decides how close two unconnected vertices can be in the
 *   embedding before feeling a repulsive force. It should be positive. Typically, 0.01 is a good
 *   number.
 * \param epochs Number of iterations of the main stochastic gradient descent loop on the
 *   cross-entropy. Usually, 500 epochs can be used if the graph is the graph is small
 *   (less than 50000 edges), 50 epochs are used for larger graphs.
 * \param sampling_prob The fraction of vertices moved at each iteration of the stochastic gradient
 *   descent (epoch). At fixed number of epochs, a higher fraction makes the algorithm slower.
 *   Vice versa, a too low number will converge very slowly, possibly too slowly.
 *
 * \return Error code.
 */
igraph_error_t igraph_layout_umap_3d(const igraph_t *graph,
                                     igraph_matrix_t *res,
                                     igraph_bool_t use_seed,
                                     const igraph_vector_t *distances,
                                     igraph_real_t min_dist,
                                     igraph_integer_t epochs,
                                     igraph_real_t sampling_prob) {
    return igraph_i_layout_umap(graph, res, use_seed,
            distances, min_dist, epochs, sampling_prob, 2);
}
