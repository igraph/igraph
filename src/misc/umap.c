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

#include "igraph_umap.h"

#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_lapack.h"
#include "igraph_layout.h"
#include "igraph_matrix.h"
#include "igraph_random.h"


/* rho is just the size of the distance from each vertex and its closest neighbor */
/* sigma is the the decay from each vertex, depends on its rho and the rest of its neighbor
 * distances */

/* Find sigma for this vertex by binary search */
static igraph_error_t igraph_i_umap_find_sigma(const igraph_t *graph,
        const igraph_vector_t *distances, igraph_integer_t i, const igraph_vector_int_t *eids,
        igraph_real_t rho, igraph_real_t *sigma_p, igraph_real_t target) {

    igraph_real_t sigma = 1;
    igraph_real_t sum;
    igraph_real_t tol = 0.01;
    igraph_integer_t maxiter = 100;
    igraph_integer_t no_of_neis = igraph_vector_int_size(eids);
    igraph_integer_t eid;
    igraph_real_t step = sigma;
    int seen_max = 0;

    /* Binary search */
    for(int iter = 0; iter < maxiter; iter++) {
        sum = 0;
        for(int j = 0; j < no_of_neis; j++) {
            eid = VECTOR(*eids)[j];
            sum += exp(-(VECTOR(*distances)[eid] - rho) / sigma);
        }

#ifdef UMAP_DEBUG
        printf("SIGMA function (i = %" IGRAPH_PRId ", no_of_neis = %" IGRAPH_PRId ")- sum: %f, "
               "target: %f, rho: %f, sigma: %f\n", i, no_of_neis, sum, target, rho, sigma);
#endif

        /* TODO: this seems fine, but is probably a little off for some corner cases */
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

    return (IGRAPH_SUCCESS);
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
        no_of_neis = igraph_vector_size(umap_weights);
        for (int j = 0; j < no_of_edges; j++) {
            VECTOR(*umap_weights)[j] = 1;
        }
        return (IGRAPH_SUCCESS);
    }
    /* alright, the graph is weighted */

    /* Initialize vectors and matrices */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eids, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&weight_seen, no_of_edges);

    /* Iterate over vertices x, like in the paper */
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        /* Edges into this vertex */
        igraph_incident(graph, &eids, i, IGRAPH_ALL);
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
            sigma_target = log(no_of_neis) / log(2);
            IGRAPH_CHECK(igraph_i_umap_find_sigma(graph, distances, i, &eids, rho, &sigma,
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
            printf("distance: %f\n", VECTOR(*distances)[eid]);
            printf("weight: %f\n", weight);
#endif
            VECTOR(*umap_weights)[eid] = weight;
            VECTOR(weight_seen)[eid] += 1;
        }

    }

    igraph_vector_int_destroy(&weight_seen);
    igraph_vector_int_destroy(&eids);
    IGRAPH_FINALLY_CLEAN(2);

    return (IGRAPH_SUCCESS);
}

/* Helper function to compute a and b parameters (smoothing probability metric in embedding space) */
static igraph_error_t igraph_i_umap_get_ab_residuals(igraph_vector_t *residuals,
        igraph_real_t *squared_sum_res, igraph_integer_t nr_points, igraph_real_t a,
        igraph_real_t b, const igraph_vector_t *powb, const igraph_vector_t *x, igraph_real_t min_dist)
{
    igraph_real_t tmp;

    *squared_sum_res = 0;
    for (int i = 0; i < nr_points; i++) {
        VECTOR(*powb)[i] = powf(VECTOR(*x)[i], 2 * b);
        tmp = 1 / (1 + a * VECTOR(*powb)[i]);
        tmp -= VECTOR(*x)[i] <= min_dist ? 1 : exp(-(VECTOR(*x)[i] - min_dist));
        VECTOR(*residuals)[i] = tmp;
        *squared_sum_res += tmp * tmp;
    }
    return (IGRAPH_SUCCESS);
}

#ifndef UMAP_DEBUG
static
#endif
igraph_error_t igraph_i_umap_fit_ab(igraph_real_t min_dist, igraph_real_t *a_p, igraph_real_t *b_p)
{
    /*We're fitting a and b, such that
     * (1 + a*d^2b)^-1
     * has the minimum least-squares error compared to
     * the fuzzy ball with min-dist ball size and
     * plain e^-d fuzzyness
     * */

    /* Grid points */
    igraph_vector_t x;
     /* Make a lattice from 0 to this distance */
    igraph_integer_t nr_points = 100;
    igraph_real_t end_point = min_dist*30; /* Need to sample decently around the kink I guess? */
    /* Initial values: a^-2b is the point where f(x) = 0.5, so around min_dist
     * Therefore, initial a is 1/sqrt(min_dist) */
    igraph_real_t b = 1.;
    igraph_real_t a = 1. / sqrt(min_dist);
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
    for (int i = 0; i < nr_points; i++) {
        VECTOR(x)[i] = (end_point / nr_points) * i + 0.001; /* added a 0.001 to prevent NaNs */
    }

#ifdef UMAP_DEBUG
    printf("start fit_ab\n");
#endif
    for (int iter = 0; iter < maxiter; iter++) {
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
        for (int i = 0; i < nr_points; i++) {
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
        for (int i = 0; i < nr_points; i++) {
            for (int j1 = 0; j1 < 2; j1++) {
                for (int j2 = 0; j2 < 2; j2++) {
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
        printf("start line search, SSR before delta: %f, current SSR:, %f\n", squared_sum_res_old,
                squared_sum_res);
#endif
        for (int k = 0; k < 30; k++) {
            /* Try new parameters */
            da /= 2.0;
            db /= 2.0;
            squared_sum_res_tmp = squared_sum_res;
            IGRAPH_CHECK(igraph_i_umap_get_ab_residuals(&residuals, &squared_sum_res, nr_points,
                        a + da, b + db, &powb, &x, min_dist));

            /* Compare and if we are going back uphill, undo last step and break */
#ifdef UMAP_DEBUG
            printf("during line search, k = %d, old SSR:, %f, new SSR (half a,b):, %f\n", k,
                    squared_sum_res_tmp, squared_sum_res);
#endif
            if (squared_sum_res > squared_sum_res_tmp - tol) {
                da *= 2;
                db *= 2;
                break;
            }
        }
#ifdef UMAP_DEBUG
        printf("end of line search and iteration, squared_sum_res: %f \n\n", squared_sum_res_tmp);
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
    printf("a, b: %f %f\n", a, b);
#endif

    *a_p = a;
    *b_p = b;

    return (IGRAPH_SUCCESS);

}

/* cross-entropy */
static igraph_error_t igraph_compute_cross_entropy(const igraph_t *graph,
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
        IGRAPH_CHECK(igraph_edge(graph, eid, &from, &to));
        /* Find distance in layout space */
        xd = (MATRIX(*layout, from, 0) - MATRIX(*layout, to, 0));
        yd = (MATRIX(*layout, from, 1) - MATRIX(*layout, to, 1));
        sqd = xd * xd + yd * yd;
        /* Find probability associated with distance using fitted Phi */
        /* NOT 2 * b since it's already squared */
        nu = 1.0 / (1 + a * powf(sqd, b));

        /* Term 1*/
        *cross_entropy -= mu * log(nu);
        /* Term 2*/
        *cross_entropy -= (1 - mu) * log(1 - nu);

        MATRIX(edge_seen, from, to) = MATRIX(edge_seen, to, from) = 1;
    }
    /* Add the entropy from the missing edges */
    for (int from = 0; from < no_of_nodes; from++) {
        for (int to = 0; to < from; to++) {
            if (MATRIX(edge_seen, from, to) > 0) {
                continue;
            }

            /* Find distance in layout space */
            xd = (MATRIX(*layout, from, 0) - MATRIX(*layout, to, 0));
            yd = (MATRIX(*layout, from, 1) - MATRIX(*layout, to, 1));
            sqd = xd * xd + yd * yd;

            /* Find probability associated with distance using fitted Phi */
            /* NOT 2 * b since it's already squared */
            nu = 1.0 / (1 + a * powf(sqd, b));

            /* Term 2*/
            *cross_entropy -= log(1 - nu);

            MATRIX(edge_seen, from, to) = MATRIX(edge_seen, to, from) = 1;
        }
    }

    igraph_matrix_destroy(&edge_seen);
    IGRAPH_FINALLY_CLEAN(1);

    return (IGRAPH_SUCCESS);
}

/*xd is difference in x direction, mu is a weight */
/* NOTE: mu is the probability of a true edge in high dimensions, not affected
 * by the embedding (in particular, xd and yd), so it's a constant for the
 * derivative/force. Same applies for the repulsion */
static igraph_error_t igraph_attract(igraph_real_t xd, igraph_real_t yd, igraph_real_t mu,
       igraph_real_t a, igraph_real_t b, igraph_real_t *force_x, igraph_real_t *force_y)
{
    igraph_real_t dsq, phi, force;

    dsq = xd * xd + yd * yd;

    phi = 1. / (1. + a * pow(dsq, b));
    force = - mu * (2 * a * b * pow(dsq, b - 1)) * phi;
    *force_x = force * xd;
    *force_y = force * yd;

#ifdef UMAP_DEBUG
    printf("force attractive: xd = %f, fx = %f\n", xd, *force_x);
    printf("force attractive: yd = %f, fy = %f\n", yd, *force_y);
#endif

    return (IGRAPH_SUCCESS);
}

/*xd is difference in x direction, mu is a weight */
static igraph_error_t igraph_repulse(igraph_real_t xd, igraph_real_t yd, igraph_real_t mu,
       igraph_real_t a, igraph_real_t b, igraph_real_t *force_x, igraph_real_t *force_y)
{
    igraph_real_t dsq, force;
    igraph_real_t min_dist = 0.01;

    dsq = fmax(min_dist * min_dist, xd * xd + yd * yd);

    /* NOTE: in practice, in negative sampling mu is always zero because we
     * *assume* the sample to be negative i.e. never a true edge */
    force = (1 - mu) * (2 * b) / dsq / (1 + a * pow(dsq, b));
    *force_x = force * xd;
    *force_y = force * yd;

#ifdef UMAP_DEBUG
    printf("force repulsive: xd = %f, fx = %f\n", xd, *force_x);
    printf("force repulsive: yd = %f, fy = %f\n", yd, *force_y);
#endif

    return (IGRAPH_SUCCESS);
}

static igraph_error_t igraph_apply_forces(const igraph_t *graph,  const igraph_vector_t *umap_weights,
       igraph_matrix_t *layout, igraph_real_t a, igraph_real_t b, igraph_real_t prob,
       igraph_real_t learning_rate, int avoid_neighbor_repulsion)
{
    igraph_integer_t no_of_nodes = igraph_matrix_nrow(layout);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t from, to;
    igraph_real_t from_x, from_y, to_x, to_y, x_diff, y_diff, fx, fy;
    /* The following is only used for small graphs, to avoid repelling your neighbors
     * For large sparse graphs, it's not necessary. For large dense graphs, you should
     * not be doing UMAP.
     * */
    igraph_vector_int_t neis;

    if (avoid_neighbor_repulsion > 0) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    }

    igraph_integer_t n_random_vertices = sqrt(no_of_nodes);

    /* iterate over a random subsample of edges */
    for (igraph_integer_t eid = 0; eid < no_of_edges; eid++) {
        if (RNG_UNIF01() > prob) {
            continue;
        }

        IGRAPH_CHECK(igraph_edge(graph, eid, &from, &to));

        /* Current coordinates of both vertices */
        from_x = MATRIX(*layout, from, 0);
        from_y = MATRIX(*layout, from, 1);
        to_x = MATRIX(*layout, to, 0) ;
        to_y = MATRIX(*layout, to, 1) ;

        /* Distance in embedding space */
        x_diff = from_x - to_x;
        y_diff = from_y - to_y;

        /* Apply attractive force since they are neighbors */
        IGRAPH_CHECK(igraph_attract(x_diff, y_diff, VECTOR(*umap_weights)[eid], a, b, &fx, &fy));
        MATRIX(*layout, from, 0) += learning_rate * fx;
        MATRIX(*layout, from, 1) += learning_rate * fy;

        /* Random other nodes are repelled from one (the first) vertex */
        for (igraph_integer_t j = 0; j < n_random_vertices; j++) {
            /* Get random neighbor, obviously you cannot repel yourself */
            to = RNG_INTEGER(0, no_of_nodes - 2);
            if (to >= from) {
                to++;
            }
            /* do not repel neighbors for small graphs, for big graphs this
             * does not matter as long as the k in knn << number of vertices */
            if (avoid_neighbor_repulsion > 0) {
                /* FIXME: This algo terribly inefficient, try to improve using
                 * adjacency lists or something but just for testing anyway */
                int skip = 0;
                igraph_incident(graph, &neis, from, IGRAPH_ALL);
                for (int k = 0; k < igraph_vector_int_size(&neis); k++) {
                    igraph_integer_t eid2 = VECTOR(neis)[k];
                    igraph_integer_t from2, to2;
                    igraph_edge(graph, eid2, &from2, &to2);
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
            to_x = MATRIX(*layout, to, 0) ;
            to_y = MATRIX(*layout, to, 1) ;
            x_diff = from_x - to_x;
            y_diff = from_y - to_y;

            /* This repels the other vertex assuming it's a negative example
             * that is no weight, no edge */
            IGRAPH_CHECK(igraph_repulse(x_diff, y_diff, 0, a, b, &fx, &fy));
            /* The repulsive force is already *away* from the other (non-neighbor) vertex */
            MATRIX(*layout, from, 0) += learning_rate * fx;
            MATRIX(*layout, from, 1) += learning_rate * fy;
        }
    }

    if (avoid_neighbor_repulsion > 0) {
        igraph_vector_int_destroy(&neis);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return (IGRAPH_SUCCESS);
}

static igraph_error_t igraph_i_umap_optimize_layout_stochastic_gradient(const igraph_t *graph,
       const igraph_vector_t *umap_weights, igraph_real_t a, igraph_real_t b,
       igraph_matrix_t *layout, igraph_integer_t epochs, igraph_real_t sampling_prob) {

    igraph_real_t learning_rate = 1;

#ifdef UMAP_DEBUG
    igraph_real_t cross_entropy, cross_entropy_old;
#endif

    /* Explicit avoidance of neighbor repulsion, only useful in small graphs
     * which are never very sparse */
    int avoid_neighbor_repulsion = 0;
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
    igraph_compute_cross_entropy(graph, umap_weights, layout, a, b, &cross_entropy);
#endif

    for (igraph_integer_t e = 0; e < epochs; e++) {
        /* Apply (stochastic) forces */
        igraph_apply_forces(graph, umap_weights, layout, a, b, sampling_prob, learning_rate, avoid_neighbor_repulsion);

#ifdef UMAP_DEBUG
        /* Recompute CE and check how it's going*/
        cross_entropy_old = cross_entropy;
        igraph_compute_cross_entropy(graph, umap_weights, layout, a, b, &cross_entropy);

        printf("Cross-entropy before shift: %f, after shift: %f\n", cross_entropy_old, cross_entropy);
#endif

         /* Adjust learning rate */
        learning_rate = 1.0 - (float)(e + 1) / epochs;
    }

    return (IGRAPH_SUCCESS);
}

/* Center layout around (0,0) at the end, just for convenience */
static igraph_error_t igraph_i_umap_center_layout(igraph_matrix_t *layout) {
    igraph_integer_t no_of_nodes = igraph_matrix_nrow(layout);
    igraph_real_t xm = 0, ym = 0;

    /* Compute center */
    xm = 0;
    ym = 0;
    for (int i = 0; i < no_of_nodes; i++) {
        xm += MATRIX(*layout, i, 0);
        ym += MATRIX(*layout, i, 1);
    }
    xm /= no_of_nodes;
    ym /= no_of_nodes;

    /* Shift vertices */
    for (int i = 0; i < no_of_nodes; i++) {
        MATRIX(*layout, i, 0) -= xm;
        MATRIX(*layout, i, 1) -= ym;
    }

    return (IGRAPH_SUCCESS);
}

static igraph_error_t igraph_i_umap_check_distances(const igraph_vector_t *distances, igraph_integer_t no_of_edges) {

    if (distances == NULL) {
        return (IGRAPH_SUCCESS);
    }

    if (igraph_vector_size(distances) != no_of_edges) {
        IGRAPH_ERROR("Distances must be the same number as the edges in the graph", IGRAPH_EINVAL);
    }

    for (int eid = 0; eid != no_of_edges; eid++) {
        if (VECTOR(*distances)[eid] < 0) {
            IGRAPH_ERROR("Distances cannot be negative", IGRAPH_EINVAL);
        }
    }

    return (IGRAPH_SUCCESS);
}

/**
 * \function igraph_layout_umap
 * \brief Layout using Uniform Manifold Approximation and Projection for Dimension Reduction
 *
 * UMAP is a mostly used to embed high-dimensional vectors in a low-dimensional space
 * (most commonly by far, 2D). The algorithm is probabilistic and introduces
 * nonlinearities, unlike e.g. PCA and similar to T-distributed Stochastic Neighbor
 * Embedding (t-SNE). Nonlinearity helps "cluster" very similar vectors together without
 * imposing a global geometry on the embedded space (e.g. a rigid rotation + compression
 * in PCA).
 *
 * However, UMAP uses a graph with distances associated to the edges as a key
 * intermediate representation of the high-dimensional space, so it is also useful as
 * a general graph layouting algorithm, hence its inclusion in igraph.
 *
 * Importantly, the edge-associated distances are derived from a similarity metric
 * between the high-dimensional vectors, often Pearson correlation:
 *
 * <code>corr(v1, v2) = v1 x v2 / [ sqrt(v1 x v1) * sqrt(v2 x v2) ]</code>
 *
 * In this case, the associated distance is usually defined as:
 *
 * <code>d(v1, v2) = 1 - corr(v1, v2)</code>
 *
 * This implementation can also work with unweighted similarity graphs, in which case
 * the distance parameter should be a null pointer and all edges beget a similarity
 * score of 1 (a distance of 0).
 *
 * While all similarity graphs are theoretically embeddable, UMAP's stochastic gradient
 * descent approach really shines when the graph is sparse. In practice, most people
 * feed a k-nearest neighbor (either computed exactly or approximated) similarity graph
 * with some additional cutoff to exclude "quasi-neighbors" that lie beyond a certain
 * distance (e.g. correlation less than 0.2).
 *
 * Therefore, if you are trying to use this function to embed high-dimensional vectors,
 * the steps are:
 *
 * 1. Compute a sparse similarity graph (either exact or approximate) from your vectors,
 *    weighted or unweighted. If unsure, compute a knn.
 * 2. If you keep the weights, convert them into distances or store them as a "weight"
 *    edge attribute and use a null pointer for the distances. If using similarity
 *    weights instead of distances, make sure they do not exceed 1.
 * 3. Feed the graph (and distances, if you have them) into this function.
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
 * Leland McInnes, John Healy, and James Melville, https://doi.org/10.21105/joss.00861,
 * https://arxiv.org/abs/1802.03426
 *
 * \param graph Pointer to the similarity graph to find a layout for (i.e. to embed). Weights in
 *   this graph object are currently ignored.
 * \param distances Pointer to a vector of "distances" between vertices. Similarity graphs for
 *   UMAP are often originally meant in terms of similarity weights (e.g. correlation between
 *   high-dimensional vectors) and converted into distances by crude dist := 1 - corr. That is
 *   fine here too. If this argument is an empty pointer (NULL), connected vertices are assumed
 *   to be very similar.
 * \param layout Pointer to the n x 2 matrix where the layout coordinates will be stored. Only 2D
 *   embedding are currently supported, as they are by far more common than any higher dimensions.
 * \param min_dist A fudge parameter that decides how close two unconnected vertices can be in the
 *   embedding before feeling a repulsive force. It should be positive. Typically, 0.01 is a good
 *   number.
 * \param epochs Number of iterations of the main stochastic gradient descent loop on the
 *   cross-entropy. Usually, 500 epochs can be used if the graph is the graph is small
 *   (less than 50k edges), 50 epochs are used for larger graphs.
 * \param sampling_prob The fraction of vertices moved at each iteration of the stochastic gradient descent (epoch). At fixed number of epochs, a higher fraction makes the algorithm slower. Vice versa, a too low number will converge very slowly, possibly too slowly. If negative of zero, a value of 0.3 is assumed.
 *
 * \return Error code.
 *
 */
igraph_error_t igraph_layout_umap(const igraph_t *graph, const igraph_vector_t *distances,
        igraph_matrix_t *layout, igraph_real_t min_dist, igraph_integer_t epochs, igraph_real_t sampling_prob) {


    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    /* probabilities of each edge being a real connection */
    igraph_vector_t umap_weights;
    /* The smoothing parameters given min_dist */
    igraph_real_t a, b;

    if (min_dist <= 0) {
        IGRAPH_ERRORF("Minimum distance should be positive, but found %f.",
                IGRAPH_EINVAL, min_dist);
    }

    if (epochs < 0) {
        IGRAPH_ERRORF("Number of epochs should be non-negative, but found %"
                IGRAPH_PRId ".",
                IGRAPH_EINVAL, epochs);
    }

    /* Default sampling_prob */
    if (sampling_prob <= 0) {
        /* FIXME: the stability seems extremely dependent on this */
        sampling_prob = 1.0 / sqrt(no_of_nodes);
    }

    /* UMAP is sometimes used on unweighted graphs, that means distances are always zero */
    IGRAPH_CHECK(igraph_i_umap_check_distances(distances, no_of_edges));

    RNG_BEGIN();
    IGRAPH_VECTOR_INIT_FINALLY(&umap_weights, no_of_edges);

    /* Make combined graph with smoothed probabilities */
    IGRAPH_CHECK(igraph_i_umap_find_prob_graph(graph, distances, &umap_weights));

    /* From now on everything lives in probability space, it does not matter whether
     * the original graph was weighted/distanced or unweighted */

    /* Skip spectral embedding for now, initialize at random */
    igraph_layout_random(graph, layout);

    /* Definition 11 */
    IGRAPH_CHECK(igraph_i_umap_fit_ab(min_dist, &a, &b));

    /* Algorithm 5 */
    IGRAPH_CHECK(igraph_i_umap_optimize_layout_stochastic_gradient(graph, &umap_weights, a, b,
                layout, epochs, sampling_prob));

    igraph_vector_destroy(&umap_weights);
    IGRAPH_FINALLY_CLEAN(1);
    RNG_END();

    /* Center layout */
    IGRAPH_CHECK(igraph_i_umap_center_layout(layout));

    return (IGRAPH_SUCCESS);
}
