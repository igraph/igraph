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

#include "igraph_matrix.h"
#include "igraph_interface.h"
#include "igraph_umap.h"
#include "igraph_constructors.h"
#include "igraph_layout.h"
#include "igraph_random.h"
#include "igraph_lapack.h"

/* rho is just the size of the distance from each vertex and its closest neighbor */
/* sigma is the the decay from each vertex, depends on its rho and the rest of its neighbor distances */

/* Find sigma for this vertex by binary search */
igraph_error_t igraph_umap_find_sigma(const igraph_t *graph, const igraph_vector_t *distances, igraph_integer_t i, igraph_vector_int_t *eids, igraph_real_t rho, igraph_real_t *sigma_p, igraph_real_t target) {

    igraph_real_t sigma = 1;
    igraph_real_t sigma_new, sum;
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

        printf("SIGMA function (i = %" IGRAPH_PRId ", no_of_neis = %" IGRAPH_PRId ")- sum: %f, target: %f, rho: %f, sigma: %f\n", i, no_of_neis, sum, target, rho, sigma);

        /* Adjust sigma FIXME: this is strictly a little optimistic? */
        if (sum < target) {
            if (seen_max == 1) {
                step /= 2;
            /* first iteration we want to increase by sigma, else double */
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
igraph_error_t igraph_umap_find_prob_graph(const igraph_t *graph, const igraph_vector_t *distances, igraph_vector_t *umap_weights) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_neis, eid;
    igraph_vector_int_t eids, weight_seen;
    igraph_real_t rho, rho_new, sigma, weight, weight_inv, sigma_target;

    /* Minimal distance to neighbor for each vertex, used to warp the scale of the embedding */
    igraph_vector_t rhos;
    /* Decay scale for each vertex, computed from its rho */
    igraph_vector_t sigmas;

    IGRAPH_VECTOR_INIT_FINALLY(&rhos, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&sigmas, no_of_nodes);

    /* Initialize vectors and matrices */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eids, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&weight_seen, igraph_vector_size(distances));

    /* Iterate over vertices x, like in the paper */
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        /* Edges into this vertex */
        igraph_incident(graph, &eids, i, IGRAPH_ALL);
        no_of_neis = igraph_vector_int_size(&eids);

        /* Find rho for this vertex, i.e. the minimal non-self distance */
        rho = VECTOR(*distances)[VECTOR(eids)[0]];
        for (igraph_integer_t j = 1; j < no_of_neis; j++) {
            rho_new = VECTOR(*distances)[VECTOR(eids)[j]];
            rho = fmin(rho, rho_new);
        }
        VECTOR(rhos)[i] = rho;

        /* Find sigma for this vertex, from its rho */
        sigma_target = log(no_of_neis) / log(2);
        IGRAPH_CHECK(igraph_umap_find_sigma(graph, distances, i, &eids, rho, &sigma, sigma_target));
        VECTOR(sigmas)[i] = sigma;

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
            weight = exp(-fmax(0, (VECTOR(*distances)[eid] - rho) / sigma));

            /* Compute the probability of either edge direction if you can */
            if (VECTOR(weight_seen)[eid] != 0) {
                weight_inv = VECTOR(*umap_weights)[eid];
                weight = weight + weight_inv - weight * weight_inv;
            }

            printf("distance: %f\n", VECTOR(*distances)[eid]);
            printf("weight: %f\n", weight);
            VECTOR(*umap_weights)[eid] = weight;
            VECTOR(weight_seen)[eid] += 1;
        }

    }

    igraph_vector_destroy(&rhos);
    igraph_vector_destroy(&sigmas);
    igraph_vector_int_destroy(&weight_seen);
    igraph_vector_int_destroy(&eids);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}


igraph_error_t igraph_umap_get_ab_residuals(igraph_vector_t *residuals, igraph_real_t *squared_sum_res, igraph_integer_t nr_points, igraph_real_t a, igraph_real_t b, igraph_vector_t *powb, igraph_vector_t *x, igraph_real_t min_dist)
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

/* FIXME this fnuction should be static after we made sure it works */
igraph_error_t igraph_umap_fit_ab(igraph_real_t min_dist, float *a_p, float *b_p)
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

    printf("start fit_ab\n");
    for (int iter = 0; iter < maxiter; iter++) {
        igraph_umap_get_ab_residuals(&residuals, &squared_sum_res, nr_points, a, b, &powb, &x, min_dist);

        /* break if good fit (conergence to truth) */
        if (squared_sum_res < tol * tol) {
            printf("convergence to zero (wow!)\n");
            break;
        }
        /* break if no change (convergence) */
        if ((iter > 0) && fabs(sqrt(squared_sum_res_old) - sqrt(squared_sum_res)) < tol) {
            printf("no-change absolute convergence\n");
            break;
        }

        /* Jacobian (first derivatives) of squared residuals at (a, b) */
        for (int i = 0; i < nr_points; i++) {
            tmp = 1 + a * VECTOR(powb)[i];
            MATRIX(jacobian, i, 0) = - 2 * VECTOR(powb)[i] / tmp / tmp;
            MATRIX(jacobian, i, 1) = MATRIX(jacobian, i, 0) * a * log(VECTOR(x)[i]) * 2;
        }

        /* At each iteration, we want to minimize the linear approximation of the sum of squared residuals:
         *
         * sum_i (Ji @ d(a,b) -r_i)^2
         *
         * Putting the first derivative to zero results in a linear system of 2 equations (for a and b):
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
        igraph_umap_get_ab_residuals(&residuals, &squared_sum_res, nr_points, a + da, b + db, &powb, &x, min_dist);

        printf("start line search, SSR before delta: %f, current SSR:, %f\n", squared_sum_res_old, squared_sum_res);
        for (int k = 0; k < 30; k++) {
            /* Try new parameters */
            da /= 2.0;
            db /= 2.0;
            squared_sum_res_tmp = squared_sum_res;
            igraph_umap_get_ab_residuals(&residuals, &squared_sum_res, nr_points, a + da, b + db, &powb, &x, min_dist);

            /* Compare and if we are going back uphill, undo last step and break */
            printf("during line search, k = %d, old SSR:, %f, new SSR (half a,b):, %f\n", k, squared_sum_res_tmp, squared_sum_res);
            if (squared_sum_res > squared_sum_res_tmp - tol) {
                da *= 2;
                db *= 2;
                break;
            }
        }
        printf("end of line search and iteration, squared_sum_res: %f \n\n", squared_sum_res_tmp);

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

    printf("a, b: %f %f\n", a, b);
    *a_p = a;
    *b_p = b;

    return IGRAPH_SUCCESS;

}

/* cross-entropy and derivatives */
static igraph_error_t igraph_compute_cross_entropy(igraph_t *umap_graph, igraph_vector_t *umap_weights, igraph_matrix_t *layout, igraph_real_t a, igraph_real_t b, igraph_real_t *cross_entropy) {

    igraph_real_t mu, nu, xd, yd, sqd;
    igraph_integer_t from, to;
    igraph_integer_t no_of_edges = igraph_ecount(umap_graph);

    /* Measure the (variable part of the) cross-entropy terms for debugging:
     * 1. - sum_edge_e mu(e) * log(nu(e))
     * 2. + sum_edge_e (1 - mu(e)) * log(1 - nu(e))
     * */
    *cross_entropy = 0;
    for (igraph_integer_t eid = 0; eid < no_of_edges; eid++) {
        mu = VECTOR(*umap_weights)[eid];

        /* Find vertices */
        IGRAPH_CHECK(igraph_edge(umap_graph, eid, &from, &to));
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
    }

    return IGRAPH_SUCCESS;
}

/*xd is difference in x direction, mu is a weight */
/* NOTE: mu is the probability of a true edge in high dimensions, not affected
 * by the embedding (in particular, xd and yd), so it's a constant for the
 * derivative/force. Same applies for the repulsion */
static igraph_error_t igraph_attract(igraph_real_t xd, igraph_real_t yd, igraph_real_t mu, igraph_real_t a, igraph_real_t b, igraph_real_t *force_x, igraph_real_t *force_y)
{
    igraph_real_t dsq;
    igraph_real_t force;
    igraph_real_t epsilon = 0.001;

    dsq = xd * xd + yd * yd;

    force = mu * (- 2 * a * b * pow(dsq, b - 1)) / (1 + dsq);
    *force_x = force * xd;
    *force_y = force * yd;

    //printf("force attractive: xd = %f, fx = %f\n", xd, *force_x);
    //printf("force attractive: yd = %f, fy = %f\n", yd, *force_y);

    return IGRAPH_SUCCESS;
}

/*xd is difference in x direction, mu is a weight */
static igraph_error_t igraph_repulse(igraph_real_t xd, igraph_real_t yd, igraph_real_t mu, igraph_real_t a, igraph_real_t b, igraph_real_t *force_x, igraph_real_t *force_y)
{
    igraph_real_t dsq;
    igraph_real_t force, phi;
    igraph_real_t epsilon = 0.001;

    dsq = xd * xd + yd * yd;

    /* NOTE: in practice, in negative sampling mu is always zero because we
     * *assume* the sample to be negative i.e. never a true edge */
    // FIXME: probably the force was already correct, but just in case
    //force = (1 - mu) * (-2 * b) / ((epsilon + dsq) * (1 + a * pow(dsq, b)));
    //phi = 1.0 / (1.0 + a * pow(dsq, b));
    //force = (1 - mu) * (2 * a * b * pow(dsq, b - 1)) * phi * phi / (1 - phi);
    force = (1 - mu) * (2 * b) / ((epsilon + dsq) * (1 + a * pow(dsq, b)));
    // FIXME: OMG, this makes the cross-entropy descent monotonic!!
    //force = 0.2;
    *force_x = force * xd;
    *force_y = force * yd;

    //printf("force repulsive: xd = %f, fx = %f\n", xd, *force_x);
    //printf("force repulsive: yd = %f, fy = %f\n", yd, *force_y);


    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_get_gradient(igraph_matrix_t *gradient, igraph_matrix_t *layout, igraph_t *umap_graph, igraph_vector_t *umap_weights, igraph_real_t a, igraph_real_t b, igraph_real_t prob)
{
    igraph_integer_t no_of_nodes = igraph_matrix_nrow(layout);
    igraph_integer_t no_of_edges = igraph_ecount(umap_graph);
    igraph_integer_t from, to;
    igraph_real_t from_x, from_y, to_x, to_y, x_diff, y_diff, fx, fy;

    // FIXME
    igraph_vector_int_t neis;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

    /* TODO: what should we use for the number of random vertices?
     * might/should change at every iteration?
     * */
    igraph_integer_t n_random_vertices = sqrt(no_of_nodes);

    /* Zero the gradient */
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        MATRIX(*gradient, i, 0) = 0;
        MATRIX(*gradient, i, 1) = 0;
    }

    /* iterate over a random subsample of edges */
    for (igraph_integer_t eid = 0; eid < no_of_edges; eid++) {
        if (RNG_UNIF01() > prob) {
            continue;
        }

        IGRAPH_CHECK(igraph_edge(umap_graph, eid, &from, &to));

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
        MATRIX(*gradient, from, 0) += fx;
        MATRIX(*gradient, from, 1) += fy;

        /* Random other nodes are repelled from one (the first) vertex */
        for (igraph_integer_t j = 0; j < n_random_vertices; j++) {
            /* Get random neighbor, obviously you cannot repel yourself */
            to = RNG_INTEGER(0, no_of_nodes - 2);
            if (to >= from) {
                to++;
            }
            /* FIXME: for now, avoid repelling neighbors */
            /* This is terribly inefficient, but just for testing anyway */
            int skip = 0;
            igraph_incident(umap_graph, &neis, from, IGRAPH_ALL);
            for (int k = 0; k < igraph_vector_int_size(&neis); k++) {
                igraph_integer_t eid2 = VECTOR(neis)[k];
                igraph_integer_t from2, to2;
                igraph_edge(umap_graph, eid2, &from2, &to2);
                if (((from2 == from) && (to2 == to)) || ((from2 == to) && (from == to2))) {
                    skip = 1;
                    break;
                }
            }
            if (skip == 1) {
                continue;
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
            MATRIX(*gradient, from, 0) += fx;
            MATRIX(*gradient, from, 1) += fy;
        }
    }

    igraph_vector_int_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_umap_optimize_layout_stochastic_gradient(igraph_t *umap_graph, igraph_vector_t *umap_weights, igraph_real_t a, igraph_real_t b,
        igraph_matrix_t *layout) {
    igraph_integer_t epochs = 50; //FIXME
    igraph_real_t learning_rate = 1;
    igraph_real_t sampling_prob = 0.2; // between 0 and 1, fraction of edges sampled for gradient at each epoch
    igraph_matrix_t gradient;
    igraph_real_t cross_entropy, cross_entropy_old;


    /* Initialize gradient */
    IGRAPH_MATRIX_INIT_FINALLY(&gradient, igraph_matrix_nrow(layout), igraph_matrix_ncol(layout));

    /* Measure the (variable part of the) cross-entropy terms for debugging:
     * 1. - sum_edge_e mu(e) * log(nu(e))
     * 2. + sum_edge_e (1 - mu(e)) * log(1 - nu(e))
     * The latter is approximated by negative sampling as:
     * 2b. + sum_random_ij 1 * log(1 - nu_ij)
     * whereby the mu = 0 because we assume there's no edge between i and j, and nu_ij
     * is basically their distance in embedding space, lensed through the probability
     * function Phi.
     * */
    igraph_compute_cross_entropy(umap_graph, umap_weights, layout, a, b, &cross_entropy);
    for (igraph_integer_t e = 0; e < epochs; e++) {
        /* Compute (stochastic) gradient */
        igraph_get_gradient(&gradient, layout, umap_graph, umap_weights, a, b, sampling_prob);

        /* Delta is the gradient times the current (decreasing) learning rate */
        igraph_matrix_scale(&gradient, learning_rate / 10);

        //printf("gradient:\n");
        //igraph_matrix_print(&gradient);
        //printf("\n");

        /* Shift the layout by the delta: the forces already have the right directions,
         * i.e. towards neighbors and away from strangers */
        igraph_matrix_add(layout, &gradient);

        /* Recompute CE and check how it's going*/
        cross_entropy_old = cross_entropy;
        igraph_compute_cross_entropy(umap_graph, umap_weights, layout, a, b, &cross_entropy);

        printf("Cross-entropy before shift: %f, after shift: %f\n", cross_entropy_old, cross_entropy);

        //printf("layout:\n");
        //igraph_matrix_print(layout);
        //printf("\n");

         /* Adjust learning rate */
        learning_rate = 1.0 - (e + 1) / epochs;
    }

    igraph_matrix_destroy(&gradient);
    IGRAPH_FINALLY_CLEAN(1);
    return (IGRAPH_SUCCESS);
}


igraph_error_t igraph_layout_umap(igraph_t *graph, igraph_vector_t *distances, igraph_matrix_t *layout) {


    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vector_t umap_weights; /* probabilities of each edge being a real connection */
    igraph_real_t min_dist = 0.01; /* This is empyrical */
    float a, b; /* The smoothing parameters given min_dist */

    if (igraph_vector_size(distances) != no_of_edges) {
        IGRAPH_ERROR("Distances must be the same number as the edges in the graph", IGRAPH_EINVAL);
    }

    RNG_BEGIN();
    IGRAPH_VECTOR_INIT_FINALLY(&umap_weights, no_of_edges);

    /* Make combined graph with smoothed probabilities */
    IGRAPH_CHECK(igraph_umap_find_prob_graph(graph, distances, &umap_weights));

    /* Skip spectral embedding for now, initialize at random */
    igraph_layout_random(graph, layout);

    /* Definition 11 */
    igraph_umap_fit_ab(min_dist, &a, &b);

    /* Algorithm 5 */
    IGRAPH_CHECK(igraph_umap_optimize_layout_stochastic_gradient(graph, &umap_weights, a, b, layout));

    igraph_vector_destroy(&umap_weights);
    IGRAPH_FINALLY_CLEAN(1);
    RNG_END();
    return (IGRAPH_SUCCESS);
}
