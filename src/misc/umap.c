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

/*open set size is just the size of the distance to the closest neighbor*/
/*the decay depends on the rest of the neighbors */

/*open set size is rho in the paper, the decay is sigma*/
static igraph_error_t igraph_umap_find_open_sets(igraph_t *knn_graph,
        igraph_vector_t *distances, igraph_vector_t *open_set_sizes,
        igraph_vector_t *open_set_decays) {
    igraph_integer_t no_of_nodes = igraph_vcount(knn_graph);
    igraph_vector_int_t eids;
    igraph_real_t l2k;
    igraph_integer_t k;
    igraph_real_t sum;

    k = igraph_vector_size(open_set_sizes);
    l2k = log(k) / log(2);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&eids, 0);
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_incident(knn_graph, &eids, i, IGRAPH_ALL);
        VECTOR(*open_set_sizes)[i] = VECTOR(*distances)[VECTOR(eids)[0]];
        sum = 0;
        for (igraph_integer_t j = 1; j < igraph_vector_int_size(&eids); j++) {
            if (VECTOR(*distances)[VECTOR(eids)[j]] < VECTOR(*open_set_sizes)[i]) {
                VECTOR(*open_set_sizes)[i] = VECTOR(*distances)[VECTOR(eids)[j]];
            }
        }
        /* TODO: according to the paper finding the decay (sigma) should be
         * done with a binary search?
         * FIXME: The current approach is just plain incorrect and leads to
         * negative decays */
        for (igraph_integer_t j = 1; j < igraph_vector_int_size(&eids); j++) {
            sum += exp(VECTOR(*open_set_sizes)[i] - VECTOR(*distances)[VECTOR(eids)[j]]);
        }
        VECTOR(*open_set_decays)[i] = log(sum / l2k);
        printf("decay: %f\n", VECTOR(*open_set_decays)[i]);
    }
    igraph_vector_int_destroy(&eids);
    IGRAPH_FINALLY_CLEAN(1);
    return (IGRAPH_SUCCESS);
}


igraph_error_t igraph_umap_decay(igraph_real_t *probability, igraph_real_t distance, igraph_real_t open_set_size, igraph_real_t open_set_decay)
{
    if (distance < open_set_size) {
        *probability = 1.0;
        return (IGRAPH_SUCCESS);
    }
    *probability = exp((open_set_size - distance) / open_set_decay);
    return (IGRAPH_SUCCESS);
}

static igraph_error_t igraph_umap_edge_weights(igraph_t *graph, igraph_vector_t *distances,
        igraph_vector_t *umap_weights, igraph_vector_t *open_set_sizes,
        igraph_vector_t *open_set_decays) {

    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_real_t weight;
    igraph_real_t weight_previous;

    igraph_vector_resize(umap_weights, igraph_vector_size(distances));
    igraph_vector_null(umap_weights);
    for (igraph_integer_t i = 0; i < no_of_edges; i++) {
        igraph_integer_t from = IGRAPH_FROM(graph, i);
        IGRAPH_CHECK(igraph_umap_decay(&weight,  VECTOR(*distances)[i],  VECTOR(*open_set_sizes)[from], VECTOR(*open_set_decays)[from]));
        weight_previous = VECTOR(*umap_weights)[i];
        weight = weight + weight_previous - weight * weight_previous;
        VECTOR(*umap_weights)[i] = weight;
    }

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_get_residuals(igraph_vector_t *residuals, igraph_real_t *squared_sum_res, igraph_integer_t nr_points, igraph_real_t a, igraph_real_t b, igraph_vector_t *powb, igraph_vector_t *x, igraph_real_t min_dist)
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
igraph_error_t igraph_fit_ab(igraph_real_t min_dist, float *a_p, float *b_p)
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
    igraph_vector_int_t ipiv; // Pivot for LAPACK
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
    IGRAPH_VECTOR_INT_INIT_FINALLY(&ipiv, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&powb, nr_points);

    printf("start fit_ab\n");
    /* Distance |x-y| (this is a lattice, there are no actual x and y) */
    for (int i = 0; i < nr_points; i++) {
        VECTOR(x)[i] = (end_point / nr_points) * i + 0.001; /* added a 0.001 to prevent NaNs */
    }

    for (int iter = 0; iter < maxiter; iter++) {
        igraph_get_residuals(&residuals, &squared_sum_res, nr_points, a, b, &powb, &x, min_dist);

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
            /* df/da * delta */
            MATRIX(jacobian, i, 0) = - 2 * VECTOR(powb)[i] / tmp / tmp;
            /* df/db * delta */
            MATRIX(jacobian, i, 1) = MATRIX(jacobian, i, 0) * a * log(2 * VECTOR(x)[i]);
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
        for (int i = 0; i < nr_points; i++) {
            for (int j1 = 0; j1 < 2; j1++) {
                for (int j2 = 0; j2 < 2; j2++) {
                    MATRIX(jTj, j1, j2) += MATRIX(jacobian, i, j1) * MATRIX(jacobian, i, j2);
                }
                MATRIX(jTr, j1, 0) += MATRIX(jacobian, i, j1) * VECTOR(residuals)[i];
            }
        }
        /* LAPACK puts solution into jTr, sometimes with row swapping (stored in ipiv) */
        igraph_lapack_dgesv(&jTj, &ipiv, &jTr, &lapack_info);

        /* This might go wrong, in which case we should fail graciously */
        if (lapack_info > 0) {
            igraph_vector_destroy(&x);
            igraph_vector_destroy(&residuals);
            igraph_matrix_destroy(&jacobian);
            igraph_matrix_destroy(&jTj);
            igraph_matrix_destroy(&jTr);
            igraph_vector_int_destroy(&ipiv);
            igraph_vector_destroy(&powb);
            IGRAPH_FINALLY_CLEAN(7);
            IGRAPH_ERROR("Singular matrix in the estimation of a and b for UMAP",  IGRAPH_EINVAL);
        }

        /* Assign deltas (LAPACK can swap rows of the solution for convenience) */
        if (VECTOR(ipiv)[0] == 1) {
            da = MATRIX(jTr, 1, 0);
            db = MATRIX(jTr, 0, 0);
        } else {
            da = MATRIX(jTr, 0, 0);
            db = MATRIX(jTr, 1, 0);
        }

        /* FIXME: perhaps wrong sign? */
        da *= -1;
        db *= -1;

        /* Improvement over GN: rough exponential line search for best delta
         * start from largest change, and keep shrinking as long as we are going down
         * */
        squared_sum_res_old = squared_sum_res;
        igraph_get_residuals(&residuals, &squared_sum_res, nr_points, a + da, b + db, &powb, &x, min_dist);

        printf("start line search, SSR before delta: %f, current SSR:, %f\n", squared_sum_res_old, squared_sum_res);
        for (int k = 0; k < 10; k++) {
            /* Try new parameters */
            da /= 2.0;
            db /= 2.0;
            squared_sum_res_tmp = squared_sum_res;
            igraph_get_residuals(&residuals, &squared_sum_res, nr_points, a + da, b + db, &powb, &x, min_dist);

            /* Compare and if we are going back uphill, undo last step and break */
            printf("during line search, k = %d, old SSR:, %f, new SSR (half a,b):, %f\n", k, squared_sum_res_old, squared_sum_res);
            if (squared_sum_res > squared_sum_res_tmp) {
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
    igraph_vector_int_destroy(&ipiv);
    igraph_vector_destroy(&powb);
    IGRAPH_FINALLY_CLEAN(7);

    printf("a, b: %f %f\n", a, b);
    *a_p = a;
    *b_p = b;

    return IGRAPH_SUCCESS;

}

/*xd is difference in x direction, w is a weight */
static igraph_error_t igraph_attract(igraph_real_t xd, igraph_real_t yd, igraph_real_t w, igraph_real_t a, igraph_real_t b, igraph_real_t *force_x, igraph_real_t *force_y)
{
    /* the paper says on page 16: "where a and b are hyper-parameters",
     * which I understood to mean that it cannot be inferred from the data.
     * Then at the top of page 21 "a and b are chosen by non-linear least squares fitting...".
     * I think it's the same a and b, so I'm surprised it's not a hyperparameter now.
     * Or doesn't the fit depend on the data? does it only depend on min-dist?
     * They're used to fit the smooth force curve to the flat ball with fuzzy falloff.
     * I also don't understand why we first make the fuzzy open balls and then use
     * the smooth curve, instead of fitting a smooth curve immediately.
     * From algorithm 5 it seems like an and b should not depend on a particular vertex,
     * but either on all the data or none at all, because the fitting should be done only once.
     * I'm guessing it should not depend on the data, but only on the min-dist, which means
     * they are actual hyperparameters. Only for fitting to curves I have to fit
     * them at certain points, or all points, which means doing some integrals
     * I guess? I'll just fit them at some random points
     * TODO:
     * we should do non-linear least squares fitting for a and b */
    igraph_real_t dsq;
    igraph_real_t force;

    dsq = xd * xd + yd * yd;
    force = (- 2 * a * b * pow(dsq, (b - 1)) * w) / (1 + dsq);
    *force_x = force * xd;
    *force_y = force * yd;
    return IGRAPH_SUCCESS;
}

/*xd is difference in x direction, w is a weight */
static igraph_error_t igraph_repulse(igraph_real_t xd, igraph_real_t yd, igraph_real_t w, igraph_real_t a, igraph_real_t b, igraph_real_t *force_x, igraph_real_t *force_y)
{
    igraph_real_t dsq;
    igraph_real_t force;
    igraph_real_t epsilon = 0.001;

    dsq = xd * xd + yd * yd;
    force = ((-2 * b) * (1 - w)) / ((epsilon + dsq) * (1 + a * pow(dsq, b)));
    *force_x = force * xd;
    *force_y = force * yd;
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_get_gradient(igraph_matrix_t *gradient, igraph_matrix_t *layout, igraph_t *umap_graph, igraph_vector_t *umap_weights, igraph_real_t a, igraph_real_t b, igraph_real_t prob)
{
    igraph_integer_t no_of_nodes = igraph_matrix_nrow(layout);
    igraph_integer_t no_of_edges = igraph_ecount(umap_graph);
    igraph_integer_t from, to;
    igraph_real_t from_x, from_y, to_x, to_y, x_diff, y_diff, fx, fy;


    /* TODO: what should we use for the number of random vertices?
     * might/should change at every iteration?
     * */
    igraph_integer_t n_random_verices = sqrt(no_of_nodes);

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
        /* FIXME: Is this correct? Or the other way? */
        MATRIX(*gradient, from, 0) += fx;
        MATRIX(*gradient, from, 1) += fy;

        /* Random other nodes are repelled from one (the first) vertex */
        for (igraph_integer_t j = 0; j < n_random_verices; j++) {
            /* Get random neighbor, obviously you cannot repel yourself */
            to = RNG_INTEGER(0, no_of_nodes - 2);
            if (to >= from) {
                to++;
            }

            /* Get layout of random neighbor and gradient in embedding */
            to_x = MATRIX(*layout, to, 0) ;
            to_y = MATRIX(*layout, to, 1) ;
            x_diff = from_x - to_x;
            y_diff = from_y - to_y;

            /* Apply repulsive force FIXME: check sign */
            /* This repels the other vertex assuming it's a negative example
             * that is no weight, no edge */
            IGRAPH_CHECK(igraph_repulse(x_diff, y_diff, 0, a, b, &fx, &fy));
            MATRIX(*gradient, from, 0) += fx;
            MATRIX(*gradient, from, 1) += fy;
        }
    }
    return IGRAPH_SUCCESS;
}


static igraph_error_t igraph_compute_cross_entropy(igraph_t *umap_graph, igraph_vector_t *umap_weights, igraph_matrix_t *layout, igraph_real_t a, igraph_real_t b, igraph_real_t *cross_entropy) {

    igraph_real_t mu, nu, dx, dy, sqdist;
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
        dx = (MATRIX(*layout, from, 0) - MATRIX(*layout, to, 0));
        dy = (MATRIX(*layout, from, 1) - MATRIX(*layout, to, 1));
        sqdist = dx * dx + dy * dy;
        /* Find probability associated with distance using fitted Phi */
        /* NOT 2 * b since it's already squared */
        nu = 1.0 / (1 + a * powf(sqdist, b));

        /* Term 1*/
        *cross_entropy -= mu * log(nu);
        /* Term 2*/
        *cross_entropy += (1 - mu) * log(1 - nu);
    }

    return IGRAPH_SUCCESS;
}


static igraph_error_t igraph_optimize_layout_stochastic_gradient(igraph_t *umap_graph, igraph_vector_t *umap_weights, igraph_real_t a, igraph_real_t b,
        igraph_matrix_t *layout) {
    igraph_integer_t epochs = 5; //FIXME
    igraph_real_t learning_rate = 1;
    igraph_real_t sampling_prob = 1.0; // between 0 and 1, fraction of edges sampled for gradient at each epoch
    igraph_matrix_t gradient;
    igraph_real_t cross_entropy, cross_entropy_old;


    /* Initialize gradient */
    IGRAPH_MATRIX_INIT_FINALLY(&gradient, igraph_matrix_nrow(layout), igraph_matrix_ncol(layout));

    /* Measure the (variable part of the) cross-entropy terms for debugging:
     * 1. - sum_edge_e mu(e) * log(nu(e))
     * 2. + sum_edge_e (1 - mu(e)) * log(1 - nu(e))
     * */
    igraph_compute_cross_entropy(umap_graph, umap_weights, layout, a, b, &cross_entropy);
    for (igraph_integer_t e = 0; e < epochs; e++) {
        /* Compute (stochastic) gradient */
        igraph_get_gradient(&gradient, layout, umap_graph, umap_weights, a, b, sampling_prob);

        //printf("gradient:\n");
        //igraph_matrix_print(&gradient);
        //printf("\n");

        /* Delta is the gradient times the current (decreasing) learning rate */
        igraph_matrix_scale(&gradient, learning_rate);
        /* FIXME: Perhaps a sign? */
        igraph_matrix_scale(&gradient, -1);

        /* Shift the layout by the delta */
        igraph_matrix_sub(layout, &gradient);

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
    igraph_vector_t open_set_sizes;
    igraph_vector_t open_set_decays;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_t umap_weights;
    igraph_real_t min_dist = 0.01; /* This is empyrical */
    float a, b; /* The smoothing parameters given min_dist */

    RNG_BEGIN();
    IGRAPH_VECTOR_INIT_FINALLY(&open_set_sizes, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&open_set_decays, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&umap_weights, 0);

    IGRAPH_CHECK(igraph_umap_find_open_sets(graph, distances, &open_set_sizes,
                &open_set_decays));
    IGRAPH_CHECK(igraph_umap_edge_weights(graph, distances, &umap_weights, &open_set_sizes,
                &open_set_decays));

    /* Skip spectral embedding for now, initialize at random */
    igraph_layout_random(graph, layout);

    /* Definition 11 */
    IGRAPH_CHECK(igraph_fit_ab(min_dist, &a, &b));
    /* Algorithm 5 */
    IGRAPH_CHECK(igraph_optimize_layout_stochastic_gradient(graph, &umap_weights, a, b, layout));

    igraph_vector_destroy(&open_set_sizes);
    igraph_vector_destroy(&open_set_decays);
    igraph_vector_destroy(&umap_weights);
    IGRAPH_FINALLY_CLEAN(3);
    RNG_END();
    return (IGRAPH_SUCCESS);
}
