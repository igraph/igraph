/*
   igraph library.
   Copyright (C) 2025  The igraph development team

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

#include "igraph_sampling.h"

/**
 * \function igraph_rng_sample_sphere_surface
 * \brief Sample points uniformly from the surface of a sphere.
 *
 * The center of the sphere is at the origin.
 *
 * \param rng The random number generator to use.
 * \param dim The dimension of the random vectors.
 * \param n The number of vectors to sample.
 * \param radius Radius of the sphere, it must be positive.
 * \param positive Whether to restrict sampling to the positive
 *    orthant.
 * \param res Pointer to an initialized matrix, the result is
 *    stored here, each column will be a sampled vector. The matrix is
 *    resized, as needed.
 * \return Error code.
 *
 * Time complexity: O(n*dim*g), where g is the time complexity of
 * generating a standard normal random number.
 *
 * \sa \ref igraph_rng_sample_sphere_volume(), \ref
 * igraph_rng_sample_dirichlet() for other similar samplers.
 */
igraph_error_t igraph_rng_sample_sphere_surface(
    igraph_rng_t* rng, igraph_int_t dim, igraph_int_t n, igraph_real_t radius,
    igraph_bool_t positive, igraph_matrix_t *res
) {
    igraph_int_t i, j;

    if (dim < 2) {
        IGRAPH_ERROR("Sphere must be at least two dimensional to sample from "
                     "surface.", IGRAPH_EINVAL);
    }
    if (n < 0) {
        IGRAPH_ERROR("Number of samples must be non-negative.", IGRAPH_EINVAL);
    }
    if (radius <= 0) {
        IGRAPH_ERROR("Sphere radius must be positive.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, dim, n));

    for (i = 0; i < n; i++) {
        igraph_real_t *col = &MATRIX(*res, 0, i);
        igraph_real_t sum = 0.0;
        for (j = 0; j < dim; j++) {
            col[j] = igraph_rng_get_normal(rng, 0, 1);
            sum += col[j] * col[j];
        }
        sum = sqrt(sum);
        for (j = 0; j < dim; j++) {
            col[j] = radius * col[j] / sum;
        }
        if (positive) {
            for (j = 0; j < dim; j++) {
                col[j] = fabs(col[j]);
            }
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_rng_sample_sphere_volume
 * \brief Sample points uniformly from the volume of a sphere.
 *
 * The center of the sphere is at the origin.
 *
 * \param rng The random number generator to use.
 * \param dim The dimension of the random vectors.
 * \param n The number of vectors to sample.
 * \param radius Radius of the sphere, it must be positive.
 * \param positive Whether to restrict sampling to the positive
 *    orthant.
 * \param res Pointer to an initialized matrix, the result is
 *    stored here, each column will be a sampled vector. The matrix is
 *    resized, as needed.
 * \return Error code.
 *
 * Time complexity: O(n*dim*g), where g is the time complexity of
 * generating a standard normal random number.
 *
 * \sa \ref igraph_rng_sample_sphere_surface(), \ref
 * igraph_rng_sample_dirichlet() for other similar samplers.
 */
igraph_error_t igraph_rng_sample_sphere_volume(
    igraph_rng_t* rng, igraph_int_t dim, igraph_int_t n, igraph_real_t radius,
    igraph_bool_t positive, igraph_matrix_t *res
) {

    igraph_int_t i, j;

    /* Arguments are checked by the following call */

    IGRAPH_CHECK(igraph_rng_sample_sphere_surface(rng, dim, n, radius, positive, res));

    for (i = 0; i < n; i++) {
        igraph_real_t *col = &MATRIX(*res, 0, i);
        igraph_real_t U = pow(igraph_rng_get_unif01(rng), 1.0 / dim);
        for (j = 0; j < dim; j++) {
            col[j] *= U;
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_rng_sample_dirichlet
 * \brief Sample points from a Dirichlet distribution.
 *
 * \param rng The random number generator to use.
 * \param n The number of vectors to sample.
 * \param alpha The parameters of the Dirichlet distribution. They
 *    must be positive. The length of this vector gives the dimension
 *    of the generated samples.
 * \param res Pointer to an initialized matrix, the result is stored
 *    here, one sample in each column. It will be resized, as needed.
 * \return Error code.
 *
 * Time complexity: O(n * dim * g), where dim is the dimension of the
 * sample vectors, set by the length of alpha, and g is the time
 * complexity of sampling from a Gamma distribution.
 *
 * \sa \ref igraph_rng_sample_sphere_surface() and
 * \ref igraph_rng_sample_sphere_volume() for other methods to sample
 * latent vectors.
 */
igraph_error_t igraph_rng_sample_dirichlet(
    igraph_rng_t* rng, igraph_int_t n, const igraph_vector_t *alpha,
    igraph_matrix_t *res
) {

    igraph_int_t len = igraph_vector_size(alpha);
    igraph_int_t i, j;
    igraph_real_t sum, num;

    if (n < 0) {
        IGRAPH_ERRORF("Number of samples should be non-negative, got %" IGRAPH_PRId ".",
                     IGRAPH_EINVAL, n);
    }

    if (len < 2) {
        IGRAPH_ERRORF("Dirichlet parameter vector too short, must "
                     "have at least two entries, got %" IGRAPH_PRId
                     ".", IGRAPH_EINVAL, len);
    }

    if (igraph_vector_min(alpha) <= 0) {
        IGRAPH_ERRORF("Dirichlet concentration parameters must be positive, got %g.",
                     IGRAPH_EINVAL, igraph_vector_min(alpha));
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, len, n));

    for (i = 0; i < n; i++) {
        for (j = 0, sum = 0.0; j < len; j++) {
            num = igraph_rng_get_gamma(rng, VECTOR(*alpha)[j], 1.0);
            sum += num;
            MATRIX(*res, j, i) = num;
        }
        for (j = 0; j < len; j++) {
            MATRIX(*res, j, i) /= sum;
        }
    }

    return IGRAPH_SUCCESS;
}
