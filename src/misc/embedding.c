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

#include "igraph_embedding.h"

#include "igraph_adjlist.h"
#include "igraph_blas.h"
#include "igraph_centrality.h"
#include "igraph_interface.h"
#include "igraph_structural.h"

typedef struct {
    const igraph_t *graph;
    const igraph_vector_t *cvec;
    const igraph_vector_t *cvec2;
    igraph_adjlist_t *outlist, *inlist;
    igraph_inclist_t *eoutlist, *einlist;
    igraph_vector_t *tmp;
    const igraph_vector_t *weights;
} igraph_i_asembedding_data_t;

/* Adjacency matrix, unweighted, undirected.
   Eigendecomposition is used */
static int igraph_i_asembeddingu(igraph_real_t *to, const igraph_real_t *from,
                          int n, void *extra) {
    igraph_i_asembedding_data_t *data = extra;
    igraph_adjlist_t *outlist = data->outlist;
    const igraph_vector_t *cvec = data->cvec;
    igraph_vector_int_t *neis;
    int i, j, nlen;

    /* to = (A+cD) from */
    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(outlist, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int nei = (long int) VECTOR(*neis)[j];
            to[i] += from[nei];
        }
        to[i] += VECTOR(*cvec)[i] * from[i];
    }

    return 0;
}

/* Adjacency matrix, weighted, undirected.
   Eigendecomposition is used. */
static int igraph_i_asembeddinguw(igraph_real_t *to, const igraph_real_t *from,
                           int n, void *extra) {
    igraph_i_asembedding_data_t *data = extra;
    igraph_inclist_t *outlist = data->eoutlist;
    const igraph_vector_t *cvec = data->cvec;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *graph = data->graph;
    igraph_vector_int_t *incs;
    int i, j, nlen;

    /* to = (A+cD) from */
    for (i = 0; i < n; i++) {
        incs = igraph_inclist_get(outlist, i);
        nlen = igraph_vector_int_size(incs);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int edge = VECTOR(*incs)[j];
            long int nei = IGRAPH_OTHER(graph, edge, i);
            igraph_real_t w = VECTOR(*weights)[edge];
            to[i] += w * from[nei];
        }
        to[i] += VECTOR(*cvec)[i] * from[i];
    }

    return 0;
}

/* Adjacency matrix, unweighted, directed. SVD. */
static int igraph_i_asembedding(igraph_real_t *to, const igraph_real_t *from,
                         int n, void *extra) {
    igraph_i_asembedding_data_t *data = extra;
    igraph_adjlist_t *outlist = data->outlist;
    igraph_adjlist_t *inlist = data->inlist;
    const igraph_vector_t *cvec = data->cvec;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_int_t *neis;
    int i, j, nlen;

    /* tmp = (A+cD)' from */
    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(inlist, i);
        nlen = igraph_vector_int_size(neis);
        VECTOR(*tmp)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int nei = (long int) VECTOR(*neis)[j];
            VECTOR(*tmp)[i] += from[nei];
        }
        VECTOR(*tmp)[i] += VECTOR(*cvec)[i] * from[i];
    }

    /* to = (A+cD) tmp */
    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(outlist, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int nei = (long int) VECTOR(*neis)[j];
            to[i] += VECTOR(*tmp)[nei];
        }
        to[i] += VECTOR(*cvec)[i] * VECTOR(*tmp)[i];
    }

    return 0;
}

/* Adjacency matrix, unweighted, directed. SVD, right eigenvectors */
static int igraph_i_asembedding_right(igraph_real_t *to, const igraph_real_t *from,
                               int n, void *extra) {
    igraph_i_asembedding_data_t *data = extra;
    igraph_adjlist_t *inlist = data->inlist;
    const igraph_vector_t *cvec = data->cvec;
    igraph_vector_int_t *neis;
    int i, j, nlen;

    /* to = (A+cD)' from */
    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(inlist, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int nei = (long int) VECTOR(*neis)[j];
            to[i] += from[nei];
        }
        to[i] += VECTOR(*cvec)[i] * from[i];
    }

    return 0;
}

/* Adjacency matrix, weighted, directed. SVD. */
static int igraph_i_asembeddingw(igraph_real_t *to, const igraph_real_t *from,
                          int n, void *extra) {
    igraph_i_asembedding_data_t *data = extra;
    igraph_inclist_t *outlist = data->eoutlist;
    igraph_inclist_t *inlist = data->einlist;
    const igraph_vector_t *cvec = data->cvec;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *graph = data->graph;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_int_t *incs;
    int i, j, nlen;

    /* tmp = (A+cD)' from */
    for (i = 0; i < n; i++) {
        incs = igraph_inclist_get(inlist, i);
        nlen = igraph_vector_int_size(incs);
        VECTOR(*tmp)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int edge = VECTOR(*incs)[j];
            long int nei = IGRAPH_OTHER(graph, edge, i);
            igraph_real_t w = VECTOR(*weights)[edge];
            VECTOR(*tmp)[i] += w * from[nei];
        }
        VECTOR(*tmp)[i] += VECTOR(*cvec)[i] * from[i];
    }

    /* to = (A+cD) tmp */
    for (i = 0; i < n; i++) {
        incs = igraph_inclist_get(outlist, i);
        nlen = igraph_vector_int_size(incs);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int edge = VECTOR(*incs)[j];
            long int nei = IGRAPH_OTHER(graph, edge, i);
            igraph_real_t w = VECTOR(*weights)[edge];
            to[i] += w * VECTOR(*tmp)[nei];
        }
        to[i] += VECTOR(*cvec)[i] * VECTOR(*tmp)[i];
    }

    return 0;
}

/* Adjacency matrix, weighted, directed. SVD, right eigenvectors. */
static int igraph_i_asembeddingw_right(igraph_real_t *to, const igraph_real_t *from,
                                int n, void *extra) {
    igraph_i_asembedding_data_t *data = extra;
    igraph_inclist_t *inlist = data->einlist;
    const igraph_vector_t *cvec = data->cvec;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *graph = data->graph;
    igraph_vector_int_t *incs;
    int i, j, nlen;

    /* to = (A+cD)' from */
    for (i = 0; i < n; i++) {
        incs = igraph_inclist_get(inlist, i);
        nlen = igraph_vector_int_size(incs);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int edge = VECTOR(*incs)[j];
            long int nei = IGRAPH_OTHER(graph, edge, i);
            igraph_real_t w = VECTOR(*weights)[edge];
            to[i] += w * from[nei];
        }
        to[i] += VECTOR(*cvec)[i] * from[i];
    }

    return 0;
}

/* Laplacian D-A, unweighted, undirected. Eigendecomposition. */
static int igraph_i_lsembedding_da(igraph_real_t *to, const igraph_real_t *from,
                            int n, void *extra) {
    igraph_i_asembedding_data_t *data = extra;
    igraph_adjlist_t *outlist = data->outlist;
    const igraph_vector_t *cvec = data->cvec;
    igraph_vector_int_t *neis;
    int i, j, nlen;

    /* to = (D-A) from */
    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(outlist, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int nei = (long int) VECTOR(*neis)[j];
            to[i] -= from[nei];
        }
        to[i] += VECTOR(*cvec)[i] * from[i];
    }

    return 0;
}

/* Laplacian D-A, weighted, undirected. Eigendecomposition. */
static int igraph_i_lsembedding_daw(igraph_real_t *to, const igraph_real_t *from,
                             int n, void *extra) {
    igraph_i_asembedding_data_t *data = extra;
    igraph_inclist_t *outlist = data->eoutlist;
    const igraph_vector_t *cvec = data->cvec;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *graph = data->graph;
    igraph_vector_int_t *incs;
    int i, j, nlen;

    /* to = (D-A) from */
    for (i = 0; i < n; i++) {
        incs = igraph_inclist_get(outlist, i);
        nlen = igraph_vector_int_size(incs);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int edge = VECTOR(*incs)[j];
            long int nei = IGRAPH_OTHER(graph, edge, i);
            igraph_real_t w = VECTOR(*weights)[edge];
            to[i] -= w * from[nei];
        }
        to[i] += VECTOR(*cvec)[i] * from[i];
    }

    return 0;
}

/* Laplacian DAD, unweighted, undirected. Eigendecomposition. */
static int igraph_i_lsembedding_dad(igraph_real_t *to, const igraph_real_t *from,
                             int n, void *extra) {

    igraph_i_asembedding_data_t *data = extra;
    igraph_adjlist_t *outlist = data->outlist;
    const igraph_vector_t *cvec = data->cvec;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_int_t *neis;
    int i, j, nlen;

    /* to = D^1/2 from */
    for (i = 0; i < n; i++) {
        to[i] = VECTOR(*cvec)[i] * from[i];
    }

    /* tmp = A to */
    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(outlist, i);
        nlen = igraph_vector_int_size(neis);
        VECTOR(*tmp)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int nei = (long int) VECTOR(*neis)[j];
            VECTOR(*tmp)[i] += to[nei];
        }
    }

    /* to = D tmp */
    for (i = 0; i < n; i++) {
        to[i] = VECTOR(*cvec)[i] * VECTOR(*tmp)[i];
    }

    return 0;
}

static int igraph_i_lsembedding_dadw(igraph_real_t *to, const igraph_real_t *from,
                              int n, void *extra) {

    igraph_i_asembedding_data_t *data = extra;
    igraph_inclist_t *outlist = data->eoutlist;
    const igraph_vector_t *cvec = data->cvec;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *graph = data->graph;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_int_t *incs;
    int i, j, nlen;

    /* to = D^-1/2 from */
    for (i = 0; i < n; i++) {
        to[i] = VECTOR(*cvec)[i] * from[i];
    }

    /* tmp = A' to */
    for (i = 0; i < n; i++) {
        incs = igraph_inclist_get(outlist, i);
        nlen = igraph_vector_int_size(incs);
        VECTOR(*tmp)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int edge = VECTOR(*incs)[j];
            long int nei = IGRAPH_OTHER(graph, edge, i);
            igraph_real_t w = VECTOR(*weights)[edge];
            VECTOR(*tmp)[i] += w * to[nei];
        }
    }

    /* to = D tmp */
    for (i = 0; i < n; i++) {
        to[i] = VECTOR(*cvec)[i] * VECTOR(*cvec)[i] * VECTOR(*tmp)[i];
    }

    /* tmp = A to */
    for (i = 0; i < n; i++) {
        incs = igraph_inclist_get(outlist, i);
        nlen = igraph_vector_int_size(incs);
        VECTOR(*tmp)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            long int edge = VECTOR(*incs)[j];
            long int nei = IGRAPH_OTHER(graph, edge, i);
            igraph_real_t w = VECTOR(*weights)[edge];
            VECTOR(*tmp)[i] += w * to[nei];
        }
    }

    /* to = D^-1/2 tmp */
    for (i = 0; i < n; i++) {
        to[i] = VECTOR(*cvec)[i] * VECTOR(*tmp)[i];
    }

    return 0;
}

/* Laplacian I-DAD, unweighted, undirected. Eigendecomposition. */
static int igraph_i_lsembedding_idad(igraph_real_t *to, const igraph_real_t *from,
                              int n, void *extra) {

    int i;

    igraph_i_lsembedding_dad(to, from, n, extra);
    for (i = 0; i < n; i++) {
        to[i] = from[i] - to[i];
    }

    return 0;
}

static int igraph_i_lsembedding_idadw(igraph_real_t *to, const igraph_real_t *from,
                               int n, void *extra) {
    int i;

    igraph_i_lsembedding_dadw(to, from, n, extra);
    for (i = 0; i < n; i++) {
        to[i] = from[i] - to[i];
    }

    return 0;
}

/* Laplacian OAP, unweighted, directed. SVD. */
static int igraph_i_lseembedding_oap(igraph_real_t *to, const igraph_real_t *from,
                              int n, void *extra) {

    igraph_i_asembedding_data_t *data = extra;
    igraph_adjlist_t *outlist = data->outlist;
    igraph_adjlist_t *inlist = data->inlist;
    const igraph_vector_t *deg_in = data->cvec;
    const igraph_vector_t *deg_out = data->cvec2;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_int_t *neis;
    int i, j, nlen;

    /* tmp = O' from */
    for (i = 0; i < n; i++) {
        VECTOR(*tmp)[i] = VECTOR(*deg_out)[i] * from[i];
    }

    /* to = A' tmp */
    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(inlist, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            int nei = VECTOR(*neis)[j];
            to[i] += VECTOR(*tmp)[nei];
        }
    }

    /* tmp = P' to */
    for (i = 0; i < n; i++) {
        VECTOR(*tmp)[i] = VECTOR(*deg_in)[i] * to[i];
    }

    /* to = P tmp */
    for (i = 0; i < n; i++) {
        to[i] = VECTOR(*deg_in)[i] * VECTOR(*tmp)[i];
    }

    /* tmp = A to */
    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(outlist, i);
        nlen = igraph_vector_int_size(neis);
        VECTOR(*tmp)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            int nei = VECTOR(*neis)[j];
            VECTOR(*tmp)[i] += to[nei];
        }
    }

    /* to = O tmp */
    for (i = 0; i < n; i++) {
        to[i] = VECTOR(*deg_out)[i] * VECTOR(*tmp)[i];
    }

    return 0;
}

/* Laplacian OAP, unweighted, directed. SVD, right eigenvectors. */
static int igraph_i_lseembedding_oap_right(igraph_real_t *to,
                                    const igraph_real_t *from,
                                    int n, void *extra) {
    igraph_i_asembedding_data_t *data = extra;
    igraph_adjlist_t *inlist = data->inlist;
    const igraph_vector_t *deg_in = data->cvec;
    const igraph_vector_t *deg_out = data->cvec2;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_int_t *neis;
    int i, j, nlen;

    /* to = O' from */
    for (i = 0; i < n; i++) {
        to[i] = VECTOR(*deg_out)[i] * from[i];
    }

    /* tmp = A' to */
    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(inlist, i);
        nlen = igraph_vector_int_size(neis);
        VECTOR(*tmp)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            int nei = VECTOR(*neis)[j];
            VECTOR(*tmp)[i] += to[nei];
        }
    }

    /* to = P' tmp */
    for (i = 0; i < n; i++) {
        to[i] = VECTOR(*deg_in)[i] * VECTOR(*tmp)[i];
    }

    return 0;
}

/* Laplacian OAP, weighted, directed. SVD. */
static int igraph_i_lseembedding_oapw(igraph_real_t *to, const igraph_real_t *from,
                               int n, void *extra) {

    igraph_i_asembedding_data_t *data = extra;
    igraph_inclist_t *outlist = data->eoutlist;
    igraph_inclist_t *inlist = data->einlist;
    const igraph_vector_t *deg_in = data->cvec;
    const igraph_vector_t *deg_out = data->cvec2;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *graph = data->graph;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_int_t *neis;
    int i, j, nlen;

    /* tmp = O' from */
    for (i = 0; i < n; i++) {
        VECTOR(*tmp)[i] = VECTOR(*deg_out)[i] * from[i];
    }

    /* to = A' tmp */
    for (i = 0; i < n; i++) {
        neis = igraph_inclist_get(inlist, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            int edge = VECTOR(*neis)[j];
            int nei = IGRAPH_OTHER(graph, edge, i);
            igraph_real_t w = VECTOR(*weights)[edge];
            to[i] += w * VECTOR(*tmp)[nei];
        }
    }

    /* tmp = P' to */
    for (i = 0; i < n; i++) {
        VECTOR(*tmp)[i] = VECTOR(*deg_in)[i] * to[i];
    }

    /* to = P tmp */
    for (i = 0; i < n; i++) {
        to[i] = VECTOR(*deg_in)[i] * VECTOR(*tmp)[i];
    }

    /* tmp = A to */
    for (i = 0; i < n; i++) {
        neis = igraph_inclist_get(outlist, i);
        nlen = igraph_vector_int_size(neis);
        VECTOR(*tmp)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            int edge = VECTOR(*neis)[j];
            int nei = IGRAPH_OTHER(graph, edge, i);
            igraph_real_t w = VECTOR(*weights)[edge];
            VECTOR(*tmp)[i] += w * to[nei];
        }
    }

    /* to = O tmp */
    for (i = 0; i < n; i++) {
        to[i] = VECTOR(*deg_out)[i] * VECTOR(*tmp)[i];
    }

    return 0;
}

/* Laplacian OAP, weighted, directed. SVD, right eigenvectors. */
static int igraph_i_lseembedding_oapw_right(igraph_real_t *to,
                                     const igraph_real_t *from,
                                     int n, void *extra) {
    igraph_i_asembedding_data_t *data = extra;
    igraph_inclist_t *inlist = data->einlist;
    const igraph_vector_t *deg_in = data->cvec;
    const igraph_vector_t *deg_out = data->cvec2;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *graph = data->graph;
    igraph_vector_t *tmp = data->tmp;
    igraph_vector_int_t *neis;
    int i, j, nlen;

    /* to = O' from */
    for (i = 0; i < n; i++) {
        to[i] = VECTOR(*deg_out)[i] * from[i];
    }

    /* tmp = A' to */
    for (i = 0; i < n; i++) {
        neis = igraph_inclist_get(inlist, i);
        nlen = igraph_vector_int_size(neis);
        VECTOR(*tmp)[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            int edge = VECTOR(*neis)[j];
            int nei = IGRAPH_OTHER(graph, edge, i);
            igraph_real_t w = VECTOR(*weights)[edge];
            VECTOR(*tmp)[i] += w * to[nei];
        }
    }

    /* to = P' tmp */
    for (i = 0; i < n; i++) {
        to[i] = VECTOR(*deg_in)[i] * VECTOR(*tmp)[i];
    }

    return 0;
}

static int igraph_i_spectral_embedding(const igraph_t *graph,
                                igraph_integer_t no,
                                const igraph_vector_t *weights,
                                igraph_eigen_which_position_t which,
                                igraph_bool_t scaled,
                                igraph_matrix_t *X,
                                igraph_matrix_t *Y,
                                igraph_vector_t *D,
                                const igraph_vector_t *cvec,
                                const igraph_vector_t *cvec2,
                                igraph_arpack_options_t *options,
                                igraph_arpack_function_t *callback,
                                igraph_arpack_function_t *callback_right,
                                igraph_bool_t symmetric,
                                igraph_bool_t eigen,
                                igraph_bool_t zapsmall) {

    igraph_integer_t vc = igraph_vcount(graph);
    igraph_vector_t tmp;
    igraph_adjlist_t outlist, inlist;
    igraph_inclist_t eoutlist, einlist;
    int i, j, cveclen = igraph_vector_size(cvec);
    igraph_i_asembedding_data_t data;
    igraph_vector_t tmpD;

    data.graph = graph;
    data.cvec = cvec;
    data.cvec2 = cvec2;
    data.outlist = &outlist;
    data.inlist = &inlist;
    data.eoutlist = &eoutlist;
    data.einlist = &einlist;
    data.tmp = &tmp;
    data.weights = weights;

    if (weights && igraph_vector_size(weights) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }

    if (which != IGRAPH_EIGEN_LM &&
        which != IGRAPH_EIGEN_LA &&
        which != IGRAPH_EIGEN_SA) {
        IGRAPH_ERROR("Invalid eigenvalue chosen, must be one of "
                     "`largest magnitude', `largest algebraic' or "
                     "`smallest algebraic'", IGRAPH_EINVAL);
    }

    if (no > vc) {
        IGRAPH_ERROR("Too many singular values requested", IGRAPH_EINVAL);
    }
    if (no <= 0) {
        IGRAPH_ERROR("No singular values requested", IGRAPH_EINVAL);
    }

    if (cveclen != 1 && cveclen != vc) {
        IGRAPH_ERROR("Augmentation vector size is invalid, it should be "
                     "the number of vertices or scalar", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_matrix_resize(X, vc, no));
    if (Y) {
        IGRAPH_CHECK(igraph_matrix_resize(Y, vc, no));
    }

    /* empty graph */
    if (igraph_ecount(graph) == 0) {
        igraph_matrix_null(X);
        if (Y) {
            igraph_matrix_null(Y);
        }
        return 0;
    }

    igraph_vector_init(&tmp, vc);
    IGRAPH_FINALLY(igraph_vector_destroy, &tmp);
    if (!weights) {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &outlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &outlist);
        if (!symmetric) {
            IGRAPH_CHECK(igraph_adjlist_init(graph, &inlist, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
            IGRAPH_FINALLY(igraph_adjlist_destroy, &inlist);
        }
    } else {
        IGRAPH_CHECK(igraph_inclist_init(graph, &eoutlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &eoutlist);
        if (!symmetric) {
            IGRAPH_CHECK(igraph_inclist_init(graph, &einlist, IGRAPH_IN, IGRAPH_LOOPS_ONCE));
            IGRAPH_FINALLY(igraph_inclist_destroy, &einlist);
        }
    }
    IGRAPH_VECTOR_INIT_FINALLY(&tmpD, no);

    options->n = vc;
    options->start = 0;   /* random start vector */
    options->nev = no;
    switch (which) {
    case IGRAPH_EIGEN_LM:
        options->which[0] = 'L'; options->which[1] = 'M';
        break;
    case IGRAPH_EIGEN_LA:
        options->which[0] = 'L'; options->which[1] = 'A';
        break;
    case IGRAPH_EIGEN_SA:
        options->which[0] = 'S'; options->which[1] = 'A';
        break;
    default:
        break;
    }
    options->ncv = no + 3;
    if (options->ncv > vc) {
        options->ncv = vc;
    }

    IGRAPH_CHECK(igraph_arpack_rssolve(callback, &data, options, 0, &tmpD, X));

    if (!symmetric) {
        /* calculate left eigenvalues */
        IGRAPH_CHECK(igraph_matrix_resize(Y, vc, no));
        for (i = 0; i < no; i++) {
            igraph_real_t norm;
            igraph_vector_t v;
            callback_right(&MATRIX(*Y, 0, i), &MATRIX(*X, 0, i), vc, &data);
            igraph_vector_view(&v, &MATRIX(*Y, 0, i), vc);
            norm = 1.0 / igraph_blas_dnrm2(&v);
            igraph_vector_scale(&v, norm);
        }
    } else if (Y) {
        IGRAPH_CHECK(igraph_matrix_update(Y, X));
    }

    if (zapsmall) {
        igraph_vector_zapsmall(&tmpD, 0);
        igraph_matrix_zapsmall(X, 0);
        if (Y) {
            igraph_matrix_zapsmall(Y, 0);
        }
    }

    if (D) {
        igraph_vector_update(D, &tmpD);
        if (!eigen) {
            for (i = 0; i < no; i++) {
                VECTOR(*D)[i] = sqrt(VECTOR(*D)[i]);
            }
        }
    }

    if (scaled) {
        if (eigen) {
            /* eigenvalues were calculated */
            for (i = 0; i < no; i++) {
                VECTOR(tmpD)[i] = sqrt(fabs(VECTOR(tmpD)[i]));
            }
        } else {
            /* singular values were calculated */
            for (i = 0; i < no; i++) {
                VECTOR(tmpD)[i] = sqrt(sqrt(VECTOR(tmpD)[i]));
            }
        }

        for (j = 0; j < vc; j++) {
            for (i = 0; i < no; i++) {
                MATRIX(*X, j, i) *= VECTOR(tmpD)[i];
            }
        }

        if (Y) {
            for (j = 0; j < vc; j++) {
                for (i = 0; i < no; i++) {
                    MATRIX(*Y, j, i) *= VECTOR(tmpD)[i];
                }
            }
        }
    }

    igraph_vector_destroy(&tmpD);
    if (!weights) {
        if (!symmetric) {
            igraph_adjlist_destroy(&inlist);
            IGRAPH_FINALLY_CLEAN(1);
        }
        igraph_adjlist_destroy(&outlist);
    } else {
        if (!symmetric) {
            igraph_inclist_destroy(&einlist);
            IGRAPH_FINALLY_CLEAN(1);
        }
        igraph_inclist_destroy(&eoutlist);
    }
    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

/**
 * \function igraph_adjacency_spectral_embedding
 * Adjacency spectral embedding
 *
 * Spectral decomposition of the adjacency matrices of graphs.
 * This function computes an <code>n</code>-dimensional Euclidean
 * representation of the graph based on its adjacency
 * matrix, A. This representation is computed via the singular value
 * decomposition of the adjacency matrix, A=U D V^T. In the case,
 * where the graph is a random dot product graph generated using latent
 * position vectors in R^n for each vertex, the embedding will
 * provide an estimate of these latent vectors.
 *
 * </para><para>
 * For undirected graphs, the latent positions are calculated as
 * X = U^n D^(1/2) where U^n equals to the first no columns of U, and
 * D^(1/2) is a diagonal matrix containing the square root of the selected
 * singular values on the diagonal.
 *
 * </para><para>
 * For directed graphs, the embedding is defined as the pair
 * X = U^n D^(1/2), Y = V^n D^(1/2).
 * (For undirected graphs U=V, so it is sufficient to keep one of them.)
 *
 * \param graph The input graph, can be directed or undirected.
 * \param n An integer scalar. This value is the embedding dimension of
 *        the spectral embedding. Should be smaller than the number of
 *        vertices. The largest n-dimensional non-zero
 *        singular values are used for the spectral embedding.
 * \param weights Optional edge weights. Supply a null pointer for
 *        unweighted graphs.
 * \param which Which eigenvalues (or singular values, for directed
 *        graphs) to use, possible values:
 *        \clist
 *          \cli IGRAPH_EIGEN_LM
 *          the ones with the largest magnitude
 *          \cli IGRAPH_EIGEN_LA
 *          the (algebraic) largest ones
 *          \cli IGRAPH_EIGEN_SA
 *          the (algebraic) smallest ones.
 *        \endclist
 *        For directed graphs, <code>IGRAPH_EIGEN_LM</code> and
 *        <code>IGRAPH_EIGEN_LA</code> are the same because singular
 *        values are used for the ordering instead of eigenvalues.
 * \param scaled Whether to return X and Y (if \c scaled is true), or
 *        U and V.
 * \param X Initialized matrix, the estimated latent positions are
 *        stored here.
 * \param Y Initialized matrix or a null pointer. If not a null
 *        pointer, then the second half of the latent positions are
 *        stored here. (For undirected graphs, this always equals X.)
 * \param D Initialized vector or a null pointer. If not a null
 *        pointer, then the eigenvalues (for undirected graphs) or the
 *        singular values (for directed graphs) are stored here.
 * \param cvec A numeric vector, its length is the number vertices in the
 *        graph. This vector is added to the diagonal of the adjacency
 *        matrix, before performing the SVD.
 * \param options Options to ARPACK. See \ref igraph_arpack_options_t
 *        for details. Note that the function overwrites the
 *        <code>n</code> (number of vertices), <code>nev</code> and
 *        <code>which</code> parameters and it always starts the
 *        calculation from a random start vector.
 * \return Error code.
 *
 */

int igraph_adjacency_spectral_embedding(const igraph_t *graph,
                                        igraph_integer_t n,
                                        const igraph_vector_t *weights,
                                        igraph_eigen_which_position_t which,
                                        igraph_bool_t scaled,
                                        igraph_matrix_t *X,
                                        igraph_matrix_t *Y,
                                        igraph_vector_t *D,
                                        const igraph_vector_t *cvec,
                                        igraph_arpack_options_t *options) {

    igraph_arpack_function_t *callback, *callback_right;
    igraph_bool_t directed = igraph_is_directed(graph);

    if (directed) {
        callback = weights ? igraph_i_asembeddingw : igraph_i_asembedding;
        callback_right = (weights ? igraph_i_asembeddingw_right :
                          igraph_i_asembedding_right);
    } else {
        callback = weights ? igraph_i_asembeddinguw : igraph_i_asembeddingu;
        callback_right = 0;
    }

    return igraph_i_spectral_embedding(graph, n, weights, which, scaled,
                                       X, Y, D, cvec, /* deg2=*/ 0,
                                       options, callback, callback_right,
                                       /*symmetric=*/ !directed,
                                       /*eigen=*/ !directed, /*zapsmall=*/ 1);
}

static int igraph_i_lse_und(const igraph_t *graph,
                     igraph_integer_t no,
                     const igraph_vector_t *weights,
                     igraph_eigen_which_position_t which,
                     igraph_laplacian_spectral_embedding_type_t type,
                     igraph_bool_t scaled,
                     igraph_matrix_t *X,
                     igraph_matrix_t *Y,
                     igraph_vector_t *D,
                     igraph_arpack_options_t *options) {

    igraph_arpack_function_t *callback;
    igraph_vector_t deg;

    switch (type) {
    case IGRAPH_EMBEDDING_D_A:
        callback = weights ? igraph_i_lsembedding_daw : igraph_i_lsembedding_da;
        break;
    case IGRAPH_EMBEDDING_DAD:
        callback = weights ? igraph_i_lsembedding_dadw : igraph_i_lsembedding_dad;
        break;
    case IGRAPH_EMBEDDING_I_DAD:
        callback = weights ? igraph_i_lsembedding_idadw : igraph_i_lsembedding_idad;
        break;
    default:
        IGRAPH_ERROR("Invalid Laplacian spectral embedding type",
                     IGRAPH_EINVAL);
        break;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&deg, 0);
    igraph_strength(graph, &deg, igraph_vss_all(), IGRAPH_ALL, /*loops=*/ 1,
                    weights);

    switch (type) {
    case IGRAPH_EMBEDDING_D_A:
        break;
    case IGRAPH_EMBEDDING_DAD:
    case IGRAPH_EMBEDDING_I_DAD: {
        int i, n = igraph_vector_size(&deg);
        for (i = 0; i < n; i++) {
            VECTOR(deg)[i] = 1.0 / sqrt(VECTOR(deg)[i]);
        }
    }
    break;
    default:
        break;
    }

    IGRAPH_CHECK(igraph_i_spectral_embedding(graph, no, weights, which,
                 scaled, X, Y, D, /*cvec=*/ &deg, /*deg2=*/ 0,
                 options, callback, 0, /*symmetric=*/ 1,
                 /*eigen=*/ 1, /*zapsmall=*/ 1));

    igraph_vector_destroy(&deg);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_lse_dir(const igraph_t *graph,
                     igraph_integer_t no,
                     const igraph_vector_t *weights,
                     igraph_eigen_which_position_t which,
                     igraph_laplacian_spectral_embedding_type_t type,
                     igraph_bool_t scaled,
                     igraph_matrix_t *X,
                     igraph_matrix_t *Y,
                     igraph_vector_t *D,
                     igraph_arpack_options_t *options) {

    igraph_arpack_function_t *callback =
        weights ? igraph_i_lseembedding_oapw : igraph_i_lseembedding_oap;
    igraph_arpack_function_t *callback_right =
        weights ? igraph_i_lseembedding_oapw_right :
        igraph_i_lseembedding_oap_right;
    igraph_vector_t deg_in, deg_out;
    int i, n = igraph_vcount(graph);

    if (type != IGRAPH_EMBEDDING_OAP) {
        IGRAPH_ERROR("Invalid Laplacian spectral embedding type", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&deg_in, n);
    IGRAPH_VECTOR_INIT_FINALLY(&deg_out, n);
    igraph_strength(graph, &deg_in, igraph_vss_all(), IGRAPH_IN, /*loops=*/ 1,
                    weights);
    igraph_strength(graph, &deg_out, igraph_vss_all(), IGRAPH_OUT, /*loops=*/ 1,
                    weights);

    for (i = 0; i < n; i++) {
        VECTOR(deg_in)[i] = 1.0 / sqrt(VECTOR(deg_in)[i]);
        VECTOR(deg_out)[i] = 1.0 / sqrt(VECTOR(deg_out)[i]);
    }

    IGRAPH_CHECK(igraph_i_spectral_embedding(graph, no, weights, which,
                 scaled, X, Y, D, /*cvec=*/ &deg_in,
                 /*deg2=*/ &deg_out, options, callback,
                 callback_right, /*symmetric=*/ 0, /*eigen=*/ 0,
                 /*zapsmall=*/ 1));

    igraph_vector_destroy(&deg_in);
    igraph_vector_destroy(&deg_out);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

/**
 * \function igraph_laplacian_spectral_embedding
 * Spectral embedding of the Laplacian of a graph
 *
 * This function essentially does the same as
 * \ref igraph_adjacency_spectral_embedding, but works on the Laplacian
 * of the graph, instead of the adjacency matrix.
 * \param graph The input graph.
 * \param n The number of eigenvectors (or singular vectors if the graph
 *        is directed) to use for the embedding.
 * \param weights Optional edge weights. Supply a null pointer for
 *        unweighted graphs.
 * \param which Which eigenvalues (or singular values, for directed
 *        graphs) to use, possible values:
 *        \clist
 *          \cli IGRAPH_EIGEN_LM
 *          the ones with the largest magnitude
 *          \cli IGRAPH_EIGEN_LA
 *          the (algebraic) largest ones
 *          \cli IGRAPH_EIGEN_SA
 *          the (algebraic) smallest ones.
 *        \endclist
 *        For directed graphs, <code>IGRAPH_EIGEN_LM</code> and
 *        <code>IGRAPH_EIGEN_LA</code> are the same because singular
 *        values are used for the ordering instead of eigenvalues.
 * \param type The type of the Laplacian to use. Various definitions
 *        exist for the Laplacian of a graph, and one can choose
 *        between them with this argument. Possible values:
 *        \clist
 *          \cli IGRAPH_EMBEDDING_D_A
 *           means D - A where D is the
 *           degree matrix and A is the adjacency matrix
 *          \cli IGRAPH_EMBEDDING_DAD
 *           means Di times A times Di,
 *           where Di is the inverse of the square root of the degree matrix;
 *          \cli IGRAPH_EMBEDDING_I_DAD
 *          means I - Di A Di, where I
 *          is the identity matrix.
 *        \endclist
 * \param scaled Whether to return X and Y (if \c scaled is true), or
 *        U and V.
 * \param X Initialized matrix, the estimated latent positions are
 *        stored here.
 * \param Y Initialized matrix or a null pointer. If not a null
 *        pointer, then the second half of the latent positions are
 *        stored here. (For undirected graphs, this always equals X.)
 * \param D Initialized vector or a null pointer. If not a null
 *        pointer, then the eigenvalues (for undirected graphs) or the
 *        singular values (for directed graphs) are stored here.
 * \param options Options to ARPACK. See \ref igraph_arpack_options_t
 *        for details. Note that the function overwrites the
 *        <code>n</code> (number of vertices), <code>nev</code> and
 *        <code>which</code> parameters and it always starts the
 *        calculation from a random start vector.
 * \return Error code.
 *
 * \sa \ref igraph_adjacency_spectral_embedding to embed the adjacency
 * matrix.
 */

int igraph_laplacian_spectral_embedding(const igraph_t *graph,
                                        igraph_integer_t n,
                                        const igraph_vector_t *weights,
                                        igraph_eigen_which_position_t which,
                                        igraph_laplacian_spectral_embedding_type_t type,
                                        igraph_bool_t scaled,
                                        igraph_matrix_t *X,
                                        igraph_matrix_t *Y,
                                        igraph_vector_t *D,
                                        igraph_arpack_options_t *options) {

    if (igraph_is_directed(graph)) {
        return igraph_i_lse_dir(graph, n, weights, which, type, scaled,
                                X, Y, D, options);
    } else {
        return igraph_i_lse_und(graph, n, weights, which, type, scaled,
                                X, Y, D, options);
    }
}

/**
 * \function igraph_dim_select
 * Dimensionality selection
 *
 * Dimensionality selection for singular values using
 * profile likelihood.
 *
 * </para><para>
 * The input of the function is a numeric vector which contains
 * the measure of "importance" for each dimension.
 *
 * </para><para>
 * For spectral embedding, these are the singular values of the adjacency
 * matrix. The singular values are assumed to be generated from a
 * Gaussian mixture distribution with two components that have different
 * means and same variance. The dimensionality d is chosen to
 * maximize the likelihood when the d largest singular values are
 * assigned to one component of the mixture and the rest of the singular
 * values assigned to the other component.
 *
 * </para><para>
 * This function can also be used for the general separation problem,
 * where we assume that the left and the right of the vector are coming
 * from two normal distributions, with different means, and we want
 * to know their border.
 *
 * \param sv A numeric vector, the ordered singular values.
 * \param dim The result is stored here.
 * \return Error code.
 *
 * Time complexity: O(n), n is the number of values in sv.
 *
 * \sa \ref igraph_adjacency_spectral_embedding().
 */

int igraph_dim_select(const igraph_vector_t *sv, igraph_integer_t *dim) {

    int i, n = igraph_vector_size(sv);
    igraph_real_t x, x2, sum1 = 0.0, sum2 = igraph_vector_sum(sv);
    igraph_real_t sumsq1 = 0.0, sumsq2 = 0.0; /* to be set */
    igraph_real_t oldmean1, oldmean2, mean1 = 0.0, mean2 = sum2 / n;
    igraph_real_t varsq1 = 0.0, varsq2 = 0.0; /* to be set */
    igraph_real_t var1, var2, sd, profile, max = IGRAPH_NEGINFINITY;

    if (n == 0) {
        IGRAPH_ERROR("Need at least one singular value for dimensionality "
                     "selection", IGRAPH_EINVAL);
    }

    if (n == 1) {
        *dim = 1;
        return 0;
    }

    for (i = 0; i < n; i++) {
        x = VECTOR(*sv)[i];
        sumsq2 += x * x;
        varsq2 += (mean2 - x) * (mean2 - x);
    }

    for (i = 0; i < n - 1; i++) {
        int n1 = i + 1, n2 = n - i - 1, n1m1 = n1 - 1, n2m1 = n2 - 1;
        x = VECTOR(*sv)[i]; x2 = x * x;
        sum1 += x; sum2 -= x;
        sumsq1 += x2; sumsq2 -= x2;
        oldmean1 = mean1; oldmean2 = mean2;
        mean1 = sum1 / n1; mean2 = sum2 / n2;
        varsq1 += (x - oldmean1) * (x - mean1);
        varsq2 -= (x - oldmean2) * (x - mean2);
        var1 = i == 0 ? 0 : varsq1 / n1m1;
        var2 = i == n - 2 ? 0 : varsq2 / n2m1;
        sd = sqrt(( n1m1 * var1 + n2m1 * var2) / (n - 2));
        profile = /* - n * log(2.0*M_PI)/2.0 */ /* This is redundant */
            - n * log(sd) -
            ((sumsq1 - 2 * mean1 * sum1 + n1 * mean1 * mean1) +
             (sumsq2 - 2 * mean2 * sum2 + n2 * mean2 * mean2)) / 2.0 / sd / sd;
        if (profile > max) {
            max = profile;
            *dim = n1;
        }
    }

    /* Plus the last case, all elements in one group */
    x = VECTOR(*sv)[n - 1];
    sum1 += x;
    oldmean1 = mean1;
    mean1 = sum1 / n;
    sumsq1 += x * x;
    varsq1 += (x - oldmean1) * (x - mean1);
    var1 = varsq1 / (n - 1);
    sd = sqrt(var1);
    profile = /* - n * log(2.0*M_PI)/2.0 */ /* This is redundant */
        - n * log(sd) -
        (sumsq1 - 2 * mean1 * sum1 + n * mean1 * mean1) / 2.0 / sd / sd;
    if (profile > max) {
        max = profile;
        *dim = n;
    }

    return 0;
}
