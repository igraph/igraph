/*
   IGraph library.
   Copyright (C) 2016-2022  The igraph development team <igraph@igraph.org>

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

#include "igraph_error.h"
#include "igraph_interface.h"

#include "core/interruption.h"
#include "cliques/cliquer_internal.h"
#include "cliques/cliquer/cliquer.h"

#include "config.h"

#include <limits.h>

/* We shall use this option struct for all calls to Cliquer */
static IGRAPH_THREAD_LOCAL clique_options igraph_cliquer_opt = {
    reorder_by_default, NULL, NULL, NULL, NULL, NULL, NULL, 0
};


/* Convert an igraph graph to a Cliquer graph */
static igraph_error_t igraph_to_cliquer(const igraph_t *ig, graph_t **cg) {
    igraph_integer_t vcount, ecount;
    igraph_integer_t i;

    if (igraph_is_directed(ig)) {
        IGRAPH_WARNING("Edge directions are ignored for clique calculations");
    }

    vcount = igraph_vcount(ig);
    ecount = igraph_ecount(ig);

    if (vcount > INT_MAX) {
        IGRAPH_ERROR("Graph too large for Cliquer", IGRAPH_EOVERFLOW);
    }

    *cg = graph_new((int) vcount);

    for (i = 0; i < ecount; ++i) {
        igraph_integer_t s, t;
        s = IGRAPH_FROM(ig, i);
        t = IGRAPH_TO(ig, i);
        if (s != t) {
            GRAPH_ADD_EDGE(*cg, s, t);
        }
    }

    return IGRAPH_SUCCESS;
}


/* Copy weights to a Cliquer graph */
static igraph_error_t set_weights(const igraph_vector_t *vertex_weights, graph_t *g) {
    igraph_integer_t i;

    IGRAPH_ASSERT(vertex_weights != NULL);

    if (igraph_vector_size(vertex_weights) != g->n) {
        IGRAPH_ERROR("Invalid vertex weight vector length", IGRAPH_EINVAL);
    }

    for (i = 0; i < g->n; ++i) {
        g->weights[i] = VECTOR(*vertex_weights)[i];
        if (g->weights[i] != VECTOR(*vertex_weights)[i]) {
            IGRAPH_WARNING("Only integer vertex weights are supported; weights will be truncated to their integer parts");
        }
        if (g->weights[i] <= 0) {
            IGRAPH_ERROR("Vertex weights must be positive", IGRAPH_EINVAL);
        }
    }

    return IGRAPH_SUCCESS;
}


/* Find all cliques. */

typedef struct {
    igraph_vector_int_t clique;
    igraph_vector_int_list_t* result;
} igraph_i_cliquer_cliques_user_data_t;

static igraph_error_t igraph_i_cliquer_cliques_user_data_init(
    igraph_i_cliquer_cliques_user_data_t* data,
    igraph_vector_int_list_t* result
) {
    data->result = result;
    igraph_vector_int_list_clear(result);
    return igraph_vector_int_init(&data->clique, 0);
}

static void igraph_i_cliquer_cliques_user_data_destroy(
    igraph_i_cliquer_cliques_user_data_t* data
) {
    igraph_vector_int_destroy(&data->clique);
    data->result = 0;
}

static igraph_error_t collect_cliques_callback(set_t s, graph_t *g, clique_options *opt) {
    int i;
    igraph_integer_t j;
    igraph_i_cliquer_cliques_user_data_t* data = (igraph_i_cliquer_cliques_user_data_t *) opt->user_data;

    IGRAPH_UNUSED(g);

    IGRAPH_ALLOW_INTERRUPTION();

    IGRAPH_CHECK(igraph_vector_int_resize(&data->clique, set_size(s)));

    i = -1; j = 0;
    while ((i = set_return_next(s, i)) >= 0) {
        VECTOR(data->clique)[j++] = i;
    }

    IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(data->result, &data->clique));

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_cliquer_cliques(const igraph_t *graph, igraph_vector_int_list_t *res,
                             igraph_integer_t min_size, igraph_integer_t max_size) {
    graph_t *g;
    igraph_integer_t vcount = igraph_vcount(graph);
    igraph_i_cliquer_cliques_user_data_t data;

    if (vcount == 0) {
        igraph_vector_int_list_clear(res);
        return IGRAPH_SUCCESS;
    }

    if (min_size <= 0) {
        min_size = 1;
    }
    if (max_size <= 0) {
        max_size = 0;
    }

    if (max_size > INT_MAX) {
        max_size = INT_MAX;
    }

    if (max_size > 0 && max_size < min_size) {
        IGRAPH_ERROR("max_size must not be smaller than min_size", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_i_cliquer_cliques_user_data_init(&data, res));
    IGRAPH_FINALLY(igraph_i_cliquer_cliques_user_data_destroy, &data);

    IGRAPH_CHECK(igraph_to_cliquer(graph, &g));
    IGRAPH_FINALLY(graph_free, g);

    igraph_cliquer_opt.user_data = &data;
    igraph_cliquer_opt.user_function = &collect_cliques_callback;

    IGRAPH_CHECK(clique_unweighted_find_all(g, (int) min_size, (int) max_size, /* maximal= */ FALSE, &igraph_cliquer_opt, NULL));

    graph_free(g);
    igraph_i_cliquer_cliques_user_data_destroy(&data);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/* Count cliques of each size. */

static igraph_error_t count_cliques_callback(set_t s, graph_t *g, clique_options *opt) {
    igraph_vector_t *hist;

    IGRAPH_UNUSED(g);

    IGRAPH_ALLOW_INTERRUPTION();

    hist = (igraph_vector_t *) opt->user_data;
    VECTOR(*hist)[set_size(s) - 1] += 1;

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_cliquer_histogram(const igraph_t *graph, igraph_vector_t *hist,
                               igraph_integer_t min_size, igraph_integer_t max_size) {
    graph_t *g;
    igraph_integer_t i;
    igraph_integer_t vcount = igraph_vcount(graph);

    if (vcount == 0) {
        igraph_vector_clear(hist);
        return IGRAPH_SUCCESS;
    }

    if (min_size <= 0) {
        min_size = 1;
    }
    if (max_size <= 0) {
        max_size = vcount;    /* also used for initial hist vector size, do not set to zero */
    }

    if (max_size > INT_MAX) {
        max_size = INT_MAX;
    }

    if (max_size < min_size) {
        IGRAPH_ERRORF("Maximum clique size (%" IGRAPH_PRId ") must not be "
                      "smaller than minimum clique size (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, max_size, min_size);
    }

    IGRAPH_CHECK(igraph_to_cliquer(graph, &g));
    IGRAPH_FINALLY(graph_free, g);

    IGRAPH_CHECK(igraph_vector_resize(hist, max_size));
    igraph_vector_null(hist);
    igraph_cliquer_opt.user_data = hist;
    igraph_cliquer_opt.user_function = &count_cliques_callback;

    IGRAPH_CHECK(clique_unweighted_find_all(g, (int) min_size, (int) max_size, /* maximal= */ FALSE, &igraph_cliquer_opt, NULL));

    for (i = max_size; i > 0; --i) {
        if (VECTOR(*hist)[i - 1] > 0) {
            break;
        }
    }
    IGRAPH_CHECK(igraph_vector_resize(hist, i));
    igraph_vector_resize_min(hist);

    graph_free(g);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/* Call function for each clique. */

struct callback_data {
    igraph_vector_int_t *clique;
    igraph_clique_handler_t *handler;
    void *arg;
};

static igraph_error_t callback_callback(set_t s, graph_t *g, clique_options *opt) {
    struct callback_data *cd;
    int i;
    igraph_integer_t j;
    igraph_error_t retval;

    IGRAPH_UNUSED(g);

    IGRAPH_ALLOW_INTERRUPTION();

    cd = (struct callback_data *) opt->user_data;

    IGRAPH_CHECK(igraph_vector_int_resize(cd->clique, set_size(s)));

    i = -1; j = 0;
    while ((i = set_return_next(s, i)) >= 0) {
        VECTOR(*cd->clique)[j++] = i;
    }

    retval = (*(cd->handler))(cd->clique, cd->arg);

    return retval;
}

igraph_error_t igraph_i_cliquer_callback(const igraph_t *graph,
                              igraph_integer_t min_size, igraph_integer_t max_size,
                              igraph_clique_handler_t *cliquehandler_fn, void *arg) {
    graph_t *g;
    igraph_vector_int_t current_clique;
    struct callback_data cd;
    igraph_integer_t vcount = igraph_vcount(graph);

    if (vcount == 0) {
        return IGRAPH_SUCCESS;
    }

    if (min_size <= 0) {
        min_size = 1;
    }
    if (max_size <= 0) {
        max_size = 0;
    }

    if (max_size > INT_MAX) {
        max_size = INT_MAX;
    }

    if (max_size > 0 && max_size < min_size) {
        IGRAPH_ERROR("max_size must not be smaller than min_size", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_to_cliquer(graph, &g));
    IGRAPH_FINALLY(graph_free, g);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&current_clique, min_size);

    cd.clique = &current_clique;
    cd.handler = cliquehandler_fn;
    cd.arg = arg;
    igraph_cliquer_opt.user_data = &cd;
    igraph_cliquer_opt.user_function = &callback_callback;

    IGRAPH_CHECK(clique_unweighted_find_all(g, (int) min_size, (int) max_size, /* maximal= */ FALSE, &igraph_cliquer_opt, NULL));

    igraph_vector_int_destroy(&current_clique);
    graph_free(g);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/* Find weighted cliques in given weight range. */

igraph_error_t igraph_i_weighted_cliques(const igraph_t *graph,
                              const igraph_vector_t *vertex_weights, igraph_vector_int_list_t *res,
                              igraph_real_t min_weight, igraph_real_t max_weight, igraph_bool_t maximal) {
    graph_t *g;
    igraph_integer_t vcount = igraph_vcount(graph);
    igraph_i_cliquer_cliques_user_data_t data;

    if (vcount == 0) {
        igraph_vector_int_list_clear(res);
        return IGRAPH_SUCCESS;
    }

    if (min_weight != (int) min_weight) {
        IGRAPH_WARNING("Only integer vertex weights are supported; the minimum weight will be truncated to its integer part");
        min_weight  = (int) min_weight;
    }

    if (max_weight != (int) max_weight) {
        IGRAPH_WARNING("Only integer vertex weights are supported; the maximum weight will be truncated to its integer part");
        max_weight = (int) max_weight;
    }

    if (min_weight <= 0) {
        min_weight = 1;
    }
    if (max_weight <= 0) {
        max_weight = 0;
    }

    if (max_weight > 0 && max_weight < min_weight) {
        IGRAPH_ERROR("max_weight must not be smaller than min_weight", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_i_cliquer_cliques_user_data_init(&data, res));
    IGRAPH_FINALLY(igraph_i_cliquer_cliques_user_data_destroy, &data);

    IGRAPH_CHECK(igraph_to_cliquer(graph, &g));
    IGRAPH_FINALLY(graph_free, g);

    IGRAPH_CHECK(set_weights(vertex_weights, g));

    igraph_cliquer_opt.user_data = &data;
    igraph_cliquer_opt.user_function = &collect_cliques_callback;

    IGRAPH_CHECK(clique_find_all(g, (int) min_weight, (int) max_weight, maximal, &igraph_cliquer_opt, NULL));

    graph_free(g);
    igraph_i_cliquer_cliques_user_data_destroy(&data);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/* Find largest weighted cliques. */

igraph_error_t igraph_i_largest_weighted_cliques(const igraph_t *graph,
                                      const igraph_vector_t *vertex_weights, igraph_vector_int_list_t *res) {
    graph_t *g;
    igraph_integer_t vcount = igraph_vcount(graph);
    igraph_i_cliquer_cliques_user_data_t data;

    if (vcount == 0) {
        igraph_vector_int_list_clear(res);
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_i_cliquer_cliques_user_data_init(&data, res));
    IGRAPH_FINALLY(igraph_i_cliquer_cliques_user_data_destroy, &data);

    IGRAPH_CHECK(igraph_to_cliquer(graph, &g));
    IGRAPH_FINALLY(graph_free, g);

    IGRAPH_CHECK(set_weights(vertex_weights, g));

    igraph_cliquer_opt.user_data = &data;
    igraph_cliquer_opt.user_function = &collect_cliques_callback;

    IGRAPH_CHECK(clique_find_all(g, 0, 0, FALSE, &igraph_cliquer_opt, NULL));

    graph_free(g);
    igraph_i_cliquer_cliques_user_data_destroy(&data);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/* Find weight of largest weight clique. */

static igraph_error_t check_interruption_callback(set_t s, graph_t *g, clique_options *opt) {
    IGRAPH_UNUSED(s); IGRAPH_UNUSED(g); IGRAPH_UNUSED(opt);
    IGRAPH_ALLOW_INTERRUPTION();
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_weighted_clique_number(const igraph_t *graph,
                                    const igraph_vector_t *vertex_weights, igraph_real_t *res) {
    graph_t *g;
    igraph_integer_t vcount = igraph_vcount(graph);
    int res_int;

    if (vcount == 0) {
        if (res) {
            *res = 0;
        }
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_to_cliquer(graph, &g));
    IGRAPH_FINALLY(graph_free, g);

    IGRAPH_CHECK(set_weights(vertex_weights, g));

    igraph_cliquer_opt.user_function = check_interruption_callback;

    IGRAPH_CHECK(clique_max_weight(g, &igraph_cliquer_opt, &res_int));

    graph_free(g);
    IGRAPH_FINALLY_CLEAN(1);

    if (res) {
        *res = res_int;
    }

    return IGRAPH_SUCCESS;
}
