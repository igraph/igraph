
#include "igraph_error.h"
#include "igraph_memory.h"
#include "igraph_constants.h"

#include "core/interruption.h"
#include "cliques/cliquer_internal.h"
#include "cliques/cliquer/cliquer.h"

#include "config.h"

/* Call this to allow for interruption in Cliquer callback functions */
#define CLIQUER_ALLOW_INTERRUPTION() \
    { \
        if (igraph_i_interruption_handler) \
            if (igraph_allow_interruption(NULL) != IGRAPH_SUCCESS) { \
                cliquer_interrupted = 1; \
                return FALSE; \
            } \
    }

/* Interruptable Cliquer functions must be wrapped in CLIQUER_INTERRUPTABLE when called */
#define CLIQUER_INTERRUPTABLE(x) \
    { \
        cliquer_interrupted = 0; \
        x; \
        if (cliquer_interrupted) return IGRAPH_INTERRUPTED; \
    }


/* Nonzero value signals interuption from Cliquer callback function */
static IGRAPH_THREAD_LOCAL int cliquer_interrupted;


/* For use with IGRAPH_FINALLY */
static void free_clique_list(igraph_vector_ptr_t *vp) {
    igraph_integer_t i, len;
    len = igraph_vector_ptr_size(vp);
    for (i = 0; i < len; ++i) {
        igraph_vector_destroy((igraph_vector_t *) VECTOR(*vp)[i]);
    }
    igraph_vector_ptr_free_all(vp);
}

/* We shall use this option struct for all calls to Cliquer */
static IGRAPH_THREAD_LOCAL clique_options igraph_cliquer_opt = {
    reorder_by_default, NULL, NULL, NULL, NULL, NULL, NULL, 0
};


/* Convert an igraph graph to a Cliquer graph */
static void igraph_to_cliquer(const igraph_t *ig, graph_t **cg) {
    igraph_integer_t vcount, ecount;
    int i;

    if (igraph_is_directed(ig)) {
        IGRAPH_WARNING("Edge directions are ignored for clique calculations");
    }

    vcount = igraph_vcount(ig);
    ecount = igraph_ecount(ig);

    *cg = graph_new(vcount);

    for (i = 0; i < ecount; ++i) {
        long s, t;
        s = IGRAPH_FROM(ig, i);
        t = IGRAPH_TO(ig, i);
        if (s != t) {
            GRAPH_ADD_EDGE(*cg, s, t);
        }
    }
}


/* Copy weights to a Cliquer graph */
static int set_weights(const igraph_vector_t *vertex_weights, graph_t *g) {
    int i;

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

static boolean collect_cliques_callback(set_t s, graph_t *g, clique_options *opt) {
    igraph_vector_ptr_t *list;
    igraph_vector_t *clique;
    int i, j;

    IGRAPH_UNUSED(g);

    CLIQUER_ALLOW_INTERRUPTION();

    list = (igraph_vector_ptr_t *) opt->user_data;
    clique = (igraph_vector_t *) malloc(sizeof(igraph_vector_t));
    igraph_vector_init(clique, set_size(s));

    i = -1; j = 0;
    while ((i = set_return_next(s, i)) >= 0) {
        VECTOR(*clique)[j++] = i;
    }

    igraph_vector_ptr_push_back(list, clique);

    return TRUE;
}

int igraph_i_cliquer_cliques(const igraph_t *graph, igraph_vector_ptr_t *res,
                             igraph_integer_t min_size, igraph_integer_t max_size) {
    graph_t *g;
    igraph_integer_t vcount = igraph_vcount(graph);

    if (vcount == 0) {
        igraph_vector_ptr_clear(res);
        return IGRAPH_SUCCESS;
    }

    if (min_size <= 0) {
        min_size = 1;
    }
    if (max_size <= 0) {
        max_size = 0;
    }

    if (max_size > 0 && max_size < min_size) {
        IGRAPH_ERROR("max_size must not be smaller than min_size", IGRAPH_EINVAL);
    }

    igraph_to_cliquer(graph, &g);
    IGRAPH_FINALLY(graph_free, g);

    igraph_vector_ptr_clear(res);
    igraph_cliquer_opt.user_data = res;
    igraph_cliquer_opt.user_function = &collect_cliques_callback;

    IGRAPH_FINALLY(free_clique_list, res);
    CLIQUER_INTERRUPTABLE(clique_unweighted_find_all(g, min_size, max_size, /* maximal= */ FALSE, &igraph_cliquer_opt));
    IGRAPH_FINALLY_CLEAN(1);

    graph_free(g);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/* Count cliques of each size. */

static boolean count_cliques_callback(set_t s, graph_t *g, clique_options *opt) {
    igraph_vector_t *hist;

    IGRAPH_UNUSED(g);

    CLIQUER_ALLOW_INTERRUPTION();

    hist = (igraph_vector_t *) opt->user_data;
    VECTOR(*hist)[set_size(s) - 1] += 1;

    return TRUE;
}

int igraph_i_cliquer_histogram(const igraph_t *graph, igraph_vector_t *hist,
                               igraph_integer_t min_size, igraph_integer_t max_size) {
    graph_t *g;
    int i;
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

    if (max_size < min_size) {
        IGRAPH_ERRORF("Maximum clique size (%" IGRAPH_PRId ") must not be "
                      "smaller than minimum clique size (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, max_size, min_size);
    }

    igraph_to_cliquer(graph, &g);
    IGRAPH_FINALLY(graph_free, g);

    igraph_vector_resize(hist, max_size);
    igraph_vector_null(hist);
    igraph_cliquer_opt.user_data = hist;
    igraph_cliquer_opt.user_function = &count_cliques_callback;

    CLIQUER_INTERRUPTABLE(clique_unweighted_find_all(g, min_size, max_size, /* maximal= */ FALSE, &igraph_cliquer_opt));

    for (i = max_size; i > 0; --i)
        if (VECTOR(*hist)[i - 1] > 0) {
            break;
        }
    igraph_vector_resize(hist, i);
    igraph_vector_resize_min(hist);

    graph_free(g);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/* Call function for each clique. */

struct callback_data {
    igraph_clique_handler_t *handler;
    void *arg;
};

static boolean callback_callback(set_t s, graph_t *g, clique_options *opt) {
    igraph_vector_t *clique;
    struct callback_data *cd;
    int i, j;

    IGRAPH_UNUSED(g);

    CLIQUER_ALLOW_INTERRUPTION();

    cd = (struct callback_data *) opt->user_data;

    clique = (igraph_vector_t *) malloc(sizeof(igraph_vector_t));
    igraph_vector_init(clique, set_size(s));

    i = -1; j = 0;
    while ((i = set_return_next(s, i)) >= 0) {
        VECTOR(*clique)[j++] = i;
    }

    return (*(cd->handler))(clique, cd->arg);
}

int igraph_i_cliquer_callback(const igraph_t *graph,
                              igraph_integer_t min_size, igraph_integer_t max_size,
                              igraph_clique_handler_t *cliquehandler_fn, void *arg) {
    graph_t *g;
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

    if (max_size > 0 && max_size < min_size) {
        IGRAPH_ERROR("max_size must not be smaller than min_size", IGRAPH_EINVAL);
    }

    igraph_to_cliquer(graph, &g);
    IGRAPH_FINALLY(graph_free, g);

    cd.handler = cliquehandler_fn;
    cd.arg = arg;
    igraph_cliquer_opt.user_data = &cd;
    igraph_cliquer_opt.user_function = &callback_callback;

    CLIQUER_INTERRUPTABLE(clique_unweighted_find_all(g, min_size, max_size, /* maximal= */ FALSE, &igraph_cliquer_opt));

    graph_free(g);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/* Find weighted cliques in given weight range. */

int igraph_i_weighted_cliques(const igraph_t *graph,
                              const igraph_vector_t *vertex_weights, igraph_vector_ptr_t *res,
                              igraph_real_t min_weight, igraph_real_t max_weight, igraph_bool_t maximal) {
    graph_t *g;
    igraph_integer_t vcount = igraph_vcount(graph);

    if (vcount == 0) {
        igraph_vector_ptr_clear(res);
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

    igraph_to_cliquer(graph, &g);
    IGRAPH_FINALLY(graph_free, g);

    IGRAPH_CHECK(set_weights(vertex_weights, g));

    igraph_vector_ptr_clear(res);
    igraph_cliquer_opt.user_data = res;
    igraph_cliquer_opt.user_function = &collect_cliques_callback;

    IGRAPH_FINALLY(free_clique_list, res);
    CLIQUER_INTERRUPTABLE(clique_find_all(g, min_weight, max_weight, maximal, &igraph_cliquer_opt));
    IGRAPH_FINALLY_CLEAN(1);

    graph_free(g);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/* Find largest weighted cliques. */

int igraph_i_largest_weighted_cliques(const igraph_t *graph,
                                      const igraph_vector_t *vertex_weights, igraph_vector_ptr_t *res) {
    graph_t *g;
    igraph_integer_t vcount = igraph_vcount(graph);

    if (vcount == 0) {
        igraph_vector_ptr_clear(res);
        return IGRAPH_SUCCESS;
    }

    igraph_to_cliquer(graph, &g);
    IGRAPH_FINALLY(graph_free, g);

    IGRAPH_CHECK(set_weights(vertex_weights, g));

    igraph_vector_ptr_clear(res);
    igraph_cliquer_opt.user_data = res;
    igraph_cliquer_opt.user_function = &collect_cliques_callback;

    IGRAPH_FINALLY(free_clique_list, res);
    CLIQUER_INTERRUPTABLE(clique_find_all(g, 0, 0, FALSE, &igraph_cliquer_opt));
    IGRAPH_FINALLY_CLEAN(1);

    graph_free(g);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/* Find weight of largest weight clique. */

int igraph_i_weighted_clique_number(const igraph_t *graph,
                                    const igraph_vector_t *vertex_weights, igraph_real_t *res) {
    graph_t *g;
    igraph_integer_t vcount = igraph_vcount(graph);

    if (vcount == 0) {
        *res = 0;
        return IGRAPH_SUCCESS;
    }

    igraph_to_cliquer(graph, &g);
    IGRAPH_FINALLY(graph_free, g);

    IGRAPH_CHECK(set_weights(vertex_weights, g));

    igraph_cliquer_opt.user_function = NULL;

    /* we are not using a callback function, thus this is not interruptable */
    *res = clique_max_weight(g, &igraph_cliquer_opt);

    graph_free(g);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
