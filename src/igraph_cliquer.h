#ifndef IGRAPH_CLIQUER_H
#define IGRAPH_CLIQUER_H

#include "igraph_types_internal.h"
#include "igraph_interface.h"
#include "igraph_cliques.h"

int igraph_i_cliquer_cliques(const igraph_t *graph, igraph_vector_ptr_t *res,
                             igraph_integer_t min_size, igraph_integer_t max_size);

int igraph_i_cliquer_histogram(const igraph_t *graph, igraph_vector_t *hist,
                               igraph_integer_t min_size, igraph_integer_t max_size);

int igraph_i_cliquer_callback(const igraph_t *graph,
                              igraph_integer_t min_size, igraph_integer_t max_size,
                              igraph_clique_handler_t *cliquehandler_fn, void *arg);

int igraph_i_weighted_cliques(const igraph_t *graph,
                              const igraph_vector_t *vertex_weights, igraph_vector_ptr_t *res,
                              igraph_real_t min_weight, igraph_real_t max_weight, igraph_bool_t maximal);

int igraph_i_largest_weighted_cliques(const igraph_t *graph,
                                      const igraph_vector_t *vertex_weights, igraph_vector_ptr_t *res);

int igraph_i_weighted_clique_number(const igraph_t *graph,
                                    const igraph_vector_t *vertex_weights, igraph_real_t *res);

#endif // IGRAPH_CLIQUER_H

