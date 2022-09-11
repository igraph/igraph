
#ifndef IGRAPH_CYCLES_H
#define IGRAPH_CYCLES_H

#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector_list.h"

__BEGIN_DECLS

IGRAPH_EXPORT igraph_error_t igraph_fundamental_cycles(
    const igraph_t *graph, igraph_vector_int_list_t *result,
    igraph_integer_t start_vid, igraph_integer_t bfs_cutoff,
    const igraph_vector_t *weights);

IGRAPH_EXPORT igraph_error_t igraph_minimum_cycle_basis(
    const igraph_t *graph, igraph_vector_int_list_t *result,
    igraph_integer_t bfs_cutoff, igraph_bool_t complete,
    igraph_bool_t use_cycle_order, const igraph_vector_t *weights);

typedef struct igraph_simple_cycle_search_state_t {
  igraph_adjlist_t AK;
  igraph_adjlist_t B;
  igraph_stack_t stack;
  igraph_vector_bool_t blocked;
} igraph_simple_cycle_search_state_t;

IGRAPH_EXPORT igraph_error_t igraph_simple_cycle_search_state_init(
    igraph_simple_cycle_search_state_t *state, const igraph_t *graph);

IGRAPH_EXPORT igraph_error_t igraph_simple_cycle_search_state_destroy(
    igraph_simple_cycle_search_state_t *state);

IGRAPH_EXPORT igraph_error_t
igraph_simple_cycles(const igraph_t *graph, igraph_vector_int_list_t *result,
                     igraph_integer_t bfs_cutoff);

__END_DECLS

#endif
