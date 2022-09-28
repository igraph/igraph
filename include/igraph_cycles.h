
#ifndef IGRAPH_CYCLES_H
#define IGRAPH_CYCLES_H

#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector_list.h"
#include "igraph_adjlist.h"
#include "igraph_stack.h"

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
  igraph_integer_t N;
  igraph_adjlist_t AK;
  igraph_adjlist_t B;
  igraph_stack_int_t stack;
  igraph_vector_bool_t blocked;
} igraph_simple_cycle_search_state_t;

IGRAPH_EXPORT igraph_error_t igraph_simple_cycle_search_state_init(
    igraph_simple_cycle_search_state_t *state, const igraph_t *graph);

IGRAPH_EXPORT igraph_error_t igraph_simple_cycle_search_state_destroy(
    igraph_simple_cycle_search_state_t *state);

IGRAPH_EXPORT igraph_error_t igraph_simple_cycles_search_one(
    igraph_simple_cycle_search_state_t *state, igraph_integer_t start,
    igraph_vector_int_list_t *result);

IGRAPH_EXPORT igraph_error_t igraph_simple_cycles_search_all(
    const igraph_t *graph, igraph_vector_int_list_t *result);

__END_DECLS

#endif
