
#ifndef IGRAPH_CYCLES_H
#define IGRAPH_CYCLES_H

#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector_list.h"

__BEGIN_DECLS

IGRAPH_EXPORT igraph_error_t igraph_fundamental_cycles(
        const igraph_t *graph,
        igraph_vector_int_list_t *result,
        igraph_integer_t start_vid,
        igraph_integer_t bfs_cutoff,
        const igraph_vector_t *weights);

IGRAPH_EXPORT igraph_error_t igraph_minimum_cycle_basis(
        const igraph_t *graph,
        igraph_vector_int_list_t *result,
        igraph_integer_t bfs_cutoff,
        igraph_bool_t complete,
        igraph_bool_t use_cycle_order,
        const igraph_vector_t *weights);

IGRAPH_EXPORT igraph_error_t igraph_find_cycle(
        const igraph_t *graph,
        igraph_vector_int_t *vertices,
        igraph_vector_int_t *edges,
        igraph_neimode_t mode);

struct igraph_simple_cycle_search_state_t;
typedef struct igraph_simple_cycle_search_state_t igraph_simple_cycle_search_state_t;

IGRAPH_EXPORT igraph_error_t igraph_simple_cycle_search_state_init(
    igraph_simple_cycle_search_state_t *state, const igraph_t *graph);

IGRAPH_EXPORT void igraph_simple_cycle_search_state_destroy(
    igraph_simple_cycle_search_state_t *state);

/**
 * \brief The interface for the callback function for when a cycle is found.
 * 
 * \param vertices The vertices of the current cycle.
 * \param edges The edges of the current cycle.
 * \return Error code; \c IGRAPH_SUCCESS to continue the search or
 *   \c IGRAPH_STOP to stop the search without signaling an error.
 */
typedef igraph_error_t igraph_simple_cycle_handler_t(const igraph_vector_int_t *vertices, const igraph_vector_int_t *edges, void *arg);

IGRAPH_EXPORT igraph_error_t igraph_simple_cycles_search_callback_from_one_vertex(
    igraph_simple_cycle_search_state_t *state, igraph_integer_t start, igraph_integer_t max_cycle_length,
    igraph_simple_cycle_handler_t *cycle_handler, void *arg
);

IGRAPH_EXPORT igraph_error_t igraph_simple_cycles_search_from_one_vertex(
    igraph_simple_cycle_search_state_t *state, igraph_integer_t start,
    igraph_vector_int_list_t *v_result, igraph_vector_int_list_t *e_result,
    igraph_integer_t max_cycle_length);

IGRAPH_EXPORT igraph_error_t igraph_simple_cycles_search_callback(
    const igraph_t *graph, igraph_integer_t max_cycle_length,
    igraph_simple_cycle_handler_t *cycle_handler, void *arg
);

IGRAPH_EXPORT igraph_error_t igraph_simple_cycles_search_all(
    const igraph_t *graph, igraph_vector_int_list_t *v_result,
    igraph_vector_int_list_t *e_result, igraph_integer_t max_cycle_length);

__END_DECLS

#endif
