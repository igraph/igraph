
#ifndef IGRAPH_CYCLES_H
#define IGRAPH_CYCLES_H

#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector_list.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Directed acyclic graphs                            */
/* -------------------------------------------------- */

IGRAPH_EXPORT igraph_error_t igraph_topological_sorting(
        const igraph_t *graph, igraph_vector_int_t *res, igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_is_dag(const igraph_t *graph, igraph_bool_t *res);

/* -------------------------------------------------- */
/* Cycle bases                                        */
/* -------------------------------------------------- */

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

__END_DECLS

#endif
