
#ifndef IGRAPH_CYCLES_H
#define IGRAPH_CYCLES_H

#include "igraph_decls.h"
#include "igraph_interface.h"
#include "igraph_types.h"
#include "igraph_vector_list.h"

__BEGIN_DECLS

IGRAPH_EXPORT igraph_error_t igraph_fundamental_cycles(const igraph_t *graph,
                                            igraph_integer_t start_vid,
                                            igraph_integer_t cutoff,
                                            igraph_vector_int_list_t *result);

IGRAPH_EXPORT igraph_error_t igraph_minimum_cycle_basis(const igraph_t *graph,
                                             igraph_integer_t cutoff,
                                             igraph_bool_t complete,
                                             igraph_vector_int_list_t *result);

__END_DECLS

#endif
