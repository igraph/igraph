
#ifndef IGRAPH_CYCLES_H
#define IGRAPH_CYCLES_H

#include "igraph_decls.h"
#include "igraph_interface.h"
#include "igraph_types.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

IGRAPH_EXPORT int igraph_fundamental_cycles(const igraph_t *graph,
                                            igraph_integer_t start_vid,
                                            igraph_integer_t cutoff,
                                            igraph_vector_ptr_t *result);

IGRAPH_EXPORT int igraph_minimum_cycle_basis(const igraph_t *graph,
                                             igraph_integer_t cutoff,
                                             igraph_bool_t complete,
                                             igraph_vector_ptr_t *result);

__END_DECLS

#endif
