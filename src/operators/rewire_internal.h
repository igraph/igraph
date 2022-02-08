#ifndef IGRAPH_OPERATORS_REWIRE_INTERNAL_H
#define IGRAPH_OPERATORS_REWIRE_INTERNAL_H

#include "igraph_decls.h"
#include "igraph_interface.h"

__BEGIN_DECLS

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_i_rewire(
    igraph_t *graph, igraph_integer_t n, igraph_rewiring_t mode,
    igraph_bool_t use_adjlist);

__END_DECLS

#endif
