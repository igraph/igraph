#ifndef IGRAPH_OPERATORS_REWIRE_INTERNAL_H
#define IGRAPH_OPERATORS_REWIRE_INTERNAL_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_error.h"

IGRAPH_BEGIN_C_DECLS

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_i_rewire(
    igraph_t *graph, igraph_int_t n, igraph_bool_t loops,
    igraph_bool_t use_adjlist, igraph_rewiring_stats_t *stats);

IGRAPH_END_C_DECLS

#endif
