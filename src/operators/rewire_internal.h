#ifndef IGRAPH_OPERATORS_REWIRE_INTERNAL_H
#define IGRAPH_OPERATORS_REWIRE_INTERNAL_H

#include "igraph_interface.h"

IGRAPH_PRIVATE_EXPORT int igraph_i_rewire(igraph_t *graph, igraph_integer_t n, igraph_rewiring_t mode, igraph_bool_t use_adjlist);

#endif
