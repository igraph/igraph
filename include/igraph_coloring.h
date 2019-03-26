#ifndef IGRAPH_COLORING_H
#define IGRAPH_COLORING_H

#include "igraph_datatype.h"

__BEGIN_DECLS

/**
 * \typedef igraph_coloring_greedy_t
 * Ordering heuristics for igraph_vertex_coloring_greedy
 *
 * \enumval IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS  Choose vertex with largest number of already colored neighbors.
 *
 */
typedef enum {
    IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS = 0
} igraph_coloring_greedy_t;

DECLDIR int igraph_vertex_coloring_greedy(const igraph_t *graph, igraph_vector_int_t *colors, igraph_coloring_greedy_t heuristic);

__END_DECLS

#endif /* IGRAPH_COLORING_H */
