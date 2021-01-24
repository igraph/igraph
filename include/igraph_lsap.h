
#ifndef IGRAPH_LSAP_H
#define IGRAPH_LSAP_H

#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"

__BEGIN_DECLS
igraph_error_t igraph_solve_lsap(igraph_matrix_t *c, igraph_long_t n,
                      igraph_vector_long_t *p);

__END_DECLS

#endif
