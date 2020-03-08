
#ifndef IGRAPH_LSAP_H
#define IGRAPH_LSAP_H

#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"

__BEGIN_DECLS

int igraph_solve_lsap(igraph_matrix_t *c, igraph_integer_t n,
                      igraph_vector_int_t *p);

__END_DECLS

#endif
