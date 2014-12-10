
#ifndef IGRAPH_LSAP
#define IGRAPH_LSAP

#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"

int igraph_solve_lsap(igraph_matrix_t *c, igraph_integer_t n,
		      igraph_vector_int_t *p);

#endif
