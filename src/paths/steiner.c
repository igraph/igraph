#include "igraph_paths.h"
#include "igraph_adjlist.h"
#include "igraph_types.h"
#include "igraph_memory.h"
#include "igraph_error.h"


#include <string.h>
#include <math.h>
#include <map>
#include <limits.h>



int igraph_get_all_shortest_paths(const igraph_t *graph,
                                  igraph_vector_t *steiner_terminals,
                                  igraph_neimode_t mode) {
        

igraph_integer_t no_of_vertices = (igraph_integer_t) igraph_vcount(&graph);
igraph_vector_t steiner_vertices;
igraph_matrix_t dp_cache; // dynamic programming table


igraph_vector_sort(&steiner_terminals);

// Creating a vector of steiner vertices. steiner vertices = vertices in graph - steiner terminals

igraph_vector_append(&steiner_vertices,&no_of_vertices);

for (int i = 0,j = 0;i < igraph_vector_size(&steiner_terminals); i++,j++)
{
	igraph_vector_remove(&steiner_vertices,i-j);		

}

igraph_matrix_init(&dp_cache,pow(2,igraph_vector_size(&steiner_terminals)),igraph_vector_size(&steiner_vertices));
igraph_matrix_fill(&dp_cache,INT_MAX);

// Singleton subset rows may be filled in trivially

for (int i =0; i < igraph_vector_size(&steiner_terminals); i++)
{
	for (int j=0; j < igraph_vector_size(&steiner_vertices); j++)
	{
		igraph_matrix_set(&dp_cache,VECTOR[steiner_terminals][i],VECTOR[steiner_vertices][j],MATRIX(distance,VECTOR[steiner_terminals[i],VECTOR[steiner_vertices][j]));	
	
	}	
}





                                 
}
