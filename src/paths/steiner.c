#include "igraph_paths.h"
#include "igraph_adjlist.h"
#include "igraph_types.h"
#include "igraph_memory.h"
#include "igraph_error.h"


#include <string.h>
#include <math.h>
#include <map>
#include <limits.h>



igraph_vector_t generateSubsets(igraph_vector_t* C, igraph_integer_t n,igraph_integer_t graphsize)
{
    igraph_integer_t count = pow(2, n);

    igraph_vector_t *allSubsets;
    
    igraph_integer_t subsetIndex = graphsize;
    
    // The outer for loop will run 2^n times to print all subset .
    // Here variable i will act as a binary counter

    for (int i = 0; i < count; i++)
    {
        // The inner for loop will run n times , As the maximum number of elements a set can have is n
        // This loop will generate a subset
        igraph_vector_t newSubset;
        for (int j = 0; j < n; j++)
        {
            // This if condition will check if jth bit in binary representation of  i  is set or not
            // if the value of (i & (1 << j)) is greater than 0 , include arr[j] in the current subset
            // otherwise exclude arr[j]
            if ((i & (1 << j)) > 0)
            {
                newSubset.insert(C[j]);
            }
        }

        if (newSubset.size() > 1)
        {
            allSubsets.push_back(newSubset);
            //subsetMap[newSubset] = subsetIndex;
            subsetMap.insert(make_pair(newSubset,subsetIndex));
            cout << "Adding subset with index" << subsetIndex << "\n";
            subsetIndex++;
        }
    }

/*
    for (auto const &x : subsetMap)
    {
        cout << "Subset = ";
        igraph_vector_t element = x.first;
        for (auto f : element)
        {
            cout << f << " ";
        }
	
        cout << "\n";
        std::cout << "Index = " << x.second << std::endl;
    }

    cout<<"Completed subsets";
*/

    return allSubsets;
}



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
