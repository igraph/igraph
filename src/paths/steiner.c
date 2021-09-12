#include "igraph_paths.h"
#include "igraph_adjlist.h"
#include "igraph_types.h"
#include "igraph_memory.h"
#include "igraph_error.h"


#include <string.h>
#include <math.h>
#include <map>
#include <limits.h>



igraph_vector_t* generateSubsets(igraph_vector_t* C, igraph_integer_t n,igraph_integer_t graphsize)
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


int fetchIndexofMapofSets(igraph_vector_t subset)
{
    int key;
    for(map< igraph_vector_t , int >::const_iterator it = subsetMap.begin(); it != subsetMap.end(); ++it)
    {

        if (it -> first == subset)
        {
            key = it ->second;
        }
    }

 return key;
}

int igraph_steiner_ratio (const igraph_t *graph_contracted , igraph_vector_t *q_terminals , igraph_integer_t vid , igraph_neimode_t mode )
{
	
	igraph_es_t es;
	igraph_eit_t eit;	
		
	if igraph_vector_size(&q_terminals) >= 2
	{
		
		es_vid = igraph_es_incident(&es,vid,mode);
		igraph_eit_create(&graph, es , &eit);		
		
		while(!IGRAPH_EIT_END(eit))
		{
			
			eid = IGRAPH_EIT_GET(eit);	
		}
						
				
	} 


}


int igraph_steiner_dreyfus_wagner(const igraph_t *graph,
                                  igraph_vector_t *steiner_terminals,
                                  igraph_neimode_t mode,const igraph_vector_t *weights) {
        


long int no_of_vertices = (igraph_integer_t) igraph_vcount(&graph);
long int no_of_edges = igraph_ecount(graph);


igraph_vector_t steiner_vertices;
igraph_matrix_t dp_cache; // dynamic programming table
igraph_integer_t q;
igraph_vector_t *allSubsets;
igraph_matrix_t distance;
    
IGRAPH_VECTOR_INIT_FINALLY(&steiner_vertices, 0);

IGRAPH_VECTOR_INIT_FINALLY(&allSubsets, 0);    
    
   
if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Weight vector length does not match", IGRAPH_EINVAL);
    }


distance = igraph_shortest_paths_floyd_warshall(&graph,&weights,mode)

igraph_vector_sort(&steiner_terminals);


// Creating a vector of steiner vertices. steiner vertices = vertices in graph - steiner terminals

igraph_vector_append(&steiner_vertices,&no_of_vertices);

for (int i = 0,j = 0;i < igraph_vector_size(&steiner_terminals); i++,j++)
{
	IGRAPH_CHECK(igraph_vector_remove(&steiner_vertices,i-j));		

}

igraph_matrix_init(&dp_cache,pow(2,igraph_vector_size(&steiner_terminals)),igraph_vector_size(&steiner_vertices));
igraph_matrix_fill(&dp_cache,INT_MAX);


q = VECTOR[steiner_terminals][0];

IGRAPH_CHECK(igraph_vector_remove(&steiner_terminals,0));

allSubsets = generateSubsets(&steiner_terminals,igraph_vector_size(&steiner_terminals),no_of_vertices);



// Singleton subset rows may be filled in trivially

for (int i =0; i < igraph_vector_size(&steiner_terminals); i++)
{
	for (int j=0; j < igraph_vector_size(&steiner_vertices); j++)
	{
		igraph_matrix_set(&dp_cache,VECTOR[steiner_terminals][i],VECTOR[steiner_vertices][j],MATRIX(*distance,VECTOR[steiner_terminals][i],VECTOR[steiner_vertices][j]));	
	
	}	
}


for (int i = igraph_vector_size(&steiner_terminals); i < igraph_matrix_capacity(&dp_cache); i++ )
{
	igraph_vector_t D = VECTOR[allsubsets][i];
	igraph_integer_t indexOfSubsetD
	
	if(D.size() == 1)
        {
        	indexOfSubsetD = VECTOR[D][0];
        }
        else
        {
                indexOfSubsetD = fetchIndexofMapofSets(&D);
        }
	
	
	for (int j = 0; j < igraph_vector_size(&steiner_vertices); j++ )
	{
		u = INT_MAX
		
		for (int subset_D_iterator = 0; subset_D_iterator < igraph_vector_size(&D) ; subset_D_iterator++)
		{
			igraph_vector_t E = VECTOR[D][subset_D_iterator];
			igraph_integer_t indexOfSubsetE = fetchIndexofMapofSets(&E);
			
			igraph_vector_t DMinusE = D;
			IGRAPH_CHECK(igraph_vector_remove(&DMinusE,E));
			
			igraph_integer_t indexOfSubsetDMinusE = fetchIndexofMapofSets(&DMinusE);
			
			if ( (MATRIX(*dp_cache,indexOfSubsetE,j )) + ( MATRIX (*dp_cache,indexOfSubsetDMinusE,j) < u )
			{
				u = (MATRIX (*dp_cache,indexOfSubsetE,j)) + ( MATRIX(*dp_cache,indexOfSubsetDMinusE,j);
			}
						
		}
		
		for (int i = 0; j < igraph_vector_size(&steiner_vertices) ; j++ )
		{
			MATRIX (*dp_cache,indexOfSubsetD,i) = min ( MATRIX (*dp_cache,indexOfSubsetD,i), MATRIX (*distance,i,j) + u )
		}
		
		
	}
		
				
}

for (int j = 0; j < igraph_vector_size(&steiner_vertices); j++ )
{
	igraph_integer_t u = INT_MAX;
	igraph_integer_t v = INT_MAX;
	
	for (int subset_C_iterator = 0; subset_C_iterator < igraph_vector_size(&steiner_terminals) ; subset_C_iterator++)
	{
		igraph_vector_t F = VECTOR[steiner_terminals][subset_C_iterator];
		igraph_integer_t indexOfSubsetF = fetchIndexofMapofSets(&F);
		
		igraph_vector_t CMinusF = C;
		IGRAPH_CHECK(igraph_vector_remove(&CMinusF,F));
		
		
		igraph_integer_t indexOfSubsetCMinusF = fetchIndexofMapofSets(&CMinusF);
		
		if ( ( MATRIX (*dp_cache,indexOfSubsetF,j) ) + ( MATRIX (*dp_cache,indexOfSubsetCMinusF,j) ) < u)
		{
			u = ( MATRIX (*dp_cache,indexOfSubsetF,j) ) + ( MATRIX (*dp_cache,indexOfSubsetCMinusF,j) );			
		}
		
	}
	
	if (MATRIX (*distance,q,j) + u < v)
	{
		v = MATRIX (*distance,q,j) + u;
	}
	
}

 IGRAPH_FINALLY_CLEAN(1);

                                 
}
