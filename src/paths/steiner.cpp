#include "igraph_paths.h"
#include "igraph_adjlist.h"
#include "igraph_types.h"
#include "igraph_memory.h"
#include "igraph_error.h"
#include "igraph.h"

#include <cstring>
#include <cmath>
#include <map>
#include <climits>
#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

std::map<std::set<igraph_integer_t>, igraph_integer_t> subsetMap;

void printSubsets(std::set<std::set<igraph_integer_t>> allSubsets)
{
	printf("Subsets :\n{\n");
	for (auto i = allSubsets.begin() ; i != allSubsets.end() ; ++i )
	{
		printf("\t{ ");
		for (auto j = (*i).begin() ; j != (*i).end() ; j++)
		{
			if (j != (*i).begin()) { printf(", ");}
			//printf("%lld",(long long)*j);
			std::cout << *j;
		}
		printf("}");
		if (i != (--allSubsets.end())) { printf(",\n");}
	}
	printf("\n}\n");
}
std::set<std::set<igraph_integer_t>> generateSubsets(igraph_vector_int_t steinerTerminals, igraph_integer_t n, igraph_integer_t graphsize)
{
	igraph_integer_t count = ((igraph_integer_t) 1 << n);
	std::set<std::set<igraph_integer_t>> allSubsets;
	igraph_integer_t subsetIndex = graphsize;

	// The outer for loop will run 2^n times to print all subset .
	// Here variable i will act as a binary counter

	for (igraph_integer_t i = 0; i < count; i++)
	{
		// The inner for loop will run n times , As the maximum number of elements a set can have is n
		// This loop will generate a subset
		std::set<igraph_integer_t> newSubset;
		for (igraph_integer_t j = 0; j < n; j++)
		{
			// This if condition will check if jth bit in binary representation of  i  is set or not
			// if the value of (i & (1 << j)) is greater than 0 , include arr[j] in the current subset
			// otherwise exclude arr[j]
			if ((i & ((igraph_integer_t)1 << j)) > 0)
			{
				newSubset.insert(VECTOR(steinerTerminals)[j]);
			}
		}
		
		if (newSubset.size() > 1)  
		{
			if (allSubsets.find(newSubset) == allSubsets.end())
			{
				allSubsets.insert(newSubset);
				subsetMap.insert(std::make_pair(newSubset, subsetIndex));
				subsetIndex++;
			}
		}
	}
	return allSubsets;
}

igraph_integer_t fetchIndexofMapofSets(std::set<igraph_integer_t> subset)
{
	igraph_integer_t key;
	std::map<std::set<igraph_integer_t>, igraph_integer_t>::iterator it;
	for (it = subsetMap.begin(); it != subsetMap.end(); ++it)
	{
		if (it->first == subset)
		{
			key = it->second;
		}
	}

	return key;
}

igraph_error_t igraph_steiner_dreyfus_wagner(const igraph_t *graph,const igraph_vector_int_t* steiner_terminals,
igraph_i_directed_t mode, const igraph_vector_t *weights,igraph_real_t *res)
{

	igraph_integer_t no_of_vertices = (igraph_integer_t)igraph_vcount(graph);
	igraph_integer_t no_of_edges = igraph_ecount(graph);

	igraph_vector_int_t steiner_terminals_copy; 
	igraph_matrix_t dp_cache; // dynamic programming table
	igraph_integer_t q;
	std::set<std::set<igraph_integer_t>> allSubsets;
	igraph_matrix_t distance;
	if (no_of_edges == 0 && no_of_vertices == 0) //graph is empty
	{
		*res = 0;
		return IGRAPH_SUCCESS;
	}
	if (igraph_vector_size(weights) != no_of_edges)
	{	
		IGRAPH_ERRORF("Weight vector length does not match %" IGRAPH_PRId "vec size and %" IGRAPH_PRId "edges \n",IGRAPH_EINVAL,igraph_vector_size(weights), no_of_edges);
	}
	IGRAPH_CHECK(igraph_matrix_init(&distance,no_of_vertices,no_of_vertices));
	IGRAPH_FINALLY(igraph_matrix_destroy,&distance);

	igraph_distances_johnson(graph, &distance, igraph_vss_all(), igraph_vss_all(), weights);
	// for (long int i = 0 ; i  <  no_of_vertices; i++)
	// {
	// 	for (long int j = 0 ; j  <  no_of_vertices; j++)
	// 	{
	// 		std::cout << igraph_matrix_get(&distance, i,j) << " ";
	// 	}
	// 	std::cout << std::endl;
	// }
	//printf("Johnson Works\n");
	
	IGRAPH_CHECK(igraph_vector_int_init_copy(&steiner_terminals_copy,steiner_terminals));
	IGRAPH_FINALLY(igraph_vector_int_destroy,&steiner_terminals_copy);
	igraph_vector_int_sort(&steiner_terminals_copy);
	
	// Creating a vector of steiner vertices. steiner vertices = vertices in graph - steiner terminals

	IGRAPH_CHECK(igraph_matrix_init(&dp_cache, pow(2, igraph_vector_int_size(&steiner_terminals_copy)), no_of_vertices));
	IGRAPH_FINALLY(igraph_matrix_destroy,&dp_cache);

    igraph_matrix_fill(&dp_cache, IGRAPH_INFINITY);
	for (long int i = 0 ; i  <  no_of_vertices; i++)
	{
		for (long int j = 0 ; j  <  no_of_vertices; j++)
		{
			//std::cout << igraph_matrix_get(&distance, i,j) << " ";
			igraph_matrix_set(&dp_cache,i,j,igraph_matrix_get(&distance,i,j));
			//std::cout << igraph_matrix_get(&dp_cache, i,j) << " ";
		}
		//std::cout << std::endl;
	}

//	printf("Matrix Filled\n");
	
	q = VECTOR(steiner_terminals_copy)[0];

	igraph_vector_int_remove(&steiner_terminals_copy, 0);

	allSubsets = generateSubsets(steiner_terminals_copy, igraph_vector_int_size(&steiner_terminals_copy), no_of_vertices);

	for (igraph_integer_t m = 2; m <= igraph_vector_int_size(&steiner_terminals_copy); m++)
	{
		for (igraph_integer_t i = 0; i < (igraph_integer_t)allSubsets.size(); i++)
		{
			auto it = allSubsets.begin();
			std::advance(it,i);
			std::set<igraph_integer_t> D = *it;
			igraph_integer_t indexOfSubsetD;
			indexOfSubsetD = fetchIndexofMapofSets(D);

			for (igraph_integer_t j = 0; j < no_of_vertices; j++){
				MATRIX(dp_cache,indexOfSubsetD,j) = IGRAPH_INFINITY;
			}
			
			for (igraph_integer_t j = 0; j < no_of_vertices; j++)
			{
				igraph_real_t distance1 = IGRAPH_INFINITY;
				std::set<igraph_integer_t>::iterator subset_D_iterator;

				for (subset_D_iterator = D.begin(); subset_D_iterator != D.end(); subset_D_iterator++)
				{
					igraph_integer_t E = *subset_D_iterator;
					if (E != j) {
						igraph_integer_t distanceEJ = MATRIX(distance, E, j);
						//std::cout << "Distance EJ" << distanceEJ << std::endl;
					
						std::set<igraph_integer_t> DMinusE = D;

					//igraph_vector_remove(&DMinusE,E);

						for (std::set<igraph_integer_t>::iterator iter = DMinusE.begin(); iter != DMinusE.end();)
						{
							if (*iter == E)
							{
								iter = DMinusE.erase(iter);
								break;
							}
							++iter;
						}
						
						igraph_integer_t indexOfSubsetDMinusE;
						if (DMinusE.size() == 1){
							std::set<igraph_integer_t>::iterator node = DMinusE.begin();
							indexOfSubsetDMinusE = *node;
						}
						else {
							indexOfSubsetDMinusE = fetchIndexofMapofSets(DMinusE);
						}
						
						//std::cout << "Index:" << indexOfSubsetDMinusE << std::endl;
					//std::cout << "Matrix Data Addition" << MATRIX(dp_cache, indexOfSubsetDMinusE, j) + distanceEJ;
					
						if ((distanceEJ + MATRIX(dp_cache, indexOfSubsetDMinusE, j)) < distance1)
						{
							distance1 = distanceEJ + (MATRIX(dp_cache, indexOfSubsetDMinusE, j));
								
						}
					}
					
				}
				//std::cout <<"Distance - 1"<< distance1 << std::endl;
				
				for (igraph_integer_t k = 0; k < no_of_vertices; k++)
				{
					igraph_matrix_set(&dp_cache,indexOfSubsetD,k,std::min(MATRIX(dp_cache, indexOfSubsetD, k), MATRIX(distance, k, j) + distance1));
				}
			}
		}
	}
	
	igraph_real_t distance2 = IGRAPH_INFINITY;

	for (igraph_integer_t j = 0; j < no_of_vertices; j++)
	{
		igraph_real_t distance1 = IGRAPH_INFINITY;
		for (igraph_integer_t subset_C_iterator = 0; subset_C_iterator < igraph_vector_int_size(steiner_terminals); subset_C_iterator++)
		{
			igraph_integer_t F = VECTOR(steiner_terminals_copy)[subset_C_iterator];
			igraph_integer_t distanceFJ = MATRIX(distance, F, j);

			std::set<igraph_integer_t> CMinusF;

			for (igraph_integer_t k = 0; k < igraph_vector_int_size(steiner_terminals); k++)
			{

				if (VECTOR(steiner_terminals_copy)[k] != F)
				{
					CMinusF.insert(VECTOR(steiner_terminals_copy)[k]);
				}
			}

			igraph_integer_t indexOfSubsetCMinusF = fetchIndexofMapofSets(CMinusF);

			if (distanceFJ != 0 && (distanceFJ + (MATRIX(dp_cache, indexOfSubsetCMinusF, j)) < distance1))
			{
				distance1 = distanceFJ + (MATRIX(dp_cache, indexOfSubsetCMinusF, j));
				//std::cout << "u:" << distance1 << std::endl;	
			}

		}
		
		
		if ( q != j && MATRIX(distance, q, j) + distance1 < distance2)
		{
			distance2 = MATRIX(distance, q, j) + distance1;
			//std::cout << "Distance-2" << distance2 <<std::endl;
		}
	}
	*res = distance2;
	//std::cout << u << " " << v << std::endl;
	
	igraph_matrix_destroy(&distance);
	
	igraph_vector_int_destroy(&steiner_terminals_copy);
	
	igraph_matrix_destroy(&dp_cache);

	IGRAPH_FINALLY_CLEAN(3);

	return IGRAPH_SUCCESS;
	

}
