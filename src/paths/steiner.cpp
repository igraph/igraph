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

static void printSubsets(std::set<std::set<igraph_integer_t>> allSubsets)
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
	std::map<std::set<igraph_integer_t>, igraph_integer_t>::iterator it;
	for (it = subsetMap.begin(); it != subsetMap.end(); ++it)
	{
		if (it->first == subset)
		{
			return it->second;
			
		}
	}
	return -1;
}

std::set<igraph_integer_t> fetchSetsBasedonIndex(igraph_integer_t index)
{
	std::map<std::set<igraph_integer_t>, igraph_integer_t>::iterator it;
	for (it = subsetMap.begin(); it != subsetMap.end(); ++it)
	{
		if (it->second == index)
		{
			return it->first;
		}
	}
	return std::set<igraph_integer_t>();
}

igraph_integer_t factorial ( igraph_integer_t n)
{
	igraph_integer_t answer = 1;
	for (igraph_integer_t i = 1 ; i <= n; i ++)
	{
		answer *= i;
	}
	return answer;
}

igraph_integer_t Combination (igraph_integer_t n, igraph_integer_t r)
{
	return factorial(n)/ (factorial(n-r) * factorial(r));
}


igraph_integer_t findMinimumK (igraph_matrix_t* dp_cache, igraph_integer_t indexD, igraph_integer_t q ) {

	igraph_integer_t min_col_num = -1;
	igraph_integer_t min_sum_for_col;

	for (igraph_integer_t i = 0; i < dp_cache->ncol ; i++)
	{
		if ( q != i)
		{

			if (min_col_num == -1)
			{
				min_col_num = i;
				min_sum_for_col = (igraph_matrix_get(dp_cache,q,i) + igraph_matrix_get(dp_cache,indexD,i));
			}
			else if ((igraph_matrix_get(dp_cache,q,i) + igraph_matrix_get(dp_cache,indexD,i)) < min_sum_for_col )
			{
				min_col_num = i;
				min_sum_for_col = (igraph_matrix_get(dp_cache,q,i) + igraph_matrix_get(dp_cache,indexD,i));
			}
		}
	}

	return min_col_num;
}



static igraph_error_t generate_steiner_tree_appx(const igraph_t* graph,const igraph_vector_t *weights,
                                       igraph_matrix_t* dp_cache, igraph_integer_t indexD , igraph_integer_t q, igraph_neimode_t mode,igraph_vector_int_t *vectorlist_all,igraph_vector_int_t *edgelist_all)
{
	
	// igraph_integer_t combination_value  = Combination(SetD.size(), SetD.size() -1);
	
	std::set<igraph_integer_t> C = fetchSetsBasedonIndex(indexD);

	igraph_integer_t m = q;
	int len = C.size();
	std::set<igraph_integer_t> D = C;		

	while(D.size() > 1) {
	
		indexD = fetchIndexofMapofSets(D);
		igraph_integer_t k = findMinimumK(dp_cache,indexD,m);
		//std::cout << "K,m" << k << ' ' << m << std::endl;
		igraph_vector_int_t vectorlist;
		IGRAPH_CHECK(igraph_vector_int_init(&vectorlist,1));
		IGRAPH_FINALLY(igraph_vector_int_destroy,&vectorlist);

		igraph_vector_int_t edgelist;
		IGRAPH_CHECK(igraph_vector_int_init(&edgelist,1));
		IGRAPH_FINALLY(igraph_vector_int_destroy,&edgelist);

		igraph_get_shortest_path_dijkstra(graph,&vectorlist,&edgelist,m,k,weights,IGRAPH_ALL);

		//std::cout << "EdgeList" << std::endl;
		

		igraph_vector_int_append(vectorlist_all,&vectorlist);
		igraph_vector_int_append(edgelist_all,&edgelist);

		igraph_integer_t min_E_value = IGRAPH_INTEGER_MAX;
		std::set<igraph_integer_t> min_F;
		if (D.size() > 2) {
			igraph_integer_t numElementsScan = Combination(len,D.size() - 1);
			igraph_integer_t min_value = IGRAPH_INTEGER_MAX;
			
			igraph_integer_t holder = fetchIndexofMapofSets(D);
			for (igraph_integer_t i=1; i <= numElementsScan; i++){
			
				igraph_integer_t value = igraph_matrix_get(dp_cache,holder - i,k);
				std::set<igraph_integer_t> F = fetchSetsBasedonIndex(holder - i);

				std::set<igraph_integer_t> E;
				//std::cout << "Till Here - Fetching E" << std::endl;
				std::set_difference(D.begin(),D.end(),F.begin(),F.end(),std::inserter(E, E.end()));
				
				//std::cout << "E" << *E.begin() << std::endl;

				igraph_integer_t temp_value = igraph_matrix_get(dp_cache,k,*E.begin()) + igraph_matrix_get(dp_cache,k,holder - i) ;
				
				if (temp_value < min_value) {
					min_value = temp_value;
					min_E_value = *E.begin();
					min_F = F;
				}
			}

			igraph_vector_int_t vectorlist_1;
			IGRAPH_CHECK(igraph_vector_int_init(&vectorlist_1,1));
			IGRAPH_FINALLY(igraph_vector_int_destroy,&vectorlist_1);

			igraph_vector_int_t edgelist_1;
			IGRAPH_CHECK(igraph_vector_int_init(&edgelist_1,1));
			IGRAPH_FINALLY(igraph_vector_int_destroy,&edgelist_1);

			igraph_get_shortest_path_dijkstra(graph,&vectorlist_1,&edgelist_1,k,min_E_value,weights,IGRAPH_ALL);
			
			
			igraph_vector_int_append(vectorlist_all,&vectorlist_1);
			igraph_vector_int_append(edgelist_all,&edgelist_1);

		
			igraph_vector_int_destroy(&vectorlist_1);
			igraph_vector_int_destroy(&edgelist_1);
			IGRAPH_FINALLY_CLEAN(2);

		}
		else {
			
			igraph_integer_t E1,F1;

			E1 = *D.begin();
			F1 = *next(D.begin(),1);;	
			
			//std::cout << "E1:" <<E1 << std::endl;

			igraph_vector_int_t vectorlist_1;
			IGRAPH_CHECK(igraph_vector_int_init(&vectorlist_1,1));
			IGRAPH_FINALLY(igraph_vector_int_destroy,&vectorlist_1);

			igraph_vector_int_t edgelist_1;
			IGRAPH_CHECK(igraph_vector_int_init(&edgelist_1,1));
			IGRAPH_FINALLY(igraph_vector_int_destroy,&edgelist_1);

			igraph_get_shortest_path_dijkstra(graph,&vectorlist_1,&edgelist_1,k,E1,weights,IGRAPH_ALL);

			igraph_vector_int_t vectorlist_2;
			IGRAPH_CHECK(igraph_vector_int_init(&vectorlist_2,1));
			IGRAPH_FINALLY(igraph_vector_int_destroy,&vectorlist_2);

			igraph_vector_int_t edgelist_2;
			IGRAPH_CHECK(igraph_vector_int_init(&edgelist_2,1));
			IGRAPH_FINALLY(igraph_vector_int_destroy,&edgelist_2);

			igraph_get_shortest_path_dijkstra(graph,&vectorlist_2,&edgelist_2,k,F1,weights,IGRAPH_ALL);

			igraph_vector_int_append(vectorlist_all,&vectorlist_1);
			igraph_vector_int_append(vectorlist_all,&vectorlist_2);

			igraph_vector_int_append(edgelist_all,&edgelist_1);
			igraph_vector_int_append(edgelist_all,&edgelist_2);

			igraph_vector_int_destroy(&vectorlist_2);
			igraph_vector_int_destroy(&vectorlist_1);
			
			igraph_vector_int_destroy(&edgelist_1);
			igraph_vector_int_destroy(&edgelist_2);

			IGRAPH_FINALLY_CLEAN(4);

			std::set<igraph_integer_t> min_F;
			min_F.insert(F1); 

			//std::cout << "Got till here!" << std::endl;

		}
		
		m = k;
		D = min_F;			
		
		igraph_vector_int_destroy(&vectorlist);
		igraph_vector_int_destroy(&edgelist);
		IGRAPH_FINALLY_CLEAN(2);
	}
	


	return IGRAPH_SUCCESS;

}


igraph_error_t igraph_steiner_dreyfus_wagner(const igraph_t *graph,const igraph_vector_int_t* steiner_terminals,
igraph_neimode_t mode, const igraph_vector_t *weights,igraph_real_t *res,igraph_vector_int_t *res_tree)
{
	if (mode != IGRAPH_ALL)
	{
		std::cout << "Currently this function only supports undirected graphs while the graph's mode is not undirected." <<std::endl;
		return IGRAPH_FAILURE;
	}
	
	
	igraph_integer_t no_of_vertices = (igraph_integer_t)igraph_vcount(graph);
	igraph_integer_t no_of_edges = igraph_ecount(graph);

	// if (igraph_vector_int_size(steiner_terminals) == no_of_vertices)
	// {
	// 	std::cout << "Getting Minimum Spanning Tree" << std::endl;
	// 	igraph_error_t ans = igraph_minimum_spanning_tree(graph,res,weights);
	// 	return ans;
		
	// } Needs to be moved to Phase 2 when we implement backtracking to get the edges out
	if (no_of_vertices == 0 || (no_of_vertices == 1)) //graph is empty
	{
		*res = 0;
		return IGRAPH_FAILURE;
	}

	igraph_vector_int_t steiner_terminals_copy; 
	igraph_matrix_t dp_cache; // dynamic programming table
	igraph_integer_t q;
	std::set<std::set<igraph_integer_t>> allSubsets;
	igraph_matrix_t distance;
	
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
	for (igraph_integer_t i = 0; i < no_of_vertices; i++)
	{
		if (igraph_matrix_get(&distance,i,i) != 0)
		{
			igraph_matrix_set(&distance,i,i,0);
			std::cout <<"Found Self-loop at node number " << i 
					<< ". Ignoring the self-loop in this function."<< std::endl;
		}
	}
	IGRAPH_CHECK(igraph_vector_int_init_copy(&steiner_terminals_copy,steiner_terminals));
	IGRAPH_FINALLY(igraph_vector_int_destroy,&steiner_terminals_copy);
	igraph_vector_int_sort(&steiner_terminals_copy);
	q = VECTOR(steiner_terminals_copy)[0];

	igraph_vector_int_remove(&steiner_terminals_copy, 0);
	// Creating a vector of steiner vertices. steiner vertices = vertices in graph - steiner terminals
	IGRAPH_CHECK(igraph_matrix_init(&dp_cache,no_of_vertices + pow(2, igraph_vector_int_size(&steiner_terminals_copy) - 1), no_of_vertices));
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
		}
	}
	*res = distance2;
	//std::cout << u << " " << v << std::endl;
	// for (igraph_integer_t i = 0 ; i < no_of_vertices +  pow(2, igraph_vector_int_size(&steiner_terminals_copy) - 1) ; i ++)
	// {
	// 	for (igraph_integer_t j = 0 ; j < no_of_vertices ; j ++)
	// 	{
	// 		std::cout << igraph_matrix_get(&dp_cache,i,j) << " ";
	// 	}
	// std::cout << std::endl;
	// }
	std::set<igraph_integer_t> newSet;
	for (auto i = 0 ; i < igraph_vector_int_size(&steiner_terminals_copy); i++)
	{
		newSet.insert(VECTOR(steiner_terminals_copy)[i]);
	}
	igraph_integer_t indexD = fetchIndexofMapofSets(newSet);
	
	igraph_vector_int_t vectorlist_all;
	

	IGRAPH_CHECK(igraph_vector_int_init(&vectorlist_all,1));
	

	IGRAPH_CHECK(generate_steiner_tree_appx(graph,weights,&dp_cache,indexD,q,IGRAPH_ALL,&vectorlist_all,res_tree));
	igraph_vector_int_remove(res_tree,0);

	/*
	for (igraph_integer_t i =0; i< igraph_vector_int_size(&vectorlist_all); i++) {
		std::cout << VECTOR(vectorlist_all)[i] << std::endl;
	}

	std::cout << std::endl ;

	
	for (igraph_integer_t i =0; i< igraph_vector_int_size(&edgelist_all); i++) {
		std::cout << VECTOR(edgelist_all)[i] << std::endl;
	}
	*/


	igraph_matrix_destroy(&distance);	
	igraph_vector_int_destroy(&steiner_terminals_copy);
	
	igraph_matrix_destroy(&dp_cache);
	igraph_vector_int_destroy(&vectorlist_all);
	

	IGRAPH_FINALLY_CLEAN(3);

	return IGRAPH_SUCCESS;
	

}


