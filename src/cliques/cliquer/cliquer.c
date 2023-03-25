
/*
 * This file contains the clique searching routines.
 *
 * Copyright (C) 2002 Sampo Niskanen, Patric Östergård.
 * Licensed under the GNU GPL, read the file LICENSE for details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
/*
#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>
*/

#include "cliquer.h"

#include "config.h"

/* Default cliquer options */
IGRAPH_THREAD_LOCAL clique_options clique_default_options = {
	reorder_by_default, NULL, /*clique_print_time*/ NULL, NULL, NULL, NULL, NULL, 0
};


/* Calculate d/q, rounding result upwards/downwards. */
#define DIV_UP(d,q) (((d)+(q)-1)/(q))
#define DIV_DOWN(d,q) ((int)((d)/(q)))


/* Global variables used: */
/* These must be saved and restored in re-entrance. */
static IGRAPH_THREAD_LOCAL int *clique_size;      /* c[i] == max. clique size in {0,1,...,i-1} */
static IGRAPH_THREAD_LOCAL set_t current_clique;  /* Current clique being searched. */
static IGRAPH_THREAD_LOCAL set_t best_clique;     /* Largest/heaviest clique found so far. */
/*static struct tms cputimer;*/      /* Timer for opts->time_function() */
/*static struct timeval realtimer;*/ /* Timer for opts->time_function() */
static IGRAPH_THREAD_LOCAL int clique_list_count=0;  /* No. of cliques in opts->clique_list[] */
static IGRAPH_THREAD_LOCAL int weight_multiplier=1;  /* Weights multiplied by this when passing
				  * to time_function(). */

/* List cache (contains memory blocks of size g->n * sizeof(int)) */
static IGRAPH_THREAD_LOCAL int **temp_list=NULL;
static IGRAPH_THREAD_LOCAL int temp_count=0;


/*
 * Macros for re-entrance.  ENTRANCE_SAVE() must be called immediately
 * after variable definitions, ENTRANCE_RESTORE() restores global
 * variables to original values.  entrance_level should be increased
 * and decreased accordingly.
 */
static IGRAPH_THREAD_LOCAL int entrance_level=0;  /* How many levels for entrance have occurred? */

#define ENTRANCE_SAVE() \
int *old_clique_size = clique_size;                     \
set_t old_current_clique = current_clique;              \
set_t old_best_clique = best_clique;                    \
int old_clique_list_count = clique_list_count;          \
int old_weight_multiplier = weight_multiplier;          \
int **old_temp_list = temp_list;                        \
int old_temp_count = temp_count;                        \
/*struct tms old_cputimer;                                \
struct timeval old_realtimer;                           \
memcpy(&old_cputimer,&cputimer,sizeof(struct tms));       \
memcpy(&old_realtimer,&realtimer,sizeof(struct timeval));*/

#define ENTRANCE_RESTORE() \
clique_size = old_clique_size;                          \
current_clique = old_current_clique;                    \
best_clique = old_best_clique;                          \
clique_list_count = old_clique_list_count;              \
weight_multiplier = old_weight_multiplier;              \
temp_list = old_temp_list;                              \
temp_count = old_temp_count;                            \
/*memcpy(&cputimer,&old_cputimer,sizeof(struct tms));       \
memcpy(&realtimer,&old_realtimer,sizeof(struct timeval));*/


/* Number of clock ticks per second (as returned by sysconf(_SC_CLK_TCK)) */
/*static int clocks_per_sec=0;*/




/* Recursion and helper functions */
static boolean sub_unweighted_single(int *table, int size, int min_size,
				     graph_t *g);
static igraph_error_t sub_unweighted_all(int *table, int size, int min_size, int max_size,
                                            boolean maximal, graph_t *g,
                                            clique_options *opts, CLIQUER_LARGE_INT *num_found);
static igraph_error_t sub_weighted_all(int *table, int size, int weight,
			    int current_weight, int prune_low, int prune_high,
			    int min_weight, int max_weight, boolean maximal,
			    graph_t *g, clique_options *opts, int *weight_found);


static igraph_error_t store_clique(set_t clique, graph_t *g, clique_options *opts);
static boolean is_maximal(set_t clique, graph_t *g);
static igraph_error_t false_function(set_t clique,graph_t *g,clique_options *opts);





/*****  Unweighted searches  *****/
/*
 * Unweighted searches are done separately from weighted searches because
 * some effective pruning methods can be used when the vertex weights
 * are all 1.  Single and all clique finding routines are separated,
 * because the single clique finding routine can store the found clique
 * while it is returning from the recursion, thus requiring no implicit
 * storing of the current clique.  When searching for all cliques the
 * current clique must be stored.
 */


/*
 * unweighted_clique_search_single()
 *
 * Searches for a single clique of size min_size.  Stores maximum clique
 * sizes into clique_size[].
 *
 *   table    - the order of the vertices in g to use
 *   min_size - minimum size of clique to search for.  If min_size==0,
 *              searches for a maximum clique.
 *   g        - the graph
 *   opts     - time printing options
 *
 * opts->time_function is called after each base-level recursion, if
 * non-NULL.
 *
 * Returns the size of the clique found, or 0 if min_size>0 and a clique
 * of that size was not found (or if time_function aborted the search).
 * The largest clique found is stored in current_clique.
 *
 * Note: Does NOT use opts->user_function of opts->clique_list.
 */
static int unweighted_clique_search_single(int *table, int min_size,
					   graph_t *g, clique_options *opts) {
	/*
    struct tms tms;
	struct timeval timeval;
	*/
	int i,j;
	int v,w;
	int *newtable;
	int newsize;

	v=table[0];
	clique_size[v]=1;
	set_empty(current_clique);
	SET_ADD_ELEMENT(current_clique,v);
	if (min_size==1)
		return 1;

	if (temp_count) {
		temp_count--;
		newtable=temp_list[temp_count];
	} else {
		newtable=malloc(g->n * sizeof(int));
	}
	for (i=1; i < g->n; i++) {
		w=v;
		v=table[i];

		newsize=0;
		for (j=0; j<i; j++) {
			if (GRAPH_IS_EDGE(g, v, table[j])) {
				newtable[newsize]=table[j];
				newsize++;
			}
		}

		if (sub_unweighted_single(newtable,newsize,clique_size[w],g)) {
			SET_ADD_ELEMENT(current_clique,v);
			clique_size[v]=clique_size[w]+1;
		} else {
			clique_size[v]=clique_size[w];
		}

		/*
		if (opts && opts->time_function) {
			gettimeofday(&timeval,NULL);
			times(&tms);
			if (!opts->time_function(entrance_level,
						 i+1,g->n,clique_size[v] *
						 weight_multiplier,
						 (double)(tms.tms_utime-
							  cputimer.tms_utime)/
						 clocks_per_sec,
						 timeval.tv_sec-
						 realtimer.tv_sec+
						 (double)(timeval.tv_usec-
							  realtimer.tv_usec)/
						 1000000,opts)) {
				temp_list[temp_count++]=newtable;
				return 0;
			}
		}
		*/

		if (min_size) {
			if (clique_size[v]>=min_size) {
				temp_list[temp_count++]=newtable;
				return clique_size[v];
			}
			if (clique_size[v]+g->n-i-1 < min_size) {
				temp_list[temp_count++]=newtable;
				return 0;
			}
		}
	}

	temp_list[temp_count++]=newtable;

	if (min_size)
		return 0;
	return clique_size[v];
}

/*
 * sub_unweighted_single()
 *
 * Recursion function for searching for a single clique of size min_size.
 *
 *    table    - subset of the vertices in graph
 *    size     - size of table
 *    min_size - size of clique to look for within the subgraph
 *               (decreased with every recursion)
 *    g        - the graph
 *
 * Returns TRUE if a clique of size min_size is found, FALSE otherwise.
 * If a clique of size min_size is found, it is stored in current_clique.
 *
 * clique_size[] for all values in table must be defined and correct,
 * otherwise inaccurate results may occur.
 */
static boolean sub_unweighted_single(int *table, int size, int min_size,
				     graph_t *g) {
	int i;
	int v;
	int *newtable;
	int *p1, *p2;

	/* Zero or one vertices needed anymore. */
	if (min_size <= 1) {
		if (size>0 && min_size==1) {
			set_empty(current_clique);
			SET_ADD_ELEMENT(current_clique,table[0]);
			return TRUE;
		}
		if (min_size==0) {
			set_empty(current_clique);
			return TRUE;
		}
		return FALSE;
	}
	if (size < min_size)
		return FALSE;

	/* Dynamic memory allocation with cache */
	if (temp_count) {
		temp_count--;
		newtable=temp_list[temp_count];
	} else {
		newtable=malloc(g->n * sizeof(int));
	}

	for (i = size-1; i >= 0; i--) {
		v = table[i];

		if (clique_size[v] < min_size)
			break;
		/* This is faster when compiling with gcc than placing
		 * this in the for-loop condition. */
		if (i+1 < min_size)
			break;

		/* Very ugly code, but works faster than "for (i=...)" */
		p1 = newtable;
		for (p2=table; p2 < table+i; p2++) {
			int w = *p2;
			if (GRAPH_IS_EDGE(g, v, w)) {
				*p1 = w;
				p1++;
			}
		}

		/* Avoid unnecessary loops (next size == p1-newtable) */
		if (p1-newtable < min_size-1)
			continue;
		/* Now p1-newtable >= min_size-1 >= 2-1 == 1, so we can use
		 * p1-newtable-1 safely. */
		if (clique_size[newtable[p1-newtable-1]] < min_size-1)
			continue;

		if (sub_unweighted_single(newtable,p1-newtable,
					  min_size-1,g)) {
			/* Clique found. */
			SET_ADD_ELEMENT(current_clique,v);
			temp_list[temp_count++]=newtable;
			return TRUE;
		}
	}
	temp_list[temp_count++]=newtable;
	return FALSE;
}


/*
 * unweighted_clique_search_all()
 *
 * Searches for all cliques with size at least min_size and at most
 * max_size.  Stores the cliques as opts declares.
 *
 *   table    - the order of the vertices in g to search
 *   start    - first index where the subgraph table[0], ..., table[start]
 *              might include a requested kind of clique
 *   min_size - minimum size of clique to search for.  min_size > 0 !
 *   max_size - maximum size of clique to search for.  If no upper limit
 *              is desired, use eg. INT_MAX
 *   maximal  - requires cliques to be maximal
 *   g        - the graph
 *   opts     - time printing and clique storage options
 *   num_found - number of cliques found
 *
 * Cliques found are stored as defined by opts->user_function and
 * opts->clique_list.  opts->time_function is called after each
 * base-level recursion, if non-NULL.
 *
 * clique_size[] must be defined and correct for all values of
 * table[0], ..., table[start-1].
 */
static igraph_error_t unweighted_clique_search_all(
	int *table, int start, int min_size, int max_size, boolean maximal,
	graph_t *g, clique_options *opts, CLIQUER_LARGE_INT *num_found
) {
	/*
	struct timeval timeval;
	struct tms tms;
	*/
	int i, j;
	int v;
	int *newtable;
	int newsize;
	CLIQUER_LARGE_INT r;
	CLIQUER_LARGE_INT count=0;
	igraph_error_t retval = IGRAPH_SUCCESS;

	if (temp_count) {
		temp_count--;
		newtable=temp_list[temp_count];
	} else {
		newtable=malloc(g->n * sizeof(int));
	}

	clique_list_count=0;
	set_empty(current_clique);
	for (i=start; i < g->n; i++) {
		v=table[i];
		clique_size[v]=min_size;  /* Do not prune here. */

		newsize=0;
		for (j=0; j<i; j++) {
			if (GRAPH_IS_EDGE(g,v,table[j])) {
				newtable[newsize]=table[j];
				newsize++;
			}
		}

		SET_ADD_ELEMENT(current_clique,v);
		retval=sub_unweighted_all(newtable,newsize,min_size-1,max_size-1,
				maximal,g,opts,&r);
		SET_DEL_ELEMENT(current_clique,v);
		count+=r;
		if (retval) {
			/* Abort. */
			break;
		}

#if 0
		if (opts->time_function) {
			gettimeofday(&timeval,NULL);
			times(&tms);
			if (!opts->time_function(entrance_level,
						 i+1,g->n,min_size *
						 weight_multiplier,
						 (double)(tms.tms_utime-
							  cputimer.tms_utime)/
						 clocks_per_sec,
						 timeval.tv_sec-
						 realtimer.tv_sec+
						 (double)(timeval.tv_usec-
							  realtimer.tv_usec)/
						 1000000,opts)) {
				/* Abort. */
				break;
			}
		}
#endif
	}
	temp_list[temp_count++]=newtable;

	if (num_found) {
		*num_found = count;
	}

	return retval;
}

/*
 * sub_unweighted_all()
 *
 * Recursion function for searching for all cliques of given size.
 *
 *   table    - subset of vertices of graph g
 *   size     - size of table
 *   min_size - minimum size of cliques to search for (decreased with
 *              every recursion)
 *   max_size - maximum size of cliques to search for (decreased with
 *              every recursion).  If no upper limit is desired, use
 *              eg. INT_MAX
 *   maximal  - require cliques to be maximal (passed through)
 *   g        - the graph
 *   opts     - storage options
 *   num_found - number of cliques found
 *
 * All cliques of suitable size found are stored according to opts.
 *
 * Returns the number of cliques found.  If user_function returns FALSE,
 * then the number of cliques is returned negative.
 *
 * Uses current_clique to store the currently-being-searched clique.
 * clique_size[] for all values in table must be defined and correct,
 * otherwise inaccurate results may occur.
 */
static igraph_error_t sub_unweighted_all(int *table, int size, int min_size, int max_size,
                                            boolean maximal, graph_t *g,
                                            clique_options *opts, CLIQUER_LARGE_INT *num_found) {
	igraph_error_t retval = IGRAPH_SUCCESS;
	int i;
	int v;
	int *newtable;
	int *p1, *p2;
	CLIQUER_LARGE_INT n;
	CLIQUER_LARGE_INT count=0;     /* Amount of cliques found */

	if (min_size <= 0) {
		if ((!maximal) || is_maximal(current_clique,g)) {
			/* We've found one.  Store it. */
			count++;
			retval = store_clique(current_clique, g, opts);
			if (retval) {
				*num_found = count;
				return retval == IGRAPH_STOP ? IGRAPH_SUCCESS : retval;
			}
		}
		if (max_size <= 0) {
			/* If we add another element, size will be too big. */
			*num_found = count;
			return IGRAPH_SUCCESS;
		}
	}

	if (size < min_size) {
		*num_found = count;
		return IGRAPH_SUCCESS;
	}

	/* Dynamic memory allocation with cache */
	if (temp_count) {
		temp_count--;
		newtable=temp_list[temp_count];
	} else {
		newtable=malloc(g->n * sizeof(int));
	}

	for (i=size-1; i>=0; i--) {
		v = table[i];
		if (clique_size[v] < min_size) {
			break;
		}
		if (i+1 < min_size) {
			break;
		}

		/* Very ugly code, but works faster than "for (i=...)" */
		p1 = newtable;
		for (p2=table; p2 < table+i; p2++) {
			int w = *p2;
			if (GRAPH_IS_EDGE(g, v, w)) {
				*p1 = w;
				p1++;
			}
		}

		/* Avoid unnecessary loops (next size == p1-newtable) */
		if (p1-newtable < min_size-1) {
			continue;
		}

		SET_ADD_ELEMENT(current_clique,v);
		retval = sub_unweighted_all(newtable,p1-newtable,
				     min_size-1,max_size-1,maximal,g,opts,&n);
		SET_DEL_ELEMENT(current_clique,v);
		count += n;
		if (retval || n < 0) {
			break;
		}
		count+=n;
	}
	temp_list[temp_count++]=newtable;

	*num_found = count;
	return retval;
}




/***** Weighted clique searches *****/
/*
 * Weighted clique searches can use the same recursive routine, because
 * in both cases (single/all) they have to search through all potential
 * permutations searching for heavier cliques.
 */


/*
 * weighted_clique_search_single()
 *
 * Searches for a single clique of weight at least min_weight, and at
 * most max_weight.  Stores maximum clique sizes into clique_size[]
 * (or min_weight-1, whichever is smaller).
 *
 *   table      - the order of the vertices in g to use
 *   min_weight - minimum weight of clique to search for.  If min_weight==0,
 *                then searches for a maximum weight clique
 *   max_weight - maximum weight of clique to search for.  If no upper limit
 *                is desired, use eg. INT_MAX
 *   g          - the graph
 *   opts       - time printing options
 *
 * opts->time_function is called after each base-level recursion, if
 * non-NULL.
 *
 * Returns 0 if a clique of requested weight was not found (also if
 * time_function requested an abort), otherwise returns >= 1.
 * If min_weight==0 (search for maximum-weight clique), then the return
 * value is the weight of the clique found.  The found clique is stored
 * in best_clique.
 *
 * Note: Does NOT use opts->user_function of opts->clique_list.
 */
static igraph_error_t weighted_clique_search_single(int *table, int min_weight,
					 int max_weight, graph_t *g,
					 clique_options *opts, int *result) {
	/*
	struct timeval timeval;
	struct tms tms;
	*/
	int i,j;
	int v;
	int *newtable;
	int newsize;
	int newweight;
	int search_weight;
	int min_w;
	clique_options localopts;
	igraph_error_t retval = IGRAPH_SUCCESS;

	ASSERT(result != NULL);

	if (min_weight==0)
		min_w=INT_MAX;
	else
		min_w=min_weight;


	if (min_weight==1) {
		/* min_weight==1 may cause trouble in the routine, and
		 * it's trivial to check as it's own case.
		 * We write nothing to clique_size[]. */
		for (i=0; i < g->n; i++) {
			if (g->weights[table[i]] <= max_weight) {
				set_empty(best_clique);
				SET_ADD_ELEMENT(best_clique,table[i]);
				*result = g->weights[table[i]];
				return IGRAPH_SUCCESS;
			}
		}

		*result = 0;
		return IGRAPH_SUCCESS;
	}

	localopts.time_function=NULL;
	localopts.reorder_function=NULL;
	localopts.reorder_map=NULL;
	localopts.user_function=false_function;
	localopts.user_data=NULL;
	localopts.clique_list=&best_clique;
	localopts.clique_list_length=1;
	clique_list_count=0;

	v=table[0];
	set_empty(best_clique);
	SET_ADD_ELEMENT(best_clique,v);
	search_weight=g->weights[v];
	if (min_weight && (search_weight >= min_weight)) {
		if (search_weight <= max_weight) {
			/* Found suitable clique. */
			*result = search_weight;
			return IGRAPH_SUCCESS;
		}
		search_weight=min_weight-1;
	}
	clique_size[v]=search_weight;
	set_empty(current_clique);

	if (temp_count) {
		temp_count--;
		newtable=temp_list[temp_count];
	} else {
		newtable=malloc(g->n * sizeof(int));
	}

	for (i = 1; i < g->n; i++) {
		v=table[i];

		newsize=0;
		newweight=0;
		for (j=0; j<i; j++) {
			if (GRAPH_IS_EDGE(g,v,table[j])) {
				newweight += g->weights[table[j]];
				newtable[newsize]=table[j];
				newsize++;
			}
		}


		SET_ADD_ELEMENT(current_clique,v);
		retval=sub_weighted_all(newtable,newsize,newweight,
					       g->weights[v],search_weight,
					       clique_size[table[i-1]] +
					       g->weights[v],
					       min_w,max_weight,FALSE,
					       g,&localopts, &search_weight);
		SET_DEL_ELEMENT(current_clique,v);
		if (retval || search_weight < 0) {
			break;
		}

		clique_size[v]=search_weight;

        /*
		if (opts->time_function) {
			gettimeofday(&timeval,NULL);
			times(&tms);
			if (!opts->time_function(entrance_level,
						 i+1,g->n,clique_size[v] *
						 weight_multiplier,
						 (double)(tms.tms_utime-
							  cputimer.tms_utime)/
						 clocks_per_sec,
						 timeval.tv_sec-
						 realtimer.tv_sec+
						 (double)(timeval.tv_usec-
							  realtimer.tv_usec)/
						 1000000,opts)) {
				set_free(current_clique);
				current_clique=NULL;
				break;
			}
		}
        */
	}
	temp_list[temp_count++]=newtable;
	if (min_weight && (search_weight > 0)) {
		/* Requested clique has not been found. */
		*result = 0;
	} else {
		*result = clique_size[table[i-1]];
	}

	return retval;
}


/*
 * weighted_clique_search_all()
 *
 * Searches for all cliques with weight at least min_weight and at most
 * max_weight.  Stores the cliques as opts declares.
 *
 *   table      - the order of the vertices in g to search
 *   start      - first index where the subgraph table[0], ..., table[start]
 *                might include a requested kind of clique
 *   min_weight - minimum weight of clique to search for.  min_weight > 0 !
 *   max_weight - maximum weight of clique to search for.  If no upper limit
 *                is desired, use eg. INT_MAX
 *   maximal    - search only for maximal cliques
 *   g          - the graph
 *   opts       - time printing and clique storage options
 *   num_found - number of cliques found
 *
 * Cliques found are stored as defined by opts->user_function and
 * opts->clique_list.  opts->time_function is called after each
 * base-level recursion, if non-NULL.
 *
 * clique_size[] must be defined and correct for all values of
 * table[0], ..., table[start-1].
 *
 * Returns the number of cliques stored (not necessarily number of cliques
 * in graph, if user/time_function aborts).
 */
static igraph_error_t weighted_clique_search_all(int *table, int start,
				      int min_weight, int max_weight,
				      boolean maximal, graph_t *g,
				      clique_options *opts, int* num_found) {
    /*
	struct timeval timeval;
	struct tms tms;
    */
	int i,j;
	int v;
	int *newtable;
	int newsize;
	int newweight;
	igraph_error_t retval = IGRAPH_SUCCESS;

	if (temp_count) {
		temp_count--;
		newtable=temp_list[temp_count];
	} else {
		newtable=malloc(g->n * sizeof(int));
	}

	clique_list_count=0;
	set_empty(current_clique);
	for (i=start; i < g->n; i++) {
		v=table[i];
		clique_size[v]=min_weight;   /* Do not prune here. */

		newsize=0;
		newweight=0;
		for (j=0; j<i; j++) {
			if (GRAPH_IS_EDGE(g,v,table[j])) {
				newtable[newsize]=table[j];
				newweight+=g->weights[table[j]];
				newsize++;
			}
		}

		SET_ADD_ELEMENT(current_clique,v);
		retval=sub_weighted_all(newtable,newsize,newweight,
				   g->weights[v],min_weight-1,INT_MAX,
				   min_weight,max_weight,maximal,g,opts,&j);
		SET_DEL_ELEMENT(current_clique,v);

		if (retval || j < 0) {
			/* Abort. */
			break;
		}

        /*
		if (opts->time_function) {
			gettimeofday(&timeval,NULL);
			times(&tms);
			if (!opts->time_function(entrance_level,
						 i+1,g->n,clique_size[v] *
						 weight_multiplier,
						 (double)(tms.tms_utime-
							  cputimer.tms_utime)/
						 clocks_per_sec,
						 timeval.tv_sec-
						 realtimer.tv_sec+
						 (double)(timeval.tv_usec-
							  realtimer.tv_usec)/
						 1000000,opts)) {
				set_free(current_clique);
				current_clique=NULL;
				break;
			}
		}
        */
	}
	temp_list[temp_count++]=newtable;

	if (num_found) {
		*num_found = clique_list_count;
	}

	return retval;
}

/*
 * sub_weighted_all()
 *
 * Recursion function for searching for all cliques of given weight.
 *
 *   table      - subset of vertices of graph g
 *   size       - size of table
 *   weight     - total weight of vertices in table
 *   current_weight - weight of clique found so far
 *   prune_low  - ignore all cliques with weight less or equal to this value
 *                (often heaviest clique found so far)  (passed through)
 *   prune_high - maximum weight possible for clique in this subgraph
 *                (passed through)
 *   min_size   - minimum weight of cliques to search for (passed through)
 *                Must be greater than 0.
 *   max_size   - maximum weight of cliques to search for (passed through)
 *                If no upper limit is desired, use eg. INT_MAX
 *   maximal    - search only for maximal cliques
 *   g          - the graph
 *   opts       - storage options
 *   weight_found - weight of the heaviest clique found (prune_low if a heavier
 *                clique hasn't been found); if a clique with weight at least
 *                min_size is found then min_size-1 is returned.
 *
 * All cliques of suitable weight found are stored according to opts.
 *
 * The largest clique found smaller than max_weight is stored in
 * best_clique, if non-NULL.
 *
 * Uses current_clique to store the currently-being-searched clique.
 * clique_size[] for all values in table must be defined and correct,
 * otherwise inaccurate results may occur.
 *
 * To search for a single maximum clique, use min_weight==max_weight==INT_MAX,
 * with best_clique non-NULL.  To search for a single given-weight clique,
 * use opts->clique_list and opts->user_function=false_function.  When
 * searching for all cliques, min_weight should be given the minimum weight
 * desired.
 */
static igraph_error_t sub_weighted_all(int *table, int size, int weight,
			    int current_weight, int prune_low, int prune_high,
			    int min_weight, int max_weight, boolean maximal,
			    graph_t *g, clique_options *opts, int* weight_found) {
	igraph_error_t retval = IGRAPH_SUCCESS;
	int i;
	int v,w;
	int *newtable;
	int *p1, *p2;
	int newweight;

	if (current_weight >= min_weight) {
		if ((current_weight <= max_weight) &&
		    ((!maximal) || is_maximal(current_clique,g))) {
			/* We've found one.  Store it. */
			retval = store_clique(current_clique,g,opts);
			if (retval) {
				*weight_found = -1;
				return retval == IGRAPH_STOP ? IGRAPH_SUCCESS : retval;
			}
		}
		if (current_weight >= max_weight) {
			/* Clique too heavy. */
			*weight_found = min_weight-1;
			return IGRAPH_SUCCESS;
		}
	}
	if (size <= 0) {
		/* current_weight < min_weight, prune_low < min_weight,
		 * so return value is always < min_weight. */
		if (current_weight>prune_low) {
			if (best_clique) {
				best_clique = set_copy(best_clique,current_clique);
			}
			if (current_weight < min_weight) {
				*weight_found = current_weight;
				return IGRAPH_SUCCESS;
			} else {
				*weight_found = min_weight-1;
				return IGRAPH_SUCCESS;
			}
		} else {
			*weight_found = prune_low;
			return IGRAPH_SUCCESS;
		}
	}

	/* Dynamic memory allocation with cache */
	if (temp_count) {
		temp_count--;
		newtable=temp_list[temp_count];
	} else {
		newtable=malloc(g->n * sizeof(int));
	}

	for (i = size-1; i >= 0; i--) {
		v = table[i];
		if (current_weight+clique_size[v] <= prune_low) {
			/* Dealing with subset without heavy enough clique. */
			break;
		}
		if (current_weight+weight <= prune_low) {
			/* Even if all elements are added, won't do. */
			break;
		}

		/* Very ugly code, but works faster than "for (i=...)" */
		p1 = newtable;
		newweight = 0;
		for (p2=table; p2 < table+i; p2++) {
			w = *p2;
			if (GRAPH_IS_EDGE(g, v, w)) {
				*p1 = w;
				newweight += g->weights[w];
				p1++;
			}
		}

		w=g->weights[v];
		weight-=w;
		/* Avoid a few unnecessary loops */
		if (current_weight+w+newweight <= prune_low) {
			continue;
		}

		SET_ADD_ELEMENT(current_clique,v);
		retval=sub_weighted_all(newtable,p1-newtable,
					   newweight,
					   current_weight+w,
					   prune_low,prune_high,
					   min_weight,max_weight,maximal,
					   g,opts, &prune_low);
		SET_DEL_ELEMENT(current_clique,v);
		if (retval || (prune_low<0) || (prune_low>=prune_high)) {
			/* Impossible to find larger clique. */
			break;
		}
	}
	temp_list[temp_count++]=newtable;

	*weight_found = prune_low;
	return IGRAPH_SUCCESS;
}




/***** Helper functions *****/


/*
 * store_clique()
 *
 * Stores a clique according to given user options.
 *
 *   clique - the clique to store
 *   opts   - storage options
 *
 * Returns the same igraph error code as the one returned by
 * opts->user_function(). Returns IGRAPH_SUCCESS if no callback is defined.
 */
static igraph_error_t store_clique(set_t clique, graph_t *g, clique_options *opts) {

	clique_list_count++;

	/* clique_list[] */
	if (opts->clique_list) {
		/*
		 * This has been a major source of bugs:
		 * Has clique_list_count been set to 0 before calling
		 * the recursions?
		 */
		if (clique_list_count <= 0) {
			IGRAPH_FATAL("CLIQUER INTERNAL ERROR: clique_list_count has negative value! Please report as a bug.");
		}
		if (clique_list_count <= opts->clique_list_length)
			opts->clique_list[clique_list_count-1] =
				set_copy(opts->clique_list[clique_list_count-1], clique);
	}

	/* user_function() */
	if (opts->user_function) {
		return opts->user_function(clique, g, opts);
	}

	return IGRAPH_SUCCESS;
}

/*
 * maximalize_clique()
 *
 * Adds greedily all possible vertices in g to set s to make it a maximal
 * clique.
 *
 *   s - clique of vertices to make maximal
 *   g - graph
 *
 * Note: Not very optimized (uses a simple O(n^2) routine), but is called
 *       at maximum once per clique_xxx() call, so it shouldn't matter.
 */
static void maximalize_clique(set_t s,graph_t *g) {
	int i,j;
	boolean add;

	for (i=0; i < g->n; i++) {
		add=TRUE;
		for (j=0; j < g->n; j++) {
			if (SET_CONTAINS_FAST(s,j) && !GRAPH_IS_EDGE(g,i,j)) {
				add=FALSE;
				break;
			}
		}
		if (add) {
			SET_ADD_ELEMENT(s,i);
		}
	}
	return;
}


/*
 * is_maximal()
 *
 * Check whether a clique is maximal or not.
 *
 *   clique - set of vertices in clique
 *   g      - graph
 *
 * Returns TRUE is clique is a maximal clique of g, otherwise FALSE.
 */
static boolean is_maximal(set_t clique, graph_t *g) {
	int i,j;
	int *table;
	int len;
	boolean addable;

	if (temp_count) {
		temp_count--;
		table=temp_list[temp_count];
	} else {
		table=malloc(g->n * sizeof(int));
	}

	len=0;
	for (i=0; i < g->n; i++)
		if (SET_CONTAINS_FAST(clique,i))
			table[len++]=i;

	for (i=0; i < g->n; i++) {
		addable=TRUE;
		for (j=0; j<len; j++) {
			if (!GRAPH_IS_EDGE(g,i,table[j])) {
				addable=FALSE;
				break;
			}
		}
		if (addable) {
			temp_list[temp_count++]=table;
			return FALSE;
		}
	}
	temp_list[temp_count++]=table;
	return TRUE;
}


/*
 * false_function()
 *
 * Returns IGRAPH_STOP.  Can be used as user_function.
 */
static igraph_error_t false_function(set_t clique,graph_t *g,clique_options *opts) {
	return IGRAPH_STOP;
}




/***** API-functions *****/

/*
 * clique_unweighted_max_weight()
 *
 *   g    - the graph
 *   opts - time printing options
 *   size - the size of the maximum (sized) clique in g
 *
 * Note: As we don't have an algorithm faster than actually finding
 *       a maximum clique, we use clique_unweighted_find_single().
 *       This incurs only very small overhead.
 */
igraph_error_t clique_unweighted_max_weight(graph_t *g, clique_options *opts, int* size) {
	set_t clique;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);

	IGRAPH_CHECK(clique_unweighted_find_single(g, 0, 0, FALSE,opts, &clique));

	if (size != NULL) {
		*size = clique ? set_size(clique) : 0;
	}

	if (clique) {
		set_free(clique);
	}

	return IGRAPH_SUCCESS;
}


/*
 * clique_unweighted_find_single()
 *
 * Returns a clique with size at least min_size and at most max_size.
 *
 *   g        - the graph
 *   min_size - minimum size of clique to search for.  If min_size==0,
 *              searches for maximum clique.
 *   max_size - maximum size of clique to search for.  If max_size==0, no
 *              upper limit is used.  If min_size==0, this must also be 0.
 *   maximal  - require returned clique to be maximal
 *   opts     - time printing options
 *
 * Returns the set of vertices forming the clique, or NULL if a clique
 * of requested size/maximality does not exist in the graph  (or if
 * opts->time_function() requests abort).
 *
 * The returned clique is newly allocated and can be freed by set_free().
 *
 * Note: Does NOT use opts->user_function() or opts->clique_list[].
 */
igraph_error_t clique_unweighted_find_single(graph_t *g,int min_size,int max_size,
				    boolean maximal, clique_options *opts, set_t *clique) {
	int i;
	int *table;
	set_t s;
	igraph_error_t retval = IGRAPH_SUCCESS;
	CLIQUER_LARGE_INT found;

	ENTRANCE_SAVE();
	entrance_level++;

	if (opts==NULL)
		opts=&clique_default_options;

	ASSERT(clique!=NULL);
	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);
	ASSERT(min_size>=0);
	ASSERT(max_size>=0);
	ASSERT((max_size==0) || (min_size <= max_size));
	ASSERT(!((min_size==0) && (max_size>0)));
	ASSERT((opts->reorder_function==NULL) || (opts->reorder_map==NULL));

	if ((max_size>0) && (min_size>max_size)) {
		/* state was not changed */
		entrance_level--;
		*clique = NULL;
		return IGRAPH_SUCCESS;
	}

    /*
	if (clocks_per_sec==0)
		clocks_per_sec=sysconf(_SC_CLK_TCK);
	ASSERT(clocks_per_sec>0);
    */

	/* Dynamic allocation */
	current_clique=set_new(g->n);
	clique_size=malloc(g->n * sizeof(int));
	/* table allocated later */
	temp_list=malloc((g->n+2)*sizeof(int *));
	temp_count=0;

	/* "start clock" */
    /*
	gettimeofday(&realtimer,NULL);
	times(&cputimer);
    */

	/* reorder */
	if (opts->reorder_function) {
		table=opts->reorder_function(g,FALSE);
	} else if (opts->reorder_map) {
		table=reorder_duplicate(opts->reorder_map,g->n);
	} else {
		table=reorder_ident(g->n);
	}
	ASSERT(reorder_is_bijection(table,g->n));

	if (unweighted_clique_search_single(table,min_size,g,opts) == 0) {
		set_free(current_clique);
		current_clique=NULL;
		goto cleanreturn;
	}
	if (maximal && (min_size>0)) {
		maximalize_clique(current_clique,g);

		if ((max_size > 0) && (set_size(current_clique) > max_size)) {
			clique_options localopts;

			s = set_new(g->n);
			localopts.time_function = opts->time_function;
			localopts.output = opts->output;
			localopts.user_function = false_function;
			localopts.clique_list = &s;
			localopts.clique_list_length = 1;

			for (i=0; i < g->n-1; i++)
				if (clique_size[table[i]]>=min_size)
					break;
			retval = unweighted_clique_search_all(
				table, i, min_size, max_size, maximal, g, &localopts, &found
			);
			set_free(current_clique);
			if (retval || !found) {
				current_clique=NULL;
			} else {
				current_clique=s;
			}
		}
	}

    cleanreturn:
	*clique = current_clique;

	/* Free resources */
	for (i=0; i < temp_count; i++)
		free(temp_list[i]);
	free(temp_list);
	free(table);
	free(clique_size);

	ENTRANCE_RESTORE();
	entrance_level--;

	return retval;
}


/*
 * clique_unweighted_find_all()
 *
 * Find all cliques with size at least min_size and at most max_size.
 *
 *   g        - the graph
 *   min_size - minimum size of cliques to search for.  If min_size==0,
 *              searches for maximum cliques.
 *   max_size - maximum size of cliques to search for.  If max_size==0, no
 *              upper limit is used.  If min_size==0, this must also be 0.
 *   maximal  - require cliques to be maximal cliques
 *   opts     - time printing and clique storage options
 *   num_found - the number of cliques found.  This can be less than the number
 *               of cliques in the graph iff opts->time_function() returns
 *               FALSE (request abort) or opts->user_function() returns an
 *               igraph error code
 *
 * The cliques found are stored in opts->clique_list[] and
 * opts->user_function() is called with them (if non-NULL).  The cliques
 * stored in opts->clique_list[] are newly allocated, and can be freed
 * by set_free().
 */
igraph_error_t clique_unweighted_find_all(
	graph_t *g, int min_size, int max_size, boolean maximal, clique_options *opts,
	CLIQUER_LARGE_INT *num_found
) {
	int i;
	int *table;
	CLIQUER_LARGE_INT count;
	igraph_error_t retval = IGRAPH_SUCCESS;

	ENTRANCE_SAVE();
	entrance_level++;

	if (opts==NULL)
		opts=&clique_default_options;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);
	ASSERT(min_size>=0);
	ASSERT(max_size>=0);
	ASSERT((max_size==0) || (min_size <= max_size));
	ASSERT(!((min_size==0) && (max_size>0)));
	ASSERT((opts->reorder_function==NULL) || (opts->reorder_map==NULL));

	if ((max_size>0) && (min_size>max_size)) {
		/* state was not changed */
		entrance_level--;
		if (num_found) {
			*num_found = 0;
		}
		return IGRAPH_SUCCESS;
	}

    /*
	if (clocks_per_sec==0)
		clocks_per_sec=sysconf(_SC_CLK_TCK);
	ASSERT(clocks_per_sec>0);
    */

	/* Dynamic allocation */
	current_clique=set_new(g->n);
	clique_size=malloc(g->n * sizeof(int));
	/* table allocated later */
	temp_list=malloc((g->n+2)*sizeof(int *));
	temp_count=0;

	clique_list_count=0;
	memset(clique_size,0,g->n * sizeof(int));

	/* "start clock" */
    /*
	gettimeofday(&realtimer,NULL);
	times(&cputimer);
    */

	/* reorder */
	if (opts->reorder_function) {
		table=opts->reorder_function(g,FALSE);
	} else if (opts->reorder_map) {
		table=reorder_duplicate(opts->reorder_map,g->n);
	} else {
		table=reorder_ident(g->n);
	}
	ASSERT(reorder_is_bijection(table,g->n));


	/* Search as normal until there is a chance to find a suitable
	 * clique. */
	if (unweighted_clique_search_single(table,min_size,g,opts) == 0) {
		count=0;
		goto cleanreturn;
	}

	if (min_size==0 && max_size==0) {
		min_size=max_size=clique_size[table[g->n-1]];
		maximal=FALSE;  /* No need to test, since we're searching
				 * for maximum cliques. */
	}
	if (max_size==0) {
		max_size=INT_MAX;
	}

	for (i=0; i < g->n-1; i++)
		if (clique_size[table[i]] >= min_size)
			break;

	retval = unweighted_clique_search_all(table, i, min_size, max_size, maximal, g, opts, &count);

  cleanreturn:
	/* Free resources */
	for (i=0; i<temp_count; i++)
		free(temp_list[i]);
	free(temp_list);
	free(table);
	free(clique_size);
	set_free(current_clique);

	ENTRANCE_RESTORE();
	entrance_level--;

	if (num_found) {
		*num_found = count;
	}

	return retval;
}




/*
 * clique_max_weight()
 *
 * Returns the weight of the maximum weight clique in the graph (or 0 if
 * the search was aborted).
 *
 *   g    - the graph
 *   opts - time printing options
 *
 * Note: As we don't have an algorithm faster than actually finding
 *       a maximum weight clique, we use clique_find_single().
 *       This incurs only very small overhead.
 */
igraph_error_t clique_max_weight(graph_t *g, clique_options *opts, int *weight_found) {
	set_t s;
	int weight = 0;
	igraph_error_t retval;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);

	retval = clique_find_single(g, 0, 0, FALSE, opts, &s);

	if (retval || s == NULL) {
		/* Search was aborted. */
		weight = 0;
	} else {
		weight = graph_subgraph_weight(g,s);
	}

	if (s != NULL) {
		set_free(s);
	}

	if (weight_found) {
		*weight_found = weight;
	}

	return retval;
}


/*
 * clique_find_single()
 *
 * Returns a clique with weight at least min_weight and at most max_weight.
 *
 *   g          - the graph
 *   min_weight - minimum weight of clique to search for.  If min_weight==0,
 *                searches for a maximum weight clique.
 *   max_weight - maximum weight of clique to search for.  If max_weight==0,
 *                no upper limit is used.  If min_weight==0, max_weight must
 *                also be 0.
 *   maximal    - require returned clique to be maximal
 *   opts       - time printing options
 *
 * Returns the set of vertices forming the clique, or NULL if a clique
 * of requested weight/maximality does not exist in the graph  (or if
 * opts->time_function() requests abort).
 *
 * The returned clique is newly allocated and can be freed by set_free().
 *
 * Note: Does NOT use opts->user_function() or opts->clique_list[].
 * Note: Automatically uses clique_unweighted_find_single if all vertex
 *       weights are the same.
 */
igraph_error_t clique_find_single(
	graph_t *g, int min_weight, int max_weight, boolean maximal,
	clique_options *opts, set_t *clique
) {
	int i;
	int *table;
	set_t s;
	igraph_error_t retval = IGRAPH_SUCCESS;
	int weight_found;
	int num_found;

	ENTRANCE_SAVE();
	entrance_level++;

	if (opts==NULL)
		opts=&clique_default_options;

	ASSERT(clique!=NULL);
	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);
	ASSERT(min_weight>=0);
	ASSERT(max_weight>=0);
	ASSERT((max_weight==0) || (min_weight <= max_weight));
	ASSERT(!((min_weight==0) && (max_weight>0)));
	ASSERT((opts->reorder_function==NULL) || (opts->reorder_map==NULL));

	if ((max_weight>0) && (min_weight>max_weight)) {
		/* state was not changed */
		entrance_level--;
		*clique = NULL;
		return IGRAPH_SUCCESS;
	}

    /*
	if (clocks_per_sec==0)
		clocks_per_sec=sysconf(_SC_CLK_TCK);
	ASSERT(clocks_per_sec>0);
    */

	/* Check whether we can use unweighted routines. */
	if (!graph_weighted(g)) {
		min_weight=DIV_UP(min_weight,g->weights[0]);
		if (max_weight) {
			max_weight=DIV_DOWN(max_weight,g->weights[0]);
			if (max_weight < min_weight) {
				/* state was not changed */
				entrance_level--;
				*clique = NULL;
				return IGRAPH_SUCCESS;
			}
		}

		weight_multiplier = g->weights[0];
		entrance_level--;
		retval = clique_unweighted_find_single(g, min_weight, max_weight, maximal, opts, &s);
		ENTRANCE_RESTORE();
		*clique = s;
		return retval;
	}

	/* Dynamic allocation */
	current_clique=set_new(g->n);
	best_clique=set_new(g->n);
	clique_size=malloc(g->n * sizeof(int));
	memset(clique_size, 0, g->n * sizeof(int));
	/* table allocated later */
	temp_list=malloc((g->n+2)*sizeof(int *));
	temp_count=0;

	clique_list_count=0;

	/* "start clock" */
    /*
	gettimeofday(&realtimer,NULL);
	times(&cputimer);
    */

	/* reorder */
	if (opts->reorder_function) {
		table=opts->reorder_function(g,TRUE);
	} else if (opts->reorder_map) {
		table=reorder_duplicate(opts->reorder_map,g->n);
	} else {
		table=reorder_ident(g->n);
	}
	ASSERT(reorder_is_bijection(table,g->n));

	if (max_weight==0)
		max_weight=INT_MAX;

	retval = weighted_clique_search_single(table, min_weight, max_weight, g, opts, &weight_found);

	if (retval || weight_found == 0) {
		/* Requested clique has not been found. */
		set_free(best_clique);
		best_clique=NULL;
		goto cleanreturn;
	}

	if (maximal && (min_weight>0)) {
		maximalize_clique(best_clique,g);
		if (graph_subgraph_weight(g,best_clique) > max_weight) {
			clique_options localopts;

			localopts.time_function = opts->time_function;
			localopts.output = opts->output;
			localopts.user_function = false_function;
			localopts.clique_list = &best_clique;
			localopts.clique_list_length = 1;

			for (i=0; i < g->n-1; i++)
				if ((clique_size[table[i]] >= min_weight) ||
				    (clique_size[table[i]] == 0))
					break;

			retval = weighted_clique_search_all(
				table, i, min_weight, max_weight, maximal, g, &localopts, &num_found
			);

			if (retval || !weight_found) {
				set_free(best_clique);
				best_clique=NULL;
			}
		}
	}

 cleanreturn:
	s=best_clique;

	/* Free resources */
	for (i=0; i < temp_count; i++)
		free(temp_list[i]);
	free(temp_list);
	temp_list=NULL;
	temp_count=0;
	free(table);
	set_free(current_clique);
	current_clique=NULL;
	free(clique_size);
	clique_size=NULL;

	ENTRANCE_RESTORE();
	entrance_level--;

	*clique = s;

	return retval;
}





/*
 * clique_find_all()
 *
 * Find all cliques with weight at least min_weight and at most max_weight.
 *
 *   g          - the graph
 *   min_weight - minimum weight of cliques to search for.  If min_weight==0,
 *                searches for maximum weight cliques.
 *   max_weight - maximum weight of cliques to search for.  If max_weight==0,
 *                no upper limit is used.  If min_weight==0, max_weight must
 *                also be 0.
 *   maximal    - require cliques to be maximal cliques
 *   opts       - time printing and clique storage options
 *   num_found  - the number of cliques found.  This can be less than the number
 *                of cliques in the graph iff opts->time_function() returns
 *                FALSE (request abort) or opts->user_function() returns an
 *                igraph error code
 *
 * The cliques found are stored in opts->clique_list[] and
 * opts->user_function() is called with them (if non-NULL).  The cliques
 * stored in opts->clique_list[] are newly allocated, and can be freed
 * by set_free().
 *
 * Note: Automatically uses clique_unweighted_find_all if all vertex
 *       weights are the same.
 */
igraph_error_t clique_find_all(graph_t *g, int min_weight, int max_weight,
		    boolean maximal, clique_options *opts, int *num_found) {
	int i,n;
	int *table;
	CLIQUER_LARGE_INT r;
	igraph_error_t retval = IGRAPH_SUCCESS;

	ENTRANCE_SAVE();
	entrance_level++;

	if (opts==NULL)
		opts=&clique_default_options;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);
	ASSERT(min_weight>=0);
	ASSERT(max_weight>=0);
	ASSERT((max_weight==0) || (min_weight <= max_weight));
	ASSERT(!((min_weight==0) && (max_weight>0)));
	ASSERT((opts->reorder_function==NULL) || (opts->reorder_map==NULL));

	if ((max_weight>0) && (min_weight>max_weight)) {
		/* state was not changed */
		entrance_level--;
		if (num_found) {
			*num_found = 0;
		}
		return IGRAPH_SUCCESS;
	}

    /*
	if (clocks_per_sec==0)
		clocks_per_sec=sysconf(_SC_CLK_TCK);
	ASSERT(clocks_per_sec>0);
    */

	if (!graph_weighted(g)) {
		min_weight=DIV_UP(min_weight,g->weights[0]);
		if (max_weight) {
			max_weight=DIV_DOWN(max_weight,g->weights[0]);
			if (max_weight < min_weight) {
				/* state was not changed */
				entrance_level--;
				if (num_found) {
					*num_found = 0;
				}
				return IGRAPH_SUCCESS;
			}
		}

		weight_multiplier = g->weights[0];
		entrance_level--;
		retval = clique_unweighted_find_all(g, min_weight, max_weight, maximal, opts, &r);
		ENTRANCE_RESTORE();
		if (num_found) {
			*num_found = r;
		}
		return retval;
	}

	/* Dynamic allocation */
	current_clique=set_new(g->n);
	best_clique=set_new(g->n);
	clique_size=malloc(g->n * sizeof(int));
	memset(clique_size, 0, g->n * sizeof(int));
	/* table allocated later */
	temp_list=malloc((g->n+2)*sizeof(int *));
	temp_count=0;

	/* "start clock" */
    /*
	gettimeofday(&realtimer,NULL);
	times(&cputimer);
    */

	/* reorder */
	if (opts->reorder_function) {
		table=opts->reorder_function(g,TRUE);
	} else if (opts->reorder_map) {
		table=reorder_duplicate(opts->reorder_map,g->n);
	} else {
		table=reorder_ident(g->n);
	}
	ASSERT(reorder_is_bijection(table,g->n));

	/* First phase */
	retval = weighted_clique_search_single(table, min_weight, INT_MAX, g, opts, &n);
	if (retval || n == 0) {
		/* Requested clique has not been found. */
		goto cleanreturn;
	}

	if (min_weight==0) {
		min_weight=n;
		max_weight=n;
		maximal=FALSE;  /* They're maximum cliques already. */
	}
	if (max_weight==0)
		max_weight=INT_MAX;

	for (i=0; i < g->n; i++)
		if ((clique_size[table[i]] >= min_weight) ||
		    (clique_size[table[i]] == 0))
			break;

	/* Second phase */
	retval = weighted_clique_search_all(table, i, min_weight, max_weight, maximal, g, opts, &n);

cleanreturn:
	/* Free resources */
	for (i=0; i < temp_count; i++)
		free(temp_list[i]);
	free(temp_list);
	free(table);
	set_free(current_clique);
	set_free(best_clique);
	free(clique_size);

	ENTRANCE_RESTORE();
	entrance_level--;

	if (num_found) {
		*num_found = n;
	}

	return retval;
}
















#if 0
/*
 * clique_print_time()
 *
 * Reports current running information every 0.1 seconds or when values
 * change.
 *
 *   level    - re-entrance level
 *   i        - current recursion level
 *   n        - maximum recursion level
 *   max      - weight of heaviest clique found
 *   cputime  - CPU time used in algorithm so far
 *   realtime - real time used in algorithm so far
 *   opts     - prints information to (FILE *)opts->output (or stdout if NULL)
 *
 * Returns always TRUE  (ie. never requests abort).
 */
boolean clique_print_time(int level, int i, int n, int max,
			  double cputime, double realtime,
			  clique_options *opts) {
	static float prev_time=100;
	static int prev_i=100;
	static int prev_max=100;
	static int prev_level=0;
	FILE *fp=opts->output;
	int j;

	if (fp==NULL)
		fp=stdout;

	if (ABS(prev_time-realtime)>0.1 || i==n || i<prev_i || max!=prev_max ||
	    level!=prev_level) {
		for (j=1; j<level; j++)
			fprintf(fp,"  ");
		if (realtime-prev_time < 0.01 || i<=prev_i)
			fprintf(fp,"%3d/%d (max %2d)  %2.2f s  "
				"(0.00 s/round)\n",i,n,max,
				realtime);
		else
			fprintf(fp,"%3d/%d (max %2d)  %2.2f s  "
				"(%2.2f s/round)\n",
				i,n,max,realtime,
				(realtime-prev_time)/(i-prev_i));
		prev_time=realtime;
		prev_i=i;
		prev_max=max;
		prev_level=level;
	}
	return TRUE;
}

/*
 * clique_print_time_always()
 *
 * Reports current running information.
 *
 *   level    - re-entrance level
 *   i        - current recursion level
 *   n        - maximum recursion level
 *   max      - largest clique found
 *   cputime  - CPU time used in algorithm so far
 *   realtime - real time used in algorithm so far
 *   opts     - prints information to (FILE *)opts->output (or stdout if NULL)
 *
 * Returns always TRUE  (ie. never requests abort).
 */
boolean clique_print_time_always(int level, int i, int n, int max,
				 double cputime, double realtime,
				 clique_options *opts) {
	static float prev_time=100;
	static int prev_i=100;
	FILE *fp=opts->output;
	int j;

	if (fp==NULL)
		fp=stdout;

	for (j=1; j<level; j++)
		fprintf(fp,"  ");

	if (realtime-prev_time < 0.01 || i<=prev_i)
		fprintf(fp,"%3d/%d (max %2d)  %2.2f s  (0.00 s/round)\n",
			i,n,max,realtime);
	else
		fprintf(fp,"%3d/%d (max %2d)  %2.2f s  (%2.2f s/round)\n",
			i,n,max,realtime,(realtime-prev_time)/(i-prev_i));
	prev_time=realtime;
	prev_i=i;

	return TRUE;
}
#endif
