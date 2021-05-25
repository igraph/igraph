
/*
 * This file contains the vertex reordering routines.
 *
 * Copyright (C) 2002 Sampo Niskanen, Patric Östergård.
 * Licensed under the GNU GPL, read the file LICENSE for details.
 */

#include "reorder.h"

#include <stdlib.h>

#include <limits.h>

#include <igraph_random.h>


/*
 * reorder_set()
 *
 * Reorders the set s with a function  i -> order[i].
 *
 * Note: Assumes that order is the same size as SET_MAX_SIZE(s).
 */
void reorder_set(set_t s,int *order) {
        set_t tmp;
        int i,j;
        setelement e;

        ASSERT(reorder_is_bijection(order,SET_MAX_SIZE(s)));

        tmp=set_new(SET_MAX_SIZE(s));

        for (i=0; i<(SET_MAX_SIZE(s)/ELEMENTSIZE); i++) {
                e=s[i];
                if (e==0)
                        continue;
                for (j=0; j<ELEMENTSIZE; j++) {
                        if (e&1) {
                                SET_ADD_ELEMENT(tmp,order[i*ELEMENTSIZE+j]);
                        }
                        e = e>>1;
                }
        }
        if (SET_MAX_SIZE(s)%ELEMENTSIZE) {
                e=s[i];
                for (j=0; j<(SET_MAX_SIZE(s)%ELEMENTSIZE); j++) {
                        if (e&1) {
                                SET_ADD_ELEMENT(tmp,order[i*ELEMENTSIZE+j]);
                        }
                        e = e>>1;
                }
        }
        set_copy(s,tmp);
        set_free(tmp);
        return;
}


/*
 * reorder_graph()
 *
 * Reorders the vertices in the graph with function  i -> order[i].
 *
 * Note: Assumes that order is of size g->n.
 */
void reorder_graph(graph_t *g, int *order) {
        int i;
        set_t *tmp_e;
        int *tmp_w;

        ASSERT(reorder_is_bijection(order,g->n));

        tmp_e=malloc(g->n * sizeof(set_t));
        tmp_w=malloc(g->n * sizeof(int));
        for (i=0; i<g->n; i++) {
                reorder_set(g->edges[i],order);
                tmp_e[order[i]]=g->edges[i];
                tmp_w[order[i]]=g->weights[i];
        }
        for (i=0; i<g->n; i++) {
                g->edges[i]=tmp_e[i];
                g->weights[i]=tmp_w[i];
        }
        free(tmp_e);
        free(tmp_w);
        return;
}



/*
 * reorder_duplicate()
 *
 * Returns a newly allocated duplicate of the given ordering.
 */
int *reorder_duplicate(int *order,int n) {
	int *new;

	new=malloc(n*sizeof(int));
	memcpy(new,order,n*sizeof(int));
	return new;
}

/*
 * reorder_invert()
 *
 * Inverts the given ordering so that new[old[i]]==i.
 *
 * Note: Asserts that order is a bijection.
 */
void reorder_invert(int *order,int n) {
	int *new;
	int i;

	ASSERT(reorder_is_bijection(order,n));

	new=malloc(n*sizeof(int));
	for (i=0; i<n; i++)
		new[order[i]]=i;
	for (i=0; i<n; i++)
		order[i]=new[i];
	free(new);
	return;
}

/*
 * reorder_reverse()
 *
 * Reverses the given ordering so that  new[i] == n-1 - old[i].
 */
void reorder_reverse(int *order,int n) {
	int i;

	for (i=0; i<n; i++)
		order[i] = n-1 - order[i];
	return;
}

/*
 * reorder_is_bijection
 *
 * Checks that an ordering is a bijection {0,...,n-1} -> {0,...,n-1}.
 *
 * Returns TRUE if it is a bijection, FALSE otherwise.
 */
boolean reorder_is_bijection(int *order,int n) {
	boolean *used;
	int i;

	used=calloc(n,sizeof(boolean));
	for (i=0; i<n; i++) {
		if (order[i]<0 || order[i]>=n) {
			free(used);
			return FALSE;
		}
		if (used[order[i]]) {
			free(used);
			return FALSE;
		}
		used[order[i]]=TRUE;
	}
	for (i=0; i<n; i++) {
		if (!used[i]) {
			free(used);
			return FALSE;
		}
	}
	free(used);
	return TRUE;
}

/*
 * reorder_ident()
 *
 * Returns a newly allocated identity ordering of size n, ie. order[i]==i.
 */
int *reorder_ident(int n) {
	int i;
	int *order;

	order=malloc(n*sizeof(int));
	for (i=0; i<n; i++)
		order[i]=i;
	return order;
}



/*** Reordering functions for use in clique_options ***/

/*
 * reorder_by_ident()
 *
 * Returns an identity ordering.
 */
int *reorder_by_ident(graph_t *g,boolean weighted) {
	return reorder_ident(g->n);
}

/*
 * reorder_by_reverse()
 *
 * Returns a reverse identity ordering.
 */
int *reorder_by_reverse(graph_t *g,boolean weighted) {
	int i;
	int *order;

	order=malloc(g->n * sizeof(int));
	for (i=0; i < g->n; i++)
		order[i]=g->n-i-1;
	return order;
}

/*
 * reorder_by_greedy_coloring()
 *
 * Equivalent to reorder_by_weighted_greedy_coloring or
 * reorder_by_unweighted_greedy_coloring according to the value of weighted.
 */
int *reorder_by_greedy_coloring(graph_t *g,boolean weighted) {
	if (weighted)
		return reorder_by_weighted_greedy_coloring(g,weighted);
	else
		return reorder_by_unweighted_greedy_coloring(g,weighted);
}


/*
 * reorder_by_unweighted_greedy_coloring()
 *
 * Returns an ordering for the graph g by coloring the clique one
 * color at a time, always adding the vertex of largest degree within
 * the uncolored graph, and numbering these vertices 0, 1, ...
 *
 * Experimentally efficient for use with unweighted graphs.
 */
int *reorder_by_unweighted_greedy_coloring(graph_t *g,boolean weighted) {
	int i,j,v;
	boolean *tmp_used;
	int *degree;   /* -1 for used vertices */
	int *order;
	int maxdegree,maxvertex=0;
	boolean samecolor;

	tmp_used=calloc(g->n,sizeof(boolean));
	degree=calloc(g->n,sizeof(int));
	order=calloc(g->n,sizeof(int));

	for (i=0; i < g->n; i++) {
		for (j=0; j < g->n; j++) {
			ASSERT(!((i==j) && GRAPH_IS_EDGE(g,i,j)));
			if (GRAPH_IS_EDGE(g,i,j))
				degree[i]++;
		}
	}

	v=0;
	while (v < g->n) {
		/* Reset tmp_used. */
		memset(tmp_used,0,g->n * sizeof(boolean));

		do {
			/* Find vertex to be colored. */
			maxdegree=0;
			samecolor=FALSE;
			for (i=0; i < g->n; i++) {
				if (!tmp_used[i] && degree[i] >= maxdegree) {
					maxvertex=i;
					maxdegree=degree[i];
					samecolor=TRUE;
				}
			}
			if (samecolor) {
				order[v]=maxvertex;
				degree[maxvertex]=-1;
				v++;

				/* Mark neighbors not to color with same
				 * color and update neighbor degrees. */
				for (i=0; i < g->n; i++) {
					if (GRAPH_IS_EDGE(g,maxvertex,i)) {
						tmp_used[i]=TRUE;
						degree[i]--;
					}
				}
			}
		} while (samecolor);
	}

	free(tmp_used);
	free(degree);
	return order;
}

/*
 * reorder_by_weighted_greedy_coloring()
 *
 * Returns an ordering for the graph g by coloring the clique one
 * color at a time, always adding the vertex that (in order of importance):
 *  1. has the minimum weight in the remaining graph
 *  2. has the largest sum of weights surrounding the vertex
 *
 * Experimentally efficient for use with weighted graphs.
 */
int *reorder_by_weighted_greedy_coloring(graph_t *g, boolean weighted) {
	int i,j,p=0;
	int cnt;
	int *nwt;    /* Sum of surrounding vertices' weights */
	int min_wt,max_nwt;
	boolean *used;
	int *order;

	nwt=malloc(g->n * sizeof(int));
	order=malloc(g->n * sizeof(int));
	used=calloc(g->n,sizeof(boolean));

	for (i=0; i < g->n; i++) {
		nwt[i]=0;
		for (j=0; j < g->n; j++)
			if (GRAPH_IS_EDGE(g, i, j))
				nwt[i] += g->weights[j];
	}

	for (cnt=0; cnt < g->n; cnt++) {
		min_wt=INT_MAX;
		max_nwt=-1;
		for (i=g->n-1; i>=0; i--)
			if ((!used[i]) && (g->weights[i] < min_wt))
				min_wt=g->weights[i];
		for (i=g->n-1; i>=0; i--) {
			if (used[i] || (g->weights[i] > min_wt))
				continue;
			if (nwt[i] > max_nwt) {
				max_nwt=nwt[i];
				p=i;
			}
		}
		order[cnt]=p;
		used[p]=TRUE;
		for (j=0; j < g->n; j++)
			if ((!used[j]) && (GRAPH_IS_EDGE(g, p, j)))
				nwt[j] -= g->weights[p];
	}

	free(nwt);
	free(used);

	ASSERT(reorder_is_bijection(order,g->n));

	return order;
}

/*
 * reorder_by_degree()
 *
 * Returns a reordering of the graph g so that the vertices with largest
 * degrees (most neighbors) are first.
 */
int *reorder_by_degree(graph_t *g, boolean weighted) {
	int i,j,v;
	int *degree;
	int *order;
	int maxdegree,maxvertex=0;

	degree=calloc(g->n,sizeof(int));
	order=calloc(g->n,sizeof(int));

	for (i=0; i < g->n; i++) {
		for (j=0; j < g->n; j++) {
			ASSERT(!((i==j) && GRAPH_IS_EDGE(g,i,j)));
			if (GRAPH_IS_EDGE(g,i,j))
				degree[i]++;
		}
	}

	for (v=0; v < g->n; v++) {
		maxdegree=0;
		for (i=0; i < g->n; i++) {
			if (degree[i] >= maxdegree) {
				maxvertex=i;
				maxdegree=degree[i];
			}
		}
		order[v]=maxvertex;
		degree[maxvertex]=-1;  /* used */
/*** Max. degree withing unselected graph:
		for (i=0; i < g->n; i++) {
			if (GRAPH_IS_EDGE(g,maxvertex,i))
				degree[i]--;
		}
***/
	}

	free(degree);
	return order;
}

/*
 * reorder_by_random()
 *
 * Returns a random reordering for graph g.
 * Note: Used the functions rand() and srand() to generate the random
 *       numbers.  srand() is re-initialized every time reorder_by_random()
 *       is called using the system time.
 */
int *reorder_by_random(graph_t *g, boolean weighted) {
	int i,r;
	int *new;
	boolean *used;

	new=calloc(g->n, sizeof(int));
	used=calloc(g->n, sizeof(boolean));
	for (i=0; i < g->n; i++) {
		do {
            r = igraph_rng_get_integer(igraph_rng_default(), 0, g->n - 1);
		} while (used[r]);
		new[i]=r;
		used[r]=TRUE;
	}
	free(used);
	return new;
}
