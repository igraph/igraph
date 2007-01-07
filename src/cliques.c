/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph.h"
#include "memory.h"
#include "random.h"


void igraph_i_cliques_free_res(igraph_vector_ptr_t *res) {
  long i, n;
  
  n = igraph_vector_ptr_size(res);
  for (i=0; i<n; i++) { 
    igraph_vector_destroy(VECTOR(*res)[i]);
    igraph_free(VECTOR(*res)[i]);
  }
}

/* Internal function for calculating cliques or independent vertex sets.
 * They are practically the same except that the complementer of the graph
 * should be used in the latter case.
 */
int igraph_i_cliques(const igraph_t *graph, igraph_vector_ptr_t *res,
		     igraph_integer_t min_size, igraph_integer_t max_size,
		     igraph_bool_t independent_vertices) {
  igraph_bool_t only_largest = 0, reached_largest = 0;
  igraph_integer_t no_of_nodes;
  igraph_vector_t *new_cliques=0, neis;
  igraph_real_t *member_storage, *new_member_storage, *c1, *c2, v1, v2;
  igraph_bool_t ok;
  long int i, j, k, l, m, n, clique_count, old_clique_count;
  
  no_of_nodes = igraph_vcount(graph);
  
  if (min_size <= 0) min_size = 0;
  if (max_size < 0) max_size = no_of_nodes;
  if (max_size == 0) {
    max_size = no_of_nodes; only_largest = 1;
  }

  igraph_vector_ptr_clear(res);
  
  /* Empty graphs have no cliques */
  if (no_of_nodes == 0) return 0;
  
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  
  /* Find all 1-cliques: every vertex will be a clique */
  member_storage=Calloc(no_of_nodes, igraph_real_t);
  if (member_storage==0) {
    IGRAPH_ERROR("cliques failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, member_storage);
  
  for (i=0; i<no_of_nodes; i++) {
    member_storage[i] = i;
  }
  clique_count = no_of_nodes;
  old_clique_count = 0;
  
  new_member_storage = Calloc(clique_count*(clique_count-1), igraph_real_t);
  if (new_member_storage==0) {
    IGRAPH_ERROR("cliques failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, new_member_storage);

  igraph_vector_ptr_clear(res);
  IGRAPH_FINALLY(igraph_i_cliques_free_res, res);
  
  for (i=2; i<=max_size+1; i++) {
    IGRAPH_ALLOW_INTERRUPTION();

    /* First, let's take i-1-cliques and copy them to the result vector if
     * necessary */    
    if (only_largest && clique_count==0) {
      /* So we had no more cliques. If that is the case, new_cliques (which has
       * been swapped with cliques at the end of the previous iteration)
       * contains the largest cliques. Swap them again so they will be returned. */
      c1 = new_member_storage;
      new_member_storage = member_storage;
      member_storage = c1;
      reached_largest = 1;
      clique_count = old_clique_count;
      i--;
    }

    /* If we reached the largest and collected all that needed */
    if (!only_largest && clique_count==0) {
      break;
    }

    /* If we had cliques AND (we requested cliques of this size OR
     * we requested for the largest and we have reached them) */
    if (clique_count>0 &&
        ((!only_largest && (i >= min_size+1 && i <= max_size+1))
	 || reached_largest)) {
      if (new_cliques != 0) {
	igraph_free(new_cliques);
	new_cliques=0;
	IGRAPH_FINALLY_CLEAN(1);
      }
      new_cliques = Calloc(clique_count, igraph_vector_t);
      if (new_cliques == 0) {
        IGRAPH_ERROR("cliques failed", IGRAPH_ENOMEM);
      }
      IGRAPH_FINALLY(igraph_free, new_cliques);
      
      /*
      printf("%d-clique storage vector:\n", (int)(i-1));
      for (j=0; j<clique_count*(i-1); j++) {
        printf("%d ", (int)member_storage[j]);
        if (j % (i-1) == (i-2)) printf("\n");
      }
      */
      
      for (j=0, k=0; j<clique_count; j++, k+=i-1) {
        IGRAPH_CHECK(igraph_vector_init_copy(&new_cliques[j], &member_storage[k], i-1));
      }

      for (j=0; j<clique_count; j++) {
	igraph_vector_t *tmp=Calloc(1, igraph_vector_t);
	*tmp=new_cliques[j];
        IGRAPH_CHECK(igraph_vector_ptr_push_back(res, tmp));
      }
      
      if (reached_largest) break;
    }
    
    if (i==max_size+1) break;
    
    /* Find all i-cliques, provided that i-1-cliques are in clique_candidates */
    
    /* First, allocate memory for the new clique sets. If there are n
     * i-1-cliques, the amount of memory needed will be i*n*(n-1)/2: there will
     * be at most n*(n-1)/2 cliques with i members each */
    new_member_storage = Realloc(new_member_storage,
                                 i*clique_count*(clique_count-1)/2,
                                 igraph_real_t);
    if (new_member_storage==0) {
      IGRAPH_ERROR("cliques failed", IGRAPH_ENOMEM);
    }
    m=n=0;
    
    /* Now consider all pairs of i-1-cliques and see if they can be merged */
    for (j=0; j<clique_count; j++) {
      for (k=j+1; k<clique_count; k++) {
        IGRAPH_ALLOW_INTERRUPTION();
        
        /* Since cliques are represented by their vertex indices in increasing
         * order, two cliques can be merged iff they have exactly the same
         * indices excluding one AND there is an edge between the two different
         * vertices */
        c1 = member_storage+j*(i-1);
        c2 = member_storage+k*(i-1);
        /* Find the longest prefixes of c1 and c2 that are equal */
        for (l=0; l<i-1 && c1[l] == c2[l]; l++)
          new_member_storage[m++]=c1[l];
        /* Now, if l == i-1, the two vectors are totally equal. This is a bug */
        if (l == i-1) {
          IGRAPH_WARNING("possible bug in igraph_cliques");
          m=n;
        } else {
          /* Assuming that j<k, c1[l] is always less than c2[l], since cliques
           * are ordered alphabetically. Now add c1[l] and store c2[l] in a
           * dummy variable */
          new_member_storage[m++]=c1[l];
          v1=c1[l];
          v2=c2[l];
          l++;
          /* Copy the remaining part of the two vectors. Every member pair
           * found in the remaining parts satisfies the following:
           * 1. If they are equal, they should be added.
           * 2. If they are not equal, the smaller must be equal to the
           *    one stored in the dummy variable. If not, the two vectors
           *    differ in more than one place. The larger will be stored in
           *    the dummy variable again.
           */
          ok=1;
          for (; l<i-1; l++) {
            if (c1[l] == c2[l]) {
              new_member_storage[m++]=c1[l];
              ok=0;
            } else if (ok) {
              if (c1[l] < c2[l]) {
                if (c1[l] == v1) {
                  new_member_storage[m++]=c1[l];
                  v2 = c2[l];
                } else break;
              } else {
                if (ok && c2[l] == v1) {
                  new_member_storage[m++]=c2[l];
                  v2 = c1[l];
                } else break;
              }
            } else break;
          }
          /* Now, if l != i-1, the two vectors had a difference in more than
           * one place, so the whole clique is invalid. */
          if (l != i-1) {
            /* Step back in new_member_storage */
            m=n;
          } else {
            /* v1 and v2 are the two different vertices. Check for an edge
	     * if we are looking for cliques and check for the absence of an
	     * edge if we are looking for independent vertex sets */
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, v1, IGRAPH_ALL));
	    l=igraph_vector_search(&neis, 0, v2, 0);
            if ((l && !independent_vertices) || (!l && independent_vertices)) {
              /* Found a new clique, step forward in new_member_storage */
              if (m==n || v2>new_member_storage[m-1]) {
                new_member_storage[m++]=v2;
                n=m;
              } else {
                m=n;
              }
            } else {
              m=n;
            }
          }
        }
      }
    }

    if (new_cliques != 0) {
      igraph_free(new_cliques);
      IGRAPH_FINALLY_CLEAN(1);
      new_cliques=0;
    }
    
    /* Now all i-cliques are in new_member_storage. Swap it with
     * member_storage */
    c1 = new_member_storage;
    new_member_storage = member_storage;
    member_storage = c1;
    
    /* Calculate how many cliques have we found */
    old_clique_count = clique_count;
    clique_count = n/i;
  }    

  igraph_free(new_member_storage);
  igraph_free(member_storage);
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(4); /* 3 here, +1 is igraph_i_cliques_free_res */
  
  return 0;
}

/**
 * \function igraph_cliques
 * Find all or some cliques in a graph
 *
 * Cliques are fully connected subgraphs of a graph.
 *
 * If you are only interested in the size of the largest clique in the graph,
 * use \ref igraph_clique_number() instead.
 *
 * \param graph The input graph.
 * \param res Pointer to a pointer vector, the result will be stored
 *   here, ie. \c res will contain pointers to \c igraph_vector_t
 *   objects which contain the indices of vertices involved in a clique.
 *   The pointer vector will be resized if needed but note that the
 *   objects in the pointer vector will not be freed.
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \return Error code.
 *
 * \sa \ref igraph_largest_cliques() and \ref igraph_clique_number().
 *
 * Time complexity: TODO
 */
int igraph_cliques(const igraph_t *graph, igraph_vector_ptr_t *res,
                   igraph_integer_t min_size, igraph_integer_t max_size) {
  if (max_size == 0) {
    max_size=-1;
  }
  return igraph_i_cliques(graph, res, min_size, max_size, 0);
}

/**
 * \function igraph_largest_cliques
 * \brief Finds the largest clique(s) in a graph.
 * 
 * A clique is largest (quite intuitively) if there is no other clique
 * in the graph which contains more vertices. 
 * 
 * </para><para>
 * Note that this is not neccessarily the same as a maximal clique,
 * ie. the largest cliques are always maximal but a maximal clique is
 * not always largest.
 * \param graph The input graph.
 * \param res Pointer to an initialized pointer vector, the result
 *        will be stored here. It will be resized as needed.
 * \return Error code.
 * 
 * \sa \ref igraph_cliques(), \ref igraph_maximal_cliques()
 * 
 * Time complexity: TODO.
 */

int igraph_largest_cliques(const igraph_t *graph, igraph_vector_ptr_t *res) {
  return igraph_i_cliques(graph, res, 0, 0, 0);
}

/**
 * \function igraph_independent_vertex_sets
 * Find all independent vertex sets in a graph
 *
 * A vertex set is considered independent if there are no edges between
 * them.
 *
 * </para><para>
 * If you are interested in the size of the largest independent vertex set,
 * use \ref igraph_independence_number() instead.
 *
 * \param graph The input graph.
 * \param res Pointer to a pointer vector, the result will be stored
 *   here, ie. \c res will contain pointers to \c igraph_vector_t
 *   objects which contain the indices of vertices involved in an independent
 *   vertex set. The pointer vector will be resized if needed but note that the
 *   objects in the pointer vector will not be freed.
 * \param min_size Integer giving the minimum size of the sets to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the sets to be
 *   returned. If negative or zero, no upper bound will be used. 
 * \return Error code.
 *
 * \sa \ref igraph_largest_independent_vertex_sets(), 
 * \ref igraph_independence_number().
 *
 * Time complexity: TODO
 */
int igraph_independent_vertex_sets(const igraph_t *graph,
				   igraph_vector_ptr_t *res,
				   igraph_integer_t min_size,
				   igraph_integer_t max_size) {
  if (max_size==0) {
    max_size=-1;
  }
  return igraph_i_cliques(graph, res, min_size, max_size, 1);
}

/**
 * \function igraph_largest_independent_vertex_sets
 * \brief Finds the largest independent vertex set(s) in a graph.
 * 
 * An independent vertex set is largest if there is no other
 * independent vertex set with more vertices in the graph.
 * \param graph The input graph.
 * \param res Pointer to a pointer vector, the result will be stored
 *     here. It will be resized as needed.
 * \return Error code.
 * 
 * \sa \ref igraph_independent_vertex_sets(), \ref
 * igraph_maximal_independent_vertex_sets().
 * 
 * Time complexity: TODO
 */

int igraph_largest_independent_vertex_sets(const igraph_t *graph,
					   igraph_vector_ptr_t *res) {
  return igraph_i_cliques(graph, res, 0, 0, 1);
}

typedef struct igraph_i_max_ind_vsets_data_t {
  igraph_integer_t matrix_size;
  igraph_i_adjlist_t adj_list;         /* Adjacency list of the graph */
  igraph_vector_t deg;                 /* Degrees of individual nodes */
  igraph_set_t* buckets;               /* Bucket array */
  /* The IS value for each node. Still to be explained :) */
  igraph_integer_t* IS;
  igraph_integer_t largest_set_size;   /* Size of the largest set encountered */
} igraph_i_max_ind_vsets_data_t;

int igraph_i_maximal_independent_vertex_sets_backtrack(const igraph_t *graph,
						       igraph_vector_ptr_t *res,
						       igraph_i_max_ind_vsets_data_t *clqdata,
						       igraph_integer_t level) {
  long int v1, v2, v3, c, j, k;
  igraph_vector_t *neis1, *neis2;
  igraph_bool_t f;
  igraph_integer_t j1;
  long int it_state;

  IGRAPH_ALLOW_INTERRUPTION();

  if (level >= clqdata->matrix_size-1) {
    igraph_vector_t *vec;
    igraph_integer_t size=0;
    vec = Calloc(1, igraph_vector_t);
    if (vec == 0)
      IGRAPH_ERROR("igraph_i_maximal_independent_vertex_sets failed", IGRAPH_ENOMEM);
    IGRAPH_VECTOR_INIT_FINALLY(vec, 0);
    if (res) {
      for (v1=0; v1<clqdata->matrix_size; v1++)
	if (clqdata->IS[v1] == 0) {
	  IGRAPH_CHECK(igraph_vector_push_back(vec, v1));
	}
      IGRAPH_CHECK(igraph_vector_ptr_push_back(res, vec));
      size=igraph_vector_size(vec);
    } else {
      for (v1=0, size=0; v1<clqdata->matrix_size; v1++)
	if (clqdata->IS[v1] == 0) size++;
    }
    if (size>clqdata->largest_set_size) clqdata->largest_set_size=size;
    IGRAPH_FINALLY_CLEAN(1);
  } else {
    v1 = level+1;
    /* Count the number of vertices with an index less than v1 that have
     * an IS value of zero */
    neis1 = igraph_i_adjlist_get(&clqdata->adj_list, v1);
    c = 0;
    j = 0;
    while (j<VECTOR(clqdata->deg)[v1] && (v2=VECTOR(*neis1)[j]) <= level) {
      if (clqdata->IS[v2] == 0) c++;
      j++;
    }

    if (c == 0) {
      /* If there are no such nodes... */
      j = 0;
      while (j<VECTOR(clqdata->deg)[v1] && (v2=VECTOR(*neis1)[j]) <= level) {
	clqdata->IS[v2]++;
	j++;
      }
      IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph,res,clqdata,v1));
      j = 0;
      while (j<VECTOR(clqdata->deg)[v1] && (v2=VECTOR(*neis1)[j]) <= level) {
	clqdata->IS[v2]--;
	j++;
      }
    } else {
      /* If there are such nodes, store the count in the IS value of v1 */
      clqdata->IS[v1] = c;
      IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph,res,clqdata,v1));
      clqdata->IS[v1] = 0;
      
      f=1;
      j=0;
      while (j<VECTOR(clqdata->deg)[v1] && (v2=VECTOR(*neis1)[j]) <= level) {
	if (clqdata->IS[v2] == 0) {
	  IGRAPH_CHECK(igraph_set_add(&clqdata->buckets[v1], j));
	  neis2 = igraph_i_adjlist_get(&clqdata->adj_list, v2);
	  k = 0;
	  while (k<VECTOR(clqdata->deg)[v2] && (v3=VECTOR(*neis2)[k])<=level) {
	    clqdata->IS[v3]--;
	    if (clqdata->IS[v3] == 0) f=0;
	    k++;
	  }
	}
	clqdata->IS[v2]++;
	j++;
      }

      if (f) 
	IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph,res,clqdata,v1));

      j=0;
      while (j<VECTOR(clqdata->deg)[v1] && (v2=VECTOR(*neis1)[j]) <= level) {
	clqdata->IS[v2]--;
	j++;
      }
      
      it_state=0;
      while (igraph_set_iterate(&clqdata->buckets[v1], &it_state, &j1)) {
	j=(long)j1;
	v2=VECTOR(*neis1)[j];
	neis2 = igraph_i_adjlist_get(&clqdata->adj_list, v2);
	k = 0;
	while (k<VECTOR(clqdata->deg)[v2] && (v3=VECTOR(*neis2)[k])<=level) {
	  clqdata->IS[v3]++;
	  k++;
	}
      }
      igraph_set_clear(&clqdata->buckets[v1]);
    }
  }

  return 0;
}

void igraph_i_free_set_array(igraph_set_t* array) {
  long int i = 0;
  while (igraph_set_inited(array+i)) {
    igraph_set_destroy(array+i);
    i++;
  }
  Free(array);
}

/**
 * \function igraph_maximal_independent_vertex_sets
 * \brief Find all maximal independent vertex sets of a graph
 *
 * A maximal independent vertex set is an independent vertex set which
 * can't be extended any more by adding a new vertex to it.
 *
 * </para><para>
 * The algorithm used here is based on the following paper:
 * S. Tsukiyama, M. Ide, H. Ariyoshi and I. Shirawaka. A new algorithm for
 * generating all the maximal independent sets. SIAM J Computing,
 * 6:505--517, 1977.
 *
 * </para><para>
 * The implementation was originally written by Kevin O'Neill and modified
 * by K M Briggs in the Very Nauty Graph Library. I simply re-wrote it to
 * use igraph's data structures.
 * 
 * </para><para>
 * If you are interested in the size of the largest independent vertex set,
 * use \ref igraph_independence_number() instead.
 *
 * \param graph The input graph.
 * \param res Pointer to a pointer vector, the result will be stored
 *   here, ie. \c res will contain pointers to \c igraph_vector_t
 *   objects which contain the indices of vertices involved in an independent
 *   vertex set. The pointer vector will be resized if needed but note that the
 *   objects in the pointer vector will not be freed.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_cliques(), \ref
 * igraph_independence_number()
 * 
 * Time complexity: TODO.
 */
int igraph_maximal_independent_vertex_sets(const igraph_t *graph,
					   igraph_vector_ptr_t *res) {
  igraph_i_max_ind_vsets_data_t clqdata;
  long int no_of_nodes = igraph_vcount(graph), i;

  clqdata.matrix_size=no_of_nodes;

  IGRAPH_CHECK(igraph_i_adjlist_init(graph, &clqdata.adj_list, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &clqdata.adj_list);
  igraph_i_adjlist_sort(&clqdata.adj_list);

  clqdata.IS = Calloc(no_of_nodes, igraph_integer_t);
  if (clqdata.IS == 0)
    IGRAPH_ERROR("igraph_maximal_independent_vertex_sets failed", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(igraph_free, clqdata.IS);

  IGRAPH_VECTOR_INIT_FINALLY(&clqdata.deg, no_of_nodes);
  for (i=0; i<no_of_nodes; i++)
    VECTOR(clqdata.deg)[i] = igraph_vector_size(igraph_i_adjlist_get(&clqdata.adj_list, i));

  clqdata.buckets = Calloc(no_of_nodes+1, igraph_set_t);
  if (clqdata.buckets == 0)
    IGRAPH_ERROR("igraph_maximal_independent_vertex_sets failed", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(igraph_i_free_set_array, clqdata.buckets);

  for (i=0; i<no_of_nodes; i++)
    IGRAPH_CHECK(igraph_set_init(&clqdata.buckets[i], 0));

  igraph_vector_ptr_clear(res);
  
  /* Do the show */
  clqdata.largest_set_size=0;
  IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, res, &clqdata, 0));

  /* Cleanup */
  for (i=0; i<no_of_nodes; i++) igraph_set_destroy(&clqdata.buckets[i]);
  igraph_i_adjlist_destroy(&clqdata.adj_list);
  igraph_vector_destroy(&clqdata.deg);
  igraph_free(clqdata.IS);
  igraph_free(clqdata.buckets);
  IGRAPH_FINALLY_CLEAN(4);
  return 0;
}

/**
 * \function igraph_independence_number
 * \brief Find the independence number of the graph
 *
 * The independence number of a graph is the cardinality of the largest
 * independent vertex set.
 *
 * \param graph The input graph.
 * \param no The independence number will be returned to the \c
 *   igraph_integer_t pointed by this variable.
 * \return Error code.
 * 
 * \sa \ref igraph_independent_vertex_sets().
 *
 * Time complexity: TODO.
 */
int igraph_independence_number(const igraph_t *graph, igraph_integer_t *no) {
  igraph_i_max_ind_vsets_data_t clqdata;
  long int no_of_nodes = igraph_vcount(graph), i;

  clqdata.matrix_size=no_of_nodes;

  IGRAPH_CHECK(igraph_i_adjlist_init(graph, &clqdata.adj_list, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &clqdata.adj_list);
  igraph_i_adjlist_sort(&clqdata.adj_list);

  clqdata.IS = Calloc(no_of_nodes, igraph_integer_t);
  if (clqdata.IS == 0)
    IGRAPH_ERROR("igraph_independence_number failed", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(igraph_free, clqdata.IS);

  IGRAPH_VECTOR_INIT_FINALLY(&clqdata.deg, no_of_nodes);
  for (i=0; i<no_of_nodes; i++)
    VECTOR(clqdata.deg)[i] = igraph_vector_size(igraph_i_adjlist_get(&clqdata.adj_list, i));

  clqdata.buckets = Calloc(no_of_nodes+1, igraph_set_t);
  if (clqdata.buckets == 0)
    IGRAPH_ERROR("igraph_independence_number failed", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(igraph_i_free_set_array, clqdata.buckets);

  for (i=0; i<no_of_nodes; i++)
    IGRAPH_CHECK(igraph_set_init(&clqdata.buckets[i], 0));

  /* Do the show */
  clqdata.largest_set_size=0;
  IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, 0, &clqdata, 0));
  *no = clqdata.largest_set_size;

  /* Cleanup */
  for (i=0; i<no_of_nodes; i++) igraph_set_destroy(&clqdata.buckets[i]);
  igraph_i_adjlist_destroy(&clqdata.adj_list);
  igraph_vector_destroy(&clqdata.deg);
  igraph_free(clqdata.IS);
  igraph_free(clqdata.buckets);
  IGRAPH_FINALLY_CLEAN(4);

  return 0;
}

/**
 * \function igraph_maximal_cliques
 * \brief Find all maximal cliques of a graph
 *
 * A maximal clique is a clique which
 * can't be extended any more by adding a new vertex to it. This is actually
 * implemented by looking for a maximal independent vertex set in the
 * complementer of the graph.
 * 
 * If you are only interested in the size of the largest clique in the graph,
 * use \ref igraph_clique_number() instead.
 *
 * \param graph The input graph.
 * \param res Pointer to a pointer vector, the result will be stored
 *   here, ie. \c res will contain pointers to \c igraph_vector_t
 *   objects which contain the indices of vertices involved in a clique.
 *   The pointer vector will be resized if needed but note that the
 *   objects in the pointer vector will not be freed.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_independent_vertex_sets(), \ref
 * igraph_clique_number() 
 * 
 * Time complexity: TODO.
 */
int igraph_maximal_cliques(const igraph_t *graph, igraph_vector_ptr_t *res) {
  igraph_i_max_ind_vsets_data_t clqdata;
  long int no_of_nodes = igraph_vcount(graph), i;

  clqdata.matrix_size=no_of_nodes;

  IGRAPH_CHECK(igraph_i_adjlist_init_complementer(graph, &clqdata.adj_list,
						  IGRAPH_ALL, 0));
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &clqdata.adj_list);
  igraph_i_adjlist_sort(&clqdata.adj_list);

  clqdata.IS = Calloc(no_of_nodes, igraph_integer_t);
  if (clqdata.IS == 0)
    IGRAPH_ERROR("igraph_maximal_cliques failed", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(igraph_free, clqdata.IS);

  IGRAPH_VECTOR_INIT_FINALLY(&clqdata.deg, no_of_nodes);
  for (i=0; i<no_of_nodes; i++)
    VECTOR(clqdata.deg)[i] = igraph_vector_size(igraph_i_adjlist_get(&clqdata.adj_list, i));

  clqdata.buckets = Calloc(no_of_nodes+1, igraph_set_t);
  if (clqdata.buckets == 0)
    IGRAPH_ERROR("igraph_maximal_cliques failed", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(igraph_i_free_set_array, clqdata.buckets);

  for (i=0; i<no_of_nodes; i++)
    IGRAPH_CHECK(igraph_set_init(&clqdata.buckets[i], 0));

  igraph_vector_ptr_clear(res);
  
  /* Do the show */
  clqdata.largest_set_size=0;
  IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, res, &clqdata, 0));

  /* Cleanup */
  for (i=0; i<no_of_nodes; i++) igraph_set_destroy(&clqdata.buckets[i]);
  igraph_i_adjlist_destroy(&clqdata.adj_list);
  igraph_vector_destroy(&clqdata.deg);
  igraph_free(clqdata.IS);
  igraph_free(clqdata.buckets);
  IGRAPH_FINALLY_CLEAN(4);
  return 0;
}

/**
 * \function igraph_clique_number
 * \brief Find the clique number of the graph
 *
 * The clique number of a graph is the size of the largest clique.
 *
 * \param graph The input graph.
 * \param no The clique number will be returned to the \c igraph_integer_t
 *   pointed by this variable.
 * \return Error code.
 * 
 * \sa \ref igraph_cliques(), \ref igraph_largest_cliques().
 * 
 * Time complexity: TODO.
 */
int igraph_clique_number(const igraph_t *graph, igraph_integer_t *no) {
  igraph_i_max_ind_vsets_data_t clqdata;
  long int no_of_nodes = igraph_vcount(graph), i;

  clqdata.matrix_size=no_of_nodes;

  IGRAPH_CHECK(igraph_i_adjlist_init_complementer(graph, &clqdata.adj_list,
						  IGRAPH_ALL, 0));
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &clqdata.adj_list);
  igraph_i_adjlist_sort(&clqdata.adj_list);

  clqdata.IS = Calloc(no_of_nodes, igraph_integer_t);
  if (clqdata.IS == 0)
    IGRAPH_ERROR("igraph_clique_number failed", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(igraph_free, clqdata.IS);

  IGRAPH_VECTOR_INIT_FINALLY(&clqdata.deg, no_of_nodes);
  for (i=0; i<no_of_nodes; i++)
    VECTOR(clqdata.deg)[i] = igraph_vector_size(igraph_i_adjlist_get(&clqdata.adj_list, i));

  clqdata.buckets = Calloc(no_of_nodes+1, igraph_set_t);
  if (clqdata.buckets == 0)
    IGRAPH_ERROR("igraph_clique_number failed", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(igraph_i_free_set_array, clqdata.buckets);

  for (i=0; i<no_of_nodes; i++)
    IGRAPH_CHECK(igraph_set_init(&clqdata.buckets[i], 0));

  /* Do the show */
  clqdata.largest_set_size=0;
  IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, 0, &clqdata, 0));
  *no = clqdata.largest_set_size;

  /* Cleanup */
  for (i=0; i<no_of_nodes; i++) igraph_set_destroy(&clqdata.buckets[i]);
  igraph_i_adjlist_destroy(&clqdata.adj_list);
  igraph_vector_destroy(&clqdata.deg);
  igraph_free(clqdata.IS);
  igraph_free(clqdata.buckets);
  IGRAPH_FINALLY_CLEAN(4);

  return 0;
}

