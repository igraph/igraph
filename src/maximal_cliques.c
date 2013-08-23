/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA
   
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

//#define DEBUG

#include "igraph_cliques.h"
#include "igraph_constants.h"
#include "igraph_interface.h"
#include "igraph_community.h"
#include "igraph_adjlist.h"
#include "igraph_interrupt_internal.h"
#include "igraph_memory.h"

#define PRINT_PX do {						       \
    int j;							       \
    printf("PX=");						       \
    for (j=0; j<PS; j++) {					       \
      printf("%i ", VECTOR(*PX)[j]);				       \
    }								       \
    printf("( ");						       \
    for (; j<=PE; j++) {					       \
      printf("%i ", VECTOR(*PX)[j]);				       \
    }								       \
    printf("| ");						       \
    for (; j<=XE; j++) {					       \
      printf("%i ", VECTOR(*PX)[j]);				       \
    }								       \
    printf(") ");						       \
    for (; j<igraph_vector_int_size(PX); j++) {			       \
      printf("%i ", VECTOR(*PX)[j]);				       \
    }								       \
    printf("\n");						       \
  } while (0);

#define PRINT_PX1 do {						       \
    int j;							       \
    printf("PX=");						       \
    for (j=0; j<PS; j++) {					       \
      printf("%i ", VECTOR(*PX)[j]);				       \
    }								       \
    printf("( ");						       \
    for (; j<=*PE; j++) {					       \
      printf("%i ", VECTOR(*PX)[j]);				       \
    }								       \
    printf("| ");						       \
    for (; j<=XE; j++) {					       \
      printf("%i ", VECTOR(*PX)[j]);				       \
    }								       \
    printf(") ");						       \
    for (; j<igraph_vector_int_size(PX); j++) {			       \
      printf("%i ", VECTOR(*PX)[j]);				       \
    }								       \
    printf("\n");						       \
  } while (0)

int igraph_i_maximal_cliques_reorder_adjlists(
			      const igraph_vector_int_t *PX,
			      int PS, int PE, int XS, int XE,
			      const igraph_vector_int_t *pos,
			      igraph_adjlist_t *adjlist) {
  int j;

  for (j=PS; j<=XE; j++) {
    int av=VECTOR(*PX)[j];
    igraph_vector_t *avneis=igraph_adjlist_get(adjlist, av);
    int pp=0, k, avlen=igraph_vector_size(avneis);

    for (k=0; k<avlen; k++) {
      int avnei=VECTOR(*avneis)[k];
      int avneipos=VECTOR(*pos)[avnei]-1;
      if (avneipos >= PS && avneipos <= PE) { 
	if (pp != k) { 
	  VECTOR(*avneis)[k]=VECTOR(*avneis)[pp];
	  VECTOR(*avneis)[pp]=avnei;
	}
	pp++;
      }
    }
  }
  return 0;
}

int igraph_i_maximal_cliques_reorder_adjlists_down(
				   const igraph_vector_int_t *PX,
				   int PS, int PE, int XS, int XE,
				   const igraph_vector_int_t *pos,
				   igraph_adjlist_t *adjlist,
				   int oldPS, int oldXE) {
  int i;
  for (i=PS; i<=XE; i++) {
    int av=VECTOR(*PX)[i];
    igraph_vector_t *avneis=igraph_adjlist_get(adjlist, av);
    int pp=0, k, avlen=igraph_vector_size(avneis);

    for (k=0; k<avlen; k++) {
      int avnei=VECTOR(*avneis)[k];
      int avneipos=VECTOR(*pos)[avnei]-1;
      if (avneipos < oldPS || avneipos > oldXE) { break; }
      if (avneipos >= PS && avneipos <= PE) {
	if (pp != k) {
	  VECTOR(*avneis)[k]=VECTOR(*avneis)[pp];
	  VECTOR(*avneis)[pp]=avnei;
	}
	pp++;
      }
    }
  }

  return 0;
}

int igraph_i_maximal_cliques_reorder_adjlists_vneis(
				    const igraph_vector_int_t *PX,
				    int PS, int PE, int XS, int XE,
				    const igraph_vector_int_t *pos,
				    igraph_adjlist_t *adjlist, int v) {

  igraph_vector_t *vneis=igraph_adjlist_get(adjlist, v);
  int i, vn=igraph_vector_size(vneis);

#ifdef DEBUG
  printf("Update, v=%i\n", v);
  PRINT_PX;
#endif

  for (i=0; i<vn; i++) {
    int av=VECTOR(*vneis)[i];
    //    int avpos=VECTOR(*pos)[av]-1;
    igraph_vector_t *avneis;
    int j, n, vpos;
    //    if (avpos < PS || avpos > XE) { break; }
    avneis=igraph_adjlist_get(adjlist, av);
    n=igraph_vector_size(avneis);
#ifdef DEBUG
    printf("before (%i): ", av); igraph_vector_print(avneis);
#endif
    /* look for 'v' */
    for (j=0; j<n; j++) {
      int neinei=VECTOR(*avneis)[j];
      if (neinei==v) { vpos=j; break; }
    }
    /* look for end of P vertices */
    for (j++; j<n; j++) {
      int neinei=VECTOR(*avneis)[j];
      int neineipos=VECTOR(*pos)[neinei]-1;
      if (neineipos < PS || neineipos > PE) { break; }
    }
    if (j-1 != vpos) {
      VECTOR(*avneis)[vpos] = VECTOR(*avneis)[j-1];
      VECTOR(*avneis)[j-1] = v;
    }
#ifdef DEBUG
    printf("swap %i-%i: ", vpos, j-1); igraph_vector_print(avneis);
#endif
  }

  return 0;
}

int igraph_i_maximal_cliques_check_order(const igraph_vector_int_t *PX,
					 int PS, int PE, int XS, int XE,
					 const igraph_vector_int_t *pos,
					 const igraph_adjlist_t *adjlist) {
  int i;
  for (i=PS; i<=XE; i++) {
    int v=VECTOR(*PX)[i];
    igraph_vector_t *neis=igraph_adjlist_get(adjlist, v);
    int x=0, j, n=igraph_vector_size(neis);
    for (j=0; j<n; j++) {
      int nei=VECTOR(*neis)[j];
      int neipos=VECTOR(*pos)[nei]-1;
      if (x==0) {
	if (neipos < PS || neipos > PE) { x=1; }
      } else {
	if (neipos >= PS && neipos <= PE) {
	  PRINT_PX;
	  printf("v: %i\n", v); igraph_vector_print(neis);
	  IGRAPH_ERROR("Adjlist ordering error", IGRAPH_EINTERNAL);
	}
      }
    }
  }

  return 0;
}

int igraph_i_maximal_cliques_select_pivot(const igraph_vector_int_t *PX, 
					  int PS, int PE, int XS, int XE,
					  const igraph_vector_int_t *pos,
					  const igraph_adjlist_t *adjlist,
					  int *pivot,
					  igraph_vector_int_t *nextv) {
  igraph_vector_t *pivotvectneis;
  int pivotvectlen, j, usize=-1;

  /* Choose a pivotvect */
  for (j = PS; j <= XE; j++) {
    int ucand=VECTOR(*PX)[j];
    int ucandsize=0;
    igraph_vector_t *ucandneis=igraph_adjlist_get(adjlist, ucand);
    int k, ucanddeg=igraph_vector_size(ucandneis);
    for (k=0; k<ucanddeg; k++) {
      int nei=VECTOR(*ucandneis)[k];
      int neipos=VECTOR(*pos)[nei]-1;
      if (PS <= neipos && neipos <= PE) {
	ucandsize++;		/* in P */
      } else {
	break;
      }
    }
    if (ucandsize > usize) { *pivot=ucand; usize=ucandsize; }
  }

  igraph_vector_int_push_back(nextv, -1);
  pivotvectneis=igraph_adjlist_get(adjlist, *pivot);
  pivotvectlen=igraph_vector_size(pivotvectneis);

  for (j=PS; j <= PE; j++) {
    int vcand=VECTOR(*PX)[j];
    igraph_bool_t nei=0;
    int k=0;
    for (k=0; k < pivotvectlen; k++) {
      int unv=VECTOR(*pivotvectneis)[k];
      int unvpos=VECTOR(*pos)[unv]-1;
      if (unvpos < PS || unvpos > PE) { break; }
      if (unv == vcand) { nei=1; break; }
    }
    if (!nei) { igraph_vector_int_push_back(nextv, vcand); }
  }

  return 0;
}

#define SWAP(p1,p2) do {			\
    int v1=VECTOR(*PX)[p1];			\
    int v2=VECTOR(*PX)[p2];			\
    VECTOR(*PX)[p1] = v2;			\
    VECTOR(*PX)[p2] = v1;			\
    VECTOR(*pos)[v1] = (p2)+1;			\
    VECTOR(*pos)[v2] = (p1)+1;			\
  } while (0)

int igraph_i_maximal_cliques_down(igraph_vector_int_t *PX,
				  int PS, int PE, int XS, int XE,
				  igraph_vector_int_t *pos, 
				  igraph_adjlist_t *adjlist, int mynextv, 
				  igraph_vector_int_t *R, 
				  int *newPS, int *newXE) {

#ifdef DEBUG
  printf("next v: %i\n", mynextv);
#endif  

  igraph_vector_t *vneis=igraph_adjlist_get(adjlist, mynextv);
  int j, vneislen=igraph_vector_size(vneis);

  *newPS=PE+1; *newXE=XS-1;
  for (j=0; j<vneislen; j++) {
    int vnei=VECTOR(*vneis)[j];
    int vneipos=VECTOR(*pos)[vnei]-1;
    if (vneipos < PS || vneipos > PE) { break; }
    (*newPS)--;
    SWAP(vneipos, *newPS);
  }

  /* Continue with the neighbors of mynextv that are in X */
  for (j=0; j<vneislen; j++) {
    int vnei=VECTOR(*vneis)[j];
    int vneipos=VECTOR(*pos)[vnei]-1;
    if (vneipos < XS || vneipos > XE) { continue; }
    (*newXE)++;
    SWAP(vneipos, *newXE); 
  }

  igraph_i_maximal_cliques_reorder_adjlists_down(PX, *newPS, PE, XS,
						 *newXE, pos, adjlist,
						 PS, XE);

  igraph_vector_int_push_back(R, mynextv);
  
  return 0;
}

#undef SWAP

int igraph_i_maximal_cliques_PX(igraph_vector_int_t *PX, int PS, int *PE, 
				int *XS, int XE, igraph_vector_int_t *pos,
				igraph_adjlist_t *adjlist, int v, 
				igraph_vector_int_t *H) {

#ifdef DEBUG
  printf("%i P->X\n", v);
#endif

  int vpos=VECTOR(*pos)[v]-1;
  int tmp=VECTOR(*PX)[*PE];
  VECTOR(*PX)[vpos]=tmp;
  VECTOR(*PX)[*PE]=v;
  VECTOR(*pos)[v]=(*PE)+1;
  VECTOR(*pos)[tmp]=vpos+1;
  (*PE)--; (*XS)--;
  igraph_vector_int_push_back(H, v);

#ifdef DEBUG
  PRINT_PX1;
#endif

  igraph_i_maximal_cliques_reorder_adjlists_vneis(PX, PS, *PE, *XS, XE,
  						  pos, adjlist, v);

  /* igraph_i_maximal_cliques_reorder_adjlists_down(PX, PS, *PE, *XS, XE, */
  /* 						 pos, adjlist, PS, XE); */

  /* igraph_i_maximal_cliques_check_order(PX, PS, *PE, *XS, XE, pos, */
  /* 				       adjlist); */
  return 0;
}

int igraph_i_maximal_cliques_up(igraph_vector_int_t *PX, int PS, int PE, 
				int XS, int XE, igraph_vector_int_t *pos, 
				igraph_adjlist_t *adjlist, 
				igraph_vector_int_t *R, 
				igraph_vector_int_t *H) {
  int vv;
  igraph_vector_int_pop_back(R);

#ifdef DEBUG
  printf("Up, X->P: ");
#endif

  while ((vv=igraph_vector_int_pop_back(H)) != -1) {
    int vvpos=VECTOR(*pos)[vv]-1;
    int tmp=VECTOR(*PX)[XS];
    VECTOR(*PX)[XS]=vv;
    VECTOR(*PX)[vvpos]=tmp;
    VECTOR(*pos)[vv]=XS+1;
    VECTOR(*pos)[tmp]=vvpos+1;
    PE++; XS++;
#ifdef DEBUG
    printf("%i ", vv);
#endif
  }
  
#ifdef DEBUG
  printf("\n");
#endif

  igraph_i_maximal_cliques_reorder_adjlists(PX, PS, PE, XS, XE, pos,
					    adjlist);

  return 0;
}

int igraph_i_maximal_cliques_bk(igraph_vector_int_t *PX, int PS, int PE, 
				int XS, int XE, igraph_vector_int_t *R, 
				igraph_vector_int_t *pos,
				igraph_adjlist_t *adjlist, 
				igraph_vector_ptr_t *res, 
				igraph_vector_int_t *nextv,
				igraph_vector_int_t *H,
				int min_size, int max_size) {

#ifdef DEBUG
  printf("<<<<\n");
  PRINT_PX;
#endif

  igraph_vector_int_push_back(H, -1); /* boundary */
  
  if (PS > PE && XS > XE) {
    /* Found a maximum clique, report it */
    int clsize=igraph_vector_int_size(R);
    if (min_size <= clsize && (clsize <= max_size || max_size <= 0)) {
      igraph_vector_t *cl=igraph_Calloc(1, igraph_vector_t);
      int j;
#ifdef DEBUG
      printf("clique: "); igraph_vector_int_print(R);
#endif
      if (!cl) {
	IGRAPH_ERROR("Cannot list maximal cliques", IGRAPH_ENOMEM);
      }
      igraph_vector_ptr_push_back(res, cl);
      igraph_vector_init(cl, clsize);
      for (j=0; j<clsize; j++) { VECTOR(*cl)[j] = VECTOR(*R)[j]; }
    }
  } else {
    /* Select a pivot element */
    int pivot, mynextv;
    igraph_i_maximal_cliques_select_pivot(PX, PS, PE, XS, XE, pos,
					  adjlist, &pivot, nextv);
#ifdef DEBUG
    printf("pivot: %i\n", pivot);
#endif
    while ((mynextv=igraph_vector_int_pop_back(nextv)) != -1) {
      int newPS, newXE;

      /* Going down, prepare */
      igraph_i_maximal_cliques_down(PX, PS, PE, XS, XE, pos, adjlist,
				    mynextv, R, &newPS, &newXE);
      /* Recursive call */
      igraph_i_maximal_cliques_bk(PX, newPS, PE, XS, newXE, R, pos, 
				  adjlist, res, nextv, H,
				  min_size, max_size);

      /* igraph_i_maximal_cliques_check_order(PX, PS, PE, XS, XE, pos, */
      /* 					   adjlist); */

#ifdef DEBUG
      printf("Restored: ");
      PRINT_PX;
#endif
      /* Putting v from P to X */
      igraph_i_maximal_cliques_PX(PX, PS, &PE, &XS, XE, pos, adjlist, 
				  mynextv, H);
    }
  }
  
  /* Putting back vertices from X to P, see notes in H */
  igraph_i_maximal_cliques_up(PX, PS, PE, XS, XE, pos, adjlist, R, H);

#ifdef DEBUG
  printf(">>>>\n");
#endif

  return 0;
}

/**
 * \function igraph_maximal_cliques
 * \brief Find all maximal cliques of a graph
 *
 * </para><para>
 * A maximal clique is a clique which can't be extended any more by
 * adding a new vertex to it.
 *
 * </para><para>
 * If you are only interested in the size of the largest clique in the
 * graph, use \ref igraph_clique_number() instead.
 *
 * </para><para>
 * The current implementation uses a modified Bron-Kerbosch
 * algorithm to find the maximal cliques, see: David Eppstein,
 * Maarten Löffler, Darren Strash: Listing All Maximal Cliques in
 * Sparse Graphs in Near-Optimal Time. Algorithms and Computation,
 * Lecture Notes in Computer Science Volume 6506, 2010, pp 403-414.
 *
 * </para><para>The implementation of this function changed between
 * igraph 0.5 and 0.6 and also between 0.6 and 0.7, so the order of
 * the cliques and the order of vertices within the cliques will
 * almost surely be different between these three versions.
 *
 * \param graph The input graph.
 * \param res Pointer to a pointer vector, the result will be stored
 *   here, ie. \c res will contain pointers to \c igraph_vector_t
 *   objects which contain the indices of vertices involved in a clique.
 *   The pointer vector will be resized if needed but note that the
 *   objects in the pointer vector will not be freed. Note that vertices
 *   of a clique may be returned in arbitrary order.
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_independent_vertex_sets(), \ref
 * igraph_clique_number() 
 * 
 * Time complexity: O(d(n-d)3^(d/3)) worst case, d is the degeneracy
 * of the graph, this is typically small for sparse graphs.
 * 
 * \example examples/simple/igraph_maximal_cliques.c
 */

int igraph_maximal_cliques(const igraph_t *graph, igraph_vector_ptr_t *res,
			   igraph_integer_t min_size, 
			   igraph_integer_t max_size) {

  /* Implementation details. TODO */

  igraph_vector_int_t PX, R, H, pos, nextv;
  igraph_vector_t coreness, order;
  igraph_vector_int_t rank;	/* TODO: this is not needed */
  int i, no_of_nodes=igraph_vcount(graph);
  igraph_adjlist_t adjlist, fulladjlist;

  if (igraph_is_directed(graph)) {
    IGRAPH_WARNING("Edge directions are ignored for maximal clique "
		   "calculation");
  }

  igraph_vector_init(&order, no_of_nodes);
  IGRAPH_FINALLY(igraph_vector_destroy, &order);
  igraph_vector_int_init(&rank, no_of_nodes);
  IGRAPH_FINALLY(igraph_vector_int_destroy, &rank);
  igraph_vector_init(&coreness, no_of_nodes);
  igraph_coreness(graph, &coreness, /*mode=*/ IGRAPH_ALL);
  igraph_vector_qsort_ind(&coreness, &order, /*descending=*/ 0);
  for (i=0; i<no_of_nodes; i++) {
    int v=VECTOR(order)[i];
    VECTOR(rank)[v] = i;
  }

#ifdef DEBUG
  printf("coreness: "); igraph_vector_print(&coreness);
  printf("order:    "); igraph_vector_print(&order);
  printf("rank:     "); igraph_vector_int_print(&rank);
#endif

  igraph_vector_destroy(&coreness);
  IGRAPH_FINALLY_CLEAN(1);
  
  igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL);
  IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
  igraph_adjlist_simplify(&adjlist);
  igraph_adjlist_init(graph, &fulladjlist, IGRAPH_ALL);
  IGRAPH_FINALLY(igraph_adjlist_destroy, &fulladjlist);
  igraph_adjlist_simplify(&fulladjlist);
  igraph_vector_int_init(&PX, 20);
  IGRAPH_FINALLY(igraph_vector_int_destroy, &PX);
  igraph_vector_int_init(&R,  20);
  IGRAPH_FINALLY(igraph_vector_int_destroy, &R);
  igraph_vector_int_init(&H, 100);
  IGRAPH_FINALLY(igraph_vector_int_destroy, &H);
  igraph_vector_int_init(&pos, no_of_nodes);
  IGRAPH_FINALLY(igraph_vector_int_destroy, &pos);
  igraph_vector_int_init(&nextv, 100);
  IGRAPH_FINALLY(igraph_vector_int_destroy, &nextv);
  
  igraph_vector_ptr_clear(res);

  for (i=0; i<no_of_nodes; i++) {
    int v=VECTOR(order)[i];
    int vrank=VECTOR(rank)[v];
    igraph_vector_t *vneis=igraph_adjlist_get(&fulladjlist, v);
    int vdeg=igraph_vector_size(vneis);
    int Pptr=0, Xptr=vdeg-1, PS=0, PE, XS, XE=vdeg-1;
    int j;

#ifdef DEBUG
    printf("----------- vertex %i\n", v);
#endif
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    igraph_vector_int_resize(&PX, vdeg);
    igraph_vector_int_resize(&R , 1);
    igraph_vector_int_resize(&H , 1);
    igraph_vector_int_null(&pos); /* TODO: makes it quadratic? */
    igraph_vector_int_resize(&nextv, 1);

    VECTOR(H)[0] = -1;		/* marks the end of the recursion */
    VECTOR(nextv)[0] = -1;

    /* ================================================================*/
    /* P <- G(v[i]) intersect { v[i+1], ..., v[n-1] }
       X <- G(v[i]) intersect { v[0], ..., v[i-1] } */
    
    VECTOR(R)[0] = v;
    for (j=0; j<vdeg; j++) {
      int vx=VECTOR(*vneis)[j];
      if (VECTOR(rank)[vx] > vrank) {
	VECTOR(PX)[Pptr] = vx;
	VECTOR(pos)[vx] = Pptr+1;
	Pptr++;
      } else if (VECTOR(rank)[vx] < vrank) {
	VECTOR(PX)[Xptr] = vx;
	VECTOR(pos)[vx] = Xptr+1;
	Xptr--;
      }
    }

    PE = Pptr-1; XS = Xptr+1;	/* end of P, start of X in PX */

    /* Create an adjacency list that is specific to the 
       v vertex. It only contains 'v' and its neighbors. Moreover, we 
       only deal with the vertices in P and X (and R). */
    igraph_vector_update(igraph_adjlist_get(&adjlist, v), 
			 igraph_adjlist_get(&fulladjlist, v));
    for (j=0; j<=vdeg-1; j++) {
      int vv=VECTOR(PX)[j];
      igraph_vector_t *fadj=igraph_adjlist_get(&fulladjlist, vv);
      igraph_vector_t *radj=igraph_adjlist_get(&adjlist, vv);
      int k, fn=igraph_vector_size(fadj);
      igraph_vector_clear(radj);
      for (k=0; k<fn; k++) {
	int nei=VECTOR(*fadj)[k];
	int neipos=VECTOR(pos)[nei]-1;
	if (neipos >= PS && neipos <= XE) {
	  igraph_vector_push_back(radj, nei);
	}
      }
    }

    /* Reorder the adjacency lists, according to P and X. */
    igraph_i_maximal_cliques_reorder_adjlists(&PX, PS, PE, XS, XE, &pos,
					      &adjlist);

    /* igraph_i_maximal_cliques_check_order(&PX, PS, PE, XS, XE, &pos, */
    /* 					 &adjlist); */

    igraph_i_maximal_cliques_bk(&PX, PS, PE, XS, XE, &R, &pos, &adjlist,
				res, &nextv, &H, min_size, max_size);

  }

  return 0;
}
