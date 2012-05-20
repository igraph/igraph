/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_revolver.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_interrupt_internal.h"
#include "igraph_structural.h"
#include "config.h"

#include <math.h>

/* This file contains tools for non-citation evolving networks */
/* Citation networks are in evolver.c */

/***********************************************/
/* degree + degree                             */
/***********************************************/

int igraph_revolver_d_d(const igraph_t *graph,
			igraph_integer_t niter,
			const igraph_vector_t *vtime,
			const igraph_vector_t *etime,
			igraph_matrix_t *kernel,
			igraph_matrix_t *sd,
			igraph_matrix_t *norm,
			igraph_matrix_t *cites,
			igraph_matrix_t *expected,
			igraph_real_t *logprob,
			igraph_real_t *lognull,
			const igraph_matrix_t *debug,
			igraph_vector_ptr_t *debugres) {
  
  igraph_integer_t no_of_events, vnoev, enoev;
  igraph_vector_t st;
  long int i;
  igraph_integer_t maxdegree;
  igraph_vector_t vtimeidx, etimeidx;
  igraph_lazy_inclist_t inclist;

  if (igraph_vector_size(vtime) != igraph_vcount(graph)) {
    IGRAPH_ERROR("Invalid vtime length", IGRAPH_EINVAL);
  }
  if (igraph_vector_size(etime) != igraph_ecount(graph)) {
    IGRAPH_ERROR("Invalid etime length", IGRAPH_EINVAL);
  }
  
  vnoev=igraph_vector_max(vtime)+1;
  enoev=igraph_vector_max(etime)+1;
  no_of_events= vnoev > enoev ? vnoev : enoev;

  IGRAPH_VECTOR_INIT_FINALLY(&st, no_of_events);
  for (i=0; i<no_of_events; i++) {
    VECTOR(st)[i]=1;
  }
  
  IGRAPH_CHECK(igraph_maxdegree(graph, &maxdegree, igraph_vss_all(),
				IGRAPH_ALL, IGRAPH_LOOPS));
  
  IGRAPH_VECTOR_INIT_FINALLY(&vtimeidx, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&etimeidx, 0);
  IGRAPH_CHECK(igraph_vector_order1(vtime, &vtimeidx, no_of_events));
  IGRAPH_CHECK(igraph_vector_order1(etime, &etimeidx, no_of_events));
  
  IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

  IGRAPH_PROGRESS("Revolver d-d", 0, NULL);
  for (i=0; i<niter; i++) {
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    if (i+1 != niter) {		/* not the last iteration */
      /* measure */
      IGRAPH_CHECK(igraph_revolver_mes_d_d(graph, &inclist, 
					   kernel, 0 /*sd*/, 0 /*norm*/,
					   0/*cites*/, 0/*debug*/, 0 /*debugres*/,
					   &st, vtime, &vtimeidx, etime, 
					   &etimeidx, no_of_events,
					   maxdegree));
      /* normalize */
      igraph_matrix_scale(kernel, 1/igraph_matrix_sum(kernel));
      
      /* update st */
      IGRAPH_CHECK(igraph_revolver_st_d_d(graph, &inclist, 
					  &st, kernel, vtime, &vtimeidx,
					  etime, &etimeidx,
					  no_of_events));
      
    } else {
      /* measure */
      IGRAPH_CHECK(igraph_revolver_mes_d_d(graph, &inclist,
					   kernel, sd, norm, cites, 
					   debug, debugres, &st, vtime, &vtimeidx,
					   etime, &etimeidx,
					   no_of_events, maxdegree));
      
      /* normalize */
      igraph_matrix_scale(kernel, 1/igraph_matrix_sum(kernel));
      
      /* update st */
      IGRAPH_CHECK(igraph_revolver_st_d_d(graph, &inclist,
					  &st, kernel, vtime, &vtimeidx,
					  etime, &etimeidx,
					  no_of_events));
      
      /* expected number of citations */
      if (expected) {
	IGRAPH_CHECK(igraph_revolver_exp_d_d(graph, &inclist,
					     expected, kernel, &st,
					     vtime, &vtimeidx, etime, &etimeidx,
					     no_of_events, 
					     maxdegree));
      }
      
      /* error calculation */
      if (logprob || lognull) {
	IGRAPH_CHECK(igraph_revolver_error_d_d(graph, &inclist,
					       kernel, &st,
					       vtime, &vtimeidx, 
					       etime, &etimeidx, no_of_events,
					       maxdegree, logprob, lognull));
      }
    }

    IGRAPH_PROGRESS("Revolver d-d", 100.0*(i+1)/niter, NULL);
  }

  igraph_lazy_inclist_destroy(&inclist);
  igraph_vector_destroy(&etimeidx);
  igraph_vector_destroy(&vtimeidx);
  igraph_vector_destroy(&st);
  IGRAPH_FINALLY_CLEAN(4);
  
  return 0;
}

#define NTKK(xidx, yidx) \
   ((xidx)==(yidx) ? (VECTOR(ntk)[(xidx)]*(VECTOR(ntk)[(xidx)]-1))/2-MATRIX(ntkk,(xidx),(yidx)) : VECTOR(ntk)[(xidx)]*VECTOR(ntk)[(yidx)]-MATRIX(ntkk,(xidx),(yidx)))

/* int print_ntkk(igraph_matrix_t *ntkk, igraph_vector_long_t *ntk) { */
/*   long int i, j, r=igraph_matrix_nrow(ntkk), c=igraph_matrix_ncol(ntkk); */
/*   for (i=0; i<r; i++) { */
/*     for (j=0; j<c; j++) { */
/*       long int val=(i==j) ? */
/* 	(VECTOR(*ntk)[i]*(VECTOR(*ntk)[i]-1)/2-MATRIX(*ntkk,i,j)) : */
/* 	(VECTOR(*ntk)[i]*VECTOR(*ntk)[j]-MATRIX(*ntkk,i,j)); */
/*       fprintf(stderr, "%li ", val); */
/*     } */
/*     fprintf(stderr, "\n"); */
/*   } */
/*   fprintf(stderr, "*************\n"); */
/*   return 0; */
/* } */

int igraph_revolver_mes_d_d(const igraph_t *graph, 
			    igraph_lazy_inclist_t *inclist,
			    igraph_matrix_t *kernel,
			    igraph_matrix_t *sd,
			    igraph_matrix_t *norm,
			    igraph_matrix_t *cites,
			    const igraph_matrix_t *debug,
			    igraph_vector_ptr_t *debugres,
			    const igraph_vector_t *st,
			    const igraph_vector_t *vtime,
			    const igraph_vector_t *vtimeidx,
			    const igraph_vector_t *etime,
			    const igraph_vector_t *etimeidx,
			    igraph_integer_t pno_of_events,
			    igraph_integer_t pmaxdegree) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  long int no_of_events=pno_of_events;
  long int maxdegree=pmaxdegree;
  
  igraph_vector_long_t degree;
  igraph_vector_char_t added;	/* is this edge already in the network? */
    
  igraph_matrix_t v_normfact, *normfact, v_notnull, *notnull;
  igraph_matrix_t ch;
  
  igraph_vector_long_t ntk;	/* # of type x vertices */
  igraph_matrix_t ntkk;	        /* # of connections between type x1, x2 vert. */
  
  igraph_vector_t *adjedges;

  long int timestep, i;
  long int nptr=0, eptr=0;
  long int nptr_save, eptr_save, eptr_new;
  
  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(&igraph_vector_long_destroy, &degree);

  IGRAPH_CHECK(igraph_vector_char_init(&added, no_of_edges));
  IGRAPH_FINALLY(igraph_vector_char_destroy, &added);

  IGRAPH_CHECK(igraph_vector_long_init(&ntk, maxdegree+1));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &ntk);
  IGRAPH_MATRIX_INIT_FINALLY(&ntkk, maxdegree+1, maxdegree+1);
  IGRAPH_MATRIX_INIT_FINALLY(&ch, maxdegree+1, maxdegree+1);
  
  if (norm) {
    normfact=norm;
    IGRAPH_CHECK(igraph_matrix_resize(normfact, maxdegree+1, maxdegree+1));
    igraph_matrix_null(normfact);
  } else {
    normfact=&v_normfact;
    IGRAPH_MATRIX_INIT_FINALLY(normfact, maxdegree+1, maxdegree+1);
  }
  
  if (cites) {
    notnull=cites;
    IGRAPH_CHECK(igraph_matrix_resize(notnull, maxdegree+1, maxdegree+1));
    igraph_matrix_null(notnull);
  } else {
    notnull=&v_notnull;
    IGRAPH_MATRIX_INIT_FINALLY(notnull, maxdegree+1, maxdegree+1);
  }
  
  IGRAPH_CHECK(igraph_matrix_resize(kernel, maxdegree+1, maxdegree+1));
  igraph_matrix_null(kernel);
  if (sd) {
    IGRAPH_CHECK(igraph_matrix_resize(sd, maxdegree+1, maxdegree+1));
    igraph_matrix_null(sd);
  }
  
  for (timestep=0; timestep<no_of_events; timestep++) {
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    /* Add the vertices in the first */
    nptr_save=nptr;
    while (nptr < no_of_nodes && 
	   VECTOR(*vtime)[(long int)VECTOR(*vtimeidx)[nptr]]==timestep) {
      nptr++;
    }
    VECTOR(ntk)[0] += (nptr-nptr_save);

    /* Update ch accordingly, could be done later as well */
    if (VECTOR(ntk)[0] == nptr-nptr_save && nptr!=nptr_save) {
      if (nptr-nptr_save >= 2) {
	MATRIX(ch, 0, 0) = eptr;
      }
      for (i=1; i<maxdegree+1; i++) {
	if (NTKK(0,i) == (nptr-nptr_save)*VECTOR(ntk)[i]) {
	  MATRIX(ch, 0, i) = MATRIX(ch, i, 0) = eptr;
	}
      }
    }    

    /* Estimate Akk */
    eptr_save=eptr;
    while (eptr < no_of_edges &&
	   VECTOR(*etime)[(long int)VECTOR(*etimeidx)[eptr] ] == timestep) {
      long int edge=VECTOR(*etimeidx)[eptr];
      long int from=IGRAPH_FROM(graph, edge), to=IGRAPH_TO(graph, edge);
      long int xidx=VECTOR(degree)[from];
      long int yidx=VECTOR(degree)[to];
      double xk, oldakk;
      
      MATRIX(*notnull, xidx, yidx) += 1;
      MATRIX(*notnull, yidx, xidx) = MATRIX(*notnull, xidx, yidx);
      
      xk=VECTOR(*st)[timestep]/NTKK(xidx, yidx);
      oldakk=MATRIX(*kernel, xidx, yidx);
      MATRIX(*kernel, xidx, yidx) +=  (xk-oldakk)/MATRIX(*notnull, xidx, yidx);
      MATRIX(*kernel, yidx, xidx) = MATRIX(*kernel, xidx, yidx);
      if (sd) {
	MATRIX(*sd, xidx, yidx) += (xk-oldakk)*(xk-MATRIX(*kernel, xidx, yidx));
	MATRIX(*sd, yidx, xidx) = MATRIX(*sd, xidx, yidx);
      }
      /* TODO: debug */

      eptr++;
    }
    
    /* Update ntkk, ntk, ch, normfact, add the edges */
    eptr_new=eptr;
    eptr=eptr_save;
    while (eptr < no_of_edges && 
	   VECTOR(*etime)[(long int) VECTOR(*etimeidx)[eptr] ] == timestep) {
      long int edge=VECTOR(*etimeidx)[eptr];
      long int from=IGRAPH_FROM(graph, edge);
      long int to=IGRAPH_TO(graph, edge);
      long int xidx=VECTOR(degree)[from];
      long int yidx=VECTOR(degree)[to];
      long int n;

      adjedges=igraph_lazy_inclist_get(inclist, from);
      n=igraph_vector_size(adjedges);
      for (i=0; i<n; i++) {
	long int edge=VECTOR(*adjedges)[i];
	if (VECTOR(added)[edge]) {
	  long int otherv=IGRAPH_OTHER(graph, edge, from); /* other than from */
	  long int deg=VECTOR(degree)[otherv];
	  MATRIX(ntkk, xidx, deg) -= 1;
	  MATRIX(ntkk, deg, xidx) = MATRIX(ntkk, xidx, deg);
	  if (NTKK(xidx, deg)==1) {
	    MATRIX(ch, deg, xidx) = eptr_new;
	    MATRIX(ch, xidx, deg) = MATRIX(ch, deg, xidx);
	  }
	  MATRIX(ntkk, xidx+1, deg) += 1;
	  MATRIX(ntkk, deg, xidx+1) = MATRIX(ntkk, xidx+1, deg);
	  if (NTKK(xidx+1, deg)==0) {
	    MATRIX(*normfact, xidx+1, deg) += eptr_new-MATRIX(ch, xidx+1, deg);
	    MATRIX(*normfact, deg, xidx+1) = MATRIX(*normfact, xidx+1, deg);
	  }
	}
      }
      adjedges=igraph_lazy_inclist_get(inclist, to);
      n=igraph_vector_size(adjedges);
      for (i=0; i<n; i++) {
	long int edge=VECTOR(*adjedges)[i];
	if (VECTOR(added)[edge]) {
	  long int otherv=IGRAPH_OTHER(graph, edge, to); /* other than to */
	  long int deg=VECTOR(degree)[otherv];
	  MATRIX(ntkk, yidx, deg) -= 1;	  
	  MATRIX(ntkk, deg, yidx) = MATRIX(ntkk, yidx, deg);
	  if (NTKK(yidx, deg)==1) {
	    MATRIX(ch, deg, yidx) = eptr_new;
	    MATRIX(ch, yidx, deg) = MATRIX(ch, deg, yidx);
	  }
	  MATRIX(ntkk, yidx+1, deg) += 1;
	  MATRIX(ntkk, deg, yidx+1) = MATRIX(ntkk, yidx+1, deg);
	  if (NTKK(yidx+1, deg)==0) {
	    MATRIX(*normfact, yidx+1, deg) += eptr_new-MATRIX(ch, yidx+1, deg);
	    MATRIX(*normfact, deg, yidx+1) = MATRIX(*normfact, yidx+1, deg);
	  }
	}
      }

      VECTOR(added)[edge]=1;

      MATRIX(ntkk, xidx+1, yidx+1) += 1;
      MATRIX(ntkk, yidx+1, xidx+1) = MATRIX(ntkk, xidx+1, yidx+1);      
      if (NTKK(xidx+1, yidx+1)==0) {
	MATRIX(*normfact, xidx+1, yidx+1) = eptr_new-MATRIX(ch, xidx+1, yidx+1);
	MATRIX(*normfact, yidx+1, xidx+1) = MATRIX(*normfact, xidx+1, yidx+1);
      }

      for (i=0; i<maxdegree+1; i++) {
	long int before, after;
	before=NTKK(xidx,i); 
	VECTOR(ntk)[xidx] -= 1;
	after=NTKK(xidx,i);
	VECTOR(ntk)[xidx] += 1;
	if (before > 0 && after==0) {
	  MATRIX(*normfact, xidx, i) += eptr_new-MATRIX(ch, xidx, i);
	  MATRIX(*normfact, i, xidx) = MATRIX(*normfact, xidx, i);
	}
      }
      VECTOR(ntk)[xidx]--;

      for (i=0; i<maxdegree+1; i++) {
	long int before, after;
	before=NTKK(yidx, i); 
	VECTOR(ntk)[yidx] -= 1;
	after=NTKK(yidx, i);
	VECTOR(ntk)[yidx] += 1;
	if (before > 0 && after==0) {
	  MATRIX(*normfact, yidx, i) += eptr_new-MATRIX(ch, yidx, i);
	  MATRIX(*normfact, i, yidx) = MATRIX(*normfact, yidx, i);
	}
      }
      VECTOR(ntk)[yidx]--;

      for (i=0; i<maxdegree+1; i++) {
	long int before, after;
	before=NTKK(xidx+1, i);
	VECTOR(ntk)[xidx+1] += 1;
	after=NTKK(xidx+1, i);
	VECTOR(ntk)[xidx+1] -= 1;
	if (before==0 && after > 0) {
	  MATRIX(ch, xidx+1, i) = eptr_new;
	  MATRIX(ch, i, xidx+1) = MATRIX(ch, xidx+1, i);
	}
      }
      VECTOR(ntk)[xidx+1]++;

      for (i=0; i<maxdegree+1; i++) {
	long int before, after;
	before=NTKK(yidx+1, i);
	VECTOR(ntk)[yidx+1] += 1;
	after=NTKK(yidx+1, i);
	VECTOR(ntk)[yidx+1] -= 1;
	if (before == 0 && after == 0) {
	  MATRIX(ch, yidx+1, i) = eptr_new;
	  MATRIX(ch, i, yidx+1) = MATRIX(ch, yidx+1, i);
	}
      }
      VECTOR(ntk)[yidx+1]++;
            
      VECTOR(degree)[from]++;
      VECTOR(degree)[to]++;

      eptr++;
    }
    
  }

  for (i=0; i<maxdegree+1; i++) {
    igraph_real_t oldakk;
    long int j;
    for (j=0; j<=i; j++) {
      if (NTKK(i, j) != 0) {
	MATRIX(*normfact, i, j) += (eptr-MATRIX(ch, i, j));
	MATRIX(*normfact, j, i) = MATRIX(*normfact, i, j);
      }
      if (MATRIX(*normfact, i, j)==0) {
	MATRIX(*kernel, i, j)=MATRIX(*kernel, j, i)=0;
	MATRIX(*normfact, i, j)=MATRIX(*normfact, j, i)=1;
      }
      oldakk=MATRIX(*kernel, i, j);
      MATRIX(*kernel, i, j) *= MATRIX(*notnull, i, j)/MATRIX(*normfact, i, j);
      MATRIX(*kernel, j, i) = MATRIX(*kernel, i, j);
      if (sd) {
	MATRIX(*sd, i, j) += oldakk * oldakk * MATRIX(*notnull, i, j) *
	  (1-MATRIX(*notnull, i, j)/MATRIX(*normfact, i, j));
	MATRIX(*sd, i, j) = sqrt(MATRIX(*sd, i, j)/(MATRIX(*normfact, i, j)-1));
	MATRIX(*sd, j, i) = MATRIX(*sd, i, j);
      }
    }
  }
  
  if (!cites) {
    igraph_matrix_destroy(notnull);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (!norm) {
    igraph_matrix_destroy(normfact);
    IGRAPH_FINALLY_CLEAN(1);
  }

  igraph_matrix_destroy(&ch);
  igraph_matrix_destroy(&ntkk);
  igraph_vector_long_destroy(&ntk);
  igraph_vector_char_destroy(&added);
  igraph_vector_long_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(5);
  
  return 0;
}

#undef NTKK

int igraph_revolver_st_d_d(const igraph_t *graph,
			   igraph_lazy_inclist_t *inclist,
			   igraph_vector_t *st,
			   const igraph_matrix_t *kernel,
			   const igraph_vector_t *vtime,
			   const igraph_vector_t *vtimeidx,
			   const igraph_vector_t *etime,
			   const igraph_vector_t *etimeidx,
			   igraph_integer_t pno_of_events) {

  long int no_of_events=pno_of_events;
  long int maxdegree=igraph_matrix_nrow(kernel)-1;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  long int timestep=0;
  
  igraph_vector_long_t degree;
  igraph_vector_long_t ntk;
  igraph_vector_char_t added;
  
  igraph_vector_t *adjedges;

  long int i;
  long int nptr=0, eptr=0;

  IGRAPH_CHECK(igraph_vector_long_init(&ntk, maxdegree+1));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &ntk);
  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);
  IGRAPH_CHECK(igraph_vector_char_init(&added, no_of_edges));
  IGRAPH_FINALLY(igraph_vector_char_destroy, &added);
  
  IGRAPH_CHECK(igraph_vector_resize(st, no_of_events));
  VECTOR(*st)[0]=0;
  
  for (timestep=0; timestep<no_of_events-1; timestep++) {
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    /* add the new nodes */
    while (nptr < no_of_nodes && 
	   VECTOR(*vtime)[ (long int) VECTOR(*vtimeidx)[nptr] ] == timestep) {
      for (i=0; i<maxdegree+1; i++) {
	VECTOR(*st)[timestep] += VECTOR(ntk)[i]*MATRIX(*kernel, i, 0);
      }
      VECTOR(ntk)[0]++;
      nptr++;
    }

    /* add the new edges as well, but this is for the next timestep
       already */
    VECTOR(*st)[timestep+1] = VECTOR(*st)[timestep];
    while (eptr < no_of_edges && 
	   VECTOR(*etime)[ (long int) VECTOR(*etimeidx)[eptr] ] == timestep) {
      long int edge=VECTOR(*etimeidx)[eptr];
      long int from=IGRAPH_FROM(graph, edge);
      long int to=IGRAPH_TO(graph, edge);
      long int xidx=VECTOR(degree)[from];
      long int yidx=VECTOR(degree)[to];
      igraph_real_t inc=0;
      long int n;
      
      inc -= MATRIX(*kernel, xidx, yidx);

      for (i=0; i<maxdegree+1; i++) {
	inc += VECTOR(ntk)[i] * (MATRIX(*kernel, i, xidx+1) -
				 MATRIX(*kernel, i, xidx)   +
				 MATRIX(*kernel, i, yidx+1) -
				 MATRIX(*kernel, i, yidx));
      }
      inc -= MATRIX(*kernel, xidx+1, xidx+1);
      inc -= MATRIX(*kernel, yidx+1, yidx+1);
      inc += MATRIX(*kernel, xidx, xidx);
      inc += MATRIX(*kernel, yidx, yidx);
            
      VECTOR(ntk)[xidx]--;
      VECTOR(ntk)[yidx]--;
      VECTOR(ntk)[xidx+1]++;
      VECTOR(ntk)[yidx+1]++;
      
      adjedges=igraph_lazy_inclist_get(inclist, from);
      n=igraph_vector_size(adjedges);
      for (i=0; i<n; i++) {
	long int edge=VECTOR(*adjedges)[i];
	if (VECTOR(added)[edge]) {
	  long int otherv=IGRAPH_OTHER(graph, edge, from);
	  long int deg=VECTOR(degree)[otherv];
	  inc += MATRIX(*kernel, xidx, deg);
	  inc -= MATRIX(*kernel, xidx+1, deg);
	}
      }
      adjedges=igraph_lazy_inclist_get(inclist, to);
      n=igraph_vector_size(adjedges);
      for (i=0; i<n; i++) {
	long int edge=VECTOR(*adjedges)[i];
	if (VECTOR(added)[edge]) {
	  long int otherv=IGRAPH_OTHER(graph, edge, to);
	  long int deg=VECTOR(degree)[otherv];
	  inc += MATRIX(*kernel, yidx, deg);
	  inc -= MATRIX(*kernel, yidx+1, deg);
	}
      }
      
      VECTOR(degree)[from] += 1;
      VECTOR(degree)[to] += 1;
      VECTOR(added)[edge]=1;
      
      VECTOR(*st)[timestep+1] += inc;
      
      eptr++;
    }        
  }

  igraph_vector_char_destroy(&added);
  igraph_vector_long_destroy(&degree);
  igraph_vector_long_destroy(&ntk);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}

int igraph_revolver_exp_d_d(const igraph_t *graph,
			    igraph_lazy_inclist_t *inclist,
			    igraph_matrix_t *expected,
			    const igraph_matrix_t *kernel,
			    const igraph_vector_t *st,
			    const igraph_vector_t *vtime,
			    const igraph_vector_t *vtimeidx,
			    const igraph_vector_t *etime,
			    const igraph_vector_t *etimeidx,
			    igraph_integer_t pno_of_events,
			    igraph_integer_t pmaxdegree) {
  /* TODO */
  return 0;
}

int igraph_revolver_error_d_d(const igraph_t *graph,
			      igraph_lazy_inclist_t *inclist,
			      const igraph_matrix_t *kernel,
			      const igraph_vector_t *st,
			      const igraph_vector_t *vtime,
			      const igraph_vector_t *vtimeidx,
			      const igraph_vector_t *etime,
			      const igraph_vector_t *etimeidx,
			      igraph_integer_t pno_of_events,
			      igraph_integer_t pmaxdegree, 
			      igraph_real_t *logprob,
			      igraph_real_t *lognull) {

  long int no_of_events=pno_of_events;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  
  igraph_vector_long_t degree;
    
  long int timestep, nptr=0, eptr=0, eptr_save;
  long int edges=0, vertices=0;
  
  igraph_real_t rlogprob, rlognull, *mylogprob=logprob, *mylognull=lognull;

  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);

  if (!logprob) { mylogprob=&rlogprob; }
  if (!lognull) { mylognull=&rlognull; }
  
  *mylogprob=0;
  *mylognull=0;
  
  for (timestep=0; timestep<no_of_events; timestep++) {
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    while (nptr < no_of_nodes && 
	   VECTOR(*vtime)[ (long int) VECTOR(*vtimeidx)[nptr] ] == timestep) {
      vertices++;
      nptr++;
    }

    eptr_save=eptr;
    while (eptr < no_of_edges && 
	   VECTOR(*etime)[ (long int) VECTOR(*etimeidx)[eptr] ] == timestep) {
      long int edge=VECTOR(*etimeidx)[eptr];
      long int from=IGRAPH_FROM(graph, edge);
      long int to=IGRAPH_TO(graph, edge);
      long int xidx=VECTOR(degree)[from];
      long int yidx=VECTOR(degree)[to];
      
      igraph_real_t prob=MATRIX(*kernel, xidx, yidx)/VECTOR(*st)[timestep];
      igraph_real_t nullprob=1.0/(vertices*(vertices-1)/2-eptr_save);
            
      *mylogprob += log(prob);
      *mylognull += log(nullprob);
      
      edges++;
      eptr++;
    }

    eptr=eptr_save;
    while (eptr < no_of_edges && 
	   VECTOR(*etime)[ (long int) VECTOR(*etimeidx)[eptr] ] == timestep) {
      long int edge=VECTOR(*etimeidx)[eptr];
      long int from=IGRAPH_FROM(graph, edge);
      long int to=IGRAPH_TO(graph, edge);
      VECTOR(degree)[from] += 1;
      VECTOR(degree)[to] += 1;
      eptr++;
    }
  }

  igraph_vector_long_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(1);
    
  return 0;
}

/***********************************************/
/* # of papers + # of papers                   */
/***********************************************/

int igraph_revolver_p_p(const igraph_t *graph,
			igraph_integer_t niter,
			const igraph_vector_t *vtime,
			const igraph_vector_t *etime,
			const igraph_vector_t *authors,
			const igraph_vector_t *eventsizes,
			igraph_matrix_t *kernel,
			igraph_matrix_t *sd,
			igraph_matrix_t *norm,
			igraph_matrix_t *cites,
			igraph_matrix_t *expected,
			igraph_real_t *logprob,
			igraph_real_t *lognull,
			const igraph_matrix_t *debug,
			igraph_vector_ptr_t *debugres) {
  
  igraph_integer_t no_of_events;
  igraph_vector_t st;
  long int i;
  igraph_integer_t maxpapers=0;
  igraph_vector_t vtimeidx, etimeidx;
  igraph_lazy_inclist_t inclist;
  igraph_vector_long_t papers;
  
  if (igraph_vector_size(vtime) != igraph_vcount(graph)) {
    IGRAPH_ERROR("Invalid vtime length", IGRAPH_EINVAL);
  }
  if (igraph_vector_size(etime) != igraph_ecount(graph)) {
    IGRAPH_ERROR("Invalid etime length", IGRAPH_EINVAL);
  }

  no_of_events=igraph_vector_size(eventsizes);
  
  IGRAPH_VECTOR_INIT_FINALLY(&st, no_of_events);
  for (i=0; i<no_of_events; i++) {
    VECTOR(st)[i]=1;
  }
  
  IGRAPH_CHECK(igraph_vector_long_init(&papers, igraph_vcount(graph)));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &papers);
  for (i=0; i<igraph_vector_size(authors); i++) {
    long int author=VECTOR(*authors)[i];
    VECTOR(papers)[author] += 1;
    if (VECTOR(papers)[author] > maxpapers) {
      maxpapers=VECTOR(papers)[author];
    }
  }
  igraph_vector_long_destroy(&papers);
  IGRAPH_FINALLY_CLEAN(1);
  
  IGRAPH_VECTOR_INIT_FINALLY(&vtimeidx, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&etimeidx, 0);
  IGRAPH_CHECK(igraph_vector_order1(vtime, &vtimeidx, no_of_events));
  IGRAPH_CHECK(igraph_vector_order1(etime, &etimeidx, no_of_events));
  
  IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

  IGRAPH_PROGRESS("Revolver p-p", 0, NULL);
  for (i=0; i<niter; i++) {
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    if (i+1 != niter) {		/* not the last iteration */
      /* measure */
      IGRAPH_CHECK(igraph_revolver_mes_p_p(graph, &inclist, 
					   kernel, 0 /*sd*/, 0 /*norm*/,
					   0/*cites*/, 0/*debug*/, 0 /*debugres*/,
					   &st, vtime, &vtimeidx, etime, 
					   &etimeidx, no_of_events,
					   authors, eventsizes,
					   maxpapers));
      /* normalize */
      igraph_matrix_scale(kernel, 1/igraph_matrix_sum(kernel));
      
      /* update st */
      IGRAPH_CHECK(igraph_revolver_st_p_p(graph, &inclist, 
					  &st, kernel, vtime, &vtimeidx,
					  etime, &etimeidx,
					  no_of_events, authors, 
					  eventsizes, maxpapers));
      
    } else {
      /* measure */
      IGRAPH_CHECK(igraph_revolver_mes_p_p(graph, &inclist,
					   kernel, sd, norm, cites, 
					   debug, debugres, &st, vtime, &vtimeidx,
					   etime, &etimeidx,
					   no_of_events, authors,
					   eventsizes, maxpapers));
      
      /* normalize */
      igraph_matrix_scale(kernel, 1/igraph_matrix_sum(kernel));
      
      /* update st */
      IGRAPH_CHECK(igraph_revolver_st_p_p(graph, &inclist,
					  &st, kernel, vtime, &vtimeidx,
					  etime, &etimeidx,
					  no_of_events, authors, eventsizes,
					  maxpapers));
      
      /* expected number of citations */
      if (expected) {
	IGRAPH_CHECK(igraph_revolver_exp_p_p(graph, &inclist,
					     expected, kernel, &st,
					     vtime, &vtimeidx, etime, &etimeidx,
					     no_of_events, authors, eventsizes,
					     maxpapers));
      }
      
      /* error calculation */
      if (logprob || lognull) {
	IGRAPH_CHECK(igraph_revolver_error_p_p(graph, &inclist,
					       kernel, &st,
					       vtime, &vtimeidx, 
					       etime, &etimeidx, no_of_events,
					       authors, eventsizes,
					       maxpapers, logprob, lognull));
      }
    }

    IGRAPH_PROGRESS("Revolver p-p", 100.0*(i+1)/niter, NULL);
  }

  igraph_lazy_inclist_destroy(&inclist);
  igraph_vector_destroy(&etimeidx);
  igraph_vector_destroy(&vtimeidx);
  igraph_vector_destroy(&st);
  IGRAPH_FINALLY_CLEAN(4);
  
  return 0;
}

#define NTKK(xidx, yidx) \
   ((xidx)==(yidx) ? (VECTOR(ntk)[(xidx)]*(VECTOR(ntk)[(xidx)]-1))/2-MATRIX(ntkk,(xidx),(yidx)) : VECTOR(ntk)[(xidx)]*VECTOR(ntk)[(yidx)]-MATRIX(ntkk,(xidx),(yidx)))

int igraph_revolver_mes_p_p(const igraph_t *graph,
			    igraph_lazy_inclist_t *inclist,
			    igraph_matrix_t *kernel,
			    igraph_matrix_t *sd,
			    igraph_matrix_t *norm,
			    igraph_matrix_t *cites,
			    const igraph_matrix_t *debug,
			    igraph_vector_ptr_t *debugres,
			    const igraph_vector_t *st,
			    const igraph_vector_t *vtime,
			    const igraph_vector_t *vtimeidx,
			    const igraph_vector_t *etime,
			    const igraph_vector_t *etimeidx,
			    igraph_integer_t pno_of_events,
			    const igraph_vector_t *authors,
			    const igraph_vector_t *eventsizes,
			    igraph_integer_t pmaxpapers) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  long int no_of_events=pno_of_events;
  long int maxpapers=pmaxpapers;
  
  igraph_vector_long_t papers;
  igraph_vector_char_t added;
  
  igraph_matrix_t v_normfact, *normfact, v_notnull, *notnull;
  igraph_matrix_t ch;
  
  igraph_vector_long_t ntk;
  igraph_matrix_t ntkk;
  
  igraph_vector_t *adjedges;
  
  long int timestep, i;
  long int nptr=0, eptr=0, aptr=0;
  long int nptr_save, eptr_save, eptr_new;

  IGRAPH_CHECK(igraph_vector_long_init(&papers, no_of_nodes));
  IGRAPH_FINALLY(&igraph_vector_long_destroy, &papers);

  IGRAPH_CHECK(igraph_vector_char_init(&added, no_of_edges));
  IGRAPH_FINALLY(igraph_vector_char_destroy, &added);

  IGRAPH_CHECK(igraph_vector_long_init(&ntk, maxpapers+1));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &ntk);
  IGRAPH_MATRIX_INIT_FINALLY(&ntkk, maxpapers+1, maxpapers+1);
  IGRAPH_MATRIX_INIT_FINALLY(&ch, maxpapers+1, maxpapers+1);
  
  if (norm) {
    normfact=norm;
    IGRAPH_CHECK(igraph_matrix_resize(normfact, maxpapers+1, maxpapers+1));
    igraph_matrix_null(normfact);
  } else {
    normfact=&v_normfact;
    IGRAPH_MATRIX_INIT_FINALLY(normfact, maxpapers+1, maxpapers+1);
  }
  
  if (cites) {
    notnull=cites;
    IGRAPH_CHECK(igraph_matrix_resize(notnull, maxpapers+1, maxpapers+1));
    igraph_matrix_null(notnull);
  } else {
    notnull=&v_notnull;
    IGRAPH_MATRIX_INIT_FINALLY(notnull, maxpapers+1, maxpapers+1);
  }
  
  IGRAPH_CHECK(igraph_matrix_resize(kernel, maxpapers+1, maxpapers+1));
  igraph_matrix_null(kernel);
  if (sd) {
    IGRAPH_CHECK(igraph_matrix_resize(sd, maxpapers+1, maxpapers+1));
    igraph_matrix_null(sd);
  }
  
  for (timestep=0; timestep<no_of_events; timestep++) {
    
    IGRAPH_ALLOW_INTERRUPTION();
        
    nptr_save=nptr;
    while (nptr < no_of_nodes && 
	   VECTOR(*vtime)[(long int)VECTOR(*vtimeidx)[nptr]]==timestep) {
      nptr++;
    }
    /* If it is a new author then she has no papers yet */
    VECTOR(ntk)[0] += (nptr-nptr_save);

    /* Update ch accordingly, could be done later as well */
    if (VECTOR(ntk)[0] == nptr-nptr_save && nptr!=nptr_save) {
      if (nptr-nptr_save >= 2) {
	MATRIX(ch, 0, 0) = eptr;
      }
      for (i=1; i<maxpapers+1; i++) {
	if (NTKK(0,i) == (nptr-nptr_save)*VECTOR(ntk)[i]) {
	  MATRIX(ch, 0, i) = MATRIX(ch, i, 0) = eptr;
	}
      }
    }    
    
/*     print_ntkk(&ntkk, &ntk); */

    /* Estimate Akk */
    eptr_save=eptr;
    while (eptr < no_of_edges &&
	   VECTOR(*etime)[(long int)VECTOR(*etimeidx)[eptr] ] == timestep) {
      long int edge=VECTOR(*etimeidx)[eptr];
      long int from=IGRAPH_FROM(graph, edge), to=IGRAPH_TO(graph, edge);
      long int xidx=VECTOR(papers)[from];
      long int yidx=VECTOR(papers)[to];
      double xk, oldakk;

      MATRIX(*notnull, xidx, yidx) += 1;
      MATRIX(*notnull, yidx, xidx) = MATRIX(*notnull, xidx, yidx);
      
      xk=VECTOR(*st)[timestep]/NTKK(xidx, yidx);
      oldakk=MATRIX(*kernel, xidx, yidx);
      MATRIX(*kernel, xidx, yidx) +=  (xk-oldakk)/MATRIX(*notnull, xidx, yidx);
      MATRIX(*kernel, yidx, xidx) = MATRIX(*kernel, xidx, yidx);
      if (sd) {
	MATRIX(*sd, xidx, yidx) += (xk-oldakk)*(xk-MATRIX(*kernel, xidx, yidx));
	MATRIX(*sd, yidx, xidx) = MATRIX(*sd, xidx, yidx);
      }
      /* TODO: debug */

      eptr++;
    }

    /* update ntkk, the new papers change the type of their authors */
    eptr_new=eptr;
    for (i=aptr; i<aptr+VECTOR(*eventsizes)[timestep]; i++) {
      long int aut=VECTOR(*authors)[i];
      long int pap=VECTOR(papers)[aut];
      long int j, n;

      adjedges=igraph_lazy_inclist_get(inclist, aut);
      n=igraph_vector_size(adjedges);
      for (j=0; j<n; j++) {
	long int edge=VECTOR(*adjedges)[j];
	if (VECTOR(added)[edge]) {
	  long int otherv=IGRAPH_OTHER(graph, edge, aut);
	  long int otherpap=VECTOR(papers)[otherv];
	  MATRIX(ntkk, pap, otherpap) -= 1;
	  MATRIX(ntkk, otherpap, pap) = MATRIX(ntkk, pap, otherpap);
	  if (NTKK(pap, otherpap)==1) {
	    MATRIX(ch, pap, otherpap) = eptr_new;
	    MATRIX(ch, otherpap, pap) = MATRIX(ch, pap, otherpap);
	  }
	  MATRIX(ntkk, pap+1, otherpap) += 1;
	  MATRIX(ntkk, otherpap, pap+1) = MATRIX(ntkk, pap+1, otherpap);
	  if (NTKK(pap+1, otherpap)==0) {
	    MATRIX(*normfact, pap+1, otherpap) += 
	      eptr_new-MATRIX(ch, pap+1, otherpap);
	    MATRIX(*normfact, otherpap, pap+1) = 
	      MATRIX(*normfact, pap+1, otherpap);
	  }
	}
      }

      /* update ntk too */
      for (j=0; j<maxpapers+1; j++) {
	long int before, after;
	before=NTKK(pap, j);
	VECTOR(ntk)[pap]-=1;
	after=NTKK(pap, j);
	VECTOR(ntk)[pap]+=1;
	if (before > 0 && after==0) {
	  MATRIX(*normfact, pap, j) += eptr_new-MATRIX(ch, pap, j);
	  MATRIX(*normfact, j, pap) = MATRIX(*normfact, pap, j);
	}
      }
      VECTOR(ntk)[pap]-=1;
      
      for (j=0; j<maxpapers+1; j++) { 
	long int before, after;
	before=NTKK(pap+1, j);
	VECTOR(ntk)[pap+1] += 1;
	after=NTKK(pap+1, j);
	VECTOR(ntk)[pap+1] -= 1;
	if (before == 0 && after > 0) {
	  MATRIX(ch, pap+1, j) = eptr_new;
	  MATRIX(ch, j, pap+1) = MATRIX(ch, pap+1, j);
	}
      }
      VECTOR(ntk)[pap+1]+=1;
      
      VECTOR(papers)[aut] += 1;
    }
    aptr += VECTOR(*eventsizes)[timestep];
    
    /* For every new edge, we lose one connection possibility, also add the edges*/
    eptr=eptr_save;
    while (eptr < no_of_edges && 
	   VECTOR(*etime)[(long int) VECTOR(*etimeidx)[eptr] ] == timestep) {
      long int edge=VECTOR(*etimeidx)[eptr];
      long int from=IGRAPH_FROM(graph, edge), to=IGRAPH_TO(graph, edge);
      long int xidx=VECTOR(papers)[from];
      long int yidx=VECTOR(papers)[to];
      
      MATRIX(ntkk, xidx, yidx) += 1;
      MATRIX(ntkk, yidx, xidx) = MATRIX(ntkk, xidx, yidx);
      if (NTKK(xidx, yidx)==0) {
	MATRIX(*normfact, xidx, yidx) += eptr_new-MATRIX(ch, xidx, yidx);
	MATRIX(*normfact, yidx, xidx) = MATRIX(*normfact, xidx, yidx);
      }

      VECTOR(added)[edge]=1;
      eptr++;
    }
  }    

  for (i=0; i<maxpapers+1; i++) {
    igraph_real_t oldakk;
    long int j;
    for (j=0; j<=i; j++) {
      if (NTKK(i, j) != 0) {
	MATRIX(*normfact, i, j) += (eptr-MATRIX(ch, i, j));
	MATRIX(*normfact, j, i) = MATRIX(*normfact, i, j);
      }
      if (MATRIX(*normfact, i, j)==0) {
	MATRIX(*kernel, i, j)=MATRIX(*kernel, j, i)=0;
	MATRIX(*normfact, i, j)=MATRIX(*normfact, j, i)=1;
      }
      oldakk=MATRIX(*kernel, i, j);
      MATRIX(*kernel, i, j) *= MATRIX(*notnull, i, j)/MATRIX(*normfact, i, j);
      MATRIX(*kernel, j, i) = MATRIX(*kernel, i, j);
      if (sd) {
	MATRIX(*sd, i, j) += oldakk * oldakk * MATRIX(*notnull, i, j) *
	  (1-MATRIX(*notnull, i, j)/MATRIX(*normfact, i, j));
	MATRIX(*sd, i, j) = sqrt(MATRIX(*sd, i, j)/(MATRIX(*normfact, i, j)-1));
	MATRIX(*sd, j, i) = MATRIX(*sd, i, j);
      }
    }
  }

  if (!cites) {
    igraph_matrix_destroy(notnull);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (!norm) {
    igraph_matrix_destroy(normfact);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  igraph_matrix_destroy(&ch);
  igraph_matrix_destroy(&ntkk);
  igraph_vector_long_destroy(&ntk);
  igraph_vector_char_destroy(&added);
  igraph_vector_long_destroy(&papers);
  IGRAPH_FINALLY_CLEAN(5);
  
  return 0;
}

#undef NTKK

int igraph_revolver_st_p_p(const igraph_t *graph,
			   igraph_lazy_inclist_t *inclist,
			   igraph_vector_t *st,
			   const igraph_matrix_t *kernel,
			   const igraph_vector_t *vtime,
			   const igraph_vector_t *vtimeidx,
			   const igraph_vector_t *etime,
			   const igraph_vector_t *etimeidx,
			   igraph_integer_t pno_of_events,
			   const igraph_vector_t *authors,
			   const igraph_vector_t *eventsizes,
			   igraph_integer_t pmaxpapers) {

  long int no_of_events=pno_of_events;
  long int maxpapers=igraph_matrix_nrow(kernel)-1;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  long int timestep=0;
  
  igraph_vector_long_t papers;
  igraph_vector_long_t ntk;
  igraph_vector_char_t added;
  
  igraph_vector_t *adjedges;
  
  long int i;
  long int nptr=0, eptr=0, aptr=0, nptr_save;
  
  IGRAPH_CHECK(igraph_vector_long_init(&ntk, maxpapers+1));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &ntk);
  IGRAPH_CHECK(igraph_vector_long_init(&papers, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &papers);
  IGRAPH_CHECK(igraph_vector_char_init(&added, no_of_edges));
  IGRAPH_FINALLY(igraph_vector_char_destroy, &added);
  
  IGRAPH_CHECK(igraph_vector_resize(st, no_of_events));
  VECTOR(*st)[0]=0;
  
  for (timestep=0; timestep<no_of_events-1; timestep++) {
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    /* add the new nodes */
    nptr_save=nptr;
    while (nptr < no_of_nodes &&
	   VECTOR(*vtime)[ (long int) VECTOR(*vtimeidx)[nptr] ] == timestep) {
      nptr++;
    }
    nptr_save=nptr-nptr_save;
    if (nptr_save != 0) {
      for (i=0; i<maxpapers+1; i++) {
	VECTOR(*st)[timestep] += VECTOR(ntk)[i]*MATRIX(*kernel, i, 0)*nptr_save;
      }
      VECTOR(*st)[timestep] += nptr_save*(nptr_save-1)/2 * MATRIX(*kernel, 0, 0);
      VECTOR(ntk)[0]+=nptr_save;
    }    
    
    VECTOR(*st)[timestep+1] = VECTOR(*st)[timestep];
    
    for (i=aptr; i<aptr+VECTOR(*eventsizes)[timestep]; i++) {
      long int aut=VECTOR(*authors)[i];
      long int pap=VECTOR(papers)[aut];
      long int j, n;
      
      for (j=0; j<maxpapers+1; j++) {
	VECTOR(*st)[timestep+1] += VECTOR(ntk)[j] * (MATRIX(*kernel, j, pap+1)-
						     MATRIX(*kernel, j, pap));
      }
      
      VECTOR(*st)[timestep+1] += MATRIX(*kernel, pap, pap);
      VECTOR(*st)[timestep+1] -= MATRIX(*kernel, pap+1, pap+1);

      VECTOR(ntk)[pap]--;
      VECTOR(ntk)[pap+1]++;

      adjedges=igraph_lazy_inclist_get(inclist, aut);
      n=igraph_vector_size(adjedges);
      for (j=0; j<n; j++) {
	long int edge=VECTOR(*adjedges)[j];
	if (VECTOR(added)[edge]) {
	  long int otherv=IGRAPH_OTHER(graph, edge, aut);
	  long int otherpap=VECTOR(papers)[otherv];
	  VECTOR(*st)[timestep+1] += MATRIX(*kernel, pap, otherpap);
	  VECTOR(*st)[timestep+1] -= MATRIX(*kernel, pap+1, otherpap);
	}
      }

      VECTOR(papers)[aut] += 1;      
    }
    aptr += VECTOR(*eventsizes)[timestep];
    
    while (eptr < no_of_edges && 
	   VECTOR(*etime)[ (long int) VECTOR(*etimeidx)[eptr] ] == timestep) {
      long int edge=VECTOR(*etimeidx)[eptr];
      long int from=IGRAPH_FROM(graph, edge);
      long int to=IGRAPH_TO(graph, edge);
      long int xidx=VECTOR(papers)[from];
      long int yidx=VECTOR(papers)[to];
      VECTOR(*st)[timestep+1] -= MATRIX(*kernel, xidx, yidx);
      VECTOR(added)[edge]=1;
      eptr++;
    }
    
  }
  
  igraph_vector_char_destroy(&added);
  igraph_vector_long_destroy(&papers);
  igraph_vector_long_destroy(&ntk);
  IGRAPH_FINALLY_CLEAN(3);
  
  return 0;
}

int igraph_revolver_exp_p_p(const igraph_t *graph,
			    igraph_lazy_inclist_t *inclist,
			    igraph_matrix_t *expected,
			    const igraph_matrix_t *kernel,
			    const igraph_vector_t *st,
			    const igraph_vector_t *vtime,
			    const igraph_vector_t *vtimeidx,
			    const igraph_vector_t *etime,
			    const igraph_vector_t *etimeidx,
			    igraph_integer_t pno_of_events,
			    const igraph_vector_t *authors,
			    const igraph_vector_t *eventsizes,
			    igraph_integer_t pmaxpapers) {

  /* TODO */
  return 0;
}

int igraph_revolver_error_p_p(const igraph_t *graph,
			      igraph_lazy_inclist_t *inclist,
			      const igraph_matrix_t *kernel,
			      const igraph_vector_t *st,
			      const igraph_vector_t *vtime,
			      const igraph_vector_t *vtimeidx,
			      const igraph_vector_t *etime,
			      const igraph_vector_t *etimeidx,
			      igraph_integer_t pno_of_events,
			      const igraph_vector_t *authors,
			      const igraph_vector_t *eventsizes,
			      igraph_integer_t pmaxpapers,
			      igraph_real_t *logprob,
			      igraph_real_t *lognull) {

  long int no_of_events=pno_of_events;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);

  igraph_vector_long_t papers;
  
  long int timestep, nptr=0, eptr=0, aptr=0, eptr_save;
  long int edges=0, vertices=0, i;
  
  igraph_real_t rlogprob, rlognull, *mylogprob=logprob, *mylognull=lognull;

  IGRAPH_CHECK(igraph_vector_long_init(&papers, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &papers);

  if (!logprob) { mylogprob=&rlogprob; }
  if (!lognull) { mylognull=&rlognull; }
  
  *mylogprob=0;
  *mylognull=0;
  
  for (timestep=0; timestep<no_of_events; timestep++) {
    
    IGRAPH_ALLOW_INTERRUPTION();
  
    while (nptr < no_of_nodes && 
	   VECTOR(*vtime)[ (long int) VECTOR(*vtimeidx)[nptr] ] == timestep) {
      vertices++;
      nptr++;
    }
    
    eptr_save=eptr;
    while (eptr < no_of_edges && 
	   VECTOR(*etime)[ (long int) VECTOR(*etimeidx)[eptr] ] == timestep) {
      long int edge=VECTOR(*etimeidx)[eptr];
      long int from=IGRAPH_FROM(graph, edge);
      long int to=IGRAPH_TO(graph, edge);
      long int xidx=VECTOR(papers)[from];
      long int yidx=VECTOR(papers)[to];
      
      igraph_real_t prob=MATRIX(*kernel, xidx, yidx)/VECTOR(*st)[timestep];
      igraph_real_t nullprob=1.0/(vertices*(vertices-1)/2-eptr_save);
            
      *mylogprob += log(prob);
      *mylognull += log(nullprob);
      
      edges++;
      eptr++;
    }
    
    for (i=aptr; i<aptr+VECTOR(*eventsizes)[timestep]; i++) {
      long int aut=VECTOR(*authors)[i];
      VECTOR(papers)[aut] += 1;
    }
    aptr += VECTOR(*eventsizes)[timestep];
    
  }
  
  igraph_vector_long_destroy(&papers);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}
