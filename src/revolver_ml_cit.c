/* -*- mode: C -*-  */
/* 
   IGraph library.
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
#include "igraph_memory.h"
#include "igraph_random.h"
#include "config.h"
#include "igraph_math.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_interrupt_internal.h"
#include "igraph_structural.h"
#include "igraph_nongraph.h"
#include "igraph_stack.h"

#include <math.h>

/* This data structure holds the context of the degree based
   optimizer. */

typedef struct igraph_i_revolver_ml_D_data_t {
  igraph_scalar_function_t *A;
  igraph_vector_function_t *dA;  
  const igraph_t *graph;
  long int no_of_nodes;
  igraph_vector_t A_vect;	/* Temporary storage */
  igraph_vector_ptr_t dA_vects;	/* Temporary storage */
  igraph_integer_t maxdegree;
  igraph_vector_long_t degree;
  igraph_vector_t neis;
  igraph_vector_t dS;		/* Temporary storage */
  igraph_vector_t par1;		/* More tmp storage */
  igraph_vector_t tmpgrad;      /* More... */

  igraph_vector_t lastparam;	/* The parameter values used last time */
  igraph_real_t lastf;		/* The evaluated function value  */
  igraph_vector_t lastgrad;	/* The evaluated gradient */

  const igraph_vector_t *filter;
} igraph_i_revolver_ml_D_data_t;

/* Evaluate the objective function and calculate its gradient too. */

int igraph_i_revolver_ml_D_eval(const igraph_vector_t *par,
				igraph_i_revolver_ml_D_data_t *data) {

  long int no_of_edges=0;
  igraph_real_t sum=0.0;
  long int t, i;
  int dim=igraph_vector_size(par);
  igraph_vector_t *grad=&data->lastgrad;
  igraph_real_t S=0.0;		/* sum N_i^t A(p,x) */
  
  /* Init */
  igraph_vector_long_null(&data->degree);
  igraph_vector_null(&data->dS);
  igraph_vector_null(grad);

  /* Calculate all possible A and dA values and store them in A_vect &
     dA_vects */
  for (t=0; t<=data->maxdegree; t++) {
    VECTOR(data->par1)[0] = t;
    VECTOR(data->A_vect)[t] = data->A(&data->par1, par, 0);
    data->dA(&data->par1, par, &data->tmpgrad, 0);
    for (i=0; i<dim; i++) {
      igraph_vector_t *v=VECTOR(data->dA_vects)[i];
      VECTOR(*v)[t] = VECTOR(data->tmpgrad)[i];
    }
  }
  
  for (t=0; t<data->no_of_nodes; t++) {
    long int n, nneis;

    IGRAPH_ALLOW_INTERRUPTION();
      
    IGRAPH_CHECK(igraph_neighbors(data->graph, &data->neis, t, IGRAPH_OUT));
    nneis=igraph_vector_size(&data->neis);

    if (! data->filter || VECTOR(*data->filter)[t] != 0) {
      
      /* Update sum(s) */
      for (n=0; n<nneis; n++) {
	long int to=VECTOR(data->neis)[n];
	long int x=VECTOR(data->degree)[to];
	
	sum -= log( VECTOR(data->A_vect)[x] );
	sum += log( S );
	for (i=0; i<dim; i++) {
	  igraph_vector_t *v=VECTOR(data->dA_vects)[i];
	  VECTOR(*grad)[i] -= VECTOR(*v)[x] / VECTOR(data->A_vect)[x];
	  VECTOR(*grad)[i] += VECTOR(data->dS)[i] / S;
	}
	no_of_edges++;
      }

    }

    /* Update S, data->dS */
    for (n=0; n<nneis; n++) {
      long int to=VECTOR(data->neis)[n];
      long int x=VECTOR(data->degree)[to];
      
      VECTOR(data->degree)[to] += 1;
      S += VECTOR(data->A_vect)[x+1];
      S -= VECTOR(data->A_vect)[x];
      for (i=0; i<dim; i++) {
	igraph_vector_t *v=VECTOR(data->dA_vects)[i];
	VECTOR(data->dS)[i] += VECTOR(*v)[x+1];
	VECTOR(data->dS)[i] -= VECTOR(*v)[x];
      }
    }
    
    S += VECTOR(data->A_vect)[0];
    for (i=0; i<dim; i++) {
      igraph_vector_t *v=VECTOR(data->dA_vects)[i];
      VECTOR(data->dS)[i] += VECTOR(*v)[0];
    }
  }

  igraph_vector_update(&data->lastparam, par);
  data->lastf=sum / no_of_edges;
  for (i=0; i<igraph_vector_size(&data->lastgrad); i++) {
    VECTOR(data->lastgrad)[i] /= no_of_edges;
  }

  return 0;
}

/* This function gives the value of the objective function at the 
   supplied parameter vector. Called by the optimizer.
*/

igraph_real_t igraph_i_revolver_ml_D_f(const igraph_vector_t *par,
				       const igraph_vector_t *garbage,
				       void* extra) {

  igraph_i_revolver_ml_D_data_t *data=extra;
  
  if (!igraph_vector_all_e(par, &data->lastparam)) {
    igraph_i_revolver_ml_D_eval(par, data);
  }

  return data->lastf;
}

/* This function gives the gradient of the objective function at
   the supplied parameter vector. Called by the optimizer.
*/

void igraph_i_revolver_ml_D_df(const igraph_vector_t *par,
			       const igraph_vector_t *garbage,
			       igraph_vector_t *res, void *extra) {

  igraph_i_revolver_ml_D_data_t *data=extra;

  if (!igraph_vector_all_e(par, &data->lastparam)) {
    igraph_i_revolver_ml_D_eval(par, data);
  }

  igraph_vector_update(res, &data->lastgrad);
}

void igraph_i_revolver_ml_D_free(igraph_vector_ptr_t *ptr) {
  long int i, n=igraph_vector_ptr_size(ptr);
  for (i=0; i<n; i++) {
    igraph_vector_t *v=VECTOR(*ptr)[i];
    if (v) {
      igraph_vector_destroy(v);
      igraph_free(v);
    }
    VECTOR(*ptr)[i]=0;
  }
}

/* This is the general degree-based optimizer. The form of the kernel
   function is given by A_fun, its first parameter (in the parameter vector) 
   is the degree and the others are the parameters to fit.

   It just initializes the context data structure and calls the optimizer.
 */

int igraph_revolver_ml_D(const igraph_t *graph,
			 igraph_vector_t *res,
			 igraph_real_t *Fmin,
			 igraph_real_t abstol, igraph_real_t reltol, int maxit,
			 igraph_scalar_function_t *A_fun,
			 igraph_vector_function_t *dA_fun, 
			 const igraph_vector_t *filter,
			 igraph_integer_t *fncount, igraph_integer_t *grcount) {

  igraph_i_revolver_ml_D_data_t info;
  igraph_integer_t maxdegree;
  long int no_of_nodes=igraph_vcount(graph);
  int dim=igraph_vector_size(res);
  int ret, i;

  if (filter && igraph_vector_size(filter) != no_of_nodes) {
    IGRAPH_ERROR("Invalid filter vector", IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_maxdegree(graph, &maxdegree, igraph_vss_all(),
				IGRAPH_IN, IGRAPH_LOOPS));
  
  /* Set up everything */
  info.A=A_fun;
  info.dA=dA_fun;
  info.graph=graph;
  info.no_of_nodes=no_of_nodes;
  IGRAPH_VECTOR_INIT_FINALLY(&info.A_vect, maxdegree+1);  
  IGRAPH_VECTOR_PTR_INIT_FINALLY(&info.dA_vects, dim);
  IGRAPH_FINALLY(igraph_i_revolver_ml_D_free, &info.dA_vects);
  for (i=0; i<dim; i++) {
    igraph_vector_t *v=igraph_Calloc(1, igraph_vector_t);
    if (!v) { IGRAPH_ERROR("Cannot perform ML D revolver", IGRAPH_ENOMEM); }
    IGRAPH_CHECK(igraph_vector_init(v, maxdegree+1));
    VECTOR(info.dA_vects)[i]=v;
  }
  info.maxdegree=maxdegree;
  IGRAPH_CHECK(igraph_vector_long_init(&info.degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &info.degree);
  IGRAPH_VECTOR_INIT_FINALLY(&info.neis, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&info.dS, dim);
  IGRAPH_VECTOR_INIT_FINALLY(&info.par1, dim);
  IGRAPH_VECTOR_INIT_FINALLY(&info.tmpgrad, dim);
  IGRAPH_VECTOR_INIT_FINALLY(&info.lastparam, dim);
  info.lastf=0.0;
  IGRAPH_VECTOR_INIT_FINALLY(&info.lastgrad, dim);
  info.filter=filter;

  igraph_i_revolver_ml_D_eval(res, &info);
  ret=igraph_bfgs(res, Fmin, igraph_i_revolver_ml_D_f,
		  igraph_i_revolver_ml_D_df, maxit, 1, abstol, reltol, 1,
		  &info, fncount, grcount);
  
  igraph_vector_destroy(&info.lastgrad);
  igraph_vector_destroy(&info.lastparam);
  igraph_vector_destroy(&info.tmpgrad);
  igraph_vector_destroy(&info.par1);
  igraph_vector_destroy(&info.dS);
  igraph_vector_destroy(&info.neis);
  igraph_vector_long_destroy(&info.degree);
  igraph_i_revolver_ml_D_free(&info.dA_vects);
  igraph_vector_ptr_destroy(&info.dA_vects);
  igraph_vector_destroy(&info.A_vect);
  IGRAPH_FINALLY_CLEAN(10);
  
  return ret;
}

/* These functions assemble the A(d)=d^alpha+1 kernel function and 
   calls the general degree-based optimizer.
*/

igraph_real_t igraph_i_revolver_ml_D_alpha_f(const igraph_vector_t *var, 
					     const igraph_vector_t *par,
					     void *extra) {
  igraph_real_t deg=VECTOR(*var)[0];
  igraph_real_t alpha=VECTOR(*par)[0];
  if (deg != 0) {
    return pow(deg, alpha) + 1.0;
  } else {
    return 1.0;
  }
}

void igraph_i_revolver_ml_D_alpha_df(const igraph_vector_t *var,
				     const igraph_vector_t *par,
				     igraph_vector_t *res, 
				     void *extra) {
  igraph_real_t deg=VECTOR(*var)[0];
  igraph_real_t alpha=VECTOR(*par)[0];
  if (deg != 0) {
    VECTOR(*res)[0] = log(deg) * pow(deg, alpha);
  } else {
    VECTOR(*res)[0] = 0.0;
  }
}

int igraph_revolver_ml_D_alpha(const igraph_t *graph,
			       igraph_real_t *alpha, igraph_real_t *Fmin,
			       igraph_real_t abstol, igraph_real_t reltol, 
			       int maxit, const igraph_vector_t *filter,
			       igraph_integer_t *fncount, 
			       igraph_integer_t *grcount) {
  
  igraph_vector_t res;
  int ret;

  IGRAPH_VECTOR_INIT_FINALLY(&res, 1);
  VECTOR(res)[0]=*alpha;
  
  ret=igraph_revolver_ml_D(graph, &res, Fmin, abstol, reltol, maxit,
			   igraph_i_revolver_ml_D_alpha_f,
			   igraph_i_revolver_ml_D_alpha_df,
			   filter,
			   fncount, grcount);

  *alpha=VECTOR(res)[0];
  igraph_vector_destroy(&res);
  IGRAPH_FINALLY_CLEAN(1);

  return ret;
}

/* These functions assemble the A(d)=d^alpha+a kernel function and 
   calls the general degree-based optimizer.
*/

igraph_real_t igraph_i_revolver_ml_D_alpha_a_f(const igraph_vector_t *var,
					       const igraph_vector_t *par,
					       void *extra) {
  igraph_real_t deg=VECTOR(*var)[0];
  igraph_real_t alpha=VECTOR(*par)[0];
  igraph_real_t a=VECTOR(*par)[1];
  if (deg != 0) {
    return pow(deg, alpha) + a;
  } else {
    return a;
  }  
}

void igraph_i_revolver_ml_D_alpha_a_df(const igraph_vector_t *var,
				       const igraph_vector_t *par,
				       igraph_vector_t *res,
				       void *extra) {
  igraph_real_t deg=VECTOR(*var)[0];
  igraph_real_t alpha=VECTOR(*par)[0];
  /* a not needed */
  if (deg != 0) {
    VECTOR(*res)[0] = log(deg) * pow(deg, alpha);
    VECTOR(*res)[1] = 1.0;
  } else {
    VECTOR(*res)[0] = 0.0;
    VECTOR(*res)[0] = 1.0;
  }    
}

int igraph_revolver_ml_D_alpha_a(const igraph_t *graph,
				 igraph_real_t *alpha, igraph_real_t *a,
				 igraph_real_t *Fmin,
				 igraph_real_t abstol, igraph_real_t reltol,
				 int maxit, const igraph_vector_t *filter,
				 igraph_integer_t *fncount, 
				 igraph_integer_t *grcount) {
  igraph_vector_t res;
  int ret;
  
  IGRAPH_VECTOR_INIT_FINALLY(&res, 2);
  VECTOR(res)[0] = *alpha;
  VECTOR(res)[1] = *a;
  
  ret=igraph_revolver_ml_D(graph, &res, Fmin, abstol, reltol, maxit,
			   igraph_i_revolver_ml_D_alpha_a_f,
			   igraph_i_revolver_ml_D_alpha_a_df, 
			   filter,
			   fncount, grcount);
  
  *alpha=VECTOR(res)[0];
  *a=VECTOR(res)[1];
  igraph_vector_destroy(&res);
  IGRAPH_FINALLY_CLEAN(1);
  
  return ret;
}

/*------------------------------------------------------------------*/

typedef struct igraph_i_revolver_ml_DE_data_t {
  igraph_scalar_function_t *A;
  igraph_vector_function_t *dA;  
  const igraph_t *graph;
  const igraph_vector_t *cats;
  long int no_of_nodes;
  igraph_matrix_t A_vect;	/* Temporary storage */
  igraph_vector_ptr_t dA_vects;	/* Temporary storage */
  igraph_integer_t maxdegree;
  igraph_integer_t nocats;
  igraph_vector_long_t degree;
  igraph_vector_t neis;
  igraph_vector_t dS;		/* Temporary storage */
  igraph_vector_t par1;		/* More tmp storage */
  igraph_vector_t tmpgrad;      /* More... */  

  igraph_vector_t lastparam;	/* The parameter values used last time */
  igraph_real_t lastf;		/* The evaluated function value  */
  igraph_vector_t lastgrad;	/* The evaluated gradient */

  const igraph_vector_t *filter;
} igraph_i_revolver_ml_DE_data_t;

int igraph_i_revolver_ml_DE_eval(const igraph_vector_t *par,
				 igraph_i_revolver_ml_DE_data_t *data) {
  
  igraph_real_t sum=0.0;
  long int t, i, j;
  int dim=igraph_vector_size(par);
  igraph_vector_t *grad=&data->lastgrad;
  igraph_real_t S=0.0;
  long int no_of_edges=0;
  
  /* Init */
  igraph_vector_long_null(&data->degree);
  igraph_vector_null(&data->dS);
  igraph_vector_null(grad);

  for (i=0; i<data->nocats; i++) {
    for (j=0; j<data->maxdegree+1; j++) {
      long int k;
      VECTOR(data->par1)[0]=i; VECTOR(data->par1)[1]=j;
      MATRIX(data->A_vect, i, j) = data->A(&data->par1, par, 0);
      data->dA(&data->par1, par, &data->tmpgrad, 0);
      for (k=0; k<dim; k++) {
	igraph_matrix_t *m=VECTOR(data->dA_vects)[k];
	MATRIX(*m, i, j)=VECTOR(data->tmpgrad)[k];
      }
    }
  }
  
  for (t=0; t<data->no_of_nodes; t++) {
    long int n, nneis;
    long int tcat=VECTOR(*data->cats)[t];
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    IGRAPH_CHECK(igraph_neighbors(data->graph, &data->neis, t, IGRAPH_OUT));
    nneis=igraph_vector_size(&data->neis);

    if (! data->filter || VECTOR(*data->filter)[t]) {

      /* Update sum(s) */
      for (n=0; n<nneis; n++) {
	long int to=VECTOR(data->neis)[n];
	long int x=VECTOR(*data->cats)[to];
	long int y=VECTOR(data->degree)[to];
	
/* 	CHECK_VALID(x,y); */
	sum -= log( MATRIX(data->A_vect, x, y) );
	sum += log( S );
	for (i=0; i<dim; i++) {
	  igraph_matrix_t *m=VECTOR(data->dA_vects)[i];
	  VECTOR(*grad)[i] -= MATRIX(*m, x, y) / MATRIX(data->A_vect, x, y);
	  VECTOR(*grad)[i] += VECTOR(data->dS)[i] / S;
	}
	no_of_edges++;
      }
    }
    
    /* Update D, data->dS */
    for (n=0; n<nneis; n++) {
      long int to=VECTOR(data->neis)[n];
      long int x=VECTOR(*data->cats)[to];
      long int y=VECTOR(data->degree)[to];
      
      VECTOR(data->degree)[to] += 1;
      S += MATRIX(data->A_vect, x, y+1);
      S -= MATRIX(data->A_vect, x, y);
      for (i=0; i<dim; i++) {
	igraph_matrix_t *m=VECTOR(data->dA_vects)[i];
	VECTOR(data->dS)[i] += MATRIX(*m, x, y+1);
	VECTOR(data->dS)[i] -= MATRIX(*m, x, y);
      }
    }
    /* New vertex */
    S += MATRIX(data->A_vect, tcat, 0);
    for (i=0; i<dim; i++) {
      igraph_matrix_t *m=VECTOR(data->dA_vects)[i];
      VECTOR(data->dS)[i] += MATRIX(*m, tcat, 0);
    }
    
  }
  
  igraph_vector_update(&data->lastparam, par);
  data->lastf=sum / no_of_edges;
  for (i=0; i<igraph_vector_size(&data->lastgrad); i++) {
    VECTOR(data->lastgrad)[i] /= no_of_edges;
  }

  return 0.0;
}

igraph_real_t igraph_i_revolver_ml_DE_f(const igraph_vector_t *par,
					const igraph_vector_t *garbage,
					void *extra) {

  igraph_i_revolver_ml_DE_data_t *data=extra;
  
  if (!igraph_vector_all_e(par, &data->lastparam)) {
    igraph_i_revolver_ml_DE_eval(par, data);
  }

  if (!igraph_finite(data->lastf)) {
    IGRAPH_WARNING("Target function evaluated to non-finite value.");
  }
  
  /* printf("eval ("); */
  /* for (i=0; i<igraph_vector_size(par); i++) { */
  /*   printf("%f ", VECTOR(*par)[i]); */
  /* } */
  /* printf(" ): "); */
  /* printf("%g\n", data->lastf); */
  return data->lastf;
}

void igraph_i_revolver_ml_DE_df(const igraph_vector_t *par,
				const igraph_vector_t *garbage,
				igraph_vector_t *res, void *extra) {

  igraph_i_revolver_ml_DE_data_t *data=extra;
  
  if (!igraph_vector_all_e(par, &data->lastparam)) {
    igraph_i_revolver_ml_DE_eval(par, data);
  }
  
  igraph_vector_update(res, &data->lastgrad);
  /* printf("derivative ("); */
  /* for (i=0; i<igraph_vector_size(par); i++) { */
  /*   printf("%f ", VECTOR(*par)[i]); */
  /* } */
  /* printf(" ): "); */
  /* for (i=0; i<igraph_vector_size(res); i++) { */
  /*   printf("%f ", VECTOR(*res)[i]); */
  /* } */
  /* printf("\n"); */
}

void igraph_i_revolver_ml_DE_free(igraph_vector_ptr_t *ptr) {
  long int i, n=igraph_vector_ptr_size(ptr);
  for (i=0; i<n; i++) {
    igraph_matrix_t *v=VECTOR(*ptr)[i];
    if (v) {
      igraph_matrix_destroy(v);
      igraph_free(v);
    }
    VECTOR(*ptr)[i]=0;
  }  
}

int igraph_revolver_ml_DE(const igraph_t *graph,
			  const igraph_vector_t *cats,
			  igraph_vector_t *res,
			  igraph_real_t *Fmin,
			  igraph_real_t abstol, igraph_real_t reltol, int maxit,
			  igraph_scalar_function_t *A_fun,
			  igraph_vector_function_t *dA_fun,
			  const igraph_vector_t *filter,
			  igraph_integer_t *fncount, 
			  igraph_integer_t *grcount,
			  igraph_vector_t *lastderiv) {
  
  igraph_i_revolver_ml_DE_data_t info;
  igraph_integer_t maxdegree;
  long int no_of_nodes=igraph_vcount(graph);
  int dim=igraph_vector_size(res);
  int ret, i;

  if (igraph_vector_size(cats) != no_of_nodes) {
    IGRAPH_ERROR("DE ML Revolver failed, invalid category vector size", 
		 IGRAPH_EINVAL);
  }
  
  IGRAPH_CHECK(igraph_maxdegree(graph, &maxdegree, igraph_vss_all(),
				IGRAPH_IN, IGRAPH_LOOPS));
  
  info.A=A_fun;
  info.dA=dA_fun;
  info.graph=graph;
  info.cats=cats;
  info.nocats=igraph_vector_max(cats)+1;
  info.no_of_nodes=no_of_nodes;
  IGRAPH_MATRIX_INIT_FINALLY(&info.A_vect, info.nocats, maxdegree+1);  
  IGRAPH_VECTOR_PTR_INIT_FINALLY(&info.dA_vects, dim);
  IGRAPH_FINALLY(igraph_i_revolver_ml_DE_free, &info.dA_vects);
  for (i=0; i<dim; i++) {
    igraph_matrix_t *m=igraph_Calloc(1, igraph_matrix_t);
    if (!m) { IGRAPH_ERROR("Cannot perform ML D revolver", IGRAPH_ENOMEM); }
    IGRAPH_CHECK(igraph_matrix_init(m, info.nocats, maxdegree+1));
    VECTOR(info.dA_vects)[i]=m;
  }
  info.maxdegree=maxdegree;
  IGRAPH_CHECK(igraph_vector_long_init(&info.degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &info.degree);
  IGRAPH_VECTOR_INIT_FINALLY(&info.neis, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&info.dS, dim);
  IGRAPH_VECTOR_INIT_FINALLY(&info.par1, dim);
  IGRAPH_VECTOR_INIT_FINALLY(&info.tmpgrad, dim);
  IGRAPH_VECTOR_INIT_FINALLY(&info.lastparam, dim);
  info.lastf=0.0;
  IGRAPH_VECTOR_INIT_FINALLY(&info.lastgrad, dim);  
  info.filter=filter;
  
  igraph_i_revolver_ml_DE_eval(res, &info);
  ret=igraph_bfgs(res, Fmin, igraph_i_revolver_ml_DE_f,
		  igraph_i_revolver_ml_DE_df, maxit, 1, abstol, reltol, 1, 
		  &info, fncount, grcount);

  if (lastderiv) {
    igraph_vector_update(lastderiv, &info.lastgrad);
  }
  
  igraph_vector_destroy(&info.lastgrad);
  igraph_vector_destroy(&info.lastparam);
  igraph_vector_destroy(&info.tmpgrad);
  igraph_vector_destroy(&info.par1);
  igraph_vector_destroy(&info.dS);
  igraph_vector_destroy(&info.neis);
  igraph_vector_long_destroy(&info.degree);
  igraph_i_revolver_ml_DE_free(&info.dA_vects);
  igraph_vector_ptr_destroy(&info.dA_vects);
  igraph_matrix_destroy(&info.A_vect);
  IGRAPH_FINALLY_CLEAN(10);

  return ret;
}
 
igraph_real_t igraph_i_revolver_ml_DE_alpha_a_f(const igraph_vector_t *var,
						const igraph_vector_t *par,
						void *extra) {
  long int cat=VECTOR(*var)[0];
  igraph_real_t deg=VECTOR(*var)[1];
  igraph_real_t alpha=VECTOR(*par)[0];
  igraph_real_t a=VECTOR(*par)[1];
  igraph_real_t c= cat==0 ? 1.0 : VECTOR(*par)[1+cat];
  if (deg != 0) {
    return c * (pow(deg, alpha) + a);
  } else {
    return c * a;
  }
}

void igraph_i_revolver_ml_DE_alpha_a_df(const igraph_vector_t *var,
					const igraph_vector_t *par,
					igraph_vector_t *res,
					void *extra) {
  long int cat=VECTOR(*var)[0];
  igraph_real_t deg=VECTOR(*var)[1];
  igraph_real_t alpha=VECTOR(*par)[0];
  igraph_real_t a=VECTOR(*par)[1];
  igraph_real_t c= cat==0 ? 1.0 : VECTOR(*par)[1+cat];
  igraph_vector_null(res);
  if (deg != 0) {
    igraph_real_t p=pow(deg, alpha);
    VECTOR(*res)[0] = c * log(deg) * p;
    VECTOR(*res)[1] = c;
    VECTOR(*res)[1+cat] = p+a;
  } else {
    VECTOR(*res)[0] = 0.0;
    VECTOR(*res)[1] = c;
    VECTOR(*res)[1+cat] = a;
  }
}

int igraph_revolver_ml_DE_alpha_a(const igraph_t *graph,
				  const igraph_vector_t *cats,
				  igraph_real_t *alpha, igraph_real_t *a,
				  igraph_vector_t *coeffs,
				  igraph_real_t *Fmin,
				  igraph_real_t abstol, igraph_real_t reltol,
				  int maxit, const igraph_vector_t *filter,
				  igraph_integer_t *fncount,
				  igraph_integer_t *grcount) {
  igraph_vector_t res;
  int ret, i;  
  
  IGRAPH_VECTOR_INIT_FINALLY(&res, igraph_vector_size(coeffs)+2);
  VECTOR(res)[0] = *alpha;
  VECTOR(res)[1] = *a;
  for (i=0; i<igraph_vector_size(coeffs); i++) {
    VECTOR(res)[i+2] = VECTOR(*coeffs)[i];
  }
  
  ret=igraph_revolver_ml_DE(graph, cats, &res, Fmin, abstol, reltol, maxit,
			    igraph_i_revolver_ml_DE_alpha_a_f,
			    igraph_i_revolver_ml_DE_alpha_a_df,
			    filter, fncount, grcount, 0);
  
  *alpha=VECTOR(res)[0];
  *a=VECTOR(res)[1];
  for (i=0; i<igraph_vector_size(coeffs); i++) {
    VECTOR(*coeffs)[i]=VECTOR(res)[i+2];
  }
  igraph_vector_destroy(&res);
  IGRAPH_FINALLY_CLEAN(1);

  return ret;
}

/*------------------------------------------------------------------*/

typedef struct igraph_i_revolver_ml_AD_data_t {
  igraph_scalar_function_t *A;
  igraph_vector_function_t *dA;
  const igraph_t *graph;
  long int no_of_nodes;
  igraph_matrix_t A_vect;	/* Temporary storage */
  igraph_vector_ptr_t dA_vects;	/* Temporary storage */
  igraph_matrix_bool_t A_valid;
  igraph_integer_t maxdegree;
  igraph_vector_long_t degree;
  igraph_vector_t neis;
  igraph_vector_t dS;		/* Temporary storage */
  igraph_vector_t par1;		/* More tmp storage */
  igraph_vector_t tmpgrad;      /* More... */
  int agebins;

  igraph_vector_t lastparam;	/* The parameter values used last time */
  igraph_real_t lastf;		/* The evaluated function value  */
  igraph_vector_t lastgrad;	/* The evaluated gradient */
  
  const igraph_vector_t *filter;
} igraph_i_revolver_ml_AD_data_t;

#define CHECK_VALID(X,Y) \
  if (!MATRIX(data->A_valid,(X),(Y))) { \
     int i; \
     VECTOR(data->par1)[0]=(X); VECTOR(data->par1)[1]=(Y); \
     MATRIX(data->A_vect, (X), (Y)) = data->A(&data->par1, par, 0); \
     data->dA(&data->par1, 0, &data->tmpgrad, 0); \
     for (i=0; i<dim; i++) { \
       igraph_matrix_t *m=VECTOR(data->dA_vects)[i]; \
       MATRIX(*m, (X), (Y)) = VECTOR(data->tmpgrad)[i]; \
     } \
  }

int igraph_i_revolver_ml_AD_eval(const igraph_vector_t *par,
				igraph_i_revolver_ml_AD_data_t *data) {
  igraph_real_t sum=0.0;
  long int t, i, j;
  int dim=igraph_vector_size(par);
  igraph_vector_t *grad=&data->lastgrad;
  igraph_real_t S=0.0;
  long int agebins=data->agebins;
  long int binwidth=data->no_of_nodes/agebins+1;
  long int no_of_edges=0;
  
  /* Init */
  igraph_vector_long_null(&data->degree);
  igraph_vector_null(&data->dS);
  igraph_vector_null(grad);
  igraph_matrix_bool_null(&data->A_valid);

  for (i=0; i<data->maxdegree+1; i++) {
    for (j=0; j<agebins; j++) {
      long int k;
      VECTOR(data->par1)[0]=(i); VECTOR(data->par1)[1]=(j);
      MATRIX(data->A_vect, (i), (j)) = data->A(&data->par1, par, 0);
      data->dA(&data->par1, par, &data->tmpgrad, 0);
      for (k=0; k<dim; k++) {
	igraph_matrix_t *m=VECTOR(data->dA_vects)[k];
	MATRIX(*m, (i), (j)) = VECTOR(data->tmpgrad)[k];
      }
    }
  }      

  for (t=0; t<data->no_of_nodes; t++) {
    long int n, nneis;

    IGRAPH_ALLOW_INTERRUPTION();
    
    IGRAPH_CHECK(igraph_neighbors(data->graph, &data->neis, t, IGRAPH_OUT));
    nneis=igraph_vector_size(&data->neis);

    if (! data->filter || VECTOR(*data->filter)[t]) {

      /* Update sum(s) */
      for (n=0; n<nneis; n++) {
	long int to=VECTOR(data->neis)[n];
	long int x=VECTOR(data->degree)[to];
	long int y=(t-to)/binwidth;
	
/* 	CHECK_VALID(x,y); */
	sum -= log( MATRIX(data->A_vect, x, y) );
	sum += log( S );
	for (i=0; i<dim; i++) {
	  igraph_matrix_t *m=VECTOR(data->dA_vects)[i];
	  VECTOR(*grad)[i] -= MATRIX(*m, x, y) / MATRIX(data->A_vect, x, y);
	  VECTOR(*grad)[i] += VECTOR(data->dS)[i] / S;
	}
	no_of_edges++;
      }
    }
    
    /* Update S, data->dS */
    for (n=0; n<nneis; n++) {
      long int to=VECTOR(data->neis)[n];
      long int x=VECTOR(data->degree)[to];
      long int y=(t-to)/binwidth;
      
/*       CHECK_VALID(x+1,y);	/\* (x,y) already checked *\/ */
      VECTOR(data->degree)[to] += 1;
      S += MATRIX(data->A_vect, x+1, y);
      S -= MATRIX(data->A_vect, x, y);
      for (i=0; i<dim; i++) {
	igraph_matrix_t *m=VECTOR(data->dA_vects)[i];
	VECTOR(data->dS)[i] += MATRIX(*m, x+1, y);
	VECTOR(data->dS)[i] -= MATRIX(*m, x, y);
      }
    }
    /* New vertex */
/*     CHECK_VALID(0,0); */
    S += MATRIX(data->A_vect, 0, 0);
    for (i=0; i<dim; i++) {
      igraph_matrix_t *m=VECTOR(data->dA_vects)[i];
      VECTOR(data->dS)[i] += MATRIX(*m, 0, 0);
    }
    /* Aging */
    for (j=1; t-binwidth*j+1>=0; j++) {
      long int shnode=t-binwidth*j+1;
      long int deg=VECTOR(data->degree)[shnode];
/*       CHECK_VALID(deg, j-1); */
/*       CHECK_VALID(deg, j); */
      S += MATRIX(data->A_vect, deg, j);
      S -= MATRIX(data->A_vect, deg, j-1);
      for (i=0; i<dim; i++) {
	igraph_matrix_t *m=VECTOR(data->dA_vects)[i];
	VECTOR(data->dS)[i] += MATRIX(*m, deg, j);
	VECTOR(data->dS)[i] -= MATRIX(*m, deg, j-1);
      }
    }
    
  }      
      
  igraph_vector_update(&data->lastparam, par);
  data->lastf=sum / no_of_edges;
  for (i=0; i<igraph_vector_size(&data->lastgrad); i++) {
    VECTOR(data->lastgrad)[i] /= no_of_edges;
  }

  return 0.0;
}

igraph_real_t igraph_i_revolver_ml_AD_f(const igraph_vector_t *par,
					const igraph_vector_t *garbage,
					void *extra) {

  igraph_i_revolver_ml_AD_data_t *data=extra;
  
  if (!igraph_vector_all_e(par, &data->lastparam)) {
    igraph_i_revolver_ml_AD_eval(par, data);
  }

  if (!igraph_finite(data->lastf)) {
    IGRAPH_WARNING("Target function evaluated to non-finite value.");
  }
  
  /* printf("eval ("); */
  /* for (i=0; i<igraph_vector_size(par); i++) { */
  /*   printf("%f ", VECTOR(*par)[i]); */
  /* } */
  /* printf(" ): "); */
  /* printf("%g\n", data->lastf); */
  return data->lastf;
}

void igraph_i_revolver_ml_AD_df(const igraph_vector_t *par,
				const igraph_vector_t *garbage,
				igraph_vector_t *res, void *extra) {

  igraph_i_revolver_ml_AD_data_t *data=extra;
  
  if (!igraph_vector_all_e(par, &data->lastparam)) {
    igraph_i_revolver_ml_AD_eval(par, data);
  }
  
  igraph_vector_update(res, &data->lastgrad);
  /* printf("derivative ("); */
  /* for (i=0; i<igraph_vector_size(par); i++) { */
  /*   printf("%f ", VECTOR(*par)[i]); */
  /* } */
  /* printf(" ): "); */
  /* for (i=0; i<igraph_vector_size(res); i++) { */
  /*   printf("%f ", VECTOR(*res)[i]); */
  /* } */
  /* printf("\n"); */
}

void igraph_i_revolver_ml_AD_free(igraph_vector_ptr_t *ptr) {
  long int i, n=igraph_vector_ptr_size(ptr);
  for (i=0; i<n; i++) {
    igraph_matrix_t *v=VECTOR(*ptr)[i];
    if (v) {
      igraph_matrix_destroy(v);
      igraph_free(v);
    }
    VECTOR(*ptr)[i]=0;
  }  
}

int igraph_revolver_ml_AD(const igraph_t *graph,
			  igraph_vector_t *res,
			  igraph_real_t *Fmin,
			  igraph_real_t abstol, igraph_real_t reltol, int maxit,
			  igraph_scalar_function_t *A_fun,
			  igraph_vector_function_t *dA_fun,
			  int agebins, const igraph_vector_t *filter,
			  igraph_integer_t *fncount, 
			  igraph_integer_t *grcount,
			  igraph_vector_t *lastderiv) {
  
  igraph_i_revolver_ml_AD_data_t info;
  igraph_integer_t maxdegree;
  long int no_of_nodes=igraph_vcount(graph);
  int dim=igraph_vector_size(res);
  int ret, i;
  
  IGRAPH_CHECK(igraph_maxdegree(graph, &maxdegree, igraph_vss_all(),
				IGRAPH_IN, IGRAPH_LOOPS));
  
  info.A=A_fun;
  info.dA=dA_fun;
  info.graph=graph;
  info.no_of_nodes=no_of_nodes;
  IGRAPH_MATRIX_INIT_FINALLY(&info.A_vect, maxdegree+1, agebins);  
  IGRAPH_VECTOR_PTR_INIT_FINALLY(&info.dA_vects, dim);
  IGRAPH_FINALLY(igraph_i_revolver_ml_AD_free, &info.dA_vects);
  for (i=0; i<dim; i++) {
    igraph_matrix_t *m=igraph_Calloc(1, igraph_matrix_t);
    if (!m) { IGRAPH_ERROR("Cannot perform ML D revolver", IGRAPH_ENOMEM); }
    IGRAPH_CHECK(igraph_matrix_init(m, maxdegree+1, agebins));
    VECTOR(info.dA_vects)[i]=m;
  }
  IGRAPH_CHECK(igraph_matrix_bool_init(&info.A_valid, maxdegree+1, agebins));
  IGRAPH_FINALLY(igraph_matrix_bool_destroy, &info.A_valid);
  info.maxdegree=maxdegree;
  IGRAPH_CHECK(igraph_vector_long_init(&info.degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &info.degree);
  IGRAPH_VECTOR_INIT_FINALLY(&info.neis, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&info.dS, dim);
  IGRAPH_VECTOR_INIT_FINALLY(&info.par1, dim);
  IGRAPH_VECTOR_INIT_FINALLY(&info.tmpgrad, dim);
  info.agebins=agebins;
  IGRAPH_VECTOR_INIT_FINALLY(&info.lastparam, dim);
  info.lastf=0.0;
  IGRAPH_VECTOR_INIT_FINALLY(&info.lastgrad, dim);  
  info.filter=filter;
  
  igraph_i_revolver_ml_AD_eval(res, &info);
  ret=igraph_bfgs(res, Fmin, igraph_i_revolver_ml_AD_f,
		  igraph_i_revolver_ml_AD_df, maxit, 1, abstol, reltol, 1, 
		  &info, fncount, grcount);

  if (lastderiv) {
    igraph_vector_update(lastderiv, &info.lastgrad);
  }
  
  igraph_vector_destroy(&info.lastgrad);
  igraph_vector_destroy(&info.lastparam);
  igraph_vector_destroy(&info.tmpgrad);
  igraph_vector_destroy(&info.par1);
  igraph_vector_destroy(&info.dS);
  igraph_vector_destroy(&info.neis);
  igraph_vector_long_destroy(&info.degree);
  igraph_matrix_bool_destroy(&info.A_valid);
  igraph_i_revolver_ml_AD_free(&info.dA_vects);
  igraph_vector_ptr_destroy(&info.dA_vects);
  igraph_matrix_destroy(&info.A_vect);
  IGRAPH_FINALLY_CLEAN(11);

  return ret;
}

igraph_real_t igraph_i_revolver_ml_AD_alpha_a_beta_f(const igraph_vector_t *var,
						     const igraph_vector_t *par,
						     void *extra) {
  igraph_real_t deg=VECTOR(*var)[0];
  igraph_real_t age=VECTOR(*var)[1]+1;
  igraph_real_t alpha=VECTOR(*par)[0];
  igraph_real_t a=VECTOR(*par)[1];
  igraph_real_t beta=VECTOR(*par)[2];
  return (pow(deg, alpha) + a) * pow(age, -beta);
}

void igraph_i_revolver_ml_AD_alpha_a_beta_df(const igraph_vector_t *var,
					     const igraph_vector_t *par,
					     igraph_vector_t *res,
					     void *extra) {
  igraph_real_t deg=VECTOR(*var)[0];
  igraph_real_t age=VECTOR(*var)[1]+1;
  igraph_real_t alpha=VECTOR(*par)[0];
  igraph_real_t a=VECTOR(*par)[1];
  igraph_real_t beta=VECTOR(*par)[2];
  igraph_real_t p1=pow(deg, alpha);
  igraph_real_t p2=pow(age, -beta);
  VECTOR(*res)[0]= deg == 0 ? 0.0 : p2*log(deg)*p1;
  VECTOR(*res)[1]= p2;
  VECTOR(*res)[2]= -(p1+a)*log(age)*p2;
}

int igraph_revolver_ml_AD_alpha_a_beta(const igraph_t *graph,
				       igraph_real_t *alpha, igraph_real_t *a,
				       igraph_real_t *beta, igraph_real_t *Fmin,
				       igraph_real_t abstol, igraph_real_t reltol,
				       int maxit, int agebins, 
				       const igraph_vector_t *filter,
				       igraph_integer_t *fncount,
				       igraph_integer_t *grcount) {
  igraph_vector_t res;
  int ret;
  
  IGRAPH_VECTOR_INIT_FINALLY(&res, 3);
  VECTOR(res)[0] = *alpha;
  VECTOR(res)[1] = *a;
  VECTOR(res)[2] = *beta;
  
  ret=igraph_revolver_ml_AD(graph, &res, Fmin, abstol, reltol, maxit,
			    igraph_i_revolver_ml_AD_alpha_a_beta_f,
			    igraph_i_revolver_ml_AD_alpha_a_beta_df,
			    agebins, filter, fncount, grcount, 0);
  
  *alpha=VECTOR(res)[0];
  *a=VECTOR(res)[1];
  *beta=VECTOR(res)[2];

  igraph_vector_destroy(&res);
  IGRAPH_FINALLY_CLEAN(1);  
  return ret;
}

igraph_real_t igraph_i_revolver_ml_AD_dpareto_f(const igraph_vector_t *var,
						const igraph_vector_t *par,
						void *extra) {
  igraph_real_t deg=VECTOR(*var)[0];
  igraph_real_t age=VECTOR(*var)[1]+1;
  igraph_real_t alpha=VECTOR(*par)[0];
  igraph_real_t a=VECTOR(*par)[1];
  igraph_real_t paralpha=VECTOR(*par)[2];
  igraph_real_t parbeta=VECTOR(*par)[3];
  igraph_real_t parscale=VECTOR(*par)[4];
  
  igraph_real_t res= age < parscale ? 
    (pow(deg,alpha)+a) * pow(age/parscale, parbeta-1) : 
    (pow(deg,alpha)+a) * pow(age/parscale, -paralpha-1);
    
/*   printf("eval at %f %f, %f %f %f %f %f: %f\n", deg, age, */
/* 	 alpha, a, paralpha, parbeta, parscale, res); */

  return res;
}

void igraph_i_revolver_ml_AD_dpareto_df(const igraph_vector_t *var,
					const igraph_vector_t *par,
					igraph_vector_t *res,
					void *extra) {
  igraph_real_t deg=VECTOR(*var)[0];
  igraph_real_t age=VECTOR(*var)[1]+1;
  igraph_real_t alpha=VECTOR(*par)[0];
  igraph_real_t a=VECTOR(*par)[1];
  igraph_real_t paralpha=VECTOR(*par)[2];
  igraph_real_t parbeta=VECTOR(*par)[3];
  igraph_real_t parscale=VECTOR(*par)[4];
  igraph_real_t exponent= age < parscale ? parbeta : -paralpha;
  igraph_real_t p1=pow(deg, alpha);  
  igraph_real_t p2=pow(age/parscale, exponent-1);
  VECTOR(*res)[0]= deg == 0 ? 0.0 : log(deg)*p1*p2;
  VECTOR(*res)[1]=p2;
  VECTOR(*res)[2]= age > parscale ? (p1+a)*log(age/parscale)*p2 : 0;
  VECTOR(*res)[3]= age < parscale ? (p1+a)*log(age/parscale)*p2 : 0;
  VECTOR(*res)[4]=-(p1+a)*(exponent-1)*pow(age/parscale, exponent-2)*
    age/parscale/parscale;
/*   printf("deriv at %f %f, %f %f %f %f %f: %f %f %f %f %f\n", deg, age, */
/* 	 alpha, a, paralpha, parbeta, parscale, */
/* 	 VECTOR(*res)[0], VECTOR(*res)[1], VECTOR(*res)[2], VECTOR(*res)[3], */
/* 	 VECTOR(*res)[4]); */
    
}

int igraph_revolver_ml_AD_dpareto_eval(const igraph_t *graph,
				       igraph_real_t alpha, igraph_real_t a,
				       igraph_real_t paralpha, 
				       igraph_real_t parbeta,
				       igraph_real_t parscale,
				       igraph_real_t *value,
				       igraph_vector_t *deriv,
				       int agebins,
				       const igraph_vector_t *filter) {
  igraph_vector_t res;
  int ret;
  igraph_integer_t fncount, grcount;
  
  IGRAPH_VECTOR_INIT_FINALLY(&res, 5);
  VECTOR(res)[0] = alpha;
  VECTOR(res)[1] = a;
  VECTOR(res)[2] = paralpha;
  VECTOR(res)[3] = parbeta;
  VECTOR(res)[4] = parscale;
  
  ret=igraph_revolver_ml_AD(graph, &res, value, 0, 0, 0,
			    igraph_i_revolver_ml_AD_dpareto_f,
			    igraph_i_revolver_ml_AD_dpareto_df,
			    agebins, filter, &fncount, &grcount, deriv);
  
  igraph_vector_destroy(&res);
  IGRAPH_FINALLY_CLEAN(1);
  return ret;
}

int igraph_revolver_ml_AD_dpareto(const igraph_t *graph,
				  igraph_real_t *alpha, igraph_real_t *a,
				  igraph_real_t *paralpha, igraph_real_t *parbeta,
				  igraph_real_t *parscale,
				  igraph_real_t *Fmin,
				  igraph_real_t abstol, igraph_real_t reltol,
				  int maxit, int agebins, 
				  const igraph_vector_t *filter,
				  igraph_integer_t *fncount,
				  igraph_integer_t *grcount) {
  igraph_vector_t res;
  int ret;
  
  IGRAPH_VECTOR_INIT_FINALLY(&res, 5);
  VECTOR(res)[0] = *alpha;
  VECTOR(res)[1] = *a;
  VECTOR(res)[2] = *paralpha;
  VECTOR(res)[3] = *parbeta;
  VECTOR(res)[4] = *parscale;
  
  ret=igraph_revolver_ml_AD(graph, &res, Fmin, abstol, reltol, maxit,
			    igraph_i_revolver_ml_AD_dpareto_f,
			    igraph_i_revolver_ml_AD_dpareto_df,
			    agebins, filter, fncount, grcount, 0);
  
  *alpha=VECTOR(res)[0];
  *a=VECTOR(res)[1];
  *paralpha=VECTOR(res)[2];
  *parbeta=VECTOR(res)[3];
  *parscale=VECTOR(res)[4];

  igraph_vector_destroy(&res);
  IGRAPH_FINALLY_CLEAN(1);
  return ret;
}

/*------------------------------------------------------------------*/

typedef struct igraph_i_revolver_ml_ADE_data_t {
  igraph_scalar_function_t *A;
  igraph_vector_function_t *dA;
  const igraph_t *graph;
  const igraph_vector_t *cats;
  long int no_of_nodes;
  igraph_array3_t A_vect;	/* Temporary storage */
  igraph_vector_ptr_t dA_vects;	/* Temporary storage */
  igraph_integer_t maxdegree;
  igraph_integer_t nocats;
  igraph_vector_long_t degree;
  igraph_vector_t neis;
  igraph_vector_t dS;		/* Temporary storage */
  igraph_vector_t par1;		/* More tmp storage */
  igraph_vector_t tmpgrad;      /* More... */
  int agebins;

  igraph_vector_t lastparam;	/* The parameter values used last time */
  igraph_real_t lastf;		/* The evaluated function value  */
  igraph_vector_t lastgrad;	/* The evaluated gradient */
  
  const igraph_vector_t *filter;
} igraph_i_revolver_ml_ADE_data_t;

int igraph_i_revolver_ml_ADE_eval(const igraph_vector_t *par,
				  igraph_i_revolver_ml_ADE_data_t *data) {
  igraph_real_t sum=0.;
  long int t, i, j, c;
  int dim=igraph_vector_size(par);
  igraph_vector_t *grad=&data->lastgrad;
  igraph_real_t S=0.0;
  long int agebins=data->agebins;
  long int binwidth=data->no_of_nodes/agebins+1;
  long int no_of_edges=0;
  
  /* Init */
  igraph_vector_long_null(&data->degree);
  igraph_vector_null(&data->dS);
  igraph_vector_null(grad);

  for (c=0; c<data->nocats; c++) {
    for (i=0; i<data->maxdegree+1; i++) {
      for (j=0; j<agebins; j++) {
	long int k;
	VECTOR(data->par1)[0]=c;
	VECTOR(data->par1)[1]=i;
	VECTOR(data->par1)[2]=j;
	ARRAY3(data->A_vect, c, (i), (j)) = data->A(&data->par1, par, 0);
	data->dA(&data->par1, par, &data->tmpgrad, 0);
	for (k=0; k<dim; k++) {
	  igraph_array3_t *m=VECTOR(data->dA_vects)[k];
	  ARRAY3(*m, c, i, j) = VECTOR(data->tmpgrad)[k];
	}
      }
    }
  }

  for (t=0; t<data->no_of_nodes; t++) {
    long int n, nneis, shnode;
    long int tcat=VECTOR(*data->cats)[t];

    IGRAPH_ALLOW_INTERRUPTION();
    
    IGRAPH_CHECK(igraph_neighbors(data->graph, &data->neis, t, IGRAPH_OUT));
    nneis=igraph_vector_size(&data->neis);

    if (! data->filter || VECTOR(*data->filter)[t]) {

      /* Update sum(s) */
      for (n=0; n<nneis; n++) {
	long int to=VECTOR(data->neis)[n];
	long int x=VECTOR(*data->cats)[to];
	long int y=VECTOR(data->degree)[to];
	long int z=(t-to)/binwidth;
	
/* 	CHECK_VALID(x,y); */
	sum -= log( ARRAY3(data->A_vect, x, y, z) );
	sum += log( S );
	for (i=0; i<dim; i++) {
	  igraph_array3_t *m=VECTOR(data->dA_vects)[i];
	  VECTOR(*grad)[i] -= ARRAY3(*m, x, y, z) / ARRAY3(data->A_vect, x, y, z);
	  VECTOR(*grad)[i] += VECTOR(data->dS)[i] / S;
	}
	no_of_edges++;
      }
    }
    
    /* Update S, data->dS */
    for (n=0; n<nneis; n++) {
      long int to=VECTOR(data->neis)[n];
      long int x=VECTOR(*data->cats)[to];
      long int y=VECTOR(data->degree)[to];
      long int z=(t-to)/binwidth;
      
/*       CHECK_VALID(x+1,y);	/\* (x,y) already checked *\/ */
      VECTOR(data->degree)[to] += 1;
      S += ARRAY3(data->A_vect, x, y+1, z);
      S -= ARRAY3(data->A_vect, x, y, z);
      for (i=0; i<dim; i++) {
	igraph_array3_t *m=VECTOR(data->dA_vects)[i];
	VECTOR(data->dS)[i] += ARRAY3(*m, x, y+1, z);
	VECTOR(data->dS)[i] -= ARRAY3(*m, x, y, z);
      }
    }
    /* New vertex */
/*     CHECK_VALID(0,0); */
    S += ARRAY3(data->A_vect, tcat, 0, 0);
    for (i=0; i<dim; i++) {
      igraph_array3_t *m=VECTOR(data->dA_vects)[i];
      VECTOR(data->dS)[i] += ARRAY3(*m, tcat, 0, 0);
    }
    /* Aging */
    for (j=1, shnode=t-binwidth+1; shnode>=0; j++, shnode-=binwidth) {
      long int cat=VECTOR(*data->cats)[shnode];
      long int deg=VECTOR(data->degree)[shnode];
/*       CHECK_VALID(deg, j-1); */
/*       CHECK_VALID(deg, j); */
      S += ARRAY3(data->A_vect, cat, deg, j);
      S -= ARRAY3(data->A_vect, cat, deg, j-1);
      for (i=0; i<dim; i++) {
	igraph_array3_t *m=VECTOR(data->dA_vects)[i];
	VECTOR(data->dS)[i] += ARRAY3(*m, cat, deg, j);
	VECTOR(data->dS)[i] -= ARRAY3(*m, cat, deg, j-1);
      }
    }
    
  }      
      
  igraph_vector_update(&data->lastparam, par);
  data->lastf=sum / no_of_edges;
  for (i=0; i<igraph_vector_size(&data->lastgrad); i++) {
    VECTOR(data->lastgrad)[i] /= no_of_edges;
  }

  return 0.0;
}

igraph_real_t igraph_i_revolver_ml_ADE_f(const igraph_vector_t *par,
					 const igraph_vector_t *garbage,
					 void *extra) {

  igraph_i_revolver_ml_ADE_data_t *data=extra;
  
  if (!igraph_vector_all_e(par, &data->lastparam)) {
    igraph_i_revolver_ml_ADE_eval(par, data);
  }

  if (!igraph_finite(data->lastf)) {
    IGRAPH_WARNING("Target function evaluated to non-finite value.");
  }
  
  /* printf("eval ("); */
  /* for (i=0; i<igraph_vector_size(par); i++) { */
  /*   printf("%f ", VECTOR(*par)[i]); */
  /* } */
  /* printf(" ): "); */
  /* printf("%g\n", data->lastf); */
  return data->lastf;
}

void igraph_i_revolver_ml_ADE_df(const igraph_vector_t *par,
				 const igraph_vector_t *garbage,
				 igraph_vector_t *res, void *extra) {

  igraph_i_revolver_ml_ADE_data_t *data=extra;
  
  if (!igraph_vector_all_e(par, &data->lastparam)) {
    igraph_i_revolver_ml_ADE_eval(par, data);
  }
  
  igraph_vector_update(res, &data->lastgrad);
  /* printf("derivative ("); */
  /* for (i=0; i<igraph_vector_size(par); i++) { */
  /*   printf("%f ", VECTOR(*par)[i]); */
  /* } */
  /* printf(" ): "); */
  /* for (i=0; i<igraph_vector_size(res); i++) { */
  /*   printf("%f ", VECTOR(*res)[i]); */
  /* } */
  /* printf("\n"); */
}

void igraph_i_revolver_ml_ADE_free(igraph_vector_ptr_t *ptr) {
  long int i, n=igraph_vector_ptr_size(ptr);
  for (i=0; i<n; i++) {
    igraph_array3_t *v=VECTOR(*ptr)[i];
    if (v) {
      igraph_array3_destroy(v);
      igraph_free(v);
    }
    VECTOR(*ptr)[i]=0;
  }  
}

int igraph_revolver_ml_ADE(const igraph_t *graph,
			   const igraph_vector_t *cats,
			   igraph_vector_t *res,
			   igraph_real_t *Fmin,
			   igraph_real_t abstol, igraph_real_t reltol, int maxit,
			   igraph_scalar_function_t *A_fun,
			   igraph_vector_function_t *dA_fun,
			   int agebins, const igraph_vector_t *filter,
			   igraph_integer_t *fncount, 
			   igraph_integer_t *grcount,
			   igraph_vector_t *lastderiv) {
  
  igraph_i_revolver_ml_ADE_data_t info;
  igraph_integer_t maxdegree;
  long int no_of_nodes=igraph_vcount(graph);
  int dim=igraph_vector_size(res);
  int ret, i;
  
  if (igraph_vector_size(cats) != no_of_nodes) {
    IGRAPH_ERROR("ADE ML Revolver failed: invalid category vector size", 
		 IGRAPH_ENOMEM);
  }

  IGRAPH_CHECK(igraph_maxdegree(graph, &maxdegree, igraph_vss_all(),
				IGRAPH_IN, IGRAPH_LOOPS));
  
  info.A=A_fun;
  info.dA=dA_fun;
  info.graph=graph;
  info.no_of_nodes=no_of_nodes;
  info.cats=cats;
  info.nocats=igraph_vector_max(cats)+1;
  IGRAPH_ARRAY3_INIT_FINALLY(&info.A_vect, info.nocats, maxdegree+1, agebins);  
  IGRAPH_VECTOR_PTR_INIT_FINALLY(&info.dA_vects, dim);
  IGRAPH_FINALLY(igraph_i_revolver_ml_ADE_free, &info.dA_vects);
  for (i=0; i<dim; i++) {
    igraph_array3_t *m=igraph_Calloc(1, igraph_array3_t);
    if (!m) { IGRAPH_ERROR("Cannot perform ML D revolver", IGRAPH_ENOMEM); }
    IGRAPH_CHECK(igraph_array3_init(m, info.nocats, maxdegree+1, agebins));
    VECTOR(info.dA_vects)[i]=m;
  }
  info.maxdegree=maxdegree;
  IGRAPH_CHECK(igraph_vector_long_init(&info.degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &info.degree);
  IGRAPH_VECTOR_INIT_FINALLY(&info.neis, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&info.dS, dim);
  IGRAPH_VECTOR_INIT_FINALLY(&info.par1, dim);
  IGRAPH_VECTOR_INIT_FINALLY(&info.tmpgrad, dim);
  info.agebins=agebins;
  IGRAPH_VECTOR_INIT_FINALLY(&info.lastparam, dim);
  info.lastf=0.0;
  IGRAPH_VECTOR_INIT_FINALLY(&info.lastgrad, dim);  
  info.filter=filter;
  
  igraph_i_revolver_ml_ADE_eval(res, &info);
  ret=igraph_bfgs(res, Fmin, igraph_i_revolver_ml_ADE_f,
		  igraph_i_revolver_ml_ADE_df, maxit, 1, abstol, reltol, 1, 
		  &info, fncount, grcount);

  if (lastderiv) {
    igraph_vector_update(lastderiv, &info.lastgrad);
  }
  
  igraph_vector_destroy(&info.lastgrad);
  igraph_vector_destroy(&info.lastparam);
  igraph_vector_destroy(&info.tmpgrad);
  igraph_vector_destroy(&info.par1);
  igraph_vector_destroy(&info.dS);
  igraph_vector_destroy(&info.neis);
  igraph_vector_long_destroy(&info.degree);
  igraph_i_revolver_ml_ADE_free(&info.dA_vects);
  igraph_vector_ptr_destroy(&info.dA_vects);
  igraph_array3_destroy(&info.A_vect);
  IGRAPH_FINALLY_CLEAN(10);

  return ret;
}

igraph_real_t igraph_i_revolver_ml_ADE_alpha_a_beta_f(const igraph_vector_t *var,
						      const igraph_vector_t *par,
						      void *extra) {
  long int cat=VECTOR(*var)[0];
  igraph_real_t deg=VECTOR(*var)[1];
  igraph_real_t age=VECTOR(*var)[2]+1;
  igraph_real_t alpha=VECTOR(*par)[0];
  igraph_real_t a=VECTOR(*par)[1];
  igraph_real_t beta=VECTOR(*par)[2];  
  igraph_real_t c= cat==0 ? 1.0 : VECTOR(*par)[2+cat];
  return c * (pow(deg, alpha) + a) * pow(age, -beta);
}

void igraph_i_revolver_ml_ADE_alpha_a_beta_df(const igraph_vector_t *var,
					      const igraph_vector_t *par,
					      igraph_vector_t *res,
					      void *extra) {
  long int cat=VECTOR(*var)[0];
  igraph_real_t deg=VECTOR(*var)[1];
  igraph_real_t age=VECTOR(*var)[2]+1;
  igraph_real_t alpha=VECTOR(*par)[0];
  igraph_real_t a=VECTOR(*par)[1];
  igraph_real_t beta=VECTOR(*par)[2];
  igraph_real_t c= cat==0 ? 1.0 : VECTOR(*par)[2+cat];
  igraph_real_t p1=pow(deg, alpha);
  igraph_real_t p2=pow(age, -beta);
  igraph_vector_null(res);
  VECTOR(*res)[0]= deg == 0 ? 0.0 : c * p2*log(deg)*p1;
  VECTOR(*res)[1]= c * p2;
  VECTOR(*res)[2]= c * -(p1+a)*log(age)*p2;
  VECTOR(*res)[2+cat] = (p1 + a) * p2;
}

int igraph_revolver_ml_ADE_alpha_a_beta(const igraph_t *graph,
					const igraph_vector_t *cats,
					igraph_real_t *alpha, igraph_real_t *a,
					igraph_real_t *beta, igraph_vector_t *coeffs,
					igraph_real_t *Fmin,
					igraph_real_t abstol, igraph_real_t reltol,
					int maxit, int agebins, 
					const igraph_vector_t *filter,
					igraph_integer_t *fncount,
					igraph_integer_t *grcount) {
  igraph_vector_t res;
  int ret, i;
  
  IGRAPH_VECTOR_INIT_FINALLY(&res, 3+igraph_vector_size(coeffs));
  VECTOR(res)[0] = *alpha;
  VECTOR(res)[1] = *a;
  VECTOR(res)[2] = *beta;
  for (i=0; i<igraph_vector_size(coeffs); i++) {
    VECTOR(res)[i+3] = VECTOR(*coeffs)[i];
  }
  
  ret=igraph_revolver_ml_ADE(graph, cats, &res, Fmin, abstol, reltol, maxit,
			     igraph_i_revolver_ml_ADE_alpha_a_beta_f,
			     igraph_i_revolver_ml_ADE_alpha_a_beta_df,
			     agebins, filter, fncount, grcount, 0);
  
  *alpha=VECTOR(res)[0];
  *a=VECTOR(res)[1];
  *beta=VECTOR(res)[2];
  for (i=0; i<igraph_vector_size(coeffs); i++) {
    VECTOR(*coeffs)[i]=VECTOR(res)[i+3];
  }

  igraph_vector_destroy(&res);
  IGRAPH_FINALLY_CLEAN(1);  
  return ret;
}

igraph_real_t igraph_i_revolver_ml_ADE_dpareto_f(const igraph_vector_t *var,
						 const igraph_vector_t *par,
						 void *extra) {
  long int cat=VECTOR(*var)[0];
  igraph_real_t deg=VECTOR(*var)[1];
  igraph_real_t age=VECTOR(*var)[2]+1;
  igraph_real_t alpha=VECTOR(*par)[0];
  igraph_real_t a=VECTOR(*par)[1];
  igraph_real_t paralpha=VECTOR(*par)[2];
  igraph_real_t parbeta=VECTOR(*par)[3];
  igraph_real_t parscale=VECTOR(*par)[4];
  igraph_real_t c= cat==0 ? 1.0 : VECTOR(*par)[4+cat];

  igraph_real_t res= age < parscale ? 
    c * (pow(deg,alpha)+a) * pow(age/parscale, parbeta-1) : 
    c * (pow(deg,alpha)+a) * pow(age/parscale, -paralpha-1);
    
/*   printf("eval at %f %f, %f %f %f %f %f: %f\n", deg, age, */
/* 	 alpha, a, paralpha, parbeta, parscale, res); */

  return res;
}

void igraph_i_revolver_ml_ADE_dpareto_df(const igraph_vector_t *var,
					 const igraph_vector_t *par,
					 igraph_vector_t *res,
					 void *extra) {
  long int cat=VECTOR(*var)[0];
  igraph_real_t deg=VECTOR(*var)[1];
  igraph_real_t age=VECTOR(*var)[2]+1;
  igraph_real_t alpha=VECTOR(*par)[0];
  igraph_real_t a=VECTOR(*par)[1];
  igraph_real_t paralpha=VECTOR(*par)[2];
  igraph_real_t parbeta=VECTOR(*par)[3];
  igraph_real_t parscale=VECTOR(*par)[4];
  igraph_real_t exponent= age < parscale ? parbeta : -paralpha;
  igraph_real_t p1=pow(deg, alpha);  
  igraph_real_t p2=pow(age/parscale, exponent-1);
  igraph_real_t c= cat==0 ? 1.0 : VECTOR(*par)[4+cat];
  igraph_vector_null(res);
  VECTOR(*res)[0]= deg == 0 ? 0.0 : c * log(deg)*p1*p2;
  VECTOR(*res)[1]= c * p2;
  VECTOR(*res)[2]= age > parscale ? c * (p1+a)*log(age/parscale)*p2 : 0;
  VECTOR(*res)[3]= age < parscale ? c * (p1+a)*log(age/parscale)*p2 : 0;
  VECTOR(*res)[4]= c * -(p1+a)*(exponent-1)*pow(age/parscale, exponent-2)*
    age/parscale/parscale;
  VECTOR(*res)[4+cat] = (p1+a) * p2;
/*   printf("deriv at %f %f, %f %f %f %f %f: %f %f %f %f %f\n", deg, age, */
/* 	 alpha, a, paralpha, parbeta, parscale, */
/* 	 VECTOR(*res)[0], VECTOR(*res)[1], VECTOR(*res)[2], VECTOR(*res)[3], */
/* 	 VECTOR(*res)[4]); */
    
}
  
int igraph_revolver_ml_ADE_dpareto_eval(const igraph_t *graph,
					const igraph_vector_t *cats,
					igraph_real_t alpha, igraph_real_t a,
					igraph_real_t paralpha, 
					igraph_real_t parbeta,
					igraph_real_t parscale,
					const igraph_vector_t *coeffs,
					igraph_real_t *value,
					igraph_vector_t *deriv,
					int agebins,
					const igraph_vector_t *filter) {
  igraph_vector_t res;
  int ret, i;
  igraph_integer_t fncount, grcount;
  
  IGRAPH_VECTOR_INIT_FINALLY(&res, 5+igraph_vector_size(coeffs));
  VECTOR(res)[0] = alpha;
  VECTOR(res)[1] = a;
  VECTOR(res)[2] = paralpha;
  VECTOR(res)[3] = parbeta;
  VECTOR(res)[4] = parscale;
  for (i=0; i<igraph_vector_size(coeffs); i++) {
    VECTOR(res)[5+i] = VECTOR(*coeffs)[i];
  }
  
  ret=igraph_revolver_ml_ADE(graph, cats, &res, value, 0, 0, 0,
			     igraph_i_revolver_ml_ADE_dpareto_f,
			     igraph_i_revolver_ml_ADE_dpareto_df,
			     agebins, filter, &fncount, &grcount, deriv);
  
  igraph_vector_destroy(&res);
  IGRAPH_FINALLY_CLEAN(1);
  return ret;
}

int igraph_revolver_ml_ADE_dpareto(const igraph_t *graph,
				   const igraph_vector_t *cats,
				   igraph_real_t *alpha, igraph_real_t *a,
				   igraph_real_t *paralpha, igraph_real_t *parbeta,
				   igraph_real_t *parscale, igraph_vector_t *coeffs,
				   igraph_real_t *Fmin,
				   igraph_real_t abstol, igraph_real_t reltol,
				   int maxit, int agebins, 
				   const igraph_vector_t *filter,
				   igraph_integer_t *fncount,
				   igraph_integer_t *grcount) {
  igraph_vector_t res;
  int ret, i;
  
  IGRAPH_VECTOR_INIT_FINALLY(&res, 5+igraph_vector_size(coeffs));
  VECTOR(res)[0] = *alpha;
  VECTOR(res)[1] = *a;
  VECTOR(res)[2] = *paralpha;
  VECTOR(res)[3] = *parbeta;
  VECTOR(res)[4] = *parscale;
  for (i=0; i<igraph_vector_size(coeffs); i++) {
    VECTOR(res)[5+i] = VECTOR(*coeffs)[i];
  }
  
  ret=igraph_revolver_ml_ADE(graph, cats, &res, Fmin, abstol, reltol, maxit,
			     igraph_i_revolver_ml_ADE_dpareto_f,
			     igraph_i_revolver_ml_ADE_dpareto_df,
			     agebins, filter, fncount, grcount, 0);
  
  *alpha=VECTOR(res)[0];
  *a=VECTOR(res)[1];
  *paralpha=VECTOR(res)[2];
  *parbeta=VECTOR(res)[3];
  *parscale=VECTOR(res)[4];
  for (i=0; i<igraph_vector_size(coeffs); i++) {
    VECTOR(*coeffs)[i] = VECTOR(res)[5+i];
  }

  igraph_vector_destroy(&res);
  IGRAPH_FINALLY_CLEAN(1);
  return ret;
}  

void igraph_i_revolver_ml_ADE_dpareto_evalf_free(igraph_vector_ptr_t *p) {
  long int i, n=igraph_vector_ptr_size(p);
  for (i=0; i<n; i++) {
    igraph_array3_t *A=VECTOR(*p)[i];
    if (A) { 
      igraph_array3_destroy(A);
      igraph_Free(A);
      VECTOR(*p)[i] = 0;
    }
  }
}

int igraph_revolver_ml_ADE_dpareto_evalf(const igraph_t *graph,
					 const igraph_vector_t *cats,
					 const igraph_matrix_t *par,
					 igraph_vector_t *value,
					 int agebins, 
					 const igraph_vector_t *filter) {
  
  igraph_vector_t S;
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth=no_of_nodes/agebins+1;
  long int no_of_edges=0;
  igraph_vector_long_t degree;
  igraph_vector_ptr_t A_ptr;
  igraph_vector_t neis;
  igraph_integer_t pmaxdegree;
  long int maxdegree;
  long int nocats=igraph_vector_max(cats)+1;
  long int nopar=igraph_matrix_nrow(par);
  long int c, i, j, t;

  if (filter && igraph_vector_size(filter) != no_of_nodes) {
    IGRAPH_ERROR("ML ADE dpareto evaf: invalid filter vector size", 
		 IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_maxdegree(graph, &pmaxdegree, igraph_vss_all(), 
				IGRAPH_IN, /*loops=*/1));
  maxdegree=pmaxdegree;
  
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&neis, maxdegree));
  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);

  IGRAPH_CHECK(igraph_vector_ptr_init(&A_ptr, nopar));
  IGRAPH_FINALLY(igraph_vector_ptr_destroy, &A_ptr);
  IGRAPH_FINALLY(igraph_i_revolver_ml_ADE_dpareto_evalf_free, &A_ptr);
  for (i=0; i<nopar; i++) {
    igraph_array3_t *A=igraph_Calloc(1, igraph_array3_t);
    igraph_array3_init(A, nocats, maxdegree+1, agebins);
    VECTOR(A_ptr)[i] = A;
  }

  IGRAPH_VECTOR_INIT_FINALLY(&S, nopar);
  IGRAPH_CHECK(igraph_vector_resize(value, nopar));
  igraph_vector_null(value);

  for (t=0; t<nopar; t++) {
    igraph_real_t alpha=MATRIX(*par,t,0);
    igraph_real_t a=MATRIX(*par,t,1);
    igraph_real_t paralpha=MATRIX(*par,t,2);
    igraph_real_t parbeta=MATRIX(*par,t,3);
    igraph_real_t parscale=MATRIX(*par,t,4);
    igraph_array3_t *A=VECTOR(A_ptr)[t];
    for (c=0; c<nocats; c++) {
      igraph_real_t cc= c==0 ? 1.0 : MATRIX(*par,t,c+4);
      for (i=0; i<maxdegree+1; i++) {
	igraph_real_t p1= i==0 ? a : pow(i, alpha)+a;
	for (j=0; j<agebins; j++) {
	  igraph_real_t age=j+1;
	  ARRAY3(*A, c, i, j) = age < parscale ? 
	    cc * p1 * pow(age/parscale, parbeta-1) :
	    cc * p1 * pow(age/parscale, -paralpha-1);
	}
      }
    }
  }    

  for (t=0; t<no_of_nodes; t++) {
    long int n, nneis, shnode;
    long int tcat=VECTOR(*cats)[t];
    
    igraph_neighbors(graph, &neis, t, IGRAPH_OUT);
    nneis=igraph_vector_size(&neis);
    
    if (! filter || VECTOR(*filter)[t]) {
      
      for (n=0; n<nneis; n++) {
	long int to=VECTOR(neis)[n];
	long int x=VECTOR(*cats)[to];
	long int y=VECTOR(degree)[to];
	long int z=(t-to)/binwidth;
	
	for (i=0; i<nopar; i++) {
	  igraph_array3_t *A=VECTOR(A_ptr)[i];
	  VECTOR(*value)[i] -= log(ARRAY3(*A, x, y, z));
	  VECTOR(*value)[i] += log(VECTOR(S)[i]);
	}
	no_of_edges++;
      }
    }
    
    for (n=0; n<nneis; n++) {
      long int to=VECTOR(neis)[n];
      long int x=VECTOR(*cats)[to];
      long int y=VECTOR(degree)[to];
      long int z=(t-to)/binwidth;
      
      VECTOR(degree)[to] += 1;
      for (i=0; i<nopar; i++) {
	igraph_array3_t *A=VECTOR(A_ptr)[i];
	VECTOR(S)[i] += ARRAY3(*A, x, y+1, z);
	VECTOR(S)[i] -= ARRAY3(*A, x, y, z);
      }
    }

    for (i=0; i<nopar; i++) {
      igraph_array3_t *A=VECTOR(A_ptr)[i];
      VECTOR(S)[i] += ARRAY3(*A, tcat, 0, 0);
    }
    for (j=1, shnode=t-binwidth+1; shnode>=0; j++, shnode-=binwidth) {
      long int cat=VECTOR(*cats)[shnode];
      long int deg=VECTOR(degree)[shnode];
      for (i=0; i<nopar; i++) {
	igraph_array3_t *A=VECTOR(A_ptr)[i];
	VECTOR(S)[i] += ARRAY3(*A, cat, deg, j);
	VECTOR(S)[i] -= ARRAY3(*A, cat, deg, j-1);
      }
    }
  } /* t < no_of_nodes */

  for (i=0; i<nopar; i++) {
    VECTOR(*value)[i] /= no_of_edges;
  }

  igraph_vector_destroy(&S);
  igraph_i_revolver_ml_ADE_dpareto_evalf_free(&A_ptr);
  igraph_vector_ptr_destroy(&A_ptr);
  igraph_vector_long_destroy(&degree);
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(5);
  
  return 0;
}
  
  
  
  

/*------------------------------------------------------------------*/

int igraph_revolver_ml_d(const igraph_t *graph,
			 igraph_integer_t niter,
			 igraph_vector_t *kernel,
			 igraph_vector_t *cites,
			 igraph_real_t delta,
			 const igraph_vector_t *filter,
			 igraph_real_t *logprob,
			 igraph_real_t *logmax) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_integer_t imaxdegree;
  long int maxdegree, actmaxdegree;
  long int it, t, i;
  igraph_vector_long_t ptk;
  igraph_vector_t *mycites, vmycites;
  igraph_vector_t neis;
  igraph_vector_long_t degree;
  igraph_real_t S=0, maxdelta, diff;

  igraph_vector_t vmykernel;
  igraph_vector_t *kernels[]={ kernel, &vmykernel };
  long int actkernel=0;
  igraph_vector_t *fromkernel=kernels[actkernel], 
    *tokernel=kernels[1-actkernel];

  if (filter && igraph_vector_size(filter) != no_of_nodes) { 
    IGRAPH_ERROR("ML d evolver: invalid filter vector size", IGRAPH_EINVAL);
  }
  
  IGRAPH_CHECK(igraph_maxdegree(graph, &imaxdegree, igraph_vss_all(),
				IGRAPH_IN, IGRAPH_LOOPS));
  maxdegree=imaxdegree;
  
  IGRAPH_CHECK(igraph_vector_long_init(&ptk, maxdegree+1));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &ptk);
  
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);
  IGRAPH_VECTOR_INIT_FINALLY(&vmykernel, maxdegree+1);
  
  if (cites) {
    IGRAPH_CHECK(igraph_vector_resize(cites, maxdegree+1));
    igraph_vector_null(cites);
    mycites=cites;
  } else {
    IGRAPH_VECTOR_INIT_FINALLY(&vmycites, maxdegree+1);
    mycites=&vmycites;
  }

  IGRAPH_CHECK(igraph_vector_resize(kernel, maxdegree+1));
  igraph_vector_fill(kernel, 1);

  IGRAPH_PROGRESS("ML Revolver d", 0, NULL);

  for (it=0; it<niter; it++) {
    
    igraph_vector_null(tokernel);
    igraph_vector_long_null(&ptk);
    igraph_vector_long_null(&degree);
    S=0.0;
    actmaxdegree=0;

    if (logprob) { *logprob=0.0; }
    if (logmax) { *logmax=0.0; }
  
    for (t=0; t<no_of_nodes; t++) {
      long int n, nneis;
      IGRAPH_CHECK(igraph_neighbors(graph, &neis, t, IGRAPH_OUT));
      nneis=igraph_vector_size(&neis);

      IGRAPH_ALLOW_INTERRUPTION();      

      if (!filter || VECTOR(*filter)[t] != 0) {
	/* Calculate some terms of the sum for the non-zero classes */
	if (S != 0) {
	  for (i=0; i<=actmaxdegree; i++) {
	    VECTOR(*tokernel)[i] += nneis * VECTOR(ptk)[i] / S;
	  }
	  
	  if (logprob || logmax || it==0) {
	    for (n=0; n<nneis; n++) {
	      long int to=VECTOR(neis)[n];
	      long int x=VECTOR(degree)[to];
	      if (logprob) { *logprob += log( VECTOR(*fromkernel)[x] / S ); }
	      if (logmax) { *logmax += log(1.0/t); }
	      
	      if (it==0) {
		VECTOR(*mycites)[x] += 1;
	      }
	    }
	  }
	}
      }
      
      /* Update ptk for the next time step */
      for (n=0; n<nneis; n++) {
	long int to=VECTOR(neis)[n];
	long int x=VECTOR(degree)[to];

	VECTOR(degree)[to] += 1;
	if (x==actmaxdegree) { actmaxdegree++; }

	VECTOR(ptk)[x+1] += 1;
	VECTOR(ptk)[x] -= 1;
	S += VECTOR(*fromkernel)[x+1];
	S -= VECTOR(*fromkernel)[x];
      		
      }
      VECTOR(ptk)[0] += 1;
      S += VECTOR(*fromkernel)[0];

    } /* t<no_of_nodes  */
    
    /* final step, Mk/sum */
    maxdelta=0.0;
    for (i=0; i<=maxdegree; i++) {
      if (VECTOR(*tokernel)[i] != 0) {
	VECTOR(*tokernel)[i] = VECTOR(*mycites)[i] / VECTOR(*tokernel)[i];      
	if ( (diff=abs(VECTOR(*tokernel)[i] - VECTOR(*fromkernel)[i])) > maxdelta) {
	  maxdelta=diff;
	}
      }
    }
    if (maxdelta < delta) { break; }
    
    /* Switch kernels */
    actkernel=1-actkernel;
    fromkernel=kernels[actkernel];
    tokernel=kernels[1-actkernel];

    IGRAPH_PROGRESS("ML Revolver d", 100*(it+1)/niter, NULL);

  } /* it<niter */

  /* switch kernels if needed */
  if (fromkernel != kernel) {
    igraph_vector_clear(kernel);
    igraph_vector_append(kernel, fromkernel);
  }
  VECTOR(*kernel)[maxdegree]=IGRAPH_NAN;

  if (!cites) {
    igraph_vector_destroy(&vmycites);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  igraph_vector_destroy(&vmykernel);
  igraph_vector_long_destroy(&degree);
  igraph_vector_destroy(&neis);
  igraph_vector_long_destroy(&ptk);
  IGRAPH_FINALLY_CLEAN(4);
  return 0;
}
			 
int igraph_random_permutation(igraph_vector_t *v) {
  long int i, j, n=igraph_vector_size(v);
  igraph_real_t tmp;
  for (i=0; i<n; i++) {
    VECTOR(*v)[i]=i;
  }
  for (i=n-1; i>=0; i--) {
    j=RNG_INTEGER(0,i);
    tmp=VECTOR(*v)[i];
    VECTOR(*v)[i] = VECTOR(*v)[j];
    VECTOR(*v)[j] = tmp;
  }
  return 0;
}

int igraph_revolver_ml_f(const igraph_t *graph,
			 igraph_integer_t niter,
			 igraph_vector_t *kernel,
			 igraph_vector_t *cites,
			 igraph_real_t delta,
			 igraph_real_t *logprob,
			 igraph_real_t *logmax) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int it, t;
  igraph_vector_long_t ptk;
  igraph_vector_t *mycites, vmycites;
  igraph_vector_t *neis, *neis2;
  igraph_adjlist_t outadjlist, inadjlist;
  igraph_vector_long_t marked;
  
  igraph_vector_t vmykernel;
  igraph_vector_t *kernels[]={ kernel, &vmykernel };
  long int actkernel=0;
  igraph_vector_t *fromkernel=kernels[actkernel], 
    *tokernel=kernels[1-actkernel];
  
  igraph_vector_t perm;
  
  IGRAPH_VECTOR_INIT_FINALLY(&perm, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&perm, no_of_nodes));

  IGRAPH_CHECK(igraph_vector_long_init(&ptk, 2));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &ptk);

  IGRAPH_CHECK(igraph_adjlist_init(graph, &outadjlist, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &outadjlist);
  IGRAPH_CHECK(igraph_adjlist_init(graph, &inadjlist, IGRAPH_IN));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &inadjlist);

  IGRAPH_VECTOR_INIT_FINALLY(&vmykernel, 2);
  IGRAPH_CHECK(igraph_vector_long_init(&marked, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &marked);
  
  if (cites) {
    IGRAPH_CHECK(igraph_vector_resize(cites, 2));
    igraph_vector_null(cites);
    mycites=cites;
  } else { 
    IGRAPH_VECTOR_INIT_FINALLY(&vmycites, 2);
    mycites=&vmycites;
  }
  
  IGRAPH_CHECK(igraph_vector_resize(kernel, 2));
  igraph_vector_fill(kernel, 1);
  
  IGRAPH_PROGRESS("ML revolver f", 0, NULL);

  RNG_BEGIN();
  
  for (it=0; it<niter; it++) {
    
    igraph_vector_null(tokernel);
    igraph_vector_long_null(&ptk);
    igraph_vector_long_null(&marked);

    if (logprob) { *logprob=0.0; }
    if (logmax) { *logmax=0.0; }

    for (t=0; t<no_of_nodes; t++) {
      long int nneis, e;
      neis=igraph_adjlist_get(&outadjlist, t);
      nneis=igraph_vector_size(neis);
      igraph_vector_resize(&perm, nneis);
      igraph_random_permutation(&perm);
      
      IGRAPH_ALLOW_INTERRUPTION();
      
      VECTOR(ptk)[0]=t;
      VECTOR(ptk)[1]=0;
      for (e=0; e<nneis; e++) {
	long int nneis2, j;
	long int which=VECTOR(perm)[e];
	long int to=VECTOR(*neis)[which];
	long int x= VECTOR(marked)[to] != t+1 ? 0 : 1;
	
	if (e != 0) {
	  igraph_real_t S=VECTOR(ptk)[0] * VECTOR(*fromkernel)[0] +
	    VECTOR(ptk)[1] * VECTOR(*fromkernel)[1];
	  VECTOR(*tokernel)[0] += VECTOR(ptk)[0] / S;
	  VECTOR(*tokernel)[1] += VECTOR(ptk)[1] / S;
	  
	  if (it==0) {
	    VECTOR(*mycites)[x] += 1;
	  }

	  if (logprob) { *logprob += log( VECTOR(*fromkernel)[x] / S ); }
	  if (logmax) { *logmax += log(1.0/t); }

	} else {
	  if (logprob) { *logprob += log(1.0/t); }
	  if (logmax) { *logmax += log(1.0/t); }
	}

 	VECTOR(ptk)[x] -= 1;
	VECTOR(marked)[to]=t+1;

	/* Update ptk, check the neighbors of 'to' */
	neis2=igraph_adjlist_get(&inadjlist, to);
	nneis2=igraph_vector_size(neis2);
	for (j=0; j<nneis2; j++) {
	  long int nei=VECTOR(*neis2)[j];
	  if (nei >= t) { break; }
	  if (VECTOR(marked)[nei] != t+1) {
	    VECTOR(marked)[nei] = t+1;
	    VECTOR(ptk)[0] -= 1;
	    VECTOR(ptk)[1] += 1;
	  }
	}
	neis2=igraph_adjlist_get(&outadjlist, to);
	nneis2=igraph_vector_size(neis2);
	for (j=0; j<nneis2; j++) {
	  long int nei=VECTOR(*neis2)[j];
	  if (VECTOR(marked)[nei] != t+1) {
	    VECTOR(marked)[nei] = t+1;
	    VECTOR(ptk)[0] -= 1;
	    VECTOR(ptk)[1] += 1;
	  }
	}
	
      }
    }
      
    /* Mk/sum */
    VECTOR(*tokernel)[0] = VECTOR(*mycites)[0] / VECTOR(*tokernel)[0];
    VECTOR(*tokernel)[1] = VECTOR(*mycites)[1] / VECTOR(*tokernel)[1];
    
/*     VECTOR(*tokernel)[1] /= VECTOR(*tokernel)[0]; */
/*     VECTOR(*tokernel)[0] = 1.0; */
    
    /* Switch kernels */
    actkernel=1-actkernel;
    fromkernel=kernels[actkernel];
    tokernel=kernels[1-actkernel];

    IGRAPH_PROGRESS("ML Revolver f", 100*(it+1)/niter, NULL);
    
  } /* it<niter */

  RNG_END();
  
  /* switch kernels if needed */
  if (fromkernel != kernel) {
    igraph_vector_clear(kernel);
    igraph_vector_append(kernel, fromkernel);
  }

  if (!cites) {
    igraph_vector_destroy(&vmycites);
    IGRAPH_FINALLY_CLEAN(1);
  }

  igraph_vector_long_destroy(&marked);
  igraph_vector_destroy(&vmykernel);
  igraph_adjlist_destroy(&inadjlist);
  igraph_adjlist_destroy(&outadjlist);
  igraph_vector_long_destroy(&ptk);
  igraph_vector_destroy(&perm);
  IGRAPH_FINALLY_CLEAN(6);
  return 0;
}

int igraph_revolver_ml_df(const igraph_t *graph,
			  igraph_integer_t niter,
			  igraph_matrix_t *kernel,
			  igraph_matrix_t *cites,
			  igraph_real_t delta,
			  igraph_real_t *logprob,
			  igraph_real_t *logmax) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int it, t, i;
  igraph_matrix_long_t ptk;
  igraph_matrix_t *mycites, vmycites;
  igraph_vector_t *neis, *neis2;
  igraph_adjlist_t outadjlist, inadjlist;
  igraph_vector_long_t marked;
  igraph_integer_t imaxdegree, omaxdegree;
  long int maxdegree;
  igraph_vector_long_t degree;
  igraph_real_t S1, S2, S3;
  long int actmaxdegree;
  
  igraph_matrix_t vmykernel;
  igraph_matrix_t *kernels[] = { kernel, &vmykernel };
  long int actkernel=0;
  igraph_matrix_t *fromkernel=kernels[actkernel],
    *tokernel=kernels[1-actkernel];

  igraph_stack_t stack;
  igraph_vector_t perm;
  
  IGRAPH_CHECK(igraph_maxdegree(graph, &imaxdegree, igraph_vss_all(), 
				IGRAPH_IN, IGRAPH_LOOPS));
  IGRAPH_CHECK(igraph_maxdegree(graph, &omaxdegree, igraph_vss_all(),
				IGRAPH_OUT, IGRAPH_LOOPS));
  maxdegree=imaxdegree;

  IGRAPH_VECTOR_INIT_FINALLY(&perm, omaxdegree);

  IGRAPH_CHECK(igraph_stack_init(&stack, maxdegree * maxdegree > no_of_nodes ?
				      no_of_nodes : maxdegree * maxdegree));
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);

  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);

  IGRAPH_CHECK(igraph_matrix_long_init(&ptk, 2, maxdegree+1));
  IGRAPH_FINALLY(igraph_matrix_long_destroy, &ptk);
  
  IGRAPH_CHECK(igraph_adjlist_init(graph, &outadjlist, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &outadjlist);
  IGRAPH_CHECK(igraph_adjlist_init(graph, &inadjlist, IGRAPH_IN));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &inadjlist);
  
  IGRAPH_MATRIX_INIT_FINALLY(&vmykernel, 3, maxdegree+1);
  IGRAPH_CHECK(igraph_vector_long_init(&marked, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &marked);
  
  if (cites) {
    IGRAPH_CHECK(igraph_matrix_resize(cites, 3, maxdegree+1));
    igraph_matrix_null(cites);
    mycites=cites;
  } else {
    IGRAPH_MATRIX_INIT_FINALLY(&vmycites, 3, maxdegree+1);
    mycites=&vmycites;
  }
  
  IGRAPH_CHECK(igraph_matrix_resize(kernel, 3, maxdegree+1));
  igraph_matrix_fill(kernel, 1);
  
  IGRAPH_PROGRESS("ML revolver df", 0, NULL);

  RNG_BEGIN();
  
  for (it=0; it<niter; it++) {
    
    igraph_matrix_null(tokernel);
    igraph_matrix_long_null(&ptk);
    igraph_vector_long_null(&marked);
    igraph_vector_long_null(&degree);
    S1=0; S2=0; S3=0;
    actmaxdegree=0;
    
    if (logprob) { *logprob=0.0; }
    if (logmax) { *logmax=0.0; }
    
    for (t=0; t<no_of_nodes; t++) {
      long int nneis, e;
      neis=igraph_adjlist_get(&outadjlist, t);
      nneis=igraph_vector_size(neis);
      igraph_vector_resize(&perm, nneis);
      igraph_random_permutation(&perm);
      
      IGRAPH_ALLOW_INTERRUPTION();

      /* Restore ptk */
      while (! igraph_stack_empty(&stack)) {
	long int deg=igraph_stack_pop(&stack);
	if (deg > 0) {
	  MATRIX(ptk, 0, deg-1) += 1;
	  MATRIX(ptk, 1, deg-1) = 0;
	} else {
	  MATRIX(ptk, 0, -deg-1) -= 1;
	  MATRIX(ptk, 1, -deg-1) = 0;
	}
      }
      S2=S3;
      
      for (e=0; e<nneis; e++) {
	long int nneis2, i, j;
	long int which=VECTOR(perm)[e];
	long int to=VECTOR(*neis)[which];
	long int x=VECTOR(marked)[to] != t+1 ? 0 : 1;
	long int y=VECTOR(degree)[to];
		  
	if (e == 0) {
	  /* First citation */
	  for (i=0; i<=actmaxdegree; i++) {
	    MATRIX(*tokernel, 0, i) += MATRIX(ptk, 0, i) / S1;
	  }
	  
	  if (it==0) {
	    MATRIX(*mycites, 0, y) += 1;
	  }
	  
	  if (logprob && MATRIX(*fromkernel, 0, y) != 0) { 
	    *logprob += log( MATRIX(*fromkernel, 0, y) / S1); 
	  }
	  if (logmax) { *logmax += log(1.0 / t); }
	  
	} else {
	  /* Subsequent citations */
	  for (i=0; i<=actmaxdegree; i++) {
	    MATRIX(*tokernel, 1, i) += MATRIX(ptk, 0, i) / S2;
	    MATRIX(*tokernel, 2, i) += MATRIX(ptk, 1, i) / S2;
	  }
	  
	  if (it==0) {
	    MATRIX(*mycites, x+1, y) += 1;
	  }

	  if (logprob && MATRIX(*fromkernel, x+1, y) != 0) { 
	    *logprob += log( MATRIX(*fromkernel, x+1, y) / S2 ); 
	  }
	  if (logmax) { *logmax += log(1.0/t); }
	}

	/* update ptk */
	VECTOR(marked)[to]=t+1;
	VECTOR(degree)[to] += 1;
	if (VECTOR(degree)[to] > actmaxdegree) { actmaxdegree ++; }
	MATRIX(ptk, x, y) -= 1; /* won't be cited again by this vertex */
	S1 -= MATRIX(*fromkernel, 0, y);
	S1 += MATRIX(*fromkernel, 0, y+1);
	S3 -= MATRIX(*fromkernel, 1, y);
	S3 += MATRIX(*fromkernel, 1, y+1);
	S2 -= MATRIX(*fromkernel, x+1, y);
	if (x==0) {
	  igraph_stack_push(&stack, y+2);
	} else {
	  igraph_stack_push(&stack, -y-1);
	  igraph_stack_push(&stack, y+2);
	}
	
	/* neighbors of 'to' */
	neis2=igraph_adjlist_get(&inadjlist, to);
	nneis2=igraph_vector_size(neis2);
	for (j=0; j<nneis2; j++) {
	  long int nei=VECTOR(*neis2)[j];
	  if (nei >= t) { break; }
	  if (VECTOR(marked)[nei] != t+1) {
	    long int neideg=VECTOR(degree)[nei];
	    VECTOR(marked)[nei] = t+1;
	    MATRIX(ptk, 0, neideg) -= 1;
	    MATRIX(ptk, 1, neideg) += 1;
	    S2 -= MATRIX(*fromkernel, 1, neideg) - MATRIX(*fromkernel, 2, neideg);
	    igraph_stack_push(&stack, neideg+1);
	  }
	}
	neis2=igraph_adjlist_get(&outadjlist, to);
	nneis2=igraph_vector_size(neis2);
	for (j=0; j<nneis2; j++) {
	  long int nei=VECTOR(*neis2)[j];
	  if (VECTOR(marked)[nei] != t+1) {
	    long int neideg=VECTOR(degree)[nei];
	    VECTOR(marked)[nei] = t+1;
	    MATRIX(ptk, 0, neideg) -= 1;
	    MATRIX(ptk, 1, neideg) += 1;
	    S2 -= MATRIX(*fromkernel, 1, neideg) - MATRIX(*fromkernel, 2, neideg);
	    igraph_stack_push(&stack, neideg+1);
	  }
	}
	
      }	/* e < nneis */
      
      S1 += MATRIX(*fromkernel, 0, 0);
      S3 += MATRIX(*fromkernel, 1, 0);
      MATRIX(ptk, 0, 0) += 1;
      
    } /* t < no_of_nodes */

    /* Mk/sum */
    for (i=0; i<maxdegree+1; i++) {
      if (MATRIX(*tokernel, 0, i) != 0) {
	MATRIX(*tokernel, 0, i) = MATRIX(*mycites, 0, i) / MATRIX(*tokernel, 0, i);
      }
      if (MATRIX(*tokernel, 1, i) != 0) {
	MATRIX(*tokernel, 1, i) = MATRIX(*mycites, 1, i) / MATRIX(*tokernel, 1, i);
      }
      if (MATRIX(*tokernel, 2, i) != 0) {
	MATRIX(*tokernel, 2, i) = MATRIX(*mycites, 2, i) / MATRIX(*tokernel, 2, i);
      }
    }
    
    /* Switch kernels */
    actkernel=1-actkernel;
    fromkernel=kernels[actkernel];
    tokernel=kernels[1-actkernel];
    
    IGRAPH_PROGRESS("ML Revolver df", 100*(it+1)/niter, NULL);
    
  } /* it < niter */

  RNG_END();

  /* switch kernels if needed */
  if (fromkernel != kernel) {
    igraph_matrix_update(kernel, fromkernel);
  }
  if (!cites) {
    igraph_matrix_destroy(&vmycites);
    IGRAPH_FINALLY_CLEAN(1);
  }

  igraph_vector_long_destroy(&marked);
  igraph_matrix_destroy(&vmykernel);
  igraph_adjlist_destroy(&inadjlist);
  igraph_adjlist_destroy(&outadjlist);
  igraph_matrix_long_destroy(&ptk);
  igraph_vector_long_destroy(&degree);
  igraph_stack_destroy(&stack);
  igraph_vector_destroy(&perm);
  IGRAPH_FINALLY_CLEAN(8);
  
  return 0;
}

int igraph_revolver_ml_ad(const igraph_t *graph,
			  igraph_integer_t niter,
			  igraph_matrix_t *kernel,
			  igraph_matrix_t *cites,
			  igraph_integer_t pagebins,
			  igraph_real_t delta,
			  const igraph_vector_t *filter,
			  igraph_real_t *logprob,
			  igraph_real_t *logmax) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_integer_t imaxdegree;
  long int maxdegree, actmaxdegree;
  long int it, t, i, j;
  igraph_matrix_long_t ptk;
  igraph_matrix_t *mycites, vmycites;
  igraph_vector_t neis;
  igraph_vector_long_t degree;
  igraph_real_t S=0, maxdelta, diff;
  long int agebins=pagebins;
  long int binwidth=no_of_nodes/agebins+1;
  
  igraph_matrix_t vmykernel;
  igraph_matrix_t *kernels[]= { kernel, &vmykernel };
  long int actkernel=0;
  igraph_matrix_t *fromkernel=kernels[actkernel],
    *tokernel=kernels[1-actkernel];

  if (filter && igraph_vector_size(filter) != no_of_nodes) {
    IGRAPH_ERROR("ML ad revolver: invalid filter vector size", IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_maxdegree(graph, &imaxdegree, igraph_vss_all(),
				IGRAPH_IN, IGRAPH_LOOPS));
  maxdegree=imaxdegree;
  
  IGRAPH_CHECK(igraph_matrix_long_init(&ptk, maxdegree+1, agebins));
  IGRAPH_FINALLY(igraph_matrix_long_destroy, &ptk);
  
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);
  IGRAPH_CHECK(igraph_matrix_init(&vmykernel, maxdegree+1, agebins));
  IGRAPH_FINALLY(igraph_matrix_destroy, &vmykernel);
  
  if (cites) {
    IGRAPH_CHECK(igraph_matrix_resize(cites, maxdegree+1, agebins));
    igraph_matrix_null(cites);
    mycites=cites;
  } else {
    IGRAPH_CHECK(igraph_matrix_init(&vmycites, maxdegree+1, agebins));
    IGRAPH_FINALLY(igraph_matrix_destroy, &vmycites);
    mycites=&vmycites;
  }
  
  IGRAPH_CHECK(igraph_matrix_resize(kernel, maxdegree+1, agebins));
  igraph_matrix_fill(kernel, 1);
  
  IGRAPH_PROGRESS("ML Revolver ad", 0, NULL);
  
  for (it=0; it<niter; it++) {
   
    igraph_matrix_null(tokernel);
    igraph_matrix_long_null(&ptk);
    igraph_vector_long_null(&degree);
    S=0.0;
    actmaxdegree=0;
    
    if (logprob) { *logprob=0.0; }
    if (logmax) { *logmax=0.0; }
    
    for (t=0; t<no_of_nodes; t++) {
      long int n, nneis;
      IGRAPH_CHECK(igraph_neighbors(graph, &neis, t, IGRAPH_OUT));
      nneis=igraph_vector_size(&neis);
      
      IGRAPH_ALLOW_INTERRUPTION();
      
      if (!filter || VECTOR(*filter)[t] != 0) {
      
	/* Calculate some terms of the sum for the non-zero classes */
	if (S != 0) {
	  for (i=0; i<=actmaxdegree; i++) {
	    for (j=0; j<=t/binwidth; j++) {
	      MATRIX(*tokernel, i, j) += nneis * MATRIX(ptk, i, j) / S;
	    }
	  }
	}
	
	if (logprob || logmax || it==0) {
	  for (n=0; n<nneis; n++) {
	    long int to=VECTOR(neis)[n];
	    long int x=VECTOR(degree)[to];
	    long int y=(t-to)/binwidth;
	    if (logprob) { *logprob += log( MATRIX(*fromkernel,x,y) / S ); }
	    if (logmax) { *logmax += log(1.0/t); }
	    if (it==0) {
	      MATRIX(*mycites, x, y) += 1;
	    }
	  }
	}
      }

      /* Update ptk for the next time step */
      for (n=0; n<nneis; n++) {
	long int to=VECTOR(neis)[n];
	long int x=VECTOR(degree)[to];
	long int y=(t-to)/binwidth;

	VECTOR(degree)[to] += 1;
	if (x==actmaxdegree) { actmaxdegree++; }
	
	MATRIX(ptk, x+1, y) += 1;
	MATRIX(ptk, x, y) -= 1;
	S += MATRIX(*fromkernel, x+1, y);
	S -= MATRIX(*fromkernel, x, y);
	
      }
      /* Aging */
      for (j=1; t-binwidth*j+1>=0; j++) {
	long int shnode=t-binwidth*j+1;
	long int deg=VECTOR(degree)[shnode];
	MATRIX(ptk, deg, j) += 1;
	MATRIX(ptk, deg, j-1) -= 1;
	S += MATRIX(*fromkernel, deg, j);
	S -= MATRIX(*fromkernel, deg, j-1);
      }
      MATRIX(ptk, 0, 0) += 1;
      S += MATRIX(*fromkernel, 0, 0);
      
    } /* t < no_of_nodes */
    
    /* Mk/sum */
    maxdelta=0.0;
    for (i=0; i<=maxdegree; i++) {
      for (j=0; j<agebins; j++) {
	MATRIX(*tokernel, i, j) = MATRIX(*mycites, i, j) / MATRIX(*tokernel, i, j);
	if ( (diff=abs(MATRIX(*tokernel,i,j) - MATRIX(*fromkernel,i,j))) > maxdelta){
	  maxdelta=diff;
	}
      }
    }
    if (maxdelta < delta) { break; }
    
    /* Switch kernels */
    actkernel=1-actkernel;
    fromkernel=kernels[actkernel];
    tokernel=kernels[1-actkernel];

    IGRAPH_PROGRESS("ML Revolver d", 100*(it+1)/niter, NULL);

  } /* it<niter */

  /* switch kernels if needed */
  if (fromkernel != kernel) {
    igraph_matrix_update(kernel, fromkernel);
  }

  if (!cites) {
    igraph_matrix_destroy(&vmycites);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  igraph_matrix_destroy(&vmykernel);
  igraph_vector_long_destroy(&degree);
  igraph_vector_destroy(&neis);
  igraph_matrix_long_destroy(&ptk);
  IGRAPH_FINALLY_CLEAN(4);
  
  return 0;
}

int igraph_revolver_ml_de(const igraph_t *graph,
			  igraph_integer_t niter,
			  igraph_matrix_t *kernel,
			  const igraph_vector_t *cats,
			  igraph_matrix_t *cites,
			  igraph_real_t delta,
			  const igraph_vector_t *filter,
			  igraph_real_t *logprob,
			  igraph_real_t *logmax) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_integer_t imaxdegree;
  long int maxdegree, actmaxdegree;
  long int it, t, i, j;
  igraph_matrix_long_t ptk;
  igraph_matrix_t *mycites, vmycites;
  igraph_vector_t neis;
  igraph_vector_long_t degree;
  igraph_real_t S=0, maxdelta, diff;
  long int no_cats=igraph_vector_max(cats)+1;
  
  igraph_matrix_t vmykernel;
  igraph_matrix_t *kernels[] = { kernel, &vmykernel };
  long int actkernel=0;
  igraph_matrix_t *fromkernel=kernels[actkernel],
    *tokernel=kernels[1-actkernel];

  if (filter && igraph_vector_size(filter) != no_of_nodes) {
    IGRAPH_ERROR("ML de evolver failed", IGRAPH_EINVAL);
  }
  
  IGRAPH_CHECK(igraph_maxdegree(graph, &imaxdegree, igraph_vss_all(),
				IGRAPH_IN, IGRAPH_LOOPS));
  maxdegree=imaxdegree;
  
  IGRAPH_CHECK(igraph_matrix_long_init(&ptk, no_cats, maxdegree+1));
  IGRAPH_FINALLY(igraph_matrix_long_destroy, &ptk);
  
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);
  IGRAPH_MATRIX_INIT_FINALLY(&vmykernel, no_cats, maxdegree+1);
  
  if (cites) {
    IGRAPH_CHECK(igraph_matrix_resize(cites, no_cats, maxdegree+1));
    igraph_matrix_null(cites);
    mycites=cites;
  } else {
    IGRAPH_MATRIX_INIT_FINALLY(&vmycites, no_cats, maxdegree+1);
    mycites=&vmycites;
  }
  
  IGRAPH_CHECK(igraph_matrix_resize(kernel, no_cats, maxdegree+1));
  igraph_matrix_fill(kernel, 1);
  
  for (it=0; it<niter; it++) {
    
    igraph_matrix_null(tokernel);
    igraph_matrix_long_null(&ptk);
    igraph_vector_long_null(&degree);
    S=0.0;
    actmaxdegree=0;
    
    if (logprob) { *logprob=0.0; }
    if (logmax) { *logmax=0.0; }
    
    for (t=0; t<no_of_nodes; t++) {
      long int n, nneis;
      long int fromcat=VECTOR(*cats)[t];
      IGRAPH_CHECK(igraph_neighbors(graph, &neis, t, IGRAPH_OUT));
      nneis=igraph_vector_size(&neis);
      
      IGRAPH_ALLOW_INTERRUPTION();

      if (!filter || VECTOR(*filter)[t] != 0) {
      
	if (S != 0) {
	  for (i=0; i<no_cats; i++) {
	    for (j=0; j<=actmaxdegree; j++) {
	      MATRIX(*tokernel, i, j) += nneis * MATRIX(ptk, i, j) / S;
	    }
	  }
	  
	  if (logprob || logmax || it==0) {
	    for (n=0; n<nneis; n++) {
	      long int to=VECTOR(neis)[n];
	      long int x=VECTOR(*cats)[to];
	      long int y=VECTOR(degree)[to];
	      if (logprob) { *logprob += log( MATRIX(*fromkernel,x,y) / S); }
	      if (logmax) { *logmax += log(1.0/t); }
	      if (it==0) {
		MATRIX(*mycites, x, y) += 1;
	      }
	    }
	  }
	}

      }
	
      /* Update ptk for the next time step */
      for (n=0; n<nneis; n++) {
	long int to=VECTOR(neis)[n];
	long int x=VECTOR(*cats)[to];
	long int y=VECTOR(degree)[to];
	
	VECTOR(degree)[to] += 1;
	if (y==actmaxdegree) { actmaxdegree++; }
	
	MATRIX(ptk, x, y+1) += 1;
	MATRIX(ptk, x, y) -= 1;
	S += MATRIX(*fromkernel, x, y+1);
	S -= MATRIX(*fromkernel, x, y);
	
      }
      
      MATRIX(ptk, fromcat, 0) += 1;
      S += MATRIX(*fromkernel, fromcat, 0);
      
    } /* t < no_of_nodes */
  
    /* Mk/sum */
    maxdelta=0.0;
    for (i=0; i<no_cats; i++) {
      for (j=0; j<=maxdegree; j++) {
	if (MATRIX(*tokernel, i, j) != 0) {
	  MATRIX(*tokernel, i, j) = MATRIX(*mycites, i, j) / MATRIX(*tokernel, i, j);
	  if ( (diff=abs(MATRIX(*tokernel,i,j)-MATRIX(*fromkernel,i,j))) > maxdelta) {
	    maxdelta=diff;
	  }
	}
      }
    }
    if (maxdelta < delta) { break; }
    
    /* Switch kernels */
    actkernel=1-actkernel;
    fromkernel=kernels[actkernel];
    tokernel=kernels[1-actkernel];
    
    IGRAPH_PROGRESS("ML Revolver de", 100*(it+1)/niter, NULL);
    
  } /* it<niter */

  /* switch kernels if needed */
  if (fromkernel != kernel) {
    igraph_matrix_update(kernel, fromkernel);
  }
  
  if (!cites) {
    igraph_matrix_destroy(&vmycites);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  igraph_matrix_destroy(&vmykernel);
  igraph_vector_long_destroy(&degree);
  igraph_vector_destroy(&neis);
  igraph_matrix_long_destroy(&ptk);
  IGRAPH_FINALLY_CLEAN(4);

  return 0;
}

int igraph_revolver_ml_ade(const igraph_t *graph,
			   igraph_integer_t niter,
			   igraph_array3_t *kernel,
			   const igraph_vector_t *cats,
			   igraph_array3_t *cites,
			   igraph_integer_t pagebins,
			   igraph_real_t delta,
			   const igraph_vector_t *filter,
			   igraph_real_t *logprob,
			   igraph_real_t *logmax) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_integer_t imaxdegree;
  long int maxdegree, actmaxdegree;
  long int it, t, i, j, k;
  igraph_array3_long_t ptk;
  igraph_array3_t *mycites, vmycites;
  igraph_vector_t neis;
  igraph_vector_long_t degree;
  igraph_real_t S=0, maxdelta, diff;
  long int no_cats=igraph_vector_max(cats)+1;

  long int agebins=pagebins;
  long int binwidth=no_of_nodes/agebins+1;
  
  igraph_array3_t vmykernel;
  igraph_array3_t *kernels[]= { kernel, &vmykernel };
  long int actkernel=0;
  igraph_array3_t *fromkernel=kernels[actkernel],
    *tokernel=kernels[1-actkernel];

  if (filter && igraph_vector_size(filter) != no_of_nodes) {
    IGRAPH_ERROR("ML ade revolver: invalid filter vector size", IGRAPH_EINVAL);
  }
  
  IGRAPH_CHECK(igraph_maxdegree(graph, &imaxdegree, igraph_vss_all(),
				IGRAPH_IN, IGRAPH_LOOPS));
  maxdegree=imaxdegree;
  
  IGRAPH_CHECK(igraph_array3_long_init(&ptk, no_cats, maxdegree+1, agebins));
  IGRAPH_FINALLY(igraph_array3_long_destroy, &ptk);
  
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);
  IGRAPH_CHECK(igraph_array3_init(&vmykernel, no_cats, maxdegree+1, agebins));
  IGRAPH_FINALLY(igraph_matrix_destroy, &vmykernel);
  
  if (cites) {
    IGRAPH_CHECK(igraph_array3_resize(cites, no_cats, maxdegree+1, agebins));
    igraph_array3_null(cites);
    mycites=cites;
  } else {
    IGRAPH_CHECK(igraph_array3_init(&vmycites, no_cats, maxdegree+1, agebins));
    IGRAPH_FINALLY(igraph_array3_destroy, &vmycites);
    mycites=&vmycites;
  }

  IGRAPH_CHECK(igraph_array3_resize(kernel, no_cats, maxdegree+1, agebins));
  igraph_array3_fill(kernel, 1);
  
  IGRAPH_PROGRESS("ML Revolver ade", 0, NULL);
  
  for (it=0; it<niter; it++) {
    
    igraph_array3_null(tokernel);
    igraph_array3_long_null(&ptk);
    igraph_vector_long_null(&degree);
    S=0.0;
    actmaxdegree=0;
    
    if (logprob) { *logprob=0.0; }
    if (logmax) { *logmax=0.0; }

    for (t=0; t<no_of_nodes; t++) {
      long int n, nneis;
      long int tcat=VECTOR(*cats)[t];
      IGRAPH_CHECK(igraph_neighbors(graph, &neis, t, IGRAPH_OUT));
      nneis=igraph_vector_size(&neis);
      
      IGRAPH_ALLOW_INTERRUPTION();

      if (!filter || VECTOR(*filter)[t] != 0) {
	
	/* Calculate some terms of the sum for the non-zero classes */
	if (S != 0) {
	  for (i=0; i<no_cats; i++) {
	    for (j=0; j<=actmaxdegree; j++) {
	      for (k=0; k<=t/binwidth; k++) {
		ARRAY3(*tokernel, i, j, k) += nneis * ARRAY3(ptk, i, j, k) / S;
	      }
	    }
	  }
	  
	  if (logprob || logmax || it==0) {
	    for (n=0; n<nneis; n++) {
	      long int to=VECTOR(neis)[n];
	      long int x=VECTOR(*cats)[to];
	      long int y=VECTOR(degree)[to];
	      long int z=(t-to)/binwidth;
	      if (logprob) { *logprob += log( ARRAY3(*fromkernel, x,y,z) / S); }
	      if (logmax) { *logmax += log(1.0/t); }
	      if (it==0) {
		ARRAY3(*mycites, x, y, z) += 1;
	      }
	    }
	  }
	}
      }
      
      /* update ptk for the next time step */
      for (n=0; n<nneis; n++) {
	long int to=VECTOR(neis)[n];
	long int x=VECTOR(*cats)[to];
	long int y=VECTOR(degree)[to];
	long int z=(t-to)/binwidth;
	
	VECTOR(degree)[to] += 1;
	if (y==actmaxdegree) { actmaxdegree++; }
	
	ARRAY3(ptk, x, y+1, z) += 1;
	ARRAY3(ptk, x, y, z) -= 1;
	S += ARRAY3(*fromkernel, x, y+1, z);
	S -= ARRAY3(*fromkernel, x, y, z);
	
      }
      /* Aging */
      for (j=1; t-binwidth*j+1>=0; j++) {
	long int shnode=t-binwidth*j+1;
	long int cat=VECTOR(*cats)[shnode];
	long int deg=VECTOR(degree)[shnode];
	ARRAY3(ptk, cat, deg, j) += 1;
	ARRAY3(ptk, cat, deg, j-1) -= 1;
	S += ARRAY3(*fromkernel, cat, deg, j);
	S -= ARRAY3(*fromkernel, cat, deg, j-1);
      }
      ARRAY3(ptk, tcat, 0, 0) += 1;
      S += ARRAY3(*fromkernel, tcat, 0, 0);
      
    } /* t < no_of_nodes */

    /* Mk/sum */
    maxdelta=0.0;
    for (i=0; i<no_cats; i++) {
      for (j=0; j<=maxdegree; j++) {
	for (k=0; k<agebins; k++) {
	  ARRAY3(*tokernel,i,j,k) = ARRAY3(*mycites,i,j,k) / ARRAY3(*tokernel,i,j,k);
	  if ( (diff=abs(ARRAY3(*tokernel,i,j,k)-ARRAY3(*fromkernel,i,j,k))) > 
	       maxdelta) {
	    maxdelta=diff;
	  }
	}
      }
    }
    if (maxdelta < delta) { break; }
    
    /* Switch kernels */
    actkernel=1-actkernel;
    fromkernel=kernels[actkernel];
    tokernel=kernels[1-actkernel];

    IGRAPH_PROGRESS("ML Revolver d", 100*(it+1)/niter, NULL);

  } /* it<niter */

  /* switch kernels if needed */
  if (fromkernel != kernel) {
    igraph_array3_update(kernel, fromkernel);
  }

  if (!cites) {
    igraph_array3_destroy(&vmycites);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  igraph_array3_destroy(&vmykernel);
  igraph_vector_long_destroy(&degree);
  igraph_vector_destroy(&neis);
  igraph_array3_long_destroy(&ptk);
  IGRAPH_FINALLY_CLEAN(4);
  
  return 0;
}	   

int igraph_revolver_ml_l(const igraph_t *graph,
			 igraph_integer_t niter,
			 igraph_vector_t *kernel,
			 igraph_vector_t *cites,
			 igraph_integer_t pagebins,
			 igraph_real_t delta,
			 igraph_real_t *logprob,
			 igraph_real_t *logmax) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int agebins=pagebins;
  long int binwidth=no_of_nodes/agebins+1;
  
  long int it, t, i, k;
  igraph_vector_long_t ptk;
  igraph_vector_t *mycites, vmycites;
  igraph_vector_t neis;
  igraph_real_t S=0, maxdelta, diff;
  
  igraph_vector_long_t lastcit;
  igraph_vector_t vmykernel;
  igraph_vector_t *kernels[]= { kernel, &vmykernel };
  long int actkernel=0;
  igraph_vector_t *fromkernel=kernels[actkernel],
    *tokernel=kernels[1-actkernel];
  
  IGRAPH_CHECK(igraph_vector_long_init(&ptk, agebins+1));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &ptk);
  
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_vector_long_init(&lastcit, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &lastcit);
  IGRAPH_VECTOR_INIT_FINALLY(&vmykernel, agebins+1);
  
  if (cites) {
    IGRAPH_CHECK(igraph_vector_resize(cites, agebins+1));
    igraph_vector_null(cites);
    mycites=cites;
  } else {
    IGRAPH_VECTOR_INIT_FINALLY(&vmycites, agebins+1);
    mycites=&vmycites;
  }
  
  IGRAPH_CHECK(igraph_vector_resize(kernel, agebins+1));
  igraph_vector_fill(kernel, 1);
  
  IGRAPH_PROGRESS("ML Revolver l", 0, NULL);
  
  for (it=0; it<niter; it++) {
    
    igraph_vector_null(tokernel);
    igraph_vector_long_null(&ptk);
    S=0.0;
    
    if (logprob) { *logprob=0.0; }
    if (logmax) { *logmax=0.0; }
    
    for (t=0; t<no_of_nodes; t++) {
      long int n, nneis;
      IGRAPH_CHECK(igraph_neighbors(graph, &neis, t, IGRAPH_OUT));
      nneis=igraph_vector_size(&neis);

      IGRAPH_ALLOW_INTERRUPTION();
      
      if (S != 0) {
	for (i=0; i<agebins+1; i++) {
	  VECTOR(*tokernel)[i] += nneis * VECTOR(ptk)[i] / S;
	}
	
	if (logprob || logmax) {
	  for (n=0; n<nneis; n++) {
	    long int to=VECTOR(neis)[n];
	    long int x= VECTOR(lastcit)[to] != 0 ? 
	      t+2-(long int)VECTOR(lastcit)[to]/binwidth : agebins;
	    if (logprob) { *logprob += log( VECTOR(*fromkernel)[x] / S ); }
	    if (logmax) { *logmax += log(1.0/t); }
	  }
	}
      }
      
      /* Update ptk for the last time step */
      for (n=0; n<nneis; n++) {
	long int to=VECTOR(neis)[n];
	long int x=VECTOR(lastcit)[to] != 0 ?
	  t+2-(long int)VECTOR(lastcit)[to]/binwidth : agebins;
	
	VECTOR(lastcit)[to]=t+2;
	VECTOR(ptk)[x] += 1;
	S += VECTOR(*fromkernel)[x];
      }
      VECTOR(ptk)[agebins] += 1;
      S += VECTOR(*fromkernel)[agebins];
      /* should we move some citations to an older bin? */
      for (k=1; t+1-binwidth*k+1>=0; k++) {
	long int shnode=t+1-binwidth*k+1;
	IGRAPH_CHECK(igraph_neighbors(graph, &neis, shnode, IGRAPH_OUT));
	nneis=igraph_vector_size(&neis);
	for (i=0; i<nneis; i++) {
	  long int cnode=VECTOR(neis)[i];
	  if (VECTOR(lastcit)[cnode]==shnode+1) {
	    VECTOR(ptk)[k-1] -= 1;
	    VECTOR(ptk)[k] += 1;
	    S -= VECTOR(*fromkernel)[k-1];
	    S += VECTOR(*fromkernel)[k];
	  }
	}
      }
      
    } /* t < no_of_nodes */
    
    /* Mk/sum */
    maxdelta=0.0;
    for (i=0; i<agebins+1; i++) {
      VECTOR(*tokernel)[i] = VECTOR(*mycites)[i] / VECTOR(*tokernel)[i];
      if ( (diff=abs(VECTOR(*tokernel)[i]-VECTOR(*fromkernel)[i])) > maxdelta) {
	maxdelta=diff;
      }
    }
    if (maxdelta < delta) { break; }
    
    /* Switch kernels */
    actkernel=1-actkernel;
    fromkernel=kernels[actkernel];
    tokernel=kernels[1-actkernel];
    
    IGRAPH_PROGRESS("ML Revolver l", 100*(it+1)/niter, NULL);

  } /* it < niter */
  
  /* Switch kernels if needed */
  if (fromkernel != kernel) {
    igraph_vector_update(kernel, fromkernel);
  }
  
  if (!cites) { 
    igraph_vector_destroy(&vmycites);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  igraph_vector_destroy(&vmykernel);
  igraph_vector_long_destroy(&lastcit);
  igraph_vector_destroy(&neis);
  igraph_vector_long_destroy(&ptk);
  IGRAPH_FINALLY_CLEAN(4);
  
  return 0;
}

/* -----------------------------------------------------------*/

int igraph_revolver_probs_d(const igraph_t *graph,
			    const igraph_vector_t *kernel,
			    igraph_vector_t *logprobs,
			    igraph_vector_t *logcited,
			    igraph_vector_t *logciting,
			    igraph_bool_t pntk) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_vector_long_t degree;
  igraph_vector_t neis;
  long int t;
  igraph_real_t S=0.0;
  igraph_vector_long_t ntk;
  long int ntksize=igraph_vector_size(kernel);
  
  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

  if (pntk) {
    IGRAPH_CHECK(igraph_vector_long_init(&ntk, ntksize));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &ntk);
  }

  if (logprobs) {
    IGRAPH_CHECK(igraph_vector_resize(logprobs, no_of_edges));
  }
  if (logciting) {
    IGRAPH_CHECK(igraph_vector_resize(logciting, no_of_nodes));
    igraph_vector_null(logciting);
  }
  if (logcited) {
    IGRAPH_CHECK(igraph_vector_resize(logcited, no_of_nodes));
    igraph_vector_null(logcited);
  }
  
  for (t=0; t<no_of_nodes; t++) {
    long int n, nneis;
    IGRAPH_CHECK(igraph_incident(graph, &neis, t, IGRAPH_OUT));
    nneis=igraph_vector_size(&neis);    
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    for (n=0; n<nneis; n++) {
      long int edge=VECTOR(neis)[n];
      long int to=IGRAPH_OTHER(graph, edge, t);
      long int x=VECTOR(degree)[to];
      igraph_real_t prob;
      if (pntk) {
	prob=log( VECTOR(ntk)[x] * VECTOR(*kernel)[x] / S );
      } else {
	prob=log( VECTOR(*kernel)[x] / S );
      }
      if (logprobs) {	
	VECTOR(*logprobs)[edge] = prob;
      }
      if (logcited) { 
	VECTOR(*logcited)[to] += prob;
      }
      if (logciting) {
	VECTOR(*logciting)[t] += prob;
      }
    }
    
    for (n=0; n<nneis; n++) {
      long int edge=VECTOR(neis)[n];
      long int to=IGRAPH_OTHER(graph, edge, t);
      long int x=VECTOR(degree)[to];
      
      VECTOR(degree)[to] += 1;
      if (pntk) {
	VECTOR(ntk)[x+1] += 1;
	VECTOR(ntk)[x] -= 1;
      }
      S += VECTOR(*kernel)[x+1];
      S -= VECTOR(*kernel)[x];
    }
    if (pntk) {
      VECTOR(ntk)[0] += 1;
    }
    S += VECTOR(*kernel)[0];
    
  } /* t < no_of_nodes */

  if (pntk) {
    igraph_vector_long_destroy(&ntk);
    IGRAPH_FINALLY_CLEAN(1);
  }

  igraph_vector_destroy(&neis);
  igraph_vector_long_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

int igraph_revolver_probs_ad(const igraph_t *graph,
			     const igraph_matrix_t *kernel,
			     igraph_vector_t *logprobs,
			     igraph_vector_t *logcited,
			     igraph_vector_t *logciting, 
			     igraph_bool_t pntk) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_vector_long_t degree;
  igraph_vector_t neis;
  long int t, j;
  igraph_real_t S=0.0;
  long int agebins=igraph_matrix_ncol(kernel);
  long int binwidth=no_of_nodes/agebins+1;
  igraph_matrix_long_t ntk;
  
  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  
  if (pntk) {
    IGRAPH_CHECK(igraph_matrix_long_init(&ntk, igraph_matrix_nrow(kernel),
					 igraph_matrix_ncol(kernel)));
    IGRAPH_FINALLY(igraph_matrix_long_destroy, &ntk);    
  }

  if (logprobs) {
    IGRAPH_CHECK(igraph_vector_resize(logprobs, no_of_edges));
  }
  if (logcited) {
    IGRAPH_CHECK(igraph_vector_resize(logcited, no_of_nodes));
    igraph_vector_null(logcited);
  }
  if (logciting) {
    IGRAPH_CHECK(igraph_vector_resize(logciting, no_of_nodes));
    igraph_vector_null(logciting);
  }  
  
  for (t=0; t<no_of_nodes; t++) {
    long int n, nneis;
    IGRAPH_CHECK(igraph_incident(graph, &neis, t, IGRAPH_OUT));
    nneis=igraph_vector_size(&neis);
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    for (n=0; n<nneis; n++) {
      long int edge=VECTOR(neis)[n];
      long int to=IGRAPH_OTHER(graph, edge, t);
      long int x=VECTOR(degree)[to];
      long int y=(t-to)/binwidth;
      igraph_real_t prob;
      if (pntk) {
	prob=log( MATRIX(ntk, x, y) * MATRIX(*kernel, x, y) / S );
      } else {
	prob=log( MATRIX(*kernel, x, y) / S );
      }
      if (logprobs) {
	VECTOR(*logprobs)[edge] = prob;
      }
      if (logcited) {
	VECTOR(*logcited)[to] += prob;
      }
      if (logciting) {
	VECTOR(*logciting)[t] += prob;
      }
    }
    
    for (n=0; n<nneis; n++) {
      long int edge=VECTOR(neis)[n];
      long int to=IGRAPH_OTHER(graph, edge, t);
      long int x=VECTOR(degree)[to];
      long int y=(t-to)/binwidth;
      
      VECTOR(degree)[to] += 1;
      if (pntk) {
	MATRIX(ntk, x+1, y) += 1;
	MATRIX(ntk, x, y) -= 1;
      }
      S += MATRIX(*kernel, x+1, y);
      S -= MATRIX(*kernel, x, y);
    }
    
    for (j=1; t-binwidth*j+1>=0; j++) {
      long int shnode=t-binwidth*j+1;
      long int deg=VECTOR(degree)[shnode];
      if (pntk) {
	MATRIX(ntk, deg, j) += 1;
	MATRIX(ntk, deg, j-1) -= 1;
      }
      S += MATRIX(*kernel, deg, j);
      S -= MATRIX(*kernel, deg, j-1);
    }
    if (pntk) {
      MATRIX(ntk, 0, 0) += 1;
    }
    S += MATRIX(*kernel, 0, 0);
    
  } /* t < no_of_nodes */

  if (pntk) {
    igraph_matrix_long_destroy(&ntk);
    IGRAPH_FINALLY_CLEAN(1);
  }

  igraph_vector_destroy(&neis);
  igraph_vector_long_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

int igraph_revolver_probs_de(const igraph_t *graph,
			     const igraph_matrix_t *kernel,
			     const igraph_vector_t *cats,
			     igraph_vector_t *logprobs,
			     igraph_vector_t *logcited,
			     igraph_vector_t *logciting) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_vector_long_t degree;
  igraph_vector_t neis;
  long int t;
  igraph_real_t S=0.0;
  
  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  
  if (logprobs) {
    IGRAPH_CHECK(igraph_vector_resize(logprobs, no_of_edges));
  }
  if (logcited) {
    IGRAPH_CHECK(igraph_vector_resize(logcited, no_of_nodes));
    igraph_vector_null(logcited);
  }
  if (logciting) {
    IGRAPH_CHECK(igraph_vector_resize(logciting, no_of_nodes));
    igraph_vector_null(logciting);
  }
  
  for (t=0; t<no_of_nodes; t++) {
    long int n, nneis;
    IGRAPH_CHECK(igraph_incident(graph, &neis, t, IGRAPH_OUT));
    nneis=igraph_vector_size(&neis);
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    for (n=0; n<nneis; n++) {
      long int edge=VECTOR(neis)[n];
      long int to=IGRAPH_OTHER(graph, edge, t);
      long int x=VECTOR(*cats)[to];
      long int y=VECTOR(degree)[to];
      igraph_real_t prob=log( MATRIX(*kernel, x, y) / S );
      if (logprobs) {
	VECTOR(*logprobs)[edge] = prob;
      }
      if (logcited) {
	VECTOR(*logcited)[to] += prob;
      }
      if (logciting) {
	VECTOR(*logciting)[t] += prob;
      }
    }
    
    for (n=0; n<nneis; n++) {
      long int edge=VECTOR(neis)[n];
      long int to=IGRAPH_OTHER(graph, edge, t);
      long int x=VECTOR(*cats)[to];
      long int y=VECTOR(degree)[to];
      
      VECTOR(degree)[to] += 1;
      S += MATRIX(*kernel, x, y+1);
      S -= MATRIX(*kernel, x, y);
    }
    S += MATRIX(*kernel, 0, 0);
    
  } /* t < no_of_nodes */
  
  igraph_vector_destroy(&neis);
  igraph_vector_long_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

int igraph_revolver_probs_ade(const igraph_t *graph,
			      const igraph_array3_t *kernel,
			      const igraph_vector_t *cats,
			      igraph_vector_t *logprobs,
			      igraph_vector_t *logcited,
			      igraph_vector_t *logciting) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_vector_long_t degree;
  igraph_vector_t neis;
  long int j, t;
  igraph_real_t S=0.0;
  long int agebins=igraph_array3_n(kernel,3);
  long int binwidth=no_of_nodes/agebins+1;
  
  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  
  if (logprobs) {
    IGRAPH_CHECK(igraph_vector_resize(logprobs, no_of_edges));
  }
  if (logcited) {
    IGRAPH_CHECK(igraph_vector_resize(logcited, no_of_nodes));
    igraph_vector_null(logcited);
  }
  if (logciting) {
    IGRAPH_CHECK(igraph_vector_resize(logciting, no_of_nodes));
    igraph_vector_null(logciting);
  }
  
  for (t=0; t<no_of_nodes; t++) {
    long int n, nneis;
    IGRAPH_CHECK(igraph_incident(graph, &neis, t, IGRAPH_OUT));
    nneis=igraph_vector_size(&neis);

    IGRAPH_ALLOW_INTERRUPTION();

    for (n=0; n<nneis; n++) {
      long int edge=VECTOR(neis)[n];
      long int to=IGRAPH_OTHER(graph, edge, t);
      long int x=VECTOR(*cats)[to];
      long int y=VECTOR(degree)[to];
      long int z=(t-to)/binwidth;
      igraph_real_t prob=log( ARRAY3(*kernel, x, y, z) / S );
      if (logprobs) {
	VECTOR(*logprobs)[edge] = prob;
      }
      if (logcited) {
	VECTOR(*logcited)[to] += prob;
      }
      if (logciting) {
	VECTOR(*logciting)[t] += prob;
      }
    }
    
    for (n=0; n<nneis; n++) {
      long int edge=VECTOR(neis)[n];
      long int to=IGRAPH_OTHER(graph, edge, t);
      long int x=VECTOR(*cats)[to];
      long int y=VECTOR(degree)[to];
      long int z=(t-to)/binwidth;
      
      VECTOR(degree)[to] += 1;
      S += ARRAY3(*kernel, x, y+1, z);
      S -= ARRAY3(*kernel, x, y, z);
    }
    
    for (j=1; t-binwidth*j+1>=0; j++) {
      long int shnode=t-binwidth*j+1;
      long int cat=VECTOR(*cats)[shnode];
      long int deg=VECTOR(degree)[shnode];
      S += ARRAY3(*kernel, cat, deg, j);
      S -= ARRAY3(*kernel, cat, deg, j-1);
    }
    S += ARRAY3(*kernel, (long int) VECTOR(*cats)[t], 0, 0);
    
  } /* t < no_of_nodes */

  igraph_vector_destroy(&neis);
  igraph_vector_long_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

int igraph_revolver_probs_ADE(const igraph_t *graph,
			      igraph_scalar_function_t *A_fun,
			      const igraph_matrix_t *par,
			      const igraph_vector_t *cats,
			      const igraph_vector_t *gcats,
			      int agebins,
			      igraph_vector_t *logprobs,
			      igraph_vector_t *logcited,
			      igraph_vector_t *logciting) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_vector_long_t degree;
  igraph_vector_t neis;
  igraph_vector_t S;
  igraph_vector_t gpar;
  igraph_vector_t var;
  int parlen=igraph_matrix_nrow(par);
  int no_gcats=igraph_matrix_ncol(par);
  long int t, i, j;
  
  long int binwidth=no_of_nodes/agebins+1;
  
  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&S, no_gcats);
  IGRAPH_VECTOR_INIT_FINALLY(&var, 3);
  
  if (logprobs) {
    IGRAPH_CHECK(igraph_vector_resize(logprobs, no_of_edges));
  }
  if (logcited) {
    IGRAPH_CHECK(igraph_vector_resize(logcited, no_of_nodes));
    igraph_vector_null(logcited);
  }
  if (logciting) {
    IGRAPH_CHECK(igraph_vector_resize(logciting, no_of_nodes));
    igraph_vector_null(logciting);
  }
  
  for (t=0; t<no_of_nodes; t++) {
    long int n, nneis;
    long int tcat=VECTOR(*gcats)[t];
    igraph_vector_view(&gpar, &MATRIX(*par,0,tcat), parlen);
    IGRAPH_CHECK(igraph_incident(graph, &neis, t, IGRAPH_OUT));
    nneis=igraph_vector_size(&neis);
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    for (n=0; n<nneis; n++) {
      long int edge=VECTOR(neis)[n];
      long int to=IGRAPH_OTHER(graph, edge, t);
      igraph_real_t prob;
      VECTOR(var)[0] = VECTOR(*cats)[to];
      VECTOR(var)[1] = VECTOR(degree)[to];
      VECTOR(var)[2] = (t-to)/binwidth;
      prob=log( A_fun(&var, &gpar, 0) / VECTOR(S)[tcat] );
      if (logprobs) {
	VECTOR(*logprobs)[edge] = prob;
      } 
      if (logcited) {
	VECTOR(*logcited)[to] += prob;
      }
      if (logciting) {
	VECTOR(*logciting)[t] += prob;
      }
    }
      
    for (n=0; n<nneis; n++) {
      long int edge=VECTOR(neis)[n];
      long int to=IGRAPH_OTHER(graph, edge, t);
      VECTOR(var)[0] = VECTOR(*cats)[to];
      VECTOR(var)[1] = VECTOR(degree)[to];
      VECTOR(var)[2] = (t-to)/binwidth;
      
      VECTOR(degree)[to] += 1;
      for (i=0; i<no_gcats; i++) {
	igraph_vector_view(&gpar, &MATRIX(*par,0,i), parlen);
	VECTOR(S)[i] -= A_fun(&var, &gpar, 0);
      }
      VECTOR(var)[1] += 1;
      for (i=0; i<no_gcats; i++) {
	igraph_vector_view(&gpar, &MATRIX(*par,0,i), parlen);
	VECTOR(S)[i] += A_fun(&var, &gpar, 0);
      }
    }
    
    for (j=1; t-binwidth*j+1>=0; j++) {
      long int shnode=t-binwidth*j+1;
      VECTOR(var)[0] = VECTOR(*cats)[shnode];
      VECTOR(var)[1] = VECTOR(degree)[shnode];
      VECTOR(var)[2] = j;
      for (i=0; i<no_gcats; i++) {
	igraph_vector_view(&gpar, &MATRIX(*par,0,i), parlen);
	VECTOR(S)[i] += A_fun(&var, &gpar, 0);
      }
      VECTOR(var)[2] = j-1;
      for (i=0; i<no_gcats; i++) {
	igraph_vector_view(&gpar, &MATRIX(*par,0,i), parlen);
	VECTOR(S)[i] += A_fun(&var, &gpar, 0);
      }
    }
    VECTOR(var)[0]=VECTOR(*cats)[t];
    VECTOR(var)[1]=0;
    VECTOR(var)[2]=0;
    for (i=0; i<no_gcats; i++) {
      igraph_vector_view(&gpar, &MATRIX(*par,0,i), parlen);
      VECTOR(S)[i] += A_fun(&var, &gpar, 0);
    }
  
  } /* t<no_of_nodes */
  
  igraph_vector_destroy(&var);
  igraph_vector_destroy(&S);
  igraph_vector_destroy(&neis);
  igraph_vector_long_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(4);
  return 0;
}
 
int igraph_revolver_probs_ADE_dpareto(const igraph_t *graph,
				      const igraph_matrix_t *par,
				      const igraph_vector_t *cats,
				      const igraph_vector_t *gcats,
				      int agebins,
				      igraph_vector_t *logprobs,
				      igraph_vector_t *logcited,
				      igraph_vector_t *logciting) {
  
  return igraph_revolver_probs_ADE(graph, igraph_i_revolver_ml_ADE_dpareto_f,
				   par, cats, gcats, agebins, logprobs, 
				   logcited, logciting);
}
