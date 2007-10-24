/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
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
} igraph_i_revolver_ml_D_data_t;

/* Evaluate the objective function and calculate its gradient too. */

int igraph_i_revolver_ml_D_eval(const igraph_vector_t *par,
				igraph_i_revolver_ml_D_data_t *data) {

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
  for (i=0; i<dim; i++) {
    VECTOR(data->par1)[i+1] = VECTOR(*par)[i];
  }
  for (t=0; t<=data->maxdegree; t++) {
    VECTOR(data->par1)[0] = t;
    VECTOR(data->A_vect)[t] = data->A(&data->par1, 0);
    data->dA(&data->par1, &data->tmpgrad, 0);
    for (i=0; i<dim; i++) {
      igraph_vector_t *v=VECTOR(data->dA_vects)[i];
      VECTOR(*v)[t] = VECTOR(data->tmpgrad)[i];
    }
  }
  
  for (t=0; t<data->no_of_nodes; t++) {
    long int n, nneis;
    IGRAPH_CHECK(igraph_neighbors(data->graph, &data->neis, t, IGRAPH_OUT));
    nneis=igraph_vector_size(&data->neis);
    
    IGRAPH_ALLOW_INTERRUPTION();
    
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
  data->lastf=sum;
  /* grad is already udpated */

  return 0;
}

/* This function gives the value of the objective function at the 
   supplied parameter vector. Called by the optimizer.
*/

igraph_real_t igraph_i_revolver_ml_D_f(const igraph_vector_t *par,
				       void* extra) {

  igraph_i_revolver_ml_D_data_t *data=extra;
  
  if (!igraph_vector_is_equal(par, &data->lastparam)) {
    igraph_i_revolver_ml_D_eval(par, data);
  }

  return data->lastf;
}

/* This function gievs the gradient of the objective function at
   the supplied parameter vector. Called by the optimizer.
*/

void igraph_i_revolver_ml_D_df(const igraph_vector_t *par,
			       igraph_vector_t *res, void *extra) {

  igraph_i_revolver_ml_D_data_t *data=extra;

  if (!igraph_vector_is_equal(par, &data->lastparam)) {
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
			 igraph_integer_t *fncount, igraph_integer_t *grcount) {

  igraph_i_revolver_ml_D_data_t info;
  igraph_integer_t maxdegree;
  long int no_of_nodes=igraph_vcount(graph);
  int dim=igraph_vector_size(res);
  int ret, i;

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
  IGRAPH_VECTOR_INIT_FINALLY(&info.par1, dim+1);
  IGRAPH_VECTOR_INIT_FINALLY(&info.tmpgrad, dim);
  IGRAPH_VECTOR_INIT_FINALLY(&info.lastparam, dim);
  info.lastf=0.0;
  IGRAPH_VECTOR_INIT_FINALLY(&info.lastgrad, dim);

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

igraph_real_t igraph_i_revolver_ml_D_alpha_f(const igraph_vector_t *par, 
					     void *extra) {
  igraph_real_t deg=VECTOR(*par)[0];
  igraph_real_t alpha=VECTOR(*par)[1];
  if (deg != 0) {
    return pow(deg, alpha) + 1.0;
  } else {
    return 1.0;
  }
}

void igraph_i_revolver_ml_D_alpha_df(const igraph_vector_t *par,
				     igraph_vector_t *res, 
				     void *extra) {
  igraph_real_t deg=VECTOR(*par)[0];
  igraph_real_t alpha=VECTOR(*par)[1];
  if (deg != 0) {
    VECTOR(*res)[0] = log(deg) * pow(deg, alpha);
  } else {
    VECTOR(*res)[0] = 0.0;
  }
}

int igraph_revolver_ml_D_alpha(const igraph_t *graph,
			       igraph_real_t *alpha, igraph_real_t *Fmin,
			       igraph_real_t abstol, igraph_real_t reltol, 
			       int maxit, igraph_integer_t *fncount, 
			       igraph_integer_t *grcount) {
  
  igraph_vector_t res;
  int ret;

  IGRAPH_VECTOR_INIT_FINALLY(&res, 1);
  VECTOR(res)[0]=*alpha;
  
  ret=igraph_revolver_ml_D(graph, &res, Fmin, abstol, reltol, maxit,
			   igraph_i_revolver_ml_D_alpha_f,
			   igraph_i_revolver_ml_D_alpha_df,
			   fncount, grcount);

  *alpha=VECTOR(res)[0];
  igraph_vector_destroy(&res);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

/* These functions assemble the A(d)=d^alpha+a kernel function and 
   calls the general degree-based optimizer.
*/

igraph_real_t igraph_i_revolver_ml_D_alpha_a_f(const igraph_vector_t *par,
					       void *extra) {
  igraph_real_t deg=VECTOR(*par)[0];
  igraph_real_t alpha=VECTOR(*par)[1];
  igraph_real_t a=VECTOR(*par)[2];
  if (deg != 0) {
    return pow(deg, alpha) + a;
  } else {
    return a;
  }  
}

void igraph_i_revolver_ml_D_alpha_a_df(const igraph_vector_t *par,
				       igraph_vector_t *res,
				       void *extra) {
  igraph_real_t deg=VECTOR(*par)[0];
  igraph_real_t alpha=VECTOR(*par)[1];
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
				 int maxit, igraph_integer_t *fncount, 
				 igraph_integer_t *grcount) {
  igraph_vector_t res;
  int ret;
  
  IGRAPH_VECTOR_INIT_FINALLY(&res, 2);
  VECTOR(res)[0] = *alpha;
  VECTOR(res)[1] = *a;
  
  ret=igraph_revolver_ml_D(graph, &res, Fmin, abstol, reltol, maxit,
			   igraph_i_revolver_ml_D_alpha_a_f,
			   igraph_i_revolver_ml_D_alpha_a_df, 
			   fncount, grcount);
  
  *alpha=VECTOR(res)[0];
  *a=VECTOR(res)[1];
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
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
} igraph_i_revolver_ml_AD_data_t;

#define CHECK_VALID(X,Y) \
  if (!MATRIX(data->A_valid,(X),(Y))) { \
     int i; \
     VECTOR(data->par1)[0]=(X); VECTOR(data->par1)[1]=(Y); \
     MATRIX(data->A_vect, (X), (Y)) = data->A(&data->par1, 0); \
     data->dA(&data->par1, &data->tmpgrad, 0); \
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
  
  /* Init */
  igraph_vector_long_null(&data->degree);
  igraph_vector_null(&data->dS);
  igraph_vector_null(grad);
  igraph_matrix_bool_null(&data->A_valid);

  for (i=0; i<dim; i++) {
    VECTOR(data->par1)[i+2] = VECTOR(*par)[i];
  }

  for (t=0; t<data->no_of_nodes; t++) {
    long int n, nneis;
    IGRAPH_CHECK(igraph_neighbors(data->graph, &data->neis, t, IGRAPH_OUT));
    nneis=igraph_vector_size(&data->neis);

    IGRAPH_ALLOW_INTERRUPTION();
    
    /* Update sum(s) */
    for (n=0; n<nneis; n++) {
      long int to=VECTOR(data->neis)[n];
      long int x=VECTOR(data->degree)[to];
      long int y=(t-to)/binwidth;

      CHECK_VALID(x,y);
      sum -= log( MATRIX(data->A_vect, x, y) );
      sum += log( S );
      for (i=0; i<dim; i++) {
	igraph_matrix_t *m=VECTOR(data->dA_vects)[i];
	VECTOR(*grad)[i] -= MATRIX(*m, x, y) / MATRIX(data->A_vect, x, y);
	VECTOR(*grad)[i] += VECTOR(data->dS)[i] / S;
      }
    }
    
    /* Update S, data->dS */
    for (n=0; n<nneis; n++) {
      long int to=VECTOR(data->neis)[n];
      long int x=VECTOR(data->degree)[to];
      long int y=(t-to)/binwidth;
      
      CHECK_VALID(x+1,y);	/* (x,y) already checked */
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
    CHECK_VALID(0,0);
    S += MATRIX(data->A_vect, 0, 0);
    for (i=0; i<dim; i++) {
      igraph_matrix_t *m=VECTOR(data->dA_vects)[i];
      VECTOR(data->dS)[i] += MATRIX(*m, 0, 0);
    }
    /* Aging */
    for (j=1; t-binwidth*j+1>=0; j++) {
      long int shnode=t-binwidth*j+1;
      long int deg=VECTOR(data->degree)[shnode];
      CHECK_VALID(deg, j-1);
      CHECK_VALID(deg, j);
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
  data->lastf=sum;
  /* grad is already updated */

  return 0.0;
}

igraph_real_t igraph_i_revolver_ml_AD_f(const igraph_vector_t *par,
					void *extra) {

  igraph_i_revolver_ml_AD_data_t *data=extra;
  
  if (!igraph_vector_is_equal(par, &data->lastparam)) {
    igraph_i_revolver_ml_AD_eval(par, data);
  }
  
  return data->lastf;
}

void igraph_i_revolver_ml_AD_df(const igraph_vector_t *par,
				igraph_vector_t *res, void *extra) {

  igraph_i_revolver_ml_AD_data_t *data=extra;
  
  if (!igraph_vector_is_equal(par, &data->lastparam)) {
    igraph_i_revolver_ml_AD_eval(par, data);
  }
  
  igraph_vector_update(res, &data->lastgrad);
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
			  int agebins, igraph_integer_t *fncount, 
			  igraph_integer_t *grcount) {
  
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
  IGRAPH_FINALLY(igraph_i_revolver_ml_D_free, &info.dA_vects);
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
  IGRAPH_VECTOR_INIT_FINALLY(&info.par1, dim+2);
  IGRAPH_VECTOR_INIT_FINALLY(&info.tmpgrad, dim);
  info.agebins=agebins;
  IGRAPH_VECTOR_INIT_FINALLY(&info.lastparam, dim);
  info.lastf=0.0;
  IGRAPH_VECTOR_INIT_FINALLY(&info.lastgrad, dim);  
  
  igraph_i_revolver_ml_AD_eval(res, &info);
  ret=igraph_bfgs(res, Fmin, igraph_i_revolver_ml_AD_f,
		  igraph_i_revolver_ml_AD_df, maxit, 1, abstol, reltol, 1, 
		  &info, fncount, grcount);
  
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

  return 0;
}

igraph_real_t igraph_i_revolver_ml_AD_alpha_a_beta_f(const igraph_vector_t *par,
						     void *extra) {
  igraph_real_t deg=VECTOR(*par)[0];
  igraph_real_t age=VECTOR(*par)[1]+1;
  igraph_real_t alpha=VECTOR(*par)[2];
  igraph_real_t a=VECTOR(*par)[3];
  igraph_real_t beta=VECTOR(*par)[4];
  return (pow(deg, alpha) + a) * pow(age, -beta);
}

void igraph_i_revolver_ml_AD_alpha_a_beta_df(const igraph_vector_t *par,
					    igraph_vector_t *res,
					    void *extra) {
  igraph_real_t deg=VECTOR(*par)[0];
  igraph_real_t age=VECTOR(*par)[1]+1;
  igraph_real_t alpha=VECTOR(*par)[2];
  igraph_real_t a=VECTOR(*par)[3];
  igraph_real_t beta=VECTOR(*par)[4];
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
			    agebins, fncount, grcount);
  
  *alpha=VECTOR(res)[0];
  *a=VECTOR(res)[1];
  *beta=VECTOR(res)[2];
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/*------------------------------------------------------------------*/

/* TODO: delta, logprob, lognull */

int igraph_revolver_ml_d(const igraph_t *graph,
			 igraph_integer_t niter,
			 igraph_vector_t *kernel,
			 igraph_vector_t *cites) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_integer_t imaxdegree;
  long int maxdegree, actmaxdegree;
  long int it, t, i;
  igraph_vector_long_t ptk;
  igraph_vector_t *mycites, vmycites;
  igraph_vector_t neis;
  igraph_vector_long_t degree;
  igraph_real_t S=0;

  igraph_vector_t vmykernel;
  igraph_vector_t *kernels[]={ kernel, &vmykernel };
  long int actkernel=0;
  igraph_vector_t *fromkernel=kernels[actkernel], 
    *tokernel=kernels[1-actkernel];
  
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

    for (t=0; t<no_of_nodes; t++) {
      long int n, nneis;
      IGRAPH_CHECK(igraph_neighbors(graph, &neis, t, IGRAPH_OUT));
      nneis=igraph_vector_size(&neis);

      IGRAPH_ALLOW_INTERRUPTION();      

      /* Calculate some terms of the sum for the non-zero classes */
      if (S != 0) {
	for (i=0; i<=actmaxdegree; i++) {
	  VECTOR(*tokernel)[i] += nneis * VECTOR(ptk)[i] / S;
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
      
	if (it==0) {
	  VECTOR(*mycites)[x] += 1;
	}
      }
      VECTOR(ptk)[0] += 1;
      S += VECTOR(*fromkernel)[0];

    } /* t<no_of_nodes  */
    
    /* final step, Mk/sum */
    for (i=0; i<maxdegree; i++) {
      VECTOR(*tokernel)[i] = VECTOR(*mycites)[i] / VECTOR(*tokernel)[i];      
    }
    
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
			 
