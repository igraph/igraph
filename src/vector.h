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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

/** 
 * Vector, dealing with arrays efficiently.
 * \ingroup types
 */

typedef struct TYPE(igraph_vector) {
  BASE* stor_begin;
  BASE* stor_end;
  BASE* end;
} TYPE(igraph_vector);

#ifndef IGRAPH_VECTOR_NULL
#define IGRAPH_VECTOR_NULL { 0,0,0 }
#endif
#ifndef IGRAPH_VECTOR_INIT_FINALLY
#define IGRAPH_VECTOR_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(igraph_vector_init(v, size)); \
  IGRAPH_FINALLY(igraph_vector_destroy, v); } while (0)
#endif

/*--------------------*/
/* Allocation         */
/*--------------------*/
 
int FUNCTION(igraph_vector,init)(TYPE(igraph_vector)* v, long int size);
int FUNCTION(igraph_vector,init_copy)(TYPE(igraph_vector)* v, 
				       BASE* data, long int length);
int FUNCTION(igraph_vector,init_seq)(TYPE(igraph_vector)*v, BASE from, BASE to);
int FUNCTION(igraph_vector,copy)(TYPE(igraph_vector) *to, 
				 const TYPE(igraph_vector) *from);
void FUNCTION(igraph_vector,destroy)(TYPE(igraph_vector)* v);

/*--------------------*/
/* Accessing elements */
/*--------------------*/

#ifndef VECTOR
/**
 * \ingroup vector
 * \define VECTOR
 * \brief Accessing an element of a vector.
 * 
 * Usage: 
 * \verbatim VECTOR(v)[0] \endverbatim 
 * to access the first element of the vector, you can also use this in
 * assignments, like: 
 * \verbatim VECTOR(v)[10]=5; \endverbatim
 *
 * Note that there are no range checks right now.
 * This functionality might be redefined later as a real function
 * instead of a <code>#define</code>. 
 * \param v The vector object.
 * 
 * Time complexity: O(1).
 */
#define VECTOR(v) ((v).stor_begin) 
#endif

BASE FUNCTION(igraph_vector,e)(const TYPE(igraph_vector)* v, long int pos);
BASE* FUNCTION(igraph_vector,e_ptr)(const TYPE(igraph_vector)* v, long int pos);
void FUNCTION(igraph_vector,set)(TYPE(igraph_vector)* v, long int pos, BASE value);
BASE FUNCTION(igraph_vector,tail)(const TYPE(igraph_vector) *v);

/*-----------------------*/
/* Initializing elements */
/*-----------------------*/

void FUNCTION(igraph_vector,null)(TYPE(igraph_vector)* v);
void FUNCTION(igraph_vector,fill)(TYPE(igraph_vector)* v, BASE e);

/*-----------------------*/
/* Vector views          */
/*-----------------------*/

const TYPE(igraph_vector) *FUNCTION(igraph_vector,view)(const TYPE(igraph_vector) *v,
							const BASE *data, 
							long int length);

/*-----------------------*/
/* Copying vectors       */
/*-----------------------*/

void FUNCTION(igraph_vector,copy_to)(const TYPE(igraph_vector) *v, BASE* to);
int FUNCTION(igraph_vector,update)(TYPE(igraph_vector) *to, 
				   const TYPE(igraph_vector) *from);
int FUNCTION(igraph_vector,append)(TYPE(igraph_vector) *to, 
				   const TYPE(igraph_vector) *from);
int FUNCTION(igraph_vector,swap)(TYPE(igraph_vector) *v1, TYPE(igraph_vector) *v2);

/*-----------------------*/
/* Exchanging elements   */
/*-----------------------*/

int FUNCTION(igraph_vector,swap_elements)(TYPE(igraph_vector) *v,
					  long int i, long int j);
int FUNCTION(igraph_vector,reverse)(TYPE(igraph_vector) *v);
int FUNCTION(igraph_vector,shuffle)(TYPE(igraph_vector) *v);

/*-----------------------*/
/* Vector operations     */
/*-----------------------*/

void FUNCTION(igraph_vector,add_constant)(TYPE(igraph_vector) *v, BASE plus);
void FUNCTION(igraph_vector,scale)(TYPE(igraph_vector) *v, BASE by);
int FUNCTION(igraph_vector,add)(TYPE(igraph_vector) *v1, 
				const TYPE(igraph_vector) *v2);
int FUNCTION(igraph_vector,sub)(TYPE(igraph_vector) *v1, 
				const TYPE(igraph_vector) *v2);
int FUNCTION(igraph_vector,mul)(TYPE(igraph_vector) *v1, 
				const TYPE(igraph_vector) *v2);
int FUNCTION(igraph_vector,div)(TYPE(igraph_vector) *v1, 
				const TYPE(igraph_vector) *v2);

/*------------------------------*/
/* Finding minimum and maximum  */
/*------------------------------*/

BASE FUNCTION(igraph_vector,min)(const TYPE(igraph_vector)* v);
BASE FUNCTION(igraph_vector,max)(const TYPE(igraph_vector)* v);
long int FUNCTION(igraph_vector,which_min)(const TYPE(igraph_vector)* v);
long int FUNCTION(igraph_vector,which_max)(const TYPE(igraph_vector)* v);
int FUNCTION(igraph_vector,minmax)(const TYPE(igraph_vector) *v,
				   BASE *min, BASE *max);
int FUNCTION(igraph_vector,which_minmax)(const TYPE(igraph_vector) *v,
					 long int *which_min, long int *which_max);

/*-------------------*/
/* Vector properties */
/*-------------------*/

igraph_bool_t FUNCTION(igraph_vector,empty)     (const TYPE(igraph_vector)* v);
long int FUNCTION(igraph_vector,size)      (const TYPE(igraph_vector)* v);
igraph_bool_t FUNCTION(igraph_vector,isnull)(const TYPE(igraph_vector) *v);
BASE FUNCTION(igraph_vector,sum)(const TYPE(igraph_vector) *v);
BASE FUNCTION(igraph_vector,prod)(const TYPE(igraph_vector) *v);
igraph_bool_t FUNCTION(igraph_vector,isininterval)(const TYPE(igraph_vector) *v, 
						   BASE low, BASE high);
igraph_bool_t FUNCTION(igraph_vector,any_smaller)(const TYPE(igraph_vector) *v, 
						  BASE limit);
igraph_bool_t FUNCTION(igraph_vector,is_equal)(const TYPE(igraph_vector) *lhs, 
					       const TYPE(igraph_vector) *rhs);
BASE FUNCTION(igraph_vector,maxdifference)(const TYPE(igraph_vector) *m1,
					   const TYPE(igraph_vector) *m2);

/*------------------------*/
/* Searching for elements */
/*------------------------*/

igraph_bool_t FUNCTION(igraph_vector,contains)(const TYPE(igraph_vector) *v, BASE e);
igraph_bool_t FUNCTION(igraph_vector,search)(const TYPE(igraph_vector) *v,
					     long int from, BASE what, 
					     long int *pos);
igraph_bool_t FUNCTION(igraph_vector,binsearch)(const TYPE(igraph_vector) *v, 
						BASE what, long int *pos);
igraph_bool_t FUNCTION(igraph_vector,binsearch2)(const TYPE(igraph_vector) *v,
						 BASE what);

/*------------------------*/
/* Resizing operations    */
/*------------------------*/

void FUNCTION(igraph_vector,clear)(TYPE(igraph_vector)* v);
int FUNCTION(igraph_vector,resize)(TYPE(igraph_vector)* v, long int newsize);
int FUNCTION(igraph_vector,reserve)(TYPE(igraph_vector)* v, long int size);
int FUNCTION(igraph_vector,push_back)(TYPE(igraph_vector)* v, BASE e);
BASE FUNCTION(igraph_vector,pop_back)(TYPE(igraph_vector)* v);
int FUNCTION(igraph_vector,insert)(TYPE(igraph_vector) *v, long int pos, BASE value);
void FUNCTION(igraph_vector,remove)(TYPE(igraph_vector) *v, long int elem);
void FUNCTION(igraph_vector,remove_section)(TYPE(igraph_vector) *v, 
					    long int from, long int to);

/*-----------*/
/* Sorting   */             
/*-----------*/

void FUNCTION(igraph_vector,sort)(TYPE(igraph_vector) *v);

/* ----------------------------------------------------------------------------*/
/* For internal use only, may be removed, rewritten ... */
/* ----------------------------------------------------------------------------*/

int FUNCTION(igraph_vector,init_real)(TYPE(igraph_vector)*v, int no, ...);
int FUNCTION(igraph_vector,init_int)(TYPE(igraph_vector)*v, int no, ...);
int FUNCTION(igraph_vector,init_real_end)(TYPE(igraph_vector)*v, BASE endmark, ...);
int FUNCTION(igraph_vector,init_int_end)(TYPE(igraph_vector)*v, int endmark, ...);

int FUNCTION(igraph_vector,move_interval)(TYPE(igraph_vector) *v, 
					  long int begin, long int end, long int to);
int FUNCTION(igraph_vector,move_interval2)(TYPE(igraph_vector) *v, 
					  long int begin, long int end, long int to);
void FUNCTION(igraph_vector,permdelete)(TYPE(igraph_vector) *v, 
					const igraph_vector_t *index, 
					long int nremove);
void FUNCTION(igraph_vector,remove_negidx)(TYPE(igraph_vector) *v, 
					   const igraph_vector_t *neg, 
					   long int nremove);
int FUNCTION(igraph_vector,filter_smaller)(TYPE(igraph_vector) *v, BASE elem);
int FUNCTION(igraph_vector,get_interval)(const TYPE(igraph_vector) *v, 
					 TYPE(igraph_vector) *res,
					 long int from, long int to);
int FUNCTION(igraph_vector,intersect_sorted)(const TYPE(igraph_vector) *v1,
  const TYPE(igraph_vector) *v2, TYPE(igraph_vector) *result,
  igraph_bool_t keep_multiplicity);

