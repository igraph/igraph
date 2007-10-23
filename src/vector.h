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
 
int FUNCTION(igraph_vector,init)      (TYPE(igraph_vector)* v, long int size);
int FUNCTION(igraph_vector,init_copy) (TYPE(igraph_vector)* v, BASE* data, long int length);
int FUNCTION(igraph_vector,init_seq)(TYPE(igraph_vector)*v, BASE from, BASE to);
int FUNCTION(igraph_vector,init_real)(TYPE(igraph_vector)*v, int no, ...);
int FUNCTION(igraph_vector,init_int)(TYPE(igraph_vector)*v, int no, ...);
int FUNCTION(igraph_vector,init_real_end)(TYPE(igraph_vector)*v, BASE endmark, ...);
int FUNCTION(igraph_vector,init_int_end)(TYPE(igraph_vector)*v, int endmark, ...);
const TYPE(igraph_vector) *FUNCTION(igraph_vector,view) (const TYPE(igraph_vector) *v, const BASE *data, 
			     long int length);
void FUNCTION(igraph_vector,destroy)   (TYPE(igraph_vector)* v);
int FUNCTION(igraph_vector,reserve)   (TYPE(igraph_vector)* v, long int size);
igraph_bool_t FUNCTION(igraph_vector,empty)     (const TYPE(igraph_vector)* v);
long int FUNCTION(igraph_vector,size)      (const TYPE(igraph_vector)* v);
void FUNCTION(igraph_vector,clear)     (TYPE(igraph_vector)* v);
void FUNCTION(igraph_vector,null)      (TYPE(igraph_vector)* v);
void FUNCTION(igraph_vector,fill)      (TYPE(igraph_vector)* v, BASE e);
int FUNCTION(igraph_vector,push_back) (TYPE(igraph_vector)* v, BASE e);
int FUNCTION(igraph_vector,insert)(TYPE(igraph_vector) *v, long int pos, BASE value);
BASE FUNCTION(igraph_vector,e)         (const TYPE(igraph_vector)* v, long int pos);
BASE* FUNCTION(igraph_vector,e_ptr)  (const TYPE(igraph_vector)* v, long int pos);
void FUNCTION(igraph_vector,set)       (TYPE(igraph_vector)* v, long int pos, BASE value);
BASE FUNCTION(igraph_vector,tail)(const TYPE(igraph_vector) *v);
BASE FUNCTION(igraph_vector,pop_back)(TYPE(igraph_vector)* v);
/* TODO: order* is meaningful only for interger vectors */
int FUNCTION(igraph_vector,order)(const TYPE(igraph_vector)* v, const TYPE(igraph_vector) *v2,
				  igraph_vector_t* res, BASE maxval);
int FUNCTION(igraph_vector,order1)(const TYPE(igraph_vector)* v, 
				   igraph_vector_t* res, BASE maxval);
void FUNCTION(igraph_vector,sort)(TYPE(igraph_vector) *v);
int FUNCTION(igraph_vector,resize)(TYPE(igraph_vector)* v, long int newsize);
BASE FUNCTION(igraph_vector,max)(const TYPE(igraph_vector)* v);
long int FUNCTION(igraph_vector,which_max)(const TYPE(igraph_vector)* v);
BASE FUNCTION(igraph_vector,min)(const TYPE(igraph_vector)* v);
long int FUNCTION(igraph_vector,which_min)(const TYPE(igraph_vector)* v);
void FUNCTION(igraph_vector,copy_to)(const TYPE(igraph_vector) *v, BASE* to);
int FUNCTION(igraph_vector,copy)(TYPE(igraph_vector) *to, const TYPE(igraph_vector) *from);
BASE FUNCTION(igraph_vector,sum)(const TYPE(igraph_vector) *v);
BASE FUNCTION(igraph_vector,prod)(const TYPE(igraph_vector) *v);
void FUNCTION(igraph_vector,remove_section)(TYPE(igraph_vector) *v, long int from, long int to);
int FUNCTION(igraph_vector,move_interval)(TYPE(igraph_vector) *v, long int begin, long int end, 
					  long int to);
void FUNCTION(igraph_vector,remove)(TYPE(igraph_vector) *v, long int elem);
void FUNCTION(igraph_vector,permdelete)(TYPE(igraph_vector) *v, const igraph_vector_t *index, long int nremove);
void FUNCTION(igraph_vector,remove_negidx)(TYPE(igraph_vector) *v, const igraph_vector_t *neg, long int nremove);
igraph_bool_t FUNCTION(igraph_vector,isininterval)(const TYPE(igraph_vector) *v, BASE low, BASE high);
igraph_bool_t FUNCTION(igraph_vector,any_smaller)(const TYPE(igraph_vector) *v, BASE limit);
igraph_bool_t FUNCTION(igraph_vector,is_equal)(const TYPE(igraph_vector) *lhs, const TYPE(igraph_vector) *rhs);
igraph_bool_t FUNCTION(igraph_vector,binsearch)(const TYPE(igraph_vector) *v, BASE what, long int *pos);
igraph_bool_t FUNCTION(igraph_vector,binsearch2)(const TYPE(igraph_vector) *v, BASE what);
void FUNCTION(igraph_vector,add_constant)(TYPE(igraph_vector) *v, BASE plus);
void FUNCTION(igraph_vector,scale)(TYPE(igraph_vector) *v, BASE by);
igraph_bool_t FUNCTION(igraph_vector,contains)(const TYPE(igraph_vector) *v, BASE e);
igraph_bool_t FUNCTION(igraph_vector,search)(const TYPE(igraph_vector) *v, long int from, BASE what, 
					     long int *pos);
int FUNCTION(igraph_vector,filter_smaller)(TYPE(igraph_vector) *v, BASE elem);
int FUNCTION(igraph_vector,append)(TYPE(igraph_vector) *to, const TYPE(igraph_vector) *from);
int FUNCTION(igraph_vector,get_interval)(const TYPE(igraph_vector) *v, TYPE(igraph_vector) *res,
					 long int from, long int to);
/* TODO: only for integers??? */
int FUNCTION(igraph_vector,rank)(const TYPE(igraph_vector) *v, igraph_vector_t *res, 
				 long int nodes);
BASE FUNCTION(igraph_vector,maxdifference)(const TYPE(igraph_vector) *m1,
					   const TYPE(igraph_vector) *m2);
int FUNCTION(igraph_vector,update)(TYPE(igraph_vector) *to, 
				   const TYPE(igraph_vector) *from);
