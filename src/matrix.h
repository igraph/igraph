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

typedef struct TYPE(igraph_matrix) {
  TYPE(igraph_vector) data;
  long int nrow, ncol;
} TYPE(igraph_matrix);

int FUNCTION(igraph_matrix,init)(TYPE(igraph_matrix) *m, long int nrow, long int ncol);
void FUNCTION(igraph_matrix,destroy)(TYPE(igraph_matrix) *m);
int FUNCTION(igraph_matrix,resize)(TYPE(igraph_matrix) *m, long int nrow, long int ncol);
long int FUNCTION(igraph_matrix,size)(const TYPE(igraph_matrix) *m);
long int FUNCTION(igraph_matrix,nrow)(const TYPE(igraph_matrix) *m);
long int FUNCTION(igraph_matrix,ncol)(const TYPE(igraph_matrix) *m);
int FUNCTION(igraph_matrix,copy_to)(const TYPE(igraph_matrix) *m, BASE *to);
int FUNCTION(igraph_matrix,null)(TYPE(igraph_matrix) *m);
int FUNCTION(igraph_matrix,add_cols)(TYPE(igraph_matrix) *m, long int n);
int FUNCTION(igraph_matrix,add_rows)(TYPE(igraph_matrix) *m, long int n);
int FUNCTION(igraph_matrix,remove_col)(TYPE(igraph_matrix) *m, long int col);
int FUNCTION(igraph_matrix,permdelete_rows)(TYPE(igraph_matrix) *m, long int *index, long int nremove);
int FUNCTION(igraph_matrix,delete_rows_neg)(TYPE(igraph_matrix) *m, 
					    const igraph_vector_t *neg, long int nremove);
int FUNCTION(igraph_matrix,copy)(TYPE(igraph_matrix) *to, const TYPE(igraph_matrix) *from);
igraph_real_t FUNCTION(igraph_matrix,max)(const TYPE(igraph_matrix) *m);
void FUNCTION(igraph_matrix,multiply)(TYPE(igraph_matrix) *m, BASE by);
int FUNCTION(igraph_matrix,select_rows)(const TYPE(igraph_matrix) *m, TYPE(igraph_matrix) *res, 
					const igraph_vector_t *rows);
int FUNCTION(igraph_matrix,get_col)(const TYPE(igraph_matrix) *m, TYPE(igraph_vector) *res,
				    long int index);
igraph_real_t FUNCTION(igraph_matrix,sum)(const TYPE(igraph_matrix) *m);
igraph_bool_t FUNCTION(igraph_matrix,is_equal)(const TYPE(igraph_matrix) *m1, 
					       const TYPE(igraph_matrix) *m2);
BASE FUNCTION(igraph_matrix,maxdifference)(const TYPE(igraph_matrix) *m1,
						    const TYPE(igraph_matrix) *m2);

