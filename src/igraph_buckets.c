/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
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

#include "types.h"
#include "config.h"

#include <stdio.h>

int igraph_buckets_init(igraph_buckets_t *b, long int bsize, long int size) {
  IGRAPH_VECTOR_INIT_FINALLY(&b->bptr, bsize);
  IGRAPH_VECTOR_INIT_FINALLY(&b->buckets, size);
  b->max=-1; b->no=0;
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

void igraph_buckets_destroy(igraph_buckets_t *b) {
  igraph_vector_destroy(&b->bptr);
  igraph_vector_destroy(&b->buckets);
}

long int igraph_buckets_popmax(igraph_buckets_t *b) {
  /* Precondition: there is at least a non-empty bucket */
  /* Search for the highest bucket first */
  long int max;
  while ( (max=VECTOR(b->bptr)[(long int) b->max]) == 0) {
    b->max --;
  }
  VECTOR(b->bptr)[(long int) b->max] = VECTOR(b->buckets)[max-1];
  b->no--;
  
  return max-1;
}

igraph_bool_t igraph_buckets_empty(const igraph_buckets_t *b) {
  return (b->no == 0);
}

void igraph_buckets_add(igraph_buckets_t *b, long int bucket,
			igraph_real_t elem) {
  
  VECTOR(b->buckets)[(long int) elem] = VECTOR(b->bptr)[(long int) bucket];
  VECTOR(b->bptr)[(long int) bucket] = elem+1;
  if (bucket > b->max) {
    b->max = bucket;
  }
  b->no++;
}
