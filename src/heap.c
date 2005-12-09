/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include "types.h"
#include "memory.h"
#include "random.h"
#include "error.h"

#include <assert.h>
#include <string.h> 		/* memcpy & co. */
#include <stdlib.h>

#define PARENT(x)     (((x)+1)/2-1)
#define LEFTCHILD(x)  (((x)+1)*2-1)
#define RIGHTCHILD(x) (((x)+1)*2)
  
/**
 * \ingroup heap
 * \brief Initializes an empty heap object (constructor).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int heap_init           (heap_t* h, long int alloc_size) {
  if (alloc_size <= 0 ) { alloc_size=1; }
  h->stor_begin=Calloc(alloc_size, real_t);
  if (h->stor_begin==0) {
    IGRAPH_FERROR("heap init failed", IGRAPH_ENOMEM);
  }
  h->stor_end=h->stor_begin + alloc_size;
  h->end=h->stor_begin;
  h->destroy=1;
  
  return 0;  
}

/**
 * \ingroup heap
 * \brief Initializes a heap object from an array, the heap is also
 * built of course (constructor).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int heap_init_array     (heap_t *h, real_t* data, long int len) {
  h->stor_begin=Calloc(len, real_t);
  if (h->stor_begin==0) {
    IGRAPH_FERROR("heap init from array failed", IGRAPH_ENOMEM);
  }
  h->stor_end=h->stor_begin+len;
  h->end=h->stor_end;
  h->destroy=1;

  memcpy(h->stor_begin, data, len*sizeof(real_t));

  heap_i_build (h->stor_begin, h->end-h->stor_begin, 0);
  
  return 0;
}

/**
 * \ingroup heap
 * \brief Destroys an initialized heap object.
 */

void heap_destroy        (heap_t* h) {  
  if (h->destroy) {
    if (h->stor_begin != 0) {
      Free(h->stor_begin);
      h->stor_begin=0;
    }
  }
}

/**
 * \ingroup heap
 * \brief Decides whether a heap object is empty.
 */

bool_t heap_empty          (heap_t* h) {
  assert(h != NULL);
  assert(h->stor_begin != NULL);
  return h->stor_begin == h->end;
}

/**
 * \ingroup heap
 * \brief Adds an element to a heap object.
 */

int heap_push           (heap_t* h, real_t elem) {
  int ret;
  assert(h != NULL);
  assert(h->stor_begin != NULL);
	
  /* full, allocate more storage */
  if (h->stor_end == h->end) {
    long int new_size = heap_size(h) * 2;
    if (new_size == 0) { new_size = 1; }
    IGRAPH_CHECK(ret=heap_reserve(h, new_size));
  }
	
  *(h->end) = elem;
  h->end += 1;

  /* maintain heap */
  heap_i_shift_up(h->stor_begin, heap_size(h), heap_size(h)-1);
	
  return 0;
}

/**
 * \ingroup heap
 * \brief Returns the largest element in a heap.
 */

real_t heap_max       (heap_t* h) {
  assert(h != NULL);
  assert(h->stor_begin != NULL);
  assert(h->stor_begin != h->end);
  
  return h->stor_begin[0];
}

/**
 * \ingroup heap
 * \brief Returns and removes the largest element in a heap.
 */

real_t heap_delete_max(heap_t* h) {
  real_t tmp;

  assert(h != NULL);
  assert(h->stor_begin != NULL);

  tmp=h->stor_begin[0];
  heap_i_switch(h->stor_begin, 0, heap_size(h)-1);
  h->end -= 1;
  heap_i_sink(h->stor_begin, h->end-h->stor_begin, 0);
  
  return tmp;
}

/**
 * \ingroup heap
 * \brief Gives the number of elements in a heap.
 */

long int heap_size      (heap_t* h) {
  assert(h != NULL);
  assert(h->stor_begin != NULL);
  return h->end - h->stor_begin;
}

/**
 * \ingroup heap
 * \brief Allocates more memory for a heap if needed.
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int heap_reserve        (heap_t* h, long int size) {
  long int actual_size=heap_size(h);
  real_t *tmp;
  assert(h != NULL);
  assert(h->stor_begin != NULL);
  
  if (size <= actual_size) { return 0; }
  
  tmp=Realloc(h->stor_begin, size, real_t);
  if (tmp==0) {
    IGRAPH_FERROR("heap reserve failed", IGRAPH_ENOMEM);
  }
  h->stor_begin=tmp;
  h->stor_end=h->stor_begin + size;
  h->end=h->stor_begin+actual_size;
  
  return 0;
}

/**
 * \ingroup heap
 * \brief Build a heap, this should not be called directly.
 */

void heap_i_build(real_t* arr, long int size, long int head) {

  if (RIGHTCHILD(head) < size) { 
    /* both subtrees */
    heap_i_build(arr, size, LEFTCHILD(head) );
    heap_i_build(arr, size, RIGHTCHILD(head));
    heap_i_sink(arr, size, head);
  } else if (LEFTCHILD(head) < size) {
    /* only left */
    heap_i_build(arr, size, LEFTCHILD(head));
    heap_i_sink(arr, size, head);
  } else {
    /* none */
  }
}

/**
 * \ingroup heap
 * \brief Shift an element upwards in a heap, this should not be
 * called directly.
 */

void heap_i_shift_up(real_t* arr, long int size, long int elem) {
  
  if (elem==0 || arr[elem] < arr[PARENT(elem)]) { 
    /* at the top */
  } else {
    heap_i_switch(arr, elem, PARENT(elem));
    heap_i_shift_up(arr, size, PARENT(elem));
  }
}

/**
 * \ingroup heap
 * \brief Moves an element down in a heap, this function should not be
 * called directly.
 */

void heap_i_sink(real_t* arr, long int size, long int head) {

  if (LEFTCHILD(head) >= size) { 
    /* no subtrees */
  } else if (RIGHTCHILD(head) == size ||
	     arr[LEFTCHILD(head)] >= arr[RIGHTCHILD(head)]) {
    /* sink to the left if needed */
    if (arr[head] < arr[LEFTCHILD(head)]) {
      heap_i_switch(arr, head, LEFTCHILD(head));
      heap_i_sink(arr, size, LEFTCHILD(head));
    }
  } else {
    /* sink to the right */
    if (arr[head] < arr[RIGHTCHILD(head)]) {
      heap_i_switch(arr, head, RIGHTCHILD(head));
      heap_i_sink(arr, size, RIGHTCHILD(head));
    }
  }
}

/**
 * \ingroup heap
 * \brief Switches two elements in a heap, this function should not be
 * called directly.
 */

void heap_i_switch(real_t* arr, long int e1, long int e2) {
  if (e1!=e2) {
    real_t tmp=arr[e1];
    arr[e1]=arr[e2];
    arr[e2]=tmp;
  }
}

/**
 * \ingroup indheap
 * \brief Initializes an indexed heap (constructor).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int indheap_init           (indheap_t* h, long int alloc_size) {
 if (alloc_size <= 0 ) { alloc_size=1; }
 h->stor_begin=Calloc(alloc_size, real_t);
 if (h->stor_begin==0) {
   h->index_begin=0;
   IGRAPH_FERROR("indheap init failed", IGRAPH_ENOMEM);
 }
 h->index_begin=Calloc(alloc_size, long int);
 if (h->index_begin==0) {
   Free(h->stor_begin);
   h->stor_begin=0;
   IGRAPH_FERROR("indheap init failed", IGRAPH_ENOMEM);
 }
 
 h->stor_end=h->stor_begin + alloc_size;
 h->end=h->stor_begin;
 h->destroy=1;
 
 return 0;  
}

/**
 * \ingroup indheap
 * \brief Initializes and build an indexed heap from a C array (constructor).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int indheap_init_array     (indheap_t *h, real_t* data, long int len) {
  long int i;

  h->stor_begin=Calloc(len, real_t);
  if (h->stor_begin==0) {
    h->index_begin=0;
    IGRAPH_FERROR("indheap init from array failed", IGRAPH_ENOMEM);
  }
  h->index_begin=Calloc(len, long int);
  if (h->index_begin==0) {
    Free(h->stor_begin);
    h->stor_begin=0;
    IGRAPH_FERROR("indheap init from array failed", IGRAPH_ENOMEM);
  }
  h->stor_end=h->stor_begin+len;
  h->end=h->stor_end;
  h->destroy=1;

  memcpy(h->stor_begin, data, len*sizeof(real_t));
  for (i=0; i<len; i++) {
    h->index_begin[i]=i+1;
  }

  indheap_i_build (h, 0);
  
  return 0;
}

/**
 * \ingroup indheap
 * \brief Destroys an initialized indexed heap. 
 */

void indheap_destroy        (indheap_t* h) {  
  assert(h != 0);
  if (h->destroy) {
    if (h->stor_begin != 0) {
      Free(h->stor_begin);
      h->stor_begin=0;
    }
    if (h->index_begin != 0) {
      Free(h->index_begin);
      h->index_begin=0;
    }
  }
}

/**
 * \ingroup indheap
 * \brief Checks whether a heap is empty.
 */

bool_t indheap_empty          (indheap_t* h) {
  assert(h != 0);
  assert(h->stor_begin != 0);
  return h->stor_begin == h->end;
}

/**
 * \ingroup indheap
 * \brief Adds an element to an indexed heap.
 */

int indheap_push           (indheap_t* h, real_t elem) {
  int ret;
  assert(h != 0);
  assert(h->stor_begin != 0);
	
  /* full, allocate more storage */
  if (h->stor_end == h->end) {
    long int new_size = indheap_size(h) * 2;
    if (new_size == 0) { new_size = 1; }
    IGRAPH_CHECK(indheap_reserve(h, new_size));
  }
	
  *(h->end) = elem;
  h->end += 1;
  *(h->index_begin+indheap_size(h)-1)=indheap_size(h)-1;

  /* maintain indheap */
  indheap_i_shift_up(h, indheap_size(h)-1);
	
  return 0;
}

/**
 * \ingroup indheap
 * \brief Returns the largest element in an indexed heap.
 */

real_t indheap_max       (indheap_t* h) {
  assert(h != NULL);
  assert(h->stor_begin != NULL);
  assert(h->stor_begin != h->end);
  
  return h->stor_begin[0];
}

/**
 * \ingroup indheap
 * \brief Removes the largest element from an indexed heap.
 */

real_t indheap_delete_max(indheap_t* h) {
  real_t tmp;

  assert(h != NULL);
  assert(h->stor_begin != NULL);

  tmp=h->stor_begin[0];
  indheap_i_switch(h, 0, indheap_size(h)-1);
  h->end -= 1;
  indheap_i_sink(h, 0);
  
  return tmp;
}

/**
 * \ingroup indheap
 * \brief Gives the number of elements in an indexed heap.
 */

long int indheap_size      (indheap_t* h) {
  assert(h != 0);
  assert(h->stor_begin != 0);
  return h->end - h->stor_begin;
}

/**
 * \ingroup indheap
 * \brief Reserves more memory for an indexed heap.
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int indheap_reserve        (indheap_t* h, long int size) {
  long int actual_size=indheap_size(h);
  real_t *tmp1;
  long int *tmp2;
  assert(h != 0);
  assert(h->stor_begin != 0);
  
  if (size <= actual_size) { return 0; }

  tmp1=Calloc(size, real_t);  
  if (tmp1==0) {
    IGRAPH_FERROR("indheap reserve failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, tmp1); 	/* TODO: hack */
  tmp2=Calloc(size, long int);
  if (tmp2==0) {
    IGRAPH_FERROR("indheap reserve failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, tmp2);
  memcpy(tmp1, h->stor_begin, actual_size*sizeof(real_t));
  memcpy(tmp2, h->index_begin, actual_size*sizeof(long int));
  Free(h->stor_begin);
  Free(h->index_begin);
  
  h->stor_begin=tmp1;
  h->index_begin=tmp2;
  h->stor_end=h->stor_begin + size;
  h->end=h->stor_begin+actual_size;
  
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

/**
 * \ingroup indheap
 * \brief Returns the index of the largest element in an indexed heap.
 */

long int indheap_max_index(indheap_t *h) {
  assert(h != 0);
  assert(h->stor_begin != 0);
  return h->index_begin[0];  
}

/**
 * \ingroup indheap
 * \brief Builds an indexed heap, this function should not be called
 * directly. 
 */

void indheap_i_build(indheap_t* h, long int head) {

  long int size=indheap_size(h);
  if (RIGHTCHILD(head) < size) { 
    /* both subtrees */
    indheap_i_build(h, LEFTCHILD(head) );
    indheap_i_build(h, RIGHTCHILD(head));
    indheap_i_sink(h, head);
  } else if (LEFTCHILD(head) < size) {
    /* only left */
    indheap_i_build(h, LEFTCHILD(head));
    indheap_i_sink(h, head);
  } else {
    /* none */
  }
}

/**
 * \ingroup indheap
 * \brief Moves an element up in the heap, don't call this function
 * directly. 
 */

void indheap_i_shift_up(indheap_t *h, long int elem) {
  
  if (elem==0 || h->stor_begin[elem] < h->stor_begin[PARENT(elem)]) { 
    /* at the top */
  } else {
    indheap_i_switch(h, elem, PARENT(elem));
    indheap_i_shift_up(h, PARENT(elem));
  }
}

/**
 * \ingroup indheap
 * \brief Moves an element down in the heap, don't call this function
 * directly. 
 */

void indheap_i_sink(indheap_t* h, long int head) {

  long int size=indheap_size(h);
  if (LEFTCHILD(head) >= size) { 
    /* no subtrees */
  } else if (RIGHTCHILD(head) == size ||
	     h->stor_begin[LEFTCHILD(head)]>=h->stor_begin[RIGHTCHILD(head)]) {
    /* sink to the left if needed */
    if (h->stor_begin[head] < h->stor_begin[LEFTCHILD(head)]) {
      indheap_i_switch(h, head, LEFTCHILD(head));
      indheap_i_sink(h, LEFTCHILD(head));
    }
  } else {
    /* sink to the right */
    if (h->stor_begin[head] < h->stor_begin[RIGHTCHILD(head)]) {
      indheap_i_switch(h, head, RIGHTCHILD(head));
      indheap_i_sink(h, RIGHTCHILD(head));
    }
  }
}

/**
 * \ingroup indheap
 * \brief Switches two elements in a heap, don't call this function
 * directly. 
 */

void indheap_i_switch(indheap_t* h, long int e1, long int e2) {
  if (e1!=e2) {
    real_t tmp=h->stor_begin[e1];
    h->stor_begin[e1]=h->stor_begin[e2];
    h->stor_begin[e2]=tmp;
    
    tmp=h->index_begin[e1];
    h->index_begin[e1]=h->index_begin[e2];
    h->index_begin[e2]=tmp;
  }
}


/**
 * \ingroup doubleindheap
 * \brief Initializes an empty doubly indexed heap object (constructor).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int d_indheap_init           (d_indheap_t* h, long int alloc_size) {
 if (alloc_size <= 0 ) { alloc_size=1; }
  h->stor_begin=Calloc(alloc_size, real_t);
  if (h->stor_begin==0) {
    h->index_begin=0;
    h->index2_begin=0;
    IGRAPH_FERROR("d_indheap init failed", IGRAPH_ENOMEM);
  }
  h->stor_end=h->stor_begin + alloc_size;
  h->end=h->stor_begin;
  h->destroy=1;
  h->index_begin=Calloc(alloc_size, long int);
  if (h->index_begin==0) {
    Free(h->stor_begin);
    h->stor_begin=0;
    h->index2_begin=0;
    IGRAPH_FERROR("d_indheap init failed", IGRAPH_ENOMEM);
  }
  h->index2_begin=Calloc(alloc_size, long int);
  if (h->index2_begin==0) {
    Free(h->stor_begin);
    Free(h->index_begin);
    h->stor_begin=0;
    h->index_begin=0;
    IGRAPH_FERROR("d_indheap init failed", IGRAPH_ENOMEM);
  }
  
  return 0;  
}

/**
 * \ingroup doubleindheap
 * \brief Destroys an initialized doubly indexed heap object.
 */

void d_indheap_destroy        (d_indheap_t* h) {  
  assert(h != 0);
  if (h->destroy) {
    if (h->stor_begin != 0) {
      Free(h->stor_begin);
      h->stor_begin=0;
    }
    if (h->index_begin != 0) {
      Free(h->index_begin);
      h->index_begin=0;
    }
    if (h->index2_begin != 0) {
      Free(h->index2_begin);
      h->index2_begin=0;
    }
  }
}

/**
 * \ingroup doubleindheap
 * \brief Decides whether a heap is empty.
 */

bool_t d_indheap_empty          (d_indheap_t* h) {
  assert(h != 0);
  assert(h->stor_begin != 0);
  return h->stor_begin == h->end;
}

/**
 * \ingroup doubleindheap
 * \brief Adds an element to the heap.
 */

int d_indheap_push           (d_indheap_t* h, real_t elem, 
			      long int idx, long int idx2) {
  assert(h != 0);
  assert(h->stor_begin != 0);
	
  /* full, allocate more storage */
  if (h->stor_end == h->end) {
    long int new_size = d_indheap_size(h) * 2;
    if (new_size == 0) { new_size = 1; }
    IGRAPH_CHECK(d_indheap_reserve(h, new_size));
  }
	
  *(h->end) = elem;
  h->end += 1;
  *(h->index_begin+d_indheap_size(h)-1)=idx ;
  *(h->index2_begin+d_indheap_size(h)-1)=idx2 ;

  /* maintain d_indheap */
  d_indheap_i_shift_up(h, d_indheap_size(h)-1);
	
  return 0;
}

/**
 * \ingroup doubleindheap
 * \brief Returns the largest element in the heap.
 */

real_t d_indheap_max       (d_indheap_t* h) {
  assert(h != NULL);
  assert(h->stor_begin != NULL);
  assert(h->stor_begin != h->end);
  
  return h->stor_begin[0];
}

/**
 * \ingroup doubleindheap
 * \brief Removes the largest element from the heap.
 */

real_t d_indheap_delete_max(d_indheap_t* h) {
  real_t tmp;

  assert(h != NULL);
  assert(h->stor_begin != NULL);

  tmp=h->stor_begin[0];
  d_indheap_i_switch(h, 0, d_indheap_size(h)-1);
  h->end -= 1;
  d_indheap_i_sink(h, 0);
  
  return tmp;
}

/**
 * \ingroup doubleindheap
 * \brief Gives the number of elements in the heap.
 */

long int d_indheap_size      (d_indheap_t* h) {
  assert(h != 0);
  assert(h->stor_begin != 0);
  return h->end - h->stor_begin;
}

/**
 * \ingroup doubleindheap
 * \brief Allocates memory for a heap.
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int d_indheap_reserve        (d_indheap_t* h, long int size) {
  long int actual_size=d_indheap_size(h);
  real_t *tmp1;
  long int *tmp2, *tmp3;
  assert(h != 0);
  assert(h->stor_begin != 0);
  
  if (size <= actual_size) { return 0; }

  tmp1=Calloc(size, real_t);
  if (tmp1==0) {
    IGRAPH_FERROR("d_indheap reserve failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, tmp1);	/* TODO: hack */
  tmp2=Calloc(size, long int);
  if (tmp2==0) {
    IGRAPH_FERROR("d_indheap reserve failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, tmp2);	/* TODO: hack */
  tmp3=Calloc(size, long int);
  if (tmp3==0) {
    IGRAPH_FERROR("d_indheap reserve failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, tmp3); 	/* TODO: hack */

  memcpy(tmp1, h->stor_begin, actual_size*sizeof(real_t));
  memcpy(tmp2, h->index_begin, actual_size*sizeof(long int));
  memcpy(tmp3, h->index2_begin, actual_size*sizeof(long int));
  Free(h->stor_begin);
  Free(h->index_begin);
  Free(h->index2_begin);

  h->stor_begin=tmp1;
  h->stor_end=h->stor_begin + size;
  h->end=h->stor_begin+actual_size;
  h->index_begin=tmp2; 
  h->index2_begin=tmp3;

  IGRAPH_FINALLY_CLEAN(3);
  return 0;
}

/**
 * \ingroup doubleindheap
 * \brief Gives the indices of the maximal element in the heap.
 */

void d_indheap_max_index(d_indheap_t *h, long int *idx, long int *idx2) {
  assert(h != 0);
  assert(h->stor_begin != 0);
  (*idx)=h->index_begin[0];
  (*idx2)=h->index2_begin[0];
}

/**
 * \ingroup doubleindheap
 * \brief Builds the heap, don't call it directly.
 */

void d_indheap_i_build(d_indheap_t* h, long int head) {

  long int size=d_indheap_size(h);
  if (RIGHTCHILD(head) < size) { 
    /* both subtrees */
    d_indheap_i_build(h, LEFTCHILD(head) );
    d_indheap_i_build(h, RIGHTCHILD(head));
    d_indheap_i_sink(h, head);
  } else if (LEFTCHILD(head) < size) {
    /* only left */
    d_indheap_i_build(h, LEFTCHILD(head));
    d_indheap_i_sink(h, head);
  } else {
    /* none */
  }
}

/**
 * \ingroup doubleindheap
 * \brief Moves an element up in the heap, don't call it directly.
 */

void d_indheap_i_shift_up(d_indheap_t *h, long int elem) {
  
  if (elem==0 || h->stor_begin[elem] < h->stor_begin[PARENT(elem)]) { 
    /* at the top */
  } else {
    d_indheap_i_switch(h, elem, PARENT(elem));
    d_indheap_i_shift_up(h, PARENT(elem));
  }
}

/**
 * \ingroup doubleindheap
 * \brief Moves an element down in the heap, don't call it directly.
 */

void d_indheap_i_sink(d_indheap_t* h, long int head) {

  long int size=d_indheap_size(h);
  if (LEFTCHILD(head) >= size) { 
    /* no subtrees */
  } else if (RIGHTCHILD(head) == size ||
	     h->stor_begin[LEFTCHILD(head)]>=h->stor_begin[RIGHTCHILD(head)]) {
    /* sink to the left if needed */
    if (h->stor_begin[head] < h->stor_begin[LEFTCHILD(head)]) {
      d_indheap_i_switch(h, head, LEFTCHILD(head));
      d_indheap_i_sink(h, LEFTCHILD(head));
    }
  } else {
    /* sink to the right */
    if (h->stor_begin[head] < h->stor_begin[RIGHTCHILD(head)]) {
      d_indheap_i_switch(h, head, RIGHTCHILD(head));
      d_indheap_i_sink(h, RIGHTCHILD(head));
    }
  }
}

/**
 * \ingroup doubleindheap
 * \brief Switches two elements in the heap, don't call it directly.
 */

void d_indheap_i_switch(d_indheap_t* h, long int e1, long int e2) {
  if (e1!=e2) {
    long int tmpi;
    real_t tmp=h->stor_begin[e1];
    h->stor_begin[e1]=h->stor_begin[e2];
    h->stor_begin[e2]=tmp;
    
    tmpi=h->index_begin[e1];
    h->index_begin[e1]=h->index_begin[e2];
    h->index_begin[e2]=tmpi;

    tmpi=h->index2_begin[e1];
    h->index2_begin[e1]=h->index2_begin[e2];
    h->index2_begin[e2]=tmpi;
  }
}

