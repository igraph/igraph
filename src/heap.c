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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

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

int igraph_heap_init           (igraph_heap_t* h, long int alloc_size) {
  if (alloc_size <= 0 ) { alloc_size=1; }
  h->stor_begin=Calloc(alloc_size, igraph_real_t);
  if (h->stor_begin==0) {
    IGRAPH_ERROR("heap init failed", IGRAPH_ENOMEM);
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

int igraph_heap_init_array     (igraph_heap_t *h, igraph_real_t* data, long int len) {
  h->stor_begin=Calloc(len, igraph_real_t);
  if (h->stor_begin==0) {
    IGRAPH_ERROR("heap init from array failed", IGRAPH_ENOMEM);
  }
  h->stor_end=h->stor_begin+len;
  h->end=h->stor_end;
  h->destroy=1;

  memcpy(h->stor_begin, data, len*sizeof(igraph_real_t));

  igraph_heap_i_build (h->stor_begin, h->end-h->stor_begin, 0);
  
  return 0;
}

/**
 * \ingroup heap
 * \brief Destroys an initialized heap object.
 */

void igraph_heap_destroy        (igraph_heap_t* h) {  
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

igraph_bool_t igraph_heap_empty          (igraph_heap_t* h) {
  assert(h != NULL);
  assert(h->stor_begin != NULL);
  return h->stor_begin == h->end;
}

/**
 * \ingroup heap
 * \brief Adds an element to a heap object.
 */

int igraph_heap_push           (igraph_heap_t* h, igraph_real_t elem) {
  int ret;
  assert(h != NULL);
  assert(h->stor_begin != NULL);
	
  /* full, allocate more storage */
  if (h->stor_end == h->end) {
    long int new_size = igraph_heap_size(h) * 2;
    if (new_size == 0) { new_size = 1; }
    IGRAPH_CHECK(ret=igraph_heap_reserve(h, new_size));
  }
	
  *(h->end) = elem;
  h->end += 1;

  /* maintain heap */
  igraph_heap_i_shift_up(h->stor_begin, igraph_heap_size(h), igraph_heap_size(h)-1);
	
  return 0;
}

/**
 * \ingroup heap
 * \brief Returns the largest element in a heap.
 */

igraph_real_t igraph_heap_max       (igraph_heap_t* h) {
  assert(h != NULL);
  assert(h->stor_begin != NULL);
  assert(h->stor_begin != h->end);
  
  return h->stor_begin[0];
}

/**
 * \ingroup heap
 * \brief Returns and removes the largest element in a heap.
 */

igraph_real_t igraph_heap_delete_max(igraph_heap_t* h) {
  igraph_real_t tmp;

  assert(h != NULL);
  assert(h->stor_begin != NULL);

  tmp=h->stor_begin[0];
  igraph_heap_i_switch(h->stor_begin, 0, igraph_heap_size(h)-1);
  h->end -= 1;
  igraph_heap_i_sink(h->stor_begin, h->end-h->stor_begin, 0);
  
  return tmp;
}

/**
 * \ingroup heap
 * \brief Gives the number of elements in a heap.
 */

long int igraph_heap_size      (igraph_heap_t* h) {
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

int igraph_heap_reserve        (igraph_heap_t* h, long int size) {
  long int actual_size=igraph_heap_size(h);
  igraph_real_t *tmp;
  assert(h != NULL);
  assert(h->stor_begin != NULL);
  
  if (size <= actual_size) { return 0; }
  
  tmp=Realloc(h->stor_begin, size, igraph_real_t);
  if (tmp==0) {
    IGRAPH_ERROR("heap reserve failed", IGRAPH_ENOMEM);
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

void igraph_heap_i_build(igraph_real_t* arr, long int size, long int head) {

  if (RIGHTCHILD(head) < size) { 
    /* both subtrees */
    igraph_heap_i_build(arr, size, LEFTCHILD(head) );
    igraph_heap_i_build(arr, size, RIGHTCHILD(head));
    igraph_heap_i_sink(arr, size, head);
  } else if (LEFTCHILD(head) < size) {
    /* only left */
    igraph_heap_i_build(arr, size, LEFTCHILD(head));
    igraph_heap_i_sink(arr, size, head);
  } else {
    /* none */
  }
}

/**
 * \ingroup heap
 * \brief Shift an element upwards in a heap, this should not be
 * called directly.
 */

void igraph_heap_i_shift_up(igraph_real_t* arr, long int size, long int elem) {
  
  if (elem==0 || arr[elem] < arr[PARENT(elem)]) { 
    /* at the top */
  } else {
    igraph_heap_i_switch(arr, elem, PARENT(elem));
    igraph_heap_i_shift_up(arr, size, PARENT(elem));
  }
}

/**
 * \ingroup heap
 * \brief Moves an element down in a heap, this function should not be
 * called directly.
 */

void igraph_heap_i_sink(igraph_real_t* arr, long int size, long int head) {

  if (LEFTCHILD(head) >= size) { 
    /* no subtrees */
  } else if (RIGHTCHILD(head) == size ||
	     arr[LEFTCHILD(head)] >= arr[RIGHTCHILD(head)]) {
    /* sink to the left if needed */
    if (arr[head] < arr[LEFTCHILD(head)]) {
      igraph_heap_i_switch(arr, head, LEFTCHILD(head));
      igraph_heap_i_sink(arr, size, LEFTCHILD(head));
    }
  } else {
    /* sink to the right */
    if (arr[head] < arr[RIGHTCHILD(head)]) {
      igraph_heap_i_switch(arr, head, RIGHTCHILD(head));
      igraph_heap_i_sink(arr, size, RIGHTCHILD(head));
    }
  }
}

/**
 * \ingroup heap
 * \brief Switches two elements in a heap, this function should not be
 * called directly.
 */

void igraph_heap_i_switch(igraph_real_t* arr, long int e1, long int e2) {
  if (e1!=e2) {
    igraph_real_t tmp=arr[e1];
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

int igraph_indheap_init           (igraph_indheap_t* h, long int alloc_size) {
 if (alloc_size <= 0 ) { alloc_size=1; }
 h->stor_begin=Calloc(alloc_size, igraph_real_t);
 if (h->stor_begin==0) {
   h->index_begin=0;
   IGRAPH_ERROR("indheap init failed", IGRAPH_ENOMEM);
 }
 h->index_begin=Calloc(alloc_size, long int);
 if (h->index_begin==0) {
   Free(h->stor_begin);
   h->stor_begin=0;
   IGRAPH_ERROR("indheap init failed", IGRAPH_ENOMEM);
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

int igraph_indheap_init_array     (igraph_indheap_t *h, igraph_real_t* data, long int len) {
  long int i;

  h->stor_begin=Calloc(len, igraph_real_t);
  if (h->stor_begin==0) {
    h->index_begin=0;
    IGRAPH_ERROR("indheap init from array failed", IGRAPH_ENOMEM);
  }
  h->index_begin=Calloc(len, long int);
  if (h->index_begin==0) {
    Free(h->stor_begin);
    h->stor_begin=0;
    IGRAPH_ERROR("indheap init from array failed", IGRAPH_ENOMEM);
  }
  h->stor_end=h->stor_begin+len;
  h->end=h->stor_end;
  h->destroy=1;

  memcpy(h->stor_begin, data, len*sizeof(igraph_real_t));
  for (i=0; i<len; i++) {
    h->index_begin[i]=i+1;
  }

  igraph_indheap_i_build (h, 0);
  
  return 0;
}

/**
 * \ingroup indheap
 * \brief Destroys an initialized indexed heap. 
 */

void igraph_indheap_destroy        (igraph_indheap_t* h) {  
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

igraph_bool_t igraph_indheap_empty          (igraph_indheap_t* h) {
  assert(h != 0);
  assert(h->stor_begin != 0);
  return h->stor_begin == h->end;
}

/**
 * \ingroup indheap
 * \brief Adds an element to an indexed heap.
 */

int igraph_indheap_push           (igraph_indheap_t* h, igraph_real_t elem) {
  assert(h != 0);
  assert(h->stor_begin != 0);
	
  /* full, allocate more storage */
  if (h->stor_end == h->end) {
    long int new_size = igraph_indheap_size(h) * 2;
    if (new_size == 0) { new_size = 1; }
    IGRAPH_CHECK(igraph_indheap_reserve(h, new_size));
  }
	
  *(h->end) = elem;
  h->end += 1;
  *(h->index_begin+igraph_indheap_size(h)-1)=igraph_indheap_size(h)-1;

  /* maintain indheap */
  igraph_indheap_i_shift_up(h, igraph_indheap_size(h)-1);
	
  return 0;
}

/**
 * \ingroup indheap
 * \brief Returns the largest element in an indexed heap.
 */

igraph_real_t igraph_indheap_max       (igraph_indheap_t* h) {
  assert(h != NULL);
  assert(h->stor_begin != NULL);
  assert(h->stor_begin != h->end);
  
  return h->stor_begin[0];
}

/**
 * \ingroup indheap
 * \brief Removes the largest element from an indexed heap.
 */

igraph_real_t igraph_indheap_delete_max(igraph_indheap_t* h) {
  igraph_real_t tmp;

  assert(h != NULL);
  assert(h->stor_begin != NULL);

  tmp=h->stor_begin[0];
  igraph_indheap_i_switch(h, 0, igraph_indheap_size(h)-1);
  h->end -= 1;
  igraph_indheap_i_sink(h, 0);
  
  return tmp;
}

/**
 * \ingroup indheap
 * \brief Gives the number of elements in an indexed heap.
 */

long int igraph_indheap_size      (igraph_indheap_t* h) {
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

int igraph_indheap_reserve        (igraph_indheap_t* h, long int size) {
  long int actual_size=igraph_indheap_size(h);
  igraph_real_t *tmp1;
  long int *tmp2;
  assert(h != 0);
  assert(h->stor_begin != 0);
  
  if (size <= actual_size) { return 0; }

  tmp1=Calloc(size, igraph_real_t);  
  if (tmp1==0) {
    IGRAPH_ERROR("indheap reserve failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, tmp1); 	/* TODO: hack */
  tmp2=Calloc(size, long int);
  if (tmp2==0) {
    IGRAPH_ERROR("indheap reserve failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, tmp2);
  memcpy(tmp1, h->stor_begin, actual_size*sizeof(igraph_real_t));
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

long int igraph_indheap_max_index(igraph_indheap_t *h) {
  assert(h != 0);
  assert(h->stor_begin != 0);
  return h->index_begin[0];  
}

/**
 * \ingroup indheap
 * \brief Builds an indexed heap, this function should not be called
 * directly. 
 */

void igraph_indheap_i_build(igraph_indheap_t* h, long int head) {

  long int size=igraph_indheap_size(h);
  if (RIGHTCHILD(head) < size) { 
    /* both subtrees */
    igraph_indheap_i_build(h, LEFTCHILD(head) );
    igraph_indheap_i_build(h, RIGHTCHILD(head));
    igraph_indheap_i_sink(h, head);
  } else if (LEFTCHILD(head) < size) {
    /* only left */
    igraph_indheap_i_build(h, LEFTCHILD(head));
    igraph_indheap_i_sink(h, head);
  } else {
    /* none */
  }
}

/**
 * \ingroup indheap
 * \brief Moves an element up in the heap, don't call this function
 * directly. 
 */

void igraph_indheap_i_shift_up(igraph_indheap_t *h, long int elem) {
  
  if (elem==0 || h->stor_begin[elem] < h->stor_begin[PARENT(elem)]) { 
    /* at the top */
  } else {
    igraph_indheap_i_switch(h, elem, PARENT(elem));
    igraph_indheap_i_shift_up(h, PARENT(elem));
  }
}

/**
 * \ingroup indheap
 * \brief Moves an element down in the heap, don't call this function
 * directly. 
 */

void igraph_indheap_i_sink(igraph_indheap_t* h, long int head) {

  long int size=igraph_indheap_size(h);
  if (LEFTCHILD(head) >= size) { 
    /* no subtrees */
  } else if (RIGHTCHILD(head) == size ||
	     h->stor_begin[LEFTCHILD(head)]>=h->stor_begin[RIGHTCHILD(head)]) {
    /* sink to the left if needed */
    if (h->stor_begin[head] < h->stor_begin[LEFTCHILD(head)]) {
      igraph_indheap_i_switch(h, head, LEFTCHILD(head));
      igraph_indheap_i_sink(h, LEFTCHILD(head));
    }
  } else {
    /* sink to the right */
    if (h->stor_begin[head] < h->stor_begin[RIGHTCHILD(head)]) {
      igraph_indheap_i_switch(h, head, RIGHTCHILD(head));
      igraph_indheap_i_sink(h, RIGHTCHILD(head));
    }
  }
}

/**
 * \ingroup indheap
 * \brief Switches two elements in a heap, don't call this function
 * directly. 
 */

void igraph_indheap_i_switch(igraph_indheap_t* h, long int e1, long int e2) {
  if (e1!=e2) {
    igraph_real_t tmp=h->stor_begin[e1];
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

int igraph_d_indheap_init           (igraph_d_indheap_t* h, long int alloc_size) {
 if (alloc_size <= 0 ) { alloc_size=1; }
  h->stor_begin=Calloc(alloc_size, igraph_real_t);
  if (h->stor_begin==0) {
    h->index_begin=0;
    h->index2_begin=0;
    IGRAPH_ERROR("d_indheap init failed", IGRAPH_ENOMEM);
  }
  h->stor_end=h->stor_begin + alloc_size;
  h->end=h->stor_begin;
  h->destroy=1;
  h->index_begin=Calloc(alloc_size, long int);
  if (h->index_begin==0) {
    Free(h->stor_begin);
    h->stor_begin=0;
    h->index2_begin=0;
    IGRAPH_ERROR("d_indheap init failed", IGRAPH_ENOMEM);
  }
  h->index2_begin=Calloc(alloc_size, long int);
  if (h->index2_begin==0) {
    Free(h->stor_begin);
    Free(h->index_begin);
    h->stor_begin=0;
    h->index_begin=0;
    IGRAPH_ERROR("d_indheap init failed", IGRAPH_ENOMEM);
  }
  
  return 0;  
}

/**
 * \ingroup doubleindheap
 * \brief Destroys an initialized doubly indexed heap object.
 */

void igraph_d_indheap_destroy        (igraph_d_indheap_t* h) {  
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

igraph_bool_t igraph_d_indheap_empty          (igraph_d_indheap_t* h) {
  assert(h != 0);
  assert(h->stor_begin != 0);
  return h->stor_begin == h->end;
}

/**
 * \ingroup doubleindheap
 * \brief Adds an element to the heap.
 */

int igraph_d_indheap_push           (igraph_d_indheap_t* h, igraph_real_t elem, 
			      long int idx, long int idx2) {
  assert(h != 0);
  assert(h->stor_begin != 0);
	
  /* full, allocate more storage */
  if (h->stor_end == h->end) {
    long int new_size = igraph_d_indheap_size(h) * 2;
    if (new_size == 0) { new_size = 1; }
    IGRAPH_CHECK(igraph_d_indheap_reserve(h, new_size));
  }
	
  *(h->end) = elem;
  h->end += 1;
  *(h->index_begin+igraph_d_indheap_size(h)-1)=idx ;
  *(h->index2_begin+igraph_d_indheap_size(h)-1)=idx2 ;

  /* maintain d_indheap */
  igraph_d_indheap_i_shift_up(h, igraph_d_indheap_size(h)-1);
	
  return 0;
}

/**
 * \ingroup doubleindheap
 * \brief Returns the largest element in the heap.
 */

igraph_real_t igraph_d_indheap_max       (igraph_d_indheap_t* h) {
  assert(h != NULL);
  assert(h->stor_begin != NULL);
  assert(h->stor_begin != h->end);
  
  return h->stor_begin[0];
}

/**
 * \ingroup doubleindheap
 * \brief Removes the largest element from the heap.
 */

igraph_real_t igraph_d_indheap_delete_max(igraph_d_indheap_t* h) {
  igraph_real_t tmp;

  assert(h != NULL);
  assert(h->stor_begin != NULL);

  tmp=h->stor_begin[0];
  igraph_d_indheap_i_switch(h, 0, igraph_d_indheap_size(h)-1);
  h->end -= 1;
  igraph_d_indheap_i_sink(h, 0);
  
  return tmp;
}

/**
 * \ingroup doubleindheap
 * \brief Gives the number of elements in the heap.
 */

long int igraph_d_indheap_size      (igraph_d_indheap_t* h) {
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

int igraph_d_indheap_reserve        (igraph_d_indheap_t* h, long int size) {
  long int actual_size=igraph_d_indheap_size(h);
  igraph_real_t *tmp1;
  long int *tmp2, *tmp3;
  assert(h != 0);
  assert(h->stor_begin != 0);
  
  if (size <= actual_size) { return 0; }

  tmp1=Calloc(size, igraph_real_t);
  if (tmp1==0) {
    IGRAPH_ERROR("d_indheap reserve failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, tmp1);	/* TODO: hack */
  tmp2=Calloc(size, long int);
  if (tmp2==0) {
    IGRAPH_ERROR("d_indheap reserve failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, tmp2);	/* TODO: hack */
  tmp3=Calloc(size, long int);
  if (tmp3==0) {
    IGRAPH_ERROR("d_indheap reserve failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, tmp3); 	/* TODO: hack */

  memcpy(tmp1, h->stor_begin, actual_size*sizeof(igraph_real_t));
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

void igraph_d_indheap_max_index(igraph_d_indheap_t *h, long int *idx, long int *idx2) {
  assert(h != 0);
  assert(h->stor_begin != 0);
  (*idx)=h->index_begin[0];
  (*idx2)=h->index2_begin[0];
}

/**
 * \ingroup doubleindheap
 * \brief Builds the heap, don't call it directly.
 */

void igraph_d_indheap_i_build(igraph_d_indheap_t* h, long int head) {

  long int size=igraph_d_indheap_size(h);
  if (RIGHTCHILD(head) < size) { 
    /* both subtrees */
    igraph_d_indheap_i_build(h, LEFTCHILD(head) );
    igraph_d_indheap_i_build(h, RIGHTCHILD(head));
    igraph_d_indheap_i_sink(h, head);
  } else if (LEFTCHILD(head) < size) {
    /* only left */
    igraph_d_indheap_i_build(h, LEFTCHILD(head));
    igraph_d_indheap_i_sink(h, head);
  } else {
    /* none */
  }
}

/**
 * \ingroup doubleindheap
 * \brief Moves an element up in the heap, don't call it directly.
 */

void igraph_d_indheap_i_shift_up(igraph_d_indheap_t *h, long int elem) {
  
  if (elem==0 || h->stor_begin[elem] < h->stor_begin[PARENT(elem)]) { 
    /* at the top */
  } else {
    igraph_d_indheap_i_switch(h, elem, PARENT(elem));
    igraph_d_indheap_i_shift_up(h, PARENT(elem));
  }
}

/**
 * \ingroup doubleindheap
 * \brief Moves an element down in the heap, don't call it directly.
 */

void igraph_d_indheap_i_sink(igraph_d_indheap_t* h, long int head) {

  long int size=igraph_d_indheap_size(h);
  if (LEFTCHILD(head) >= size) { 
    /* no subtrees */
  } else if (RIGHTCHILD(head) == size ||
	     h->stor_begin[LEFTCHILD(head)]>=h->stor_begin[RIGHTCHILD(head)]) {
    /* sink to the left if needed */
    if (h->stor_begin[head] < h->stor_begin[LEFTCHILD(head)]) {
      igraph_d_indheap_i_switch(h, head, LEFTCHILD(head));
      igraph_d_indheap_i_sink(h, LEFTCHILD(head));
    }
  } else {
    /* sink to the right */
    if (h->stor_begin[head] < h->stor_begin[RIGHTCHILD(head)]) {
      igraph_d_indheap_i_switch(h, head, RIGHTCHILD(head));
      igraph_d_indheap_i_sink(h, RIGHTCHILD(head));
    }
  }
}

/**
 * \ingroup doubleindheap
 * \brief Switches two elements in the heap, don't call it directly.
 */

void igraph_d_indheap_i_switch(igraph_d_indheap_t* h, long int e1, long int e2) {
  if (e1!=e2) {
    long int tmpi;
    igraph_real_t tmp=h->stor_begin[e1];
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

