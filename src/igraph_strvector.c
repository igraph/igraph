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

/**
 * \ingroup strvector
 * \brief Initializes a string vector (constructor).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int igraph_strvector_init(igraph_strvector_t *sv, long int len) {
  long int i;
  sv->data=Calloc(len, char*);
  if (sv->data==0) {
    IGRAPH_FERROR("strvector init failed", IGRAPH_ENOMEM);
  }
  for (i=0; i<len; i++) {
    sv->data[i]=Calloc(1, char);
    if (sv->data[i]==0) {
      igraph_strvector_destroy(sv);
      IGRAPH_FERROR("strvector init failed", IGRAPH_ENOMEM);
    }
    sv->data[i][0]='\0';
  }
  sv->len=len;

  return 0;
}

/**
 * \ingroup strvector
 * \brief Frees memory allocated for a string vector.
 */

void igraph_strvector_destroy(igraph_strvector_t *sv) {
  long int i;
  assert(sv != 0);
  if (sv->data != 0) {
    for (i=0; i<sv->len; i++) {
      if (sv->data[i] != 0) {
	Free(sv->data[i]);
      }
    }
    Free(sv->data);
  }
}

/**
 * \ingroup strvector
 * \brief Returns an element of a string vector.
 */

void igraph_strvector_get(igraph_strvector_t *sv, long int idx, char **value) {
  assert(sv != 0);
  assert(sv->data != 0);
  assert(sv->data[idx] != 0);
  *value = sv->data[idx];
}

/**
 * \ingroup strvector
 * \brief Sets an element of a string vector.
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int igraph_strvector_set(igraph_strvector_t *sv, long int idx, 
			 const char *value) {
  assert(sv != 0);
  assert(sv->data != 0);
  if (sv->data[idx] == 0) {
    sv->data[idx] = Calloc(strlen(value)+1, char);
    if (sv->data[idx]==0) {
      IGRAPH_FERROR("strvector set failed", IGRAPH_ENOMEM);
    }
  } else {
    char *tmp=Realloc(sv->data[idx], strlen(value)+1, char);
    if (tmp==0) { 
      IGRAPH_FERROR("strvector set failed", IGRAPH_ENOMEM);
    }
    sv->data[idx]=tmp;
  }
  strcpy(sv->data[idx], value);
  
  return 0;
}

/**
 * \ingroup strvector
 * \brief Removes a section from a string vector.
 * \todo repair realloc
 */

void igraph_strvector_remove_section(igraph_strvector_t *v, long int from, 
				    long int to) {
  long int i;
/*   char **tmp; */
  
  assert(v != 0);
  assert(v->data != 0);

  for (i=from; i<to; i++) {
    if (v->data[i] != 0) {
      Free(v->data[i]);
    }
  }
  for (i=0; i<v->len-to; i++) {
    v->data[from+i]=v->data[to+i];
  }
  
  v->len -= (to-from);

  /* try to make it smaller */
/*   tmp=Realloc(v->data, v->len, char*); */
/*   if (tmp!=0) { */
/*     v->data=tmp; */
/*   } */
}

/**
 * \ingroup strvector
 * \brief Removes a single element from a string vector.
 */

void igraph_strvector_remove(igraph_strvector_t *v, long int elem) {
  assert(v != 0);
  assert(v->data != 0);
  igraph_strvector_remove_section(v, elem, elem+1);
}

/**
 * \ingroup strvector
 * \brief Copies an interval of a string vector.
 */

void igraph_strvector_move_interval(igraph_strvector_t *v, long int begin, 
				   long int end, long int to) {
  long int i;
  assert(v != 0);
  assert(v->data != 0);
  for (i=to; i<to+end-begin; i++) {
    if (v->data[i] != 0) {
      Free(v->data[i]);
    }
  }
  for (i=0; i<end-begin; i++) {
    v->data[to+i]=v->data[begin+i];
  }
}

/**
 * \ingroup strvector
 * \brief Copies a string vector (constructor).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 * \todo why does the assert fail
 */

int igraph_strvector_copy(igraph_strvector_t *to, igraph_strvector_t *from) {
  long int i;
  char *str;
  assert(from != 0);
/*   assert(from->data != 0); */
  to->data=Calloc(from->len, char*);
  if (to->data==0) { 
    IGRAPH_FERROR("Cannot copy string vector", IGRAPH_ENOMEM);
  }
  to->len=from->len;
  
  for (i=0; i<from->len; i++) {
    int ret;
    igraph_strvector_get(from, i, &str);
    ret=igraph_strvector_set(to, i, str);
    if (ret != 0) {
      igraph_strvector_destroy(to);
      IGRAPH_FERROR("cannot copy string vector", ret);
    }
  }
  
  return 0;
}

/**
 * \ingroup strvector
 * \brief Resizes a string vector.
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int igraph_strvector_resize(igraph_strvector_t* v, long int newsize) {
  long int toadd=newsize-v->len, i, j;
  char **tmp;
  
  assert(v != 0);
  assert(v->data != 0);
  if (newsize < v->len) { 
    long int i;
    for (i=newsize; i<v->len; i++) {
      Free(v->data[i]);
    }
    /* try to give back some space */
    tmp=Realloc(v->data, newsize, char*);
    if (tmp != 0) {
      v->data=tmp;
    }
  } else if (newsize > v->len) {
    bool_t error=0;
    tmp=Realloc(v->data, newsize, char*);
    if (tmp==0) {
      IGRAPH_FERROR("cannot resize string vector", IGRAPH_ENOMEM);
    }
    v->data = tmp;
    
    for (i=0; i<toadd; i++) {
      v->data[v->len+i] = Calloc(1, char);
      if (v->data[v->len+i] == 0) {
	error=1;
	break;
      }
      v->data[v->len+i][0]='\0';
    }
    if (error) {
      /* There was an error, free everything we've allocated so far */
      for (j=0; j<i; j++) {
	if (v->data[v->len+i] != 0) {
	  Free(v->data[v->len+i]);
	}
      }
      /* Try to give back space */
      tmp=Realloc(v->data, v->len, char*);
      if (tmp != 0) {
	v->data=tmp;
      }
      IGRAPH_FERROR("canot resize string vector", IGRAPH_ENOMEM);
    }
  }
  v->len = newsize;
  
  return 0;
}

/**
 * \ingroup strvector
 * \brief Gives the size of a string vector.
 * \todo repair assert
 */

long int igraph_strvector_size(igraph_strvector_t *sv) {
  assert(sv != 0);
  assert(sv->data != 0);
  return sv->len;
}

/**
 * \ingroup strvector
 * \brief Adds an element to the back of a string vector.
 */

int igraph_strvector_add(igraph_strvector_t *v, const char *value) {
  long int s=igraph_strvector_size(v);
  char **tmp;
  assert(v != 0);
  assert(v->data != 0);
  tmp=Realloc(v->data, s+1, char*);
  if (tmp == 0) {
    IGRAPH_FERROR("cannot add string to string vector", IGRAPH_ENOMEM);
  }
  v->data=tmp;
  v->data[s]=Calloc(strlen(value)+1, char);
  if (v->data[s]==0) {
    IGRAPH_FERROR("cannot add string to string vector", IGRAPH_ENOMEM);
  }
  strcpy(v->data[s], value);
  v->len += 1;

  return 0;
}

/**
 * \ingroup strvector
 * \brief Removes elements from a string vector (for internal use)
 */

void igraph_strvector_permdelete(igraph_strvector_t *v, long int *index, 
				long int nremove) {
  long int i;
  char **tmp;
  assert(v != 0);
  assert(v->data != 0);

  for (i=0; i<igraph_strvector_size(v); i++) {
    if (index[i] != 0) {
      v->data[ index[i]-1 ] = v->data[i];
    } else {
      Free(v->data[i]);
    }
  }
  /* Try to make it shorter */
  tmp=Realloc(v->data, v->len-nremove, char*);
  if (tmp != 0) {
    v->data=tmp;
  }
  v->len -= nremove;
}

/**
 * \ingroup strvector
 * \brief Removes elements from a string vector (for internal use)
 */

void igraph_strvector_remove_negidx(igraph_strvector_t *v, vector_t *neg, 
				   long int nremove) {
  long int i, idx=0;
  char **tmp;
  assert(v != 0);
  assert(v->data != 0);
  for (i=0; i<igraph_strvector_size(v); i++) {
    if (VECTOR(*neg)[i] >= 0) {
      v->data[idx++] = v->data[i];
    } else {
      Free(v->data[i]);
    }
  }
  /* Try to give back some memory */
  tmp=Realloc(v->data, v->len-nremove, char*);
  if (tmp != 0) {
    v->data=tmp;
  }
  v->len -= nremove;
}

