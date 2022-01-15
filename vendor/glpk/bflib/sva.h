/* sva.h (sparse vector area) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2012-2013 Free Software Foundation, Inc.
*  Written by Andrew Makhorin <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#ifndef SVA_H
#define SVA_H

/***********************************************************************
*  Sparse Vector Area (SVA) is a container for sparse vectors. This
*  program object is used mainly on computing factorization, where the
*  sparse vectors are rows and columns of sparse matrices.
*
*  The SVA storage is a set of locations numbered 1, 2, ..., size,
*  where size is the size of SVA, which is the total number of
*  locations currently allocated. Each location is identified by its
*  pointer p, 1 <= p <= size, and is the pair (ind[p], val[p]), where
*  ind[p] and val[p] are, respectively, the index and value fields used
*  to store the index and numeric value of a particular vector element.
*
*  Each sparse vector is identified by its reference number k,
*  1 <= k <= n, where n is the total number of vectors currently stored
*  in SVA, and defined by the triplet (ptr[k], len[k], cap[k]), where:
*  ptr[k] is a pointer to the first location of the vector; len[k] is
*  the vector length, which is the number of its non-zero elements,
*  len[k] >= 0; and cap[k] is the capacity of the vector, which is the
*  total number of adjacent locations allocated to that vector,
*  cap[k] >= len[k]. Thus, non-zero elements of k-th vector are stored
*  in locations ptr[k], ptr[k]+1, ..., ptr[k]+len[k]-1, and locations
*  ptr[k]+len[k], ptr[k]+len[k]+1, ..., ptr[k]+cap[k]-1 are reserved.
*
*  The SVA storage is divided into three parts as follows:
*
*  Locations 1, 2, ..., m_ptr-1 constitute the left (dynamic) part of
*  SVA. This part is used to store vectors, whose capacity may change.
*  Note that all vectors stored in the left part are also included in
*  a doubly linked list, where they are ordered by increasing their
*  pointers ptr[k] (this list is needed for efficient implementation
*  of the garbage collector used to defragment the left part of SVA);
*
*  Locations m_ptr, m_ptr+1, ..., r_ptr-1 are free and constitute the
*  middle (free) part of SVA.
*
*  Locations r_ptr, r_ptr+1, ..., size constitute the right (static)
*  part of SVA. This part is used to store vectors, whose capacity is
*  not changed. */

typedef struct SVA SVA;

struct SVA
{     /* sparse vector area */
      int n_max;
      /* maximal value of n (enlarged automatically) */
      int n;
      /* number of currently allocated vectors, 0 <= n <= n_max */
      int *ptr; /* int ptr[1+n_max]; */
      /* ptr[0] is not used;
       * ptr[k], 1 <= i <= n, is pointer to first location of k-th
       * vector in the arrays ind and val */
      int *len; /* int len[1+n_max]; */
      /* len[0] is not used;
       * len[k], 1 <= k <= n, is length of k-th vector, len[k] >= 0 */
      int *cap; /* int cap[1+n_max]; */
      /* cap[0] is not used;
       * cap[k], 1 <= k <= n, is capacity of k-th vector (the number
       * of adjacent locations allocated to it), cap[k] >= len[k] */
      /* NOTE: if cap[k] = 0, then ptr[k] = 0 and len[k] = 0 */
      int size;
      /* total number of locations in SVA */
      int m_ptr, r_ptr;
      /* partitioning pointers that define the left, middle, and right
       * parts of SVA (see above); 1 <= m_ptr <= r_ptr <= size+1 */
      int head;
      /* number of first (leftmost) vector in the linked list */
      int tail;
      /* number of last (rightmost) vector in the linked list */
      int *prev; /* int prev[1+n_max]; */
      /* prev[0] is not used;
       * prev[k] is number of vector which precedes k-th vector in the
       * linked list;
       * prev[k] < 0 means that k-th vector is not in the list */
      int *next; /* int next[1+n_max]; */
      /* next[0] is not used;
       * next[k] is number of vector which succedes k-th vector in the
       * linked list;
       * next[k] < 0 means that k-th vector is not in the list */
      /* NOTE: only vectors having non-zero capacity and stored in the
       *       left part of SVA are included in this linked list */
      int *ind; /* int ind[1+size]; */
      /* ind[0] is not used;
       * ind[p], 1 <= p <= size, is index field of location p */
      double *val; /* double val[1+size]; */
      /* val[0] is not used;
       * val[p], 1 <= p <= size, is value field of location p */
#if 1
      int talky;
      /* option to enable talky mode */
#endif
};

#define sva_create_area _glp_sva_create_area
SVA *sva_create_area(int n_max, int size);
/* create sparse vector area (SVA) */

#define sva_alloc_vecs _glp_sva_alloc_vecs
int sva_alloc_vecs(SVA *sva, int nnn);
/* allocate new vectors in SVA */

#define sva_resize_area _glp_sva_resize_area
void sva_resize_area(SVA *sva, int delta);
/* change size of SVA storage */

#define sva_defrag_area _glp_sva_defrag_area
void sva_defrag_area(SVA *sva);
/* defragment left part of SVA */

#define sva_more_space _glp_sva_more_space
void sva_more_space(SVA *sva, int m_size);
/* increase size of middle (free) part of SVA */

#define sva_enlarge_cap _glp_sva_enlarge_cap
void sva_enlarge_cap(SVA *sva, int k, int new_cap, int skip);
/* enlarge capacity of specified vector */

#define sva_reserve_cap _glp_sva_reserve_cap
void sva_reserve_cap(SVA *sva, int k, int new_cap);
/* reserve locations for specified vector */

#define sva_make_static _glp_sva_make_static
void sva_make_static(SVA *sva, int k);
/* relocate specified vector to right part of SVA */

#define sva_check_area _glp_sva_check_area
void sva_check_area(SVA *sva);
/* check sparse vector area (SVA) */

#define sva_delete_area _glp_sva_delete_area
void sva_delete_area(SVA *sva);
/* delete sparse vector area (SVA) */

#endif

/* eof */
