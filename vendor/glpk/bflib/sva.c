/* sva.c (sparse vector area) */

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

#include "env.h"
#include "sva.h"

/***********************************************************************
*  sva_create_area - create sparse vector area (SVA)
*
*  This routine creates the sparse vector area (SVA), which initially
*  is empty.
*
*  The parameter n_max specifies the initial number of vectors that can
*  be allocated in the SVA, n_max > 0.
*
*  The parameter size specifies the initial number of free locations in
*  the SVA, size > 0.
*
*  On exit the routine returns a pointer to the SVA created. */

SVA *sva_create_area(int n_max, int size)
{     SVA *sva;
      xassert(0 < n_max && n_max < INT_MAX);
      xassert(0 < size && size < INT_MAX);
      sva = talloc(1, SVA);
      sva->n_max = n_max;
      sva->n = 0;
      sva->ptr = talloc(1+n_max, int);
      sva->len = talloc(1+n_max, int);
      sva->cap = talloc(1+n_max, int);
      sva->size = size;
      sva->m_ptr = 1;
      sva->r_ptr = size+1;
      sva->head = sva->tail = 0;
      sva->prev = talloc(1+n_max, int);
      sva->next = talloc(1+n_max, int);
      sva->ind = talloc(1+size, int);
      sva->val = talloc(1+size, double);
      sva->talky = 0;
      return sva;
}

/***********************************************************************
*  sva_alloc_vecs - allocate new vectors in SVA
*
*  This routine allocates nnn new empty vectors, nnn > 0, in the sparse
*  vector area (SVA).
*
*  The new vectors are assigned reference numbers k, k+1, ..., k+nnn-1,
*  where k is a reference number assigned to the very first new vector,
*  which is returned by the routine on exit. */

int sva_alloc_vecs(SVA *sva, int nnn)
{     int n = sva->n;
      int n_max = sva->n_max;
      int *ptr = sva->ptr;
      int *len = sva->len;
      int *cap = sva->cap;
      int *prev = sva->prev;
      int *next = sva->next;
      int k, new_n;
#if 1
      if (sva->talky)
         xprintf("sva_alloc_vecs: nnn = %d\n", nnn);
#endif
      xassert(nnn > 0);
      /* determine new number of vectors in SVA */
      new_n = n + nnn;
      xassert(new_n > n);
      if (n_max < new_n)
      {  /* enlarge the SVA arrays */
         while (n_max < new_n)
         {  n_max += n_max;
            xassert(n_max > 0);
         }
         sva->n_max = n_max;
         sva->ptr = ptr = trealloc(ptr, 1+n_max, int);
         sva->len = len = trealloc(len, 1+n_max, int);
         sva->cap = cap = trealloc(cap, 1+n_max, int);
         sva->prev = prev = trealloc(prev, 1+n_max, int);
         sva->next = next = trealloc(next, 1+n_max, int);
      }
      /* initialize new vectors */
      sva->n = new_n;
      for (k = n+1; k <= new_n; k++)
      {  ptr[k] = len[k] = cap[k] = 0;
         prev[k] = next[k] = -1;
      }
#if 1
      if (sva->talky)
         xprintf("now sva->n_max = %d, sva->n = %d\n",
            sva->n_max, sva->n);
#endif
      /* return reference number of very first new vector */
      return n+1;
}

/***********************************************************************
*  sva_resize_area - change size of SVA storage
*
*  This routine increases or decrases the size of the SVA storage by
*  reallocating it.
*
*  The parameter delta specifies the number of location by which the
*  current size of the SVA storage should be increased (if delta > 0)
*  or decreased (if delta < 0). Note that if delta is negative, it
*  should not be less than the current size of the middle part.
*
*  As a result of this operation the size of the middle part of SVA is
*  increased/decreased by delta locations.
*
*  NOTE: This operation changes ptr[k] for all vectors stored in the
*        right part of SVA. */

void sva_resize_area(SVA *sva, int delta)
{     int n = sva->n;
      int *ptr = sva->ptr;
      int size = sva->size;
      int m_ptr = sva->m_ptr;
      int r_ptr = sva->r_ptr;
      int k, r_size;
#if 1
      if (sva->talky)
         xprintf("sva_resize_area: delta = %d\n", delta);
#endif
      xassert(delta != 0);
      /* determine size of the right part, in locations */
      r_size = size - r_ptr + 1;
      /* relocate the right part in case of negative delta */
      if (delta < 0)
      {  xassert(delta >= m_ptr - r_ptr);
         sva->r_ptr += delta;
         memmove(&sva->ind[sva->r_ptr], &sva->ind[r_ptr],
            r_size * sizeof(int));
         memmove(&sva->val[sva->r_ptr], &sva->val[r_ptr],
            r_size * sizeof(double));
      }
      /* reallocate the storage arrays */
      xassert(delta < INT_MAX - sva->size);
      sva->size += delta;
      sva->ind = trealloc(sva->ind, 1+sva->size, int);
      sva->val = trealloc(sva->val, 1+sva->size, double);
      /* relocate the right part in case of positive delta */
      if (delta > 0)
      {  sva->r_ptr += delta;
         memmove(&sva->ind[sva->r_ptr], &sva->ind[r_ptr],
            r_size * sizeof(int));
         memmove(&sva->val[sva->r_ptr], &sva->val[r_ptr],
            r_size * sizeof(double));
      }
      /* update pointers to vectors stored in the right part */
      for (k = 1; k <= n; k++)
      {  if (ptr[k] >= r_ptr)
            ptr[k] += delta;
      }
#if 1
      if (sva->talky)
         xprintf("now sva->size = %d\n", sva->size);
#endif
      return;
}

/***********************************************************************
*  sva_defrag_area - defragment left part of SVA
*
*  This routine performs "garbage" collection to defragment the left
*  part of SVA.
*
*  NOTE: This operation may change ptr[k] and cap[k] for all vectors
*        stored in the left part of SVA. */

void sva_defrag_area(SVA *sva)
{     int *ptr = sva->ptr;
      int *len = sva->len;
      int *cap = sva->cap;
      int *prev = sva->prev;
      int *next = sva->next;
      int *ind = sva->ind;
      double *val = sva->val;
      int k, next_k, ptr_k, len_k, m_ptr, head, tail;
#if 1
      if (sva->talky)
      {  xprintf("sva_defrag_area:\n");
         xprintf("before defragmenting = %d %d %d\n", sva->m_ptr - 1,
            sva->r_ptr - sva->m_ptr, sva->size + 1 - sva->r_ptr);
      }
#endif
      m_ptr = 1;
      head = tail = 0;
      /* walk through the linked list of vectors stored in the left
       * part of SVA */
      for (k = sva->head; k != 0; k = next_k)
      {  /* save number of next vector in the list */
         next_k = next[k];
         /* determine length of k-th vector */
         len_k = len[k];
         if (len_k == 0)
         {  /* k-th vector is empty; remove it from the left part */
            ptr[k] = cap[k] = 0;
            prev[k] = next[k] = -1;
         }
         else
         {  /* determine pointer to first location of k-th vector */
            ptr_k = ptr[k];
            xassert(m_ptr <= ptr_k);
            /* relocate k-th vector to the beginning of the left part,
             * if necessary */
            if (m_ptr < ptr_k)
            {  memmove(&ind[m_ptr], &ind[ptr_k],
                  len_k * sizeof(int));
               memmove(&val[m_ptr], &val[ptr_k],
                  len_k * sizeof(double));
               ptr[k] = m_ptr;
            }
            /* remove unused locations from k-th vector */
            cap[k] = len_k;
            /* the left part of SVA has been enlarged */
            m_ptr += len_k;
            /* add k-th vector to the end of the new linked list */
            prev[k] = tail;
            next[k] = 0;
            if (head == 0)
               head = k;
            else
               next[tail] = k;
            tail = k;
         }
      }
      /* set new pointer to the middle part of SVA */
      xassert(m_ptr <= sva->r_ptr);
      sva->m_ptr = m_ptr;
      /* set new head and tail of the linked list */
      sva->head = head;
      sva->tail = tail;
#if 1
      if (sva->talky)
         xprintf("after defragmenting = %d %d %d\n", sva->m_ptr - 1,
            sva->r_ptr - sva->m_ptr, sva->size + 1 - sva->r_ptr);
#endif
      return;
}

/***********************************************************************
*  sva_more_space - increase size of middle (free) part of SVA
*
*  This routine increases the size of the middle (free) part of the
*  sparse vector area (SVA).
*
*  The parameter m_size specifies the minimal size, in locations, of
*  the middle part to be provided. This new size should be greater than
*  the current size of the middle part.
*
*  First, the routine defragments the left part of SVA. Then, if the
*  size of the left part has not sufficiently increased, the routine
*  increases the total size of the SVA storage by reallocating it. */

void sva_more_space(SVA *sva, int m_size)
{     int size, delta;
#if 1
      if (sva->talky)
         xprintf("sva_more_space: m_size = %d\n", m_size);
#endif
      xassert(m_size > sva->r_ptr - sva->m_ptr);
      /* defragment the left part */
      sva_defrag_area(sva);
      /* set, heuristically, the minimal size of the middle part to be
       * not less than the size of the defragmented left part */
      if (m_size < sva->m_ptr - 1)
         m_size = sva->m_ptr - 1;
      /* if there is still not enough room, increase the total size of
       * the SVA storage */
      if (sva->r_ptr - sva->m_ptr < m_size)
      {  size = sva->size; /* new sva size */
         for (;;)
         {  delta = size - sva->size;
            if (sva->r_ptr - sva->m_ptr + delta >= m_size)
               break;
            size += size;
            xassert(size > 0);
         }
         sva_resize_area(sva, delta);
         xassert(sva->r_ptr - sva->m_ptr >= m_size);
      }
      return;
}

/***********************************************************************
*  sva_enlarge_cap - enlarge capacity of specified vector
*
*  This routine enlarges the current capacity of the specified vector
*  by relocating its content.
*
*  The parameter k specifies the reference number of the vector whose
*  capacity should be enlarged, 1 <= k <= n. This vector should either
*  have zero capacity or be stored in the left (dynamic) part of SVA.
*
*  The parameter new_cap specifies the new capacity of the vector,
*  in locations. This new capacity should be greater than the current
*  capacity of the vector.
*
*  The parameter skip is a flag. If this flag is set, the routine does
*  *not* copy numerical values of elements of the vector on relocating
*  its content, i.e. only element indices are copied.
*
*  NOTE: On entry to the routine the middle part of SVA should have at
*        least new_cap free locations. */

void sva_enlarge_cap(SVA *sva, int k, int new_cap, int skip)
{     int *ptr = sva->ptr;
      int *len = sva->len;
      int *cap = sva->cap;
      int *prev = sva->prev;
      int *next = sva->next;
      int *ind = sva->ind;
      double *val = sva->val;
      xassert(1 <= k && k <= sva->n);
      xassert(new_cap > cap[k]);
      /* there should be at least new_cap free locations */
      xassert(sva->r_ptr - sva->m_ptr >= new_cap);
      /* relocate the vector */
      if (cap[k] == 0)
      {  /* the vector is empty */
         xassert(ptr[k] == 0);
         xassert(len[k] == 0);
      }
      else
      {  /* the vector has non-zero capacity */
         xassert(ptr[k] + len[k] <= sva->m_ptr);
         /* copy the current vector content to the beginning of the
          * middle part */
         if (len[k] > 0)
         {  memcpy(&ind[sva->m_ptr], &ind[ptr[k]],
               len[k] * sizeof(int));
            if (!skip)
               memcpy(&val[sva->m_ptr], &val[ptr[k]],
                  len[k] * sizeof(double));
         }
         /* remove the vector from the linked list */
         if (prev[k] == 0)
            sva->head = next[k];
         else
         {  /* preceding vector exists; increase its capacity */
            cap[prev[k]] += cap[k];
            next[prev[k]] = next[k];
         }
         if (next[k] == 0)
            sva->tail = prev[k];
         else
            prev[next[k]] = prev[k];
      }
      /* set new pointer and capacity of the vector */
      ptr[k] = sva->m_ptr;
      cap[k] = new_cap;
      /* add the vector to the end of the linked list */
      prev[k] = sva->tail;
      next[k] = 0;
      if (sva->head == 0)
         sva->head = k;
      else
         next[sva->tail] = k;
      sva->tail = k;
      /* new_cap free locations have been consumed */
      sva->m_ptr += new_cap;
      xassert(sva->m_ptr <= sva->r_ptr);
      return;
}

/***********************************************************************
*  sva_reserve_cap - reserve locations for specified vector
*
*  This routine reserves locations for the specified vector in the
*  right (static) part of SVA.
*
*  The parameter k specifies the reference number of the vector (this
*  vector should have zero capacity), 1 <= k <= n.
*
*  The parameter new_cap specifies a non-zero capacity of the vector,
*  in locations.
*
*  NOTE: On entry to the routine the middle part of SVA should have at
*        least new_cap free locations. */

void sva_reserve_cap(SVA *sva, int k, int new_cap)
{     int *ptr = sva->ptr;
      int *len = sva->len;
      int *cap = sva->cap;
      xassert(1 <= k && k <= sva->n);
      xassert(new_cap > 0);
      xassert(ptr[k] == 0 && len[k] == 0 && cap[k] == 0);
      /* there should be at least new_cap free locations */
      xassert(sva->r_ptr - sva->m_ptr >= new_cap);
      /* set the pointer and capacity of the vector */
      ptr[k] = sva->r_ptr - new_cap;
      cap[k] = new_cap;
      /* new_cap free locations have been consumed */
      sva->r_ptr -= new_cap;
      return;
}

/***********************************************************************
*  sva_make_static - relocate specified vector to right part of SVA
*
*  Assuming that the specified vector is stored in the left (dynamic)
*  part of SVA, this routine makes the vector static by relocating its
*  content to the right (static) part of SVA. However, if the specified
*  vector has zero capacity, the routine does nothing.
*
*  The parameter k specifies the reference number of the vector to be
*  relocated, 1 <= k <= n.
*
*  NOTE: On entry to the routine the middle part of SVA should have at
*        least len[k] free locations, where len[k] is the length of the
*        vector to be relocated. */

void sva_make_static(SVA *sva, int k)
{     int *ptr = sva->ptr;
      int *len = sva->len;
      int *cap = sva->cap;
      int *prev = sva->prev;
      int *next = sva->next;
      int *ind = sva->ind;
      double *val = sva->val;
      int ptr_k, len_k;
      xassert(1 <= k && k <= sva->n);
      /* if the vector has zero capacity, do nothing */
      if (cap[k] == 0)
      {  xassert(ptr[k] == 0);
         xassert(len[k] == 0);
         goto done;
      }
      /* there should be at least len[k] free locations */
      len_k = len[k];
      xassert(sva->r_ptr - sva->m_ptr >= len_k);
      /* remove the vector from the linked list */
      if (prev[k] == 0)
         sva->head = next[k];
      else
      {  /* preceding vector exists; increase its capacity */
         cap[prev[k]] += cap[k];
         next[prev[k]] = next[k];
      }
      if (next[k] == 0)
         sva->tail = prev[k];
      else
         prev[next[k]] = prev[k];
      /* if the vector has zero length, make it empty */
      if (len_k == 0)
      {  ptr[k] = cap[k] = 0;
         goto done;
      }
      /* copy the vector content to the beginning of the right part */
      ptr_k = sva->r_ptr - len_k;
      memcpy(&ind[ptr_k], &ind[ptr[k]], len_k * sizeof(int));
      memcpy(&val[ptr_k], &val[ptr[k]], len_k * sizeof(double));
      /* set new pointer and capacity of the vector */
      ptr[k] = ptr_k;
      cap[k] = len_k;
      /* len[k] free locations have been consumed */
      sva->r_ptr -= len_k;
done: return;
}

/***********************************************************************
*  sva_check_area - check sparse vector area (SVA)
*
*  This routine checks the SVA data structures for correctness.
*
*  NOTE: For testing/debugging only. */

void sva_check_area(SVA *sva)
{     int n_max = sva->n_max;
      int n = sva->n;
      int *ptr = sva->ptr;
      int *len = sva->len;
      int *cap = sva->cap;
      int size = sva->size;
      int m_ptr = sva->m_ptr;
      int r_ptr = sva->r_ptr;
      int head = sva->head;
      int tail = sva->tail;
      int *prev = sva->prev;
      int *next = sva->next;
      int k;
#if 0 /* 16/II-2004; SVA may be empty */
      xassert(1 <= n && n <= n_max);
#else
      xassert(0 <= n && n <= n_max);
#endif
      xassert(1 <= m_ptr && m_ptr <= r_ptr && r_ptr <= size+1);
      /* all vectors included the linked list should have non-zero
       * capacity and be stored in the left part */
      for (k = head; k != 0; k = next[k])
      {  xassert(1 <= k && k <= n);
         xassert(cap[k] > 0);
         xassert(0 <= len[k] && len[k] <= cap[k]);
         if (prev[k] == 0)
            xassert(k == head);
         else
         {  xassert(1 <= prev[k] && prev[k] <= n);
            xassert(next[prev[k]] == k);
         }
         if (next[k] == 0)
         {  xassert(k == tail);
            xassert(ptr[k] + cap[k] <= m_ptr);
         }
         else
         {  xassert(1 <= next[k] && next[k] <= n);
            xassert(prev[next[k]] == k);
            xassert(ptr[k] + cap[k] <= ptr[next[k]]);
         }
         cap[k] = -cap[k];
      }
      /* all other vectors should either have zero capacity or be
       * stored in the right part */
      for (k = 1; k <= n; k++)
      {  if (cap[k] < 0)
         {  /* k-th vector is stored in the left part */
            cap[k] = -cap[k];
         }
         else if (cap[k] == 0)
         {  /* k-th vector has zero capacity */
            xassert(ptr[k] == 0);
            xassert(len[k] == 0);
         }
         else /* cap[k] > 0 */
         {  /* k-th vector is stored in the right part */
            xassert(0 <= len[k] && len[k] <= cap[k]);
            xassert(r_ptr <= ptr[k] && ptr[k] + cap[k] <= size+1);
         }
      }
      return;
}

/***********************************************************************
*  sva_delete_area - delete sparse vector area (SVA)
*
*  This routine deletes the sparse vector area (SVA) freeing all the
*  memory allocated to it. */

void sva_delete_area(SVA *sva)
{     tfree(sva->ptr);
      tfree(sva->len);
      tfree(sva->cap);
      tfree(sva->prev);
      tfree(sva->next);
      tfree(sva->ind);
      tfree(sva->val);
      tfree(sva);
      return;
}

/* eof */
