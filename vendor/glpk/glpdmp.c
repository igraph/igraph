/* glpdmp.c (dynamic memory pool) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
*  2009, 2010 Andrew Makhorin, Department for Applied Informatics,
*  Moscow Aviation Institute, Moscow, Russia. All rights reserved.
*  E-mail: <mao@gnu.org>.
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

#include "glpdmp.h"

#if 1 /* 29/VIII-2008 */
/* some processors need data to be properly aligned; the macro
   align_datasize enlarges the specified size of a data item to provide
   a proper alignment of immediately following data */

#define align_datasize(size) ((((size) + 7) / 8) * 8)
/* 8 bytes is sufficient in both 32- and 64-bit environments */
#endif

#ifdef GLP_DEBUG
struct info
{     DMP *pool;
      int size;
};
#endif

/***********************************************************************
*  NAME
*
*  dmp_create_pool - create dynamic memory pool
*
*  SYNOPSIS
*
*  #include "glpdmp.h"
*  DMP *dmp_create_pool(void);
*
*  DESCRIPTION
*
*  The routine dmp_create_pool creates a dynamic memory pool.
*
*  RETURNS
*
*  The routine returns a pointer to the memory pool created. */

DMP *dmp_create_pool(void)
{     DMP *pool;
      int k;
#ifdef GLP_DEBUG
      xprintf("dmp_create_pool: warning: debug mode enabled\n");
#endif
      pool = xmalloc(sizeof(DMP));
#if 0
      pool->size = 0;
#endif
      for (k = 0; k <= 31; k++) pool->avail[k] = NULL;
      pool->block = NULL;
      pool->used = DMP_BLK_SIZE;
      pool->count.lo = pool->count.hi = 0;
      return pool;
}

/***********************************************************************
*  NAME
*
*  dmp_get_atom - get free atom from dynamic memory pool
*
*  SYNOPSIS
*
*  #include "glpdmp.h"
*  void *dmp_get_atom(DMP *pool, int size);
*
*  DESCRIPTION
*
*  The routine dmp_get_atom obtains a free atom (memory block) from the
*  specified memory pool.
*
*  The parameter size is the atom size, in bytes, 1 <= size <= 256.
*
*  Note that the free atom contains arbitrary data, not binary zeros.
*
*  RETURNS
*
*  The routine returns a pointer to the free atom obtained. */

void *dmp_get_atom(DMP *pool, int size)
{     void *atom;
      int k;
#ifdef GLP_DEBUG
      int orig_size = size;
#endif
      if (!(1 <= size && size <= 256))
         xerror("dmp_get_atom: size = %d; invalid atom size\n", size);
#if 0
      if (!(pool->size == 0 || pool->size == size))
         xerror("dmp_get_atom: size = %d; wrong atom size\n", size);
#endif
      /* adjust the size to provide the proper data alignment */
      size = align_datasize(size);
#ifdef GLP_DEBUG
      size += align_datasize(sizeof(struct info));
#endif
      /* adjust the size to make it multiple of 8 bytes, if needed */
      size = ((size + 7) / 8) * 8;
      /* determine the corresponding list of free cells */
      k = size / 8 - 1;
      xassert(0 <= k && k <= 31);
      /* obtain a free atom */
      if (pool->avail[k] == NULL)
      {  /* the list of free cells is empty */
         if (pool->used + size > DMP_BLK_SIZE)
         {  /* allocate a new memory block */
            void *block = xmalloc(DMP_BLK_SIZE);
            *(void **)block = pool->block;
            pool->block = block;
            pool->used = align_datasize(sizeof(void *));
         }
         /* place the atom in the current memory block */
         atom = (char *)pool->block + pool->used;
         pool->used += size;
      }
      else
      {  /* obtain the atom from the list of free cells */
         atom = pool->avail[k];
         pool->avail[k] = *(void **)atom;
      }
      memset(atom, '?', size);
      /* increase the number of atoms which are currently in use */
      pool->count.lo++;
      if (pool->count.lo == 0) pool->count.hi++;
#ifdef GLP_DEBUG
      ((struct info *)atom)->pool = pool;
      ((struct info *)atom)->size = orig_size;
      atom = (char *)atom + align_datasize(sizeof(struct info));
#endif
      return atom;
}

/***********************************************************************
*  NAME
*
*  dmp_free_atom - return atom to dynamic memory pool
*
*  SYNOPSIS
*
*  #include "glpdmp.h"
*  void dmp_free_atom(DMP *pool, void *atom, int size);
*
*  DESCRIPTION
*
*  The routine dmp_free_atom returns the specified atom (memory block)
*  to the specified memory pool, making it free.
*
*  The parameter size is the atom size, in bytes, 1 <= size <= 256.
*
*  Note that the atom can be returned only to the pool, from which it
*  was obtained, and its size must be exactly the same as on obtaining
*  it from the pool. */

void dmp_free_atom(DMP *pool, void *atom, int size)
{     int k;
      if (!(1 <= size && size <= 256))
         xerror("dmp_free_atom: size = %d; invalid atom size\n", size);
#if 0
      if (!(pool->size == 0 || pool->size == size))
         xerror("dmp_free_atom: size = %d; wrong atom size\n", size);
#endif
      if (pool->count.lo == 0 && pool->count.hi == 0)
         xerror("dmp_free_atom: pool allocation error\n");
#ifdef GLP_DEBUG
      atom = (char *)atom - align_datasize(sizeof(struct info));
      xassert(((struct info *)atom)->pool == pool);
      xassert(((struct info *)atom)->size == size);
#endif
      /* adjust the size to provide the proper data alignment */
      size = align_datasize(size);
#ifdef GLP_DEBUG
      size += align_datasize(sizeof(struct info));
#endif
      /* adjust the size to make it multiple of 8 bytes, if needed */
      size = ((size + 7) / 8) * 8;
      /* determine the corresponding list of free cells */
      k = size / 8 - 1;
      xassert(0 <= k && k <= 31);
      /* return the atom to the list of free cells */
      *(void **)atom = pool->avail[k];
      pool->avail[k] = atom;
      /* decrease the number of atoms which are currently in use */
      pool->count.lo--;
      if (pool->count.lo == 0xFFFFFFFF) pool->count.hi--;
      return;
}

/***********************************************************************
*  NAME
*
*  dmp_in_use - determine how many atoms are still in use
*
*  SYNOPSIS
*
*  #include "glpdmp.h"
*  glp_long dmp_in_use(DMP *pool);
*
*  DESCRIPTION
*
*  The routine dmp_in_use determines how many atoms allocated from the
*  specified memory pool with the routine dmp_get_atom are still in use,
*  i.e. not returned to the pool with the routine dmp_free_atom.
*
*  RETURNS
*
*  The routine returns the number of atoms which are still in use. */

glp_long dmp_in_use(DMP *pool)
{     return
         pool->count;
}

/***********************************************************************
*  NAME
*
*  dmp_delete_pool - delete dynamic memory pool
*
*  SYNOPSIS
*
*  #include "glpdmp.h"
*  void dmp_delete_pool(DMP *pool);
*
*  DESCRIPTION
*
*  The routine dmp_delete_pool deletes the specified dynamic memory
*  pool and frees all the memory allocated to this object. */

void dmp_delete_pool(DMP *pool)
{     while (pool->block != NULL)
      {  void *block = pool->block;
         pool->block = *(void **)block;
         xfree(block);
      }
      xfree(pool);
      return;
}

/* eof */
