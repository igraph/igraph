/* glpenv05.c (memory allocation) */

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

#include "glpapi.h"

/* some processors need data to be properly aligned; the macro
   align_datasize enlarges the specified size of a data item to provide
   a proper alignment of immediately following data */

#define align_datasize(size) ((((size) + 15) / 16) * 16)
/* 16 bytes is sufficient in both 32- and 64-bit environments
   (8 bytes is not sufficient in 64-bit environment due to jmp_buf) */

/***********************************************************************
*  NAME
*
*  glp_malloc - allocate memory block
*
*  SYNOPSIS
*
*  void *glp_malloc(int size);
*
*  DESCRIPTION
*
*  The routine glp_malloc allocates a memory block of size bytes long.
*
*  Note that being allocated the memory block contains arbitrary data
*  (not binary zeros).
*
*  RETURNS
*
*  The routine glp_malloc returns a pointer to the allocated block.
*  To free this block the routine glp_free (not free!) must be used. */

void *glp_malloc(int size)
{     ENV *env = get_env_ptr();
      MEM *desc;
      int size_of_desc = align_datasize(sizeof(MEM));
      if (size < 1 || size > INT_MAX - size_of_desc)
         xerror("glp_malloc: size = %d; invalid parameter\n", size);
      size += size_of_desc;
      if (xlcmp(xlset(size),
          xlsub(env->mem_limit, env->mem_total)) > 0)
         xerror("glp_malloc: memory limit exceeded\n");
      if (env->mem_count == INT_MAX)
         xerror("glp_malloc: too many memory blocks allocated\n");
      desc = malloc(size);
      if (desc == NULL)
         xerror("glp_malloc: no memory available\n");
      memset(desc, '?', size);
      desc->flag = MEM_MAGIC;
      desc->size = size;
      desc->prev = NULL;
      desc->next = env->mem_ptr;
      if (desc->next != NULL) desc->next->prev = desc;
      env->mem_ptr = desc;
      env->mem_count++;
      if (env->mem_cpeak < env->mem_count)
         env->mem_cpeak = env->mem_count;
      env->mem_total = xladd(env->mem_total, xlset(size));
      if (xlcmp(env->mem_tpeak, env->mem_total) < 0)
         env->mem_tpeak = env->mem_total;
      return (void *)((char *)desc + size_of_desc);
}

/***********************************************************************
*  NAME
*
*  glp_calloc - allocate memory block
*
*  SYNOPSIS
*
*  void *glp_calloc(int n, int size);
*
*  DESCRIPTION
*
*  The routine glp_calloc allocates a memory block of (n*size) bytes
*  long.
*
*  Note that being allocated the memory block contains arbitrary data
*  (not binary zeros).
*
*  RETURNS
*
*  The routine glp_calloc returns a pointer to the allocated block.
*  To free this block the routine glp_free (not free!) must be used. */

void *glp_calloc(int n, int size)
{     if (n < 1)
         xerror("glp_calloc: n = %d; invalid parameter\n", n);
      if (size < 1)
         xerror("glp_calloc: size = %d; invalid parameter\n", size);
      if (n > INT_MAX / size)
         xerror("glp_calloc: n = %d; size = %d; array too big\n", n,
            size);
      return xmalloc(n * size);
}

/***********************************************************************
*  NAME
*
*  glp_free - free memory block
*
*  SYNOPSIS
*
*  void glp_free(void *ptr);
*
*  DESCRIPTION
*
*  The routine glp_free frees a memory block pointed to by ptr, which
*  was previuosly allocated by the routine glp_malloc or glp_calloc. */

void glp_free(void *ptr)
{     ENV *env = get_env_ptr();
      MEM *desc;
      int size_of_desc = align_datasize(sizeof(MEM));
      if (ptr == NULL)
         xerror("glp_free: ptr = %p; null pointer\n", ptr);
      desc = (void *)((char *)ptr - size_of_desc);
      if (desc->flag != MEM_MAGIC)
         xerror("glp_free: ptr = %p; invalid pointer\n", ptr);
      if (env->mem_count == 0 ||
          xlcmp(env->mem_total, xlset(desc->size)) < 0)
         xerror("glp_free: memory allocation error\n");
      if (desc->prev == NULL)
         env->mem_ptr = desc->next;
      else
         desc->prev->next = desc->next;
      if (desc->next == NULL)
         ;
      else
         desc->next->prev = desc->prev;
      env->mem_count--;
      env->mem_total = xlsub(env->mem_total, xlset(desc->size));
      memset(desc, '?', size_of_desc);
      free(desc);
      return;
}

/***********************************************************************
*  NAME
*
*  glp_mem_limit - set memory usage limit
*
*  SYNOPSIS
*
*  void glp_mem_limit(int limit);
*
*  DESCRIPTION
*
*  The routine glp_mem_limit limits the amount of memory available for
*  dynamic allocation (in GLPK routines) to limit megabytes. */

void glp_mem_limit(int limit)
{     ENV *env = get_env_ptr();
      if (limit < 0)
         xerror("glp_mem_limit: limit = %d; invalid parameter\n",
            limit);
      env->mem_limit = xlmul(xlset(limit), xlset(1 << 20));
      return;
}

/***********************************************************************
*  NAME
*
*  glp_mem_usage - get memory usage information
*
*  SYNOPSIS
*
*  void glp_mem_usage(int *count, int *cpeak, glp_long *total,
*     glp_long *tpeak);
*
*  DESCRIPTION
*
*  The routine glp_mem_usage reports some information about utilization
*  of the memory by GLPK routines. Information is stored to locations
*  specified by corresponding parameters (see below). Any parameter can
*  be specified as NULL, in which case corresponding information is not
*  stored.
*
*  *count is the number of the memory blocks currently allocated by the
*  routines xmalloc and xcalloc (one call to xmalloc or xcalloc results
*  in allocating one memory block).
*
*  *cpeak is the peak value of *count reached since the initialization
*  of the GLPK library environment.
*
*  *total is the total amount, in bytes, of the memory blocks currently
*  allocated by the routines xmalloc and xcalloc.
*
*  *tpeak is the peak value of *total reached since the initialization
*  of the GLPK library envirionment. */

void glp_mem_usage(int *count, int *cpeak, glp_long *total,
      glp_long *tpeak)
{     ENV *env = get_env_ptr();
      if (count != NULL) *count = env->mem_count;
      if (cpeak != NULL) *cpeak = env->mem_cpeak;
      if (total != NULL) *total = env->mem_total;
      if (tpeak != NULL) *tpeak = env->mem_tpeak;
      return;
}

/* eof */
