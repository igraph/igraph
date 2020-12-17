/* glpdmp.h (dynamic memory pool) */

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

#ifndef GLPDMP_H
#define GLPDMP_H

#include "glpenv.h"

typedef struct DMP DMP;

#define DMP_BLK_SIZE 8000
/* size of memory blocks, in bytes, allocated for memory pools */

struct DMP
{     /* dynamic memory pool */
#if 0
      int size;
      /* size of atoms, in bytes, 1 <= size <= 256; if size = 0, atoms
         may have different sizes */
#endif
      void *avail[32];
      /* avail[k], 0 <= k <= 31, is a pointer to the first available
         (free) cell of (k+1)*8 bytes long; in the beginning of each
         free cell there is a pointer to another free cell of the same
         length */
      void *block;
      /* pointer to the most recently allocated memory block; in the
         beginning of each allocated memory block there is a pointer to
         the previously allocated memory block */
      int used;
      /* number of bytes used in the most recently allocated memory
         block */
      glp_long count;
      /* number of atoms which are currently in use */
};

#define dmp_create_pool _glp_dmp_create_pool
DMP *dmp_create_pool(void);
/* create dynamic memory pool */

#define dmp_get_atom _glp_dmp_get_atom
void *dmp_get_atom(DMP *pool, int size);
/* get free atom from dynamic memory pool */

#define dmp_free_atom _glp_dmp_free_atom
void dmp_free_atom(DMP *pool, void *atom, int size);
/* return atom to dynamic memory pool */

#define dmp_in_use _glp_dmp_in_use
glp_long dmp_in_use(DMP *pool);
/* determine how many atoms are still in use */

#define dmp_delete_pool _glp_dmp_delete_pool
void dmp_delete_pool(DMP *pool);
/* delete dynamic memory pool */

#endif

/* eof */
