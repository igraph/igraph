/* dmp.h (dynamic memory pool) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2000-2013 Free Software Foundation, Inc.
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

#ifndef DMP_H
#define DMP_H

#include "stdc.h"

typedef struct DMP DMP;

#define dmp_debug _glp_dmp_debug
extern int dmp_debug;
/* debug mode flag */

#define dmp_create_pool _glp_dmp_create_pool
DMP *dmp_create_pool(void);
/* create dynamic memory pool */

#define dmp_talloc(pool, type) \
      ((type *)dmp_get_atom(pool, sizeof(type)))

#define dmp_get_atom _glp_dmp_get_atom
void *dmp_get_atom(DMP *pool, int size);
/* get free atom from dynamic memory pool */

#define dmp_tfree(pool, atom) \
      dmp_free_atom(pool, atom, sizeof(*(atom)))

#define dmp_free_atom _glp_dmp_free_atom
void dmp_free_atom(DMP *pool, void *atom, int size);
/* return atom to dynamic memory pool */

#define dmp_in_use _glp_dmp_in_use
size_t dmp_in_use(DMP *pool);
/* determine how many atoms are still in use */

#define dmp_delete_pool _glp_dmp_delete_pool
void dmp_delete_pool(DMP *pool);
/* delete dynamic memory pool */

#endif

/* eof */
