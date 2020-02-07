/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#ifndef IGRAPH_GLPK_SUPPORT_H
#define IGRAPH_GLPK_SUPPORT_H

#include "config.h"

/* Note: only files calling the GLPK routines directly need to
   include this header.
*/

#ifdef HAVE_GLPK

#include <glpk.h>

int igraph_i_glpk_check(int retval, const char* message);
void igraph_i_glpk_interruption_hook(glp_tree *tree, void *info);
#define IGRAPH_GLPK_CHECK(func, message) do {\
        int igraph_i_ret = igraph_i_glpk_check(func, message); \
        if (IGRAPH_UNLIKELY(igraph_i_ret != 0)) {\
            return igraph_i_ret; \
        } } while (0)

#endif

#endif
