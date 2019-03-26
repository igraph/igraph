/* -*- mode: C -*-  */
/* vim:set ts=2 sw=2 sts=2 et: */
/* 
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "config.h"

#ifdef HAVE_GLPK

#include "igraph_types.h"
#include "igraph_error.h"
#include "igraph_interrupt_internal.h"
#include <glpk.h>
#include <memory.h>
#include <stdio.h>

void igraph_i_glpk_interruption_hook(glp_tree *tree, void *info) {
  IGRAPH_UNUSED(tree);
  IGRAPH_UNUSED(info);
  IGRAPH_ALLOW_INTERRUPTION_NORETURN();
}

int igraph_i_glpk_check(int retval, const char* message) {
  char* code = "none";
  char message_and_code[4096];

  if (retval == IGRAPH_SUCCESS)
    return IGRAPH_SUCCESS;

  /* handle errors */
#define HANDLE_CODE(c) case c: code = #c; retval = IGRAPH_##c; break;
#define HANDLE_CODE2(c) case c: code = #c; retval = IGRAPH_FAILURE; break;
  switch (retval) {
    HANDLE_CODE(GLP_EBOUND);
    HANDLE_CODE(GLP_EROOT);
    HANDLE_CODE(GLP_ENOPFS);
    HANDLE_CODE(GLP_ENODFS);
    HANDLE_CODE(GLP_EFAIL);
    HANDLE_CODE(GLP_EMIPGAP);
    HANDLE_CODE(GLP_ETMLIM);
    HANDLE_CODE(GLP_ESTOP);

    HANDLE_CODE2(GLP_EBADB);
    HANDLE_CODE2(GLP_ESING);
    HANDLE_CODE2(GLP_ECOND);
    HANDLE_CODE2(GLP_EOBJLL);
    HANDLE_CODE2(GLP_EOBJUL);
    HANDLE_CODE2(GLP_EITLIM);

    default:
      IGRAPH_ERROR("unknown GLPK error", IGRAPH_FAILURE);
  }
#undef HANDLE_CODE

  sprintf(message_and_code, "%s (%s)", message, code);
  IGRAPH_ERROR(message_and_code, retval);
}

#endif

#ifdef USING_R

int igraph_glpk_dummy() {
  return 'b' + 'a' + 's' + 's' + 'z' + 'a' + 't' + 'o' + 'k' +
    'm' + 'e' + 'g';
}

#endif
