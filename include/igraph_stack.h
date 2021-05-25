/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef IGRAPH_STACK_H
#define IGRAPH_STACK_H

#include "igraph_decls.h"
#include "igraph_types.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Plain stack                                        */
/* -------------------------------------------------- */

#define BASE_IGRAPH_REAL
#include "igraph_pmt.h"
#include "igraph_stack_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_IGRAPH_REAL

#define BASE_LONG
#include "igraph_pmt.h"
#include "igraph_stack_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_LONG

#define BASE_INT
#include "igraph_pmt.h"
#include "igraph_stack_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_INT

#define BASE_CHAR
#include "igraph_pmt.h"
#include "igraph_stack_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_CHAR

#define BASE_BOOL
#include "igraph_pmt.h"
#include "igraph_stack_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_BOOL

#define BASE_PTR
#include "igraph_pmt.h"
#include "igraph_stack_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_PTR

#define IGRAPH_STACK_NULL { 0,0,0 }

IGRAPH_EXPORT void igraph_stack_ptr_free_all(igraph_stack_ptr_t* s);
IGRAPH_EXPORT void igraph_stack_ptr_destroy_all(igraph_stack_ptr_t* s);

__END_DECLS

#endif
