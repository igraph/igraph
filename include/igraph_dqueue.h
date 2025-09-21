/*
   igraph library.
   Copyright (C) 2009-2025  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef IGRAPH_DQUEUE_H
#define IGRAPH_DQUEUE_H

#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"

IGRAPH_BEGIN_C_DECLS

/* -------------------------------------------------- */
/* double ended queue, very useful                    */
/* -------------------------------------------------- */

#define BASE_IGRAPH_REAL
#include "igraph_pmt.h"
#include "igraph_dqueue_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_IGRAPH_REAL

#define BASE_CHAR
#include "igraph_pmt.h"
#include "igraph_dqueue_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_CHAR

#define BASE_BOOL
#include "igraph_pmt.h"
#include "igraph_dqueue_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_BOOL

#define BASE_INT
#include "igraph_pmt.h"
#include "igraph_dqueue_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_INT

#define IGRAPH_DQUEUE_NULL { 0,0,0,0 }
#define IGRAPH_DQUEUE_INIT_FINALLY(q, capacity) \
    do { IGRAPH_CHECK(igraph_dqueue_init(q, capacity)); \
        IGRAPH_FINALLY(igraph_dqueue_destroy, q); } while (0)
#define IGRAPH_DQUEUE_INT_INIT_FINALLY(q, capacity) \
    do { IGRAPH_CHECK(igraph_dqueue_int_init(q, capacity)); \
        IGRAPH_FINALLY(igraph_dqueue_int_destroy, q); } while (0)

IGRAPH_END_C_DECLS

#endif
