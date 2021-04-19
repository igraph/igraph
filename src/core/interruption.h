/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_INTERRUPT_INTERNAL_H
#define IGRAPH_INTERRUPT_INTERNAL_H

#include "igraph_decls.h"
#include "igraph_interrupt.h"
#include "config.h"

__BEGIN_DECLS

extern IGRAPH_THREAD_LOCAL igraph_interruption_handler_t
*igraph_i_interruption_handler;

/**
 * \define IGRAPH_ALLOW_INTERRUPTION
 * \brief
 *
 * This macro should be called when interruption is allowed.  It calls
 * \ref igraph_allow_interruption() with the proper parameters and if that returns
 * anything but \c IGRAPH_SUCCESS then
 * the macro returns the "calling" function as well, with the proper
 * error code (\c IGRAPH_INTERRUPTED).
 */

#define IGRAPH_ALLOW_INTERRUPTION() \
    do { \
        if (igraph_i_interruption_handler) { if (igraph_allow_interruption(NULL) != IGRAPH_SUCCESS) return IGRAPH_INTERRUPTED; \
        } } while (0)

#define IGRAPH_ALLOW_INTERRUPTION_NORETURN() \
    do { \
        if (igraph_i_interruption_handler) { igraph_allow_interruption(NULL); } \
    } while (0)

__END_DECLS

#endif

