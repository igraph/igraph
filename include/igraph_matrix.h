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

#ifndef IGRAPH_MATRIX_H
#define IGRAPH_MATRIX_H

#include "igraph_decls.h"
#include "igraph_vector.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Matrix, very similar to vector                     */
/* -------------------------------------------------- */

#define BASE_IGRAPH_REAL
#include "igraph_pmt.h"
#include "igraph_matrix_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_IGRAPH_REAL

#define BASE_INT
#include "igraph_pmt.h"
#include "igraph_matrix_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_INT

#define BASE_LONG
#include "igraph_pmt.h"
#include "igraph_matrix_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_LONG

#define BASE_CHAR
#include "igraph_pmt.h"
#include "igraph_matrix_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_CHAR

#define BASE_BOOL
#include "igraph_pmt.h"
#include "igraph_matrix_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_BOOL

#define BASE_COMPLEX
#include "igraph_pmt.h"
#include "igraph_matrix_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_COMPLEX

#define IGRAPH_MATRIX_NULL { IGRAPH_VECTOR_NULL, 0, 0 }
#define IGRAPH_MATRIX_INIT_FINALLY(m, nr, nc) \
    do { IGRAPH_CHECK(igraph_matrix_init(m, nr, nc)); \
        IGRAPH_FINALLY(igraph_matrix_destroy, m); } while (0)

/**
 * \ingroup matrix
 * \define MATRIX
 * \brief Accessing an element of a matrix.
 *
 * Note that there are no range checks right now.
 * This functionality might be redefined as a proper function later.
 * \param m The matrix object.
 * \param i The index of the row, starting with zero.
 * \param j The index of the column, starting with zero.
 *
 * Time complexity: O(1).
 */
#define MATRIX(m,i,j) ((m).data.stor_begin[(m).nrow*(j)+(i)])

IGRAPH_EXPORT igraph_bool_t igraph_matrix_all_e_tol(const igraph_matrix_t *lhs,
                                                    const igraph_matrix_t *rhs,
                                                    igraph_real_t tol);

IGRAPH_EXPORT int igraph_matrix_zapsmall(igraph_matrix_t *m, igraph_real_t tol);

__END_DECLS

#endif
