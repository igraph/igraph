/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include "types.h"
#include "memory.h"
#include "random.h"
#include "error.h"

#include <assert.h>
#include <string.h> 		/* memcpy & co. */
#include <stdlib.h>

/**
 * \ingroup matrix
 * \brief Initializes a matrix.
 * 
 * Every matrix needs to be initialized before using it, this is done
 * by calling this function. A matrix has to be destroyed if it is not
 * needed any more, see matrix_destroy().
 * @param m Pointer to a not yet initialized matrix object to be
 *        initialized. 
 * @param nrow The number of rows in the matrix.
 * @param ncol The number of columns in the matrix.
 * @return Error code.
 *
 * Time complexity: ususally <code>O(n)</code>, <code>n</code> is the
 * number of elements in the matrix.
 */
int matrix_init(matrix_t *m, long int nrow, long int ncol) {
  int ret1;
  ret1=vector_init(&m->data, nrow*ncol);
  m->nrow=nrow;
  m->ncol=ncol;
  return ret1;
}

/** 
 * \ingroup matrix
 * \brief Destroys a matrix object.
 * 
 * This function frees all the memory allocated for a matrix
 * object. The destroyed object needs to be reinitialized before using
 * it again.
 * @param m The matrix to destroy.
 * @return Error code.
 *
 * Time complexity: operating system dependent.
 */ 

void matrix_destroy(matrix_t *m) {
  vector_destroy(&m->data);
}

/**
 * \ingroup matrix
 * \brief Resizes a matrix.
 *
 * This function resizes a matrix by adding more elements to it.
 * The matrix contains arbitrary data after resizing it.
 * Ie. after calling this function you cannot expect that element
 * <code>(i,j)</code> in the matrix remains the same as before. 
 * @param m Pointer to an already initialized matrix object.
 * @param nrow The number of rows in the resized matrix.
 * @param ncol The number of columns in the resized matrix.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code> if the matrix gets smaller,
 * usually <code>O(n)</code> if it gets larger, <code>n</code> is the
 * number of elements in the resized matrix.
 */

int matrix_resize(matrix_t *m, long int nrow, long int ncol) {
  vector_resize(&m->data, nrow*ncol);
  m->nrow=nrow;
  m->ncol=ncol;
  return 0;
}

/**
 * \ingroup matrix
 * \brief The number of elements in a matrix.
 * 
 * @param m Pointer to an initialized matrix object.
 * @return The size of the matrix.
 *
 * Time complexity: <code>O(1)</code>.
 */

long int matrix_size(matrix_t *m) {
  return (m->nrow) * (m->ncol);
}

/**
 * \ingroup matrix
 * \brief The number of rows in a matrix.
 * 
 * @param m Pointer to an initialized matrix object.
 * @return The number of rows in the matrix.
 * 
 * Time complexity: <code>O(1)</code>.
 */

long int matrix_nrow(matrix_t *m) {
  return m->nrow;
}

/**
 * \ingroup matrix
 * \brief The number of columns in a matrix.
 * 
 * @param m Pointer to an initialized matrix object.
 * @return The number of columns in the matrix.
 * 
 * Time complexity: <code>O(1)</code>.
 */

long int matrix_ncol(matrix_t *m) {
  return m->ncol;
}

/** 
 * \ingroup matrix
 * \brief Copies a matrix to a regular C array.
 *
 * The matrix is copied columnwise, as this is the format most
 * programs and languages use.
 * The C array should be of sufficient size, there are (of course) not
 * range checks done.
 * @param m Pointer to an initialized matrix object.
 * @param to Pointer to a C array, the place to copy the data to.
 * @return Error code.
 *
 * Time complexity: <code>O(n)</code>, <code>n</code> is the number of 
 * elements in the matrix.
 */

int matrix_copy_to(matrix_t *m, real_t *to) {
  vector_copy_to(&m->data, to);
  return 0;
}

/** 
 * \ingroup matrix
 * \brief Sets all element in a matrix to zero.
 * 
 * @param m Pointer to an initialized matrix object.
 * @return Error code.
 * 
 * Time complexity: <code>O(n)</code>, <code>n</code> is the number of 
 * elements in the matrix.
 */

int matrix_null(matrix_t *m) {
  vector_null(&m->data);
  return 0;
}

/**
 * \ingroup matrix
 * \brief Adds columns to a matrix.
 */

int matrix_add_cols(matrix_t *m, long int n) {
  matrix_resize(m, m->nrow, m->ncol+n);
  return 0;
}

/**
 * \ingroup matrix
 * \brief Adds rows to a matrix.
 */

int matrix_add_rows(matrix_t *m, long int n) {
  long int i;
  vector_resize(&m->data, (m->ncol)*(m->nrow+n));
  for (i=m->ncol-1; i>=0; i--) {
    vector_move_interval(&m->data, (m->nrow)*i, (m->nrow)*(i+1),
			 (m->nrow+n)*i);
  }
  m->nrow += n;
  return 0;
}

/**
 * \ingroup matrix
 * \brief Removes a column from a matrix.
 */

int matrix_remove_col(matrix_t *m, long int col) {
  vector_remove_section(&m->data, (m->nrow)*col, (m->nrow)*(col+1));
  m->ncol--;
  return 0;
}

/**
 * \ingroup matrix
 * \brief Removes columns from a matrix (for internal use).
 */

int matrix_permdelete_rows(matrix_t *m, long int *index, long int nremove) {
  long int i, j;
  for (i=0; i<m->ncol; i++) {
    for (j=0; j<m->nrow; j++) {
      if (index[j] != 0) {
	MATRIX(*m, index[j]-1, i) = MATRIX(*m, j, i);
      }
    }
  }
  matrix_resize(m, m->nrow-nremove, m->ncol);

  return 0;
}

/**
 * \ingroup matrix
 * \brief Removes columns from a matrix (for internal use).
 */

int matrix_delete_rows_neg(matrix_t *m, vector_t *neg, long int nremove) {
  long int i, j, idx=0;
  for (i=0; i<m->ncol; i++) {
    for (j=0; j<m->nrow; j++) {
      if (VECTOR(*neg)[j] >= 0) {
	MATRIX(*m, idx++, i) = MATRIX(*m, j, i);
      } 
    }
    idx=0;
  }
  matrix_resize(m, m->nrow-nremove, m->ncol);

  return 0;
}

/**
 * \ingroup matrix
 * \brief Copies a matrix.
 */

int matrix_copy(matrix_t *to, matrix_t *from) {
  to->nrow = from->nrow;
  to->ncol = from->ncol;
  return vector_copy(&to->data, &from->data);
}

