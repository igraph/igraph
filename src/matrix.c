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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "types.h"
#include "memory.h"
#include "random.h"
#include "error.h"

#include <assert.h>
#include <string.h> 		/* memcpy & co. */
#include <stdlib.h>

/**
 * \section igraph_matrix_constructor_and_destructor Matrix constructors and
 * destructors
 */

/**
 * \ingroup matrix
 * \function igraph_matrix_init
 * \brief Initializes a matrix.
 * 
 * </para><para>
 * Every matrix needs to be initialized before using it, this is done
 * by calling this function. A matrix has to be destroyed if it is not
 * needed any more, see \ref igraph_matrix_destroy().
 * \param m Pointer to a not yet initialized matrix object to be
 *        initialized. 
 * \param nrow The number of rows in the matrix.
 * \param ncol The number of columns in the matrix.
 * \return Error code.
 *
 * Time complexity: usually O(n), 
 * n is the
 * number of elements in the matrix.
 */

int igraph_matrix_init(igraph_matrix_t *m, long int nrow, long int ncol) {
  int ret1;
  ret1=igraph_vector_init(&m->data, nrow*ncol);
  m->nrow=nrow;
  m->ncol=ncol;
  return ret1;
}

/** 
 * \ingroup matrix
 * \function igraph_matrix_destroy
 * \brief Destroys a matrix object.
 * 
 * </para><para>
 * This function frees all the memory allocated for a matrix
 * object. The destroyed object needs to be reinitialized before using
 * it again.
 * \param m The matrix to destroy.
 * 
 * Time complexity: operating system dependent.
 */ 

void igraph_matrix_destroy(igraph_matrix_t *m) {
  igraph_vector_destroy(&m->data);
}

/**
 * \section igraph_matrix_accessing_elements Accessing elements of a matrix
 */

/**
 * \ingroup matrix
 * \function igraph_matrix_resize
 * \brief Resizes a matrix.
 *
 * </para><para>
 * This function resizes a matrix by adding more elements to it.
 * The matrix contains arbitrary data after resizing it.
 * Ie. after calling this function you cannot expect that element
 * (i,j) in the matrix remains the
 * same as before.  
 * \param m Pointer to an already initialized matrix object.
 * \param nrow The number of rows in the resized matrix.
 * \param ncol The number of columns in the resized matrix.
 * \return Error code.
 * 
 * Time complexity: O(1) if the
 * matrix gets smaller, usually O(n)
 * if it gets larger, n is the 
 * number of elements in the resized matrix.
 */

int igraph_matrix_resize(igraph_matrix_t *m, long int nrow, long int ncol) {
  igraph_vector_resize(&m->data, nrow*ncol);
  m->nrow=nrow;
  m->ncol=ncol;
  return 0;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_size
 * \brief The number of elements in a matrix.
 * 
 * \param m Pointer to an initialized matrix object.
 * \return The size of the matrix.
 *
 * Time complexity: O(1).
 */

long int igraph_matrix_size(const igraph_matrix_t *m) {
  return (m->nrow) * (m->ncol);
}

/**
 * \ingroup matrix
 * \function igraph_matrix_nrow
 * \brief The number of rows in a matrix.
 * 
 * \param m Pointer to an initialized matrix object.
 * \return The number of rows in the matrix.
 * 
 * Time complexity: O(1).
 */

long int igraph_matrix_nrow(const igraph_matrix_t *m) {
  return m->nrow;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_ncol
 * \brief The number of columns in a matrix.
 * 
 * \param m Pointer to an initialized matrix object.
 * \return The number of columns in the matrix.
 * 
 * Time complexity: O(1).
 */

long int igraph_matrix_ncol(const igraph_matrix_t *m) {
  return m->ncol;
}

/** 
 * \ingroup matrix
 * \brief Copies a matrix to a regular C array.
 *
 * </para><para>
 * The matrix is copied columnwise, as this is the format most
 * programs and languages use.
 * The C array should be of sufficient size, there are (of course) not
 * range checks done.
 * \param m Pointer to an initialized matrix object.
 * \param to Pointer to a C array, the place to copy the data to.
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the number of 
 * elements in the matrix.
 */

int igraph_matrix_copy_to(const igraph_matrix_t *m, igraph_real_t *to) {
  igraph_vector_copy_to(&m->data, to);
  return 0;
}

/** 
 * \ingroup matrix
 * \brief Sets all element in a matrix to zero.
 * 
 * \param m Pointer to an initialized matrix object.
 * \return Error code, always returns with success.
 * 
 * Time complexity: O(n),
 * n is the number of  elements in
 * the matrix. 
 */

int igraph_matrix_null(igraph_matrix_t *m) {
  igraph_vector_null(&m->data);
  return 0;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_add_cols
 * \brief Adds columns to a matrix.
 * \param m The matrix object.
 * \param n The number of columns to add.
 * \return Error code, \c IGRAPH_ENOMEM if there is
 *   not enough memory to perform the operation.
 *
 * Time complexity: linear with the number of elements of the new,
 * resized matrix.
 */

int igraph_matrix_add_cols(igraph_matrix_t *m, long int n) {
  igraph_matrix_resize(m, m->nrow, m->ncol+n);
  return 0;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_add_rows
 * \brief Adds rows to a matrix.
 * \param m The matrix object.
 * \param n The number of rows to add.
 * \return Error code, \c IGRAPH_ENOMEM if there
 *   isn't enough memory for the operation.
 * 
 * Time complexity: linear with the number of elements of the new,
 * resized, matrix.
 */

int igraph_matrix_add_rows(igraph_matrix_t *m, long int n) {
  long int i;
  igraph_vector_resize(&m->data, (m->ncol)*(m->nrow+n));
  for (i=m->ncol-1; i>=0; i--) {
    igraph_vector_move_interval(&m->data, (m->nrow)*i, (m->nrow)*(i+1),
			 (m->nrow+n)*i);
  }
  m->nrow += n;
  return 0;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_remove_col
 * \brief Removes a column from a matrix.
 * 
 * \param m The matrix object.
 * \param col The column to remove.
 * \return Error code, always returns with success. 
 * 
 * Time complexity: linear with the number of elements of the new,
 * resized matrix.
 */

int igraph_matrix_remove_col(igraph_matrix_t *m, long int col) {
  igraph_vector_remove_section(&m->data, (m->nrow)*col, (m->nrow)*(col+1));
  m->ncol--;
  return 0;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_permdelete_rows
 * \brief Removes columns from a matrix (for internal use).
 * 
 * Time complexity: linear with the number of elements of the original
 * matrix. 
 */

int igraph_matrix_permdelete_rows(igraph_matrix_t *m, long int *index, long int nremove) {
  long int i, j;
  for (i=0; i<m->ncol; i++) {
    for (j=0; j<m->nrow; j++) {
      if (index[j] != 0) {
	MATRIX(*m, index[j]-1, i) = MATRIX(*m, j, i);
      }
    }
  }
  igraph_matrix_resize(m, m->nrow-nremove, m->ncol);

  return 0;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_delete_rows_neg
 * \brief Removes columns from a matrix (for internal use).
 * 
 * Time complexity: linear with the number of elements of the original
 * matrix. 
 */

int igraph_matrix_delete_rows_neg(igraph_matrix_t *m, igraph_vector_t *neg, long int nremove) {
  long int i, j, idx=0;
  for (i=0; i<m->ncol; i++) {
    for (j=0; j<m->nrow; j++) {
      if (VECTOR(*neg)[j] >= 0) {
	MATRIX(*m, idx++, i) = MATRIX(*m, j, i);
      } 
    }
    idx=0;
  }
  igraph_matrix_resize(m, m->nrow-nremove, m->ncol);

  return 0;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_copy
 * \brief Copies a matrix.
 *
 * </para><para>
 * Creates a matrix object by copying another one.
 * \param to Pointer to an uninitialized matrix object.
 * \param from The initialized matrix object to copy.
 * \return Error code, \c IGRAPH_ENOMEM if there
 *   isn't enough memory to allocate the new matrix.
 * 
 * Time complexity: O(n), the number
 * of elements in the matrix.
 */

int igraph_matrix_copy(igraph_matrix_t *to, const igraph_matrix_t *from) {
  to->nrow = from->nrow;
  to->ncol = from->ncol;
  return igraph_vector_copy(&to->data, &from->data);
}

/**
 * \function igraph_matrix_max
 * 
 * Returns the maximal element of a matrix.
 * \param m The matrix object.
 * \return The maximum element. For empty matrix the returned value is
 * undefined. 
 * 
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(n), the number of elements in the matrix.
 */

igraph_real_t igraph_matrix_max(const igraph_matrix_t *m) {
  return igraph_vector_max(&m->data);
}

/**
 * \function igraph_matrix_multiply
 * 
 * Multiplies each element of the matrix by a constant.
 * \param m The matrix.
 * \param by The constant.
 *
 * Added in version 0.2.</para><para>
 * 
 * Time complexity: O(n), the number of elements in the matrix.
 */

void igraph_matrix_multiply(igraph_matrix_t *m, igraph_real_t by) {
  igraph_vector_multiply(&m->data, by);
}

int igraph_matrix_select_rows(const igraph_matrix_t *m, igraph_matrix_t *res, 
			      const igraph_vector_t *rows) {
  long int norows=igraph_vector_size(rows);
  long int i, j, ncols=igraph_matrix_ncol(m);
  
  IGRAPH_CHECK(igraph_matrix_resize(res, norows, ncols));
  for (i=0; i<norows; i++) {
    for (j=0; j<ncols; j++) {
      MATRIX(*res, i, j) = MATRIX(*m, (long int)VECTOR(*rows)[i], j);
    }
  }
  
  return 0;
}

int igraph_matrix_get_col(const igraph_matrix_t *m, igraph_vector_t *res,
			  long int index) {
  long int nrow=igraph_matrix_nrow(m);

  IGRAPH_CHECK(igraph_vector_get_interval(&m->data, res, 
					  nrow*index, nrow*(index+1)));
  return 0;
}

igraph_real_t igraph_matrix_sum(const igraph_matrix_t *m) {
  return igraph_vector_sum(&m->data);
}
