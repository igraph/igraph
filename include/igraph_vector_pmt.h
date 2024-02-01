/* -*- mode: C -*-  */
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

/*--------------------*/
/* Allocation         */
/*--------------------*/

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, init)(
        TYPE(igraph_vector)* v, igraph_integer_t size);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, init_array)(
        TYPE(igraph_vector)* v, const BASE* data, igraph_integer_t length);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, init_copy)(
        TYPE(igraph_vector) *to, const TYPE(igraph_vector) *from);

#ifndef NOTORDERED
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, init_range)(TYPE(igraph_vector)*v, BASE start, BASE end);
IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t FUNCTION(igraph_vector, init_seq)(TYPE(igraph_vector)*v, BASE from, BASE to);
#endif

IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t FUNCTION(igraph_vector, copy)(
        TYPE(igraph_vector) *to, const TYPE(igraph_vector) *from);
IGRAPH_EXPORT void FUNCTION(igraph_vector, destroy)(TYPE(igraph_vector)* v);

IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_integer_t FUNCTION(igraph_vector, capacity)(const TYPE(igraph_vector)*v);

/*--------------------*/
/* Accessing elements */
/*--------------------*/

#ifndef VECTOR
/**
 * \ingroup vector
 * \define VECTOR
 * \brief Accessing an element of a vector.
 *
 * Usage:
 * \verbatim VECTOR(v)[0] \endverbatim
 * to access the first element of the vector, you can also use this in
 * assignments, like:
 * \verbatim VECTOR(v)[10]=5; \endverbatim
 *
 * Note that there are no range checks right now.
 *
 * \param v The vector object.
 *
 * Time complexity: O(1).
 */
#define VECTOR(v) ((v).stor_begin)
#endif

IGRAPH_EXPORT IGRAPH_DEPRECATED BASE FUNCTION(igraph_vector, e)(const TYPE(igraph_vector)* v, igraph_integer_t pos);
IGRAPH_EXPORT IGRAPH_DEPRECATED BASE* FUNCTION(igraph_vector, e_ptr)(const TYPE(igraph_vector)* v, igraph_integer_t pos);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE FUNCTION(igraph_vector, get)(const TYPE(igraph_vector)* v, igraph_integer_t pos);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE* FUNCTION(igraph_vector, get_ptr)(const TYPE(igraph_vector)* v, igraph_integer_t pos);
IGRAPH_EXPORT void FUNCTION(igraph_vector, set)(TYPE(igraph_vector)* v, igraph_integer_t pos, BASE value);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE FUNCTION(igraph_vector, tail)(const TYPE(igraph_vector) *v);

/*-----------------------*/
/* Initializing elements */
/*-----------------------*/

IGRAPH_EXPORT void FUNCTION(igraph_vector, null)(TYPE(igraph_vector)* v);
IGRAPH_EXPORT void FUNCTION(igraph_vector, fill)(TYPE(igraph_vector)* v, BASE e);

#ifndef NOTORDERED
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, range)(TYPE(igraph_vector)*v, BASE start, BASE end);
#endif

/*-----------------------*/
/* Vector views          */
/*-----------------------*/

IGRAPH_EXPORT const TYPE(igraph_vector) *FUNCTION(igraph_vector, view)(const TYPE(igraph_vector) *v,
                                                                       const BASE *data,
                                                                       igraph_integer_t length);

/*-----------------------*/
/* Copying vectors       */
/*-----------------------*/

IGRAPH_EXPORT void FUNCTION(igraph_vector, copy_to)(const TYPE(igraph_vector) *v, BASE* to);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, update)(TYPE(igraph_vector) *to,
                                                  const TYPE(igraph_vector) *from);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, append)(TYPE(igraph_vector) *to,
                                                  const TYPE(igraph_vector) *from);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, swap)(TYPE(igraph_vector) *v1, TYPE(igraph_vector) *v2);

/*-----------------------*/
/* Exchanging elements   */
/*-----------------------*/

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, swap_elements)(
        TYPE(igraph_vector) *v, igraph_integer_t i, igraph_integer_t j);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, reverse)(TYPE(igraph_vector) *v);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, permute)(TYPE(igraph_vector) *v,
                                                         const igraph_vector_int_t *ind);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, shuffle)(TYPE(igraph_vector) *v);

/*-----------------------*/
/* Vector operations     */
/*-----------------------*/

IGRAPH_EXPORT void FUNCTION(igraph_vector, add_constant)(TYPE(igraph_vector) *v, BASE plus);
IGRAPH_EXPORT void FUNCTION(igraph_vector, scale)(TYPE(igraph_vector) *v, BASE by);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, add)(TYPE(igraph_vector) *v1,
                                               const TYPE(igraph_vector) *v2);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, sub)(TYPE(igraph_vector) *v1,
                                               const TYPE(igraph_vector) *v2);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, mul)(TYPE(igraph_vector) *v1,
                                               const TYPE(igraph_vector) *v2);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, div)(TYPE(igraph_vector) *v1,
                                               const TYPE(igraph_vector) *v2);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, cumsum)(TYPE(igraph_vector) *to,
                                                  const TYPE(igraph_vector) *from);

#ifndef NOABS
    IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, abs)(TYPE(igraph_vector) *v);
#endif

/*------------------------------*/
/* Comparison                   */
/*------------------------------*/

IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_vector, all_e)(const TYPE(igraph_vector) *lhs,
                                                           const TYPE(igraph_vector) *rhs);
#ifndef NOTORDERED
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_vector, all_l)(const TYPE(igraph_vector) *lhs,
                                                           const TYPE(igraph_vector) *rhs);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_vector, all_g)(const TYPE(igraph_vector) *lhs,
                                                           const TYPE(igraph_vector) *rhs);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_vector, all_le)(const TYPE(igraph_vector) *lhs,
                                                            const TYPE(igraph_vector) *rhs);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_vector, all_ge)(const TYPE(igraph_vector) *lhs,
                                                            const TYPE(igraph_vector) *rhs);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE int FUNCTION(igraph_vector, lex_cmp)(
        const TYPE(igraph_vector) *lhs, const TYPE(igraph_vector) *rhs);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE int FUNCTION(igraph_vector, colex_cmp)(
        const TYPE(igraph_vector) *lhs, const TYPE(igraph_vector) *rhs);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE int FUNCTION(igraph_vector, lex_cmp_untyped)(const void *lhs, const void *rhs);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE int FUNCTION(igraph_vector, colex_cmp_untyped)(const void *lhs, const void *rhs);
#endif

/*------------------------------*/
/* Finding minimum and maximum  */
/*------------------------------*/

#ifndef NOTORDERED
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE FUNCTION(igraph_vector, min)(const TYPE(igraph_vector)* v);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE FUNCTION(igraph_vector, max)(const TYPE(igraph_vector)* v);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_integer_t FUNCTION(igraph_vector, which_min)(const TYPE(igraph_vector)* v);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_integer_t FUNCTION(igraph_vector, which_max)(const TYPE(igraph_vector)* v);
IGRAPH_EXPORT void FUNCTION(igraph_vector, minmax)(
        const TYPE(igraph_vector) *v, BASE *min, BASE *max);
IGRAPH_EXPORT void FUNCTION(igraph_vector, which_minmax)(
        const TYPE(igraph_vector) *v, igraph_integer_t *which_min, igraph_integer_t *which_max);
#endif

/*-------------------*/
/* Vector properties */
/*-------------------*/

IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_vector, empty)(const TYPE(igraph_vector)* v);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_integer_t FUNCTION(igraph_vector, size)(const TYPE(igraph_vector)* v);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_vector, isnull)(const TYPE(igraph_vector) *v);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE FUNCTION(igraph_vector, sum)(const TYPE(igraph_vector) *v);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_real_t FUNCTION(igraph_vector, sumsq)(const TYPE(igraph_vector) *v);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE FUNCTION(igraph_vector, prod)(const TYPE(igraph_vector) *v);
#ifndef NOTORDERED
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_vector, isininterval)(const TYPE(igraph_vector) *v,
                                                                  BASE low, BASE high);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_vector, any_smaller)(const TYPE(igraph_vector) *v,
                                                                 BASE limit);
#endif
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_vector, is_equal)(const TYPE(igraph_vector) *lhs,
                                                              const TYPE(igraph_vector) *rhs);
#ifndef NOTORDERED
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_real_t FUNCTION(igraph_vector, maxdifference)(const TYPE(igraph_vector) *m1,
                                                                   const TYPE(igraph_vector) *m2);
#endif

/*------------------------*/
/* Searching for elements */
/*------------------------*/

IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_vector, contains)(const TYPE(igraph_vector) *v, BASE e);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_vector, search)(
        const TYPE(igraph_vector) *v, igraph_integer_t from, BASE what, igraph_integer_t *pos);
#ifndef NOTORDERED
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_vector, binsearch_slice)(
        const TYPE(igraph_vector) *v, BASE what, igraph_integer_t *pos,
        igraph_integer_t start, igraph_integer_t end);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_vector, binsearch)(
        const TYPE(igraph_vector) *v, BASE what, igraph_integer_t *pos);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_vector, binsearch2)(
        const TYPE(igraph_vector) *v, BASE what);
#endif

/*------------------------*/
/* Resizing operations    */
/*------------------------*/

IGRAPH_EXPORT void FUNCTION(igraph_vector, clear)(TYPE(igraph_vector)* v);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, resize)(
        TYPE(igraph_vector)* v, igraph_integer_t new_size);
IGRAPH_EXPORT void FUNCTION(igraph_vector, resize_min)(TYPE(igraph_vector)*v);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, reserve)(
        TYPE(igraph_vector)* v, igraph_integer_t capacity);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, push_back)(TYPE(igraph_vector)* v, BASE e);
IGRAPH_EXPORT BASE FUNCTION(igraph_vector, pop_back)(TYPE(igraph_vector)* v);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, insert)(
        TYPE(igraph_vector) *v, igraph_integer_t pos, BASE value);
IGRAPH_EXPORT void FUNCTION(igraph_vector, remove)(
        TYPE(igraph_vector) *v, igraph_integer_t elem);
IGRAPH_EXPORT void FUNCTION(igraph_vector, remove_fast)(
        TYPE(igraph_vector) *v, igraph_integer_t elem);
IGRAPH_EXPORT void FUNCTION(igraph_vector, remove_section)(
        TYPE(igraph_vector) *v, igraph_integer_t from, igraph_integer_t to);

/*-----------*/
/* Sorting   */
/*-----------*/

#ifndef NOTORDERED

IGRAPH_EXPORT void FUNCTION(igraph_vector, sort)(TYPE(igraph_vector) *v);
IGRAPH_EXPORT void FUNCTION(igraph_vector, reverse_sort)(TYPE(igraph_vector) *v);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, qsort_ind)(
        const TYPE(igraph_vector) *v, igraph_vector_int_t *inds, igraph_order_t order);

#endif

/*-----------*/
/* Printing  */
/*-----------*/

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, print)(const TYPE(igraph_vector) *v);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, fprint)(const TYPE(igraph_vector) *v, FILE *file);

#ifdef OUT_FORMAT
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, printf)(const TYPE(igraph_vector) *v, const char *format);
#endif /* OUT_FORMAT */

/*----------------------------------------*/
/* Operations specific to complex vectors */
/*----------------------------------------*/

#ifdef BASE_COMPLEX

IGRAPH_EXPORT igraph_error_t igraph_vector_complex_real(const igraph_vector_complex_t *v,
                                                        igraph_vector_t *real);
IGRAPH_EXPORT igraph_error_t igraph_vector_complex_imag(const igraph_vector_complex_t *v,
                                                        igraph_vector_t *imag);
IGRAPH_EXPORT igraph_error_t igraph_vector_complex_realimag(const igraph_vector_complex_t *v,
                                                            igraph_vector_t *real,
                                                            igraph_vector_t *imag);
IGRAPH_EXPORT igraph_error_t igraph_vector_complex_create(igraph_vector_complex_t *v,
                                                          const igraph_vector_t *real,
                                                          const igraph_vector_t *imag);
IGRAPH_EXPORT igraph_error_t igraph_vector_complex_create_polar(igraph_vector_complex_t *v,
                                                                const igraph_vector_t *r,
                                                                const igraph_vector_t *theta);

IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t igraph_vector_complex_all_almost_e(const igraph_vector_complex_t *lhs,
                                                               const igraph_vector_complex_t *rhs,
                                                               igraph_real_t eps);

#endif /* BASE_COMPLEX */

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, init_real)(TYPE(igraph_vector)*v, int no, ...);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, init_int)(TYPE(igraph_vector)*v, int no, ...);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, init_real_end)(TYPE(igraph_vector)*v, double endmark, ...);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, init_int_end)(TYPE(igraph_vector)*v, int endmark, ...);

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, move_interval)(
        TYPE(igraph_vector) *v, igraph_integer_t begin, igraph_integer_t end,
        igraph_integer_t to);
IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t FUNCTION(igraph_vector, move_interval2)(
        TYPE(igraph_vector) *v, igraph_integer_t begin, igraph_integer_t end,
        igraph_integer_t to);
#ifndef NOTORDERED
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, filter_smaller)(TYPE(igraph_vector) *v, BASE elem);
#endif
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, get_interval)(
        const TYPE(igraph_vector) *v, TYPE(igraph_vector) *res,
        igraph_integer_t from, igraph_integer_t to);
#ifndef NOTORDERED
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, difference_sorted)(const TYPE(igraph_vector) *v1,
                                                             const TYPE(igraph_vector) *v2, TYPE(igraph_vector) *result);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, intersect_sorted)(const TYPE(igraph_vector) *v1,
                                                            const TYPE(igraph_vector) *v2, TYPE(igraph_vector) *result);
#endif
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, index)(const TYPE(igraph_vector) *v,
                                                 TYPE(igraph_vector) *newv,
                                                 const igraph_vector_int_t *idx);

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, index_int)(TYPE(igraph_vector) *v,
                                                     const igraph_vector_int_t *idx);
