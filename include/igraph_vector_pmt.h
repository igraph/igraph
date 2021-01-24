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

DECLDIR igraph_error_t FUNCTION(igraph_vector, init)(TYPE(igraph_vector)* v, igraph_long_t size);
DECLDIR igraph_error_t FUNCTION(igraph_vector, init_copy)(TYPE(igraph_vector)* v,
        const BASE* data, igraph_long_t length);
DECLDIR igraph_error_t FUNCTION(igraph_vector, init_seq)(TYPE(igraph_vector)*v, BASE from, BASE to);
DECLDIR igraph_error_t FUNCTION(igraph_vector, copy)(TYPE(igraph_vector) *to,
        const TYPE(igraph_vector) *from);
DECLDIR void FUNCTION(igraph_vector, destroy)(TYPE(igraph_vector)* v);

DECLDIR igraph_error_t FUNCTION(igraph_vector, capacity)(const TYPE(igraph_vector)*v);

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
 * This functionality might be redefined later as a real function
 * instead of a <code>#define</code>.
 * \param v The vector object.
 *
 * Time complexity: O(1).
 */
#define VECTOR(v) ((v).stor_begin)
#endif

DECLDIR BASE FUNCTION(igraph_vector, e)(const TYPE(igraph_vector)* v, igraph_long_t pos);
DECLDIR BASE* FUNCTION(igraph_vector, e_ptr)(const TYPE(igraph_vector)* v, igraph_long_t pos);
DECLDIR void FUNCTION(igraph_vector, set)(TYPE(igraph_vector)* v, igraph_long_t pos, BASE value);
DECLDIR BASE FUNCTION(igraph_vector, tail)(const TYPE(igraph_vector) *v);

/*-----------------------*/
/* Initializing elements */
/*-----------------------*/

DECLDIR void FUNCTION(igraph_vector, null)(TYPE(igraph_vector)* v);
DECLDIR void FUNCTION(igraph_vector, fill)(TYPE(igraph_vector)* v, BASE e);

/*-----------------------*/
/* Vector views          */
/*-----------------------*/

DECLDIR const TYPE(igraph_vector) *FUNCTION(igraph_vector, view)(const TYPE(igraph_vector) *v,
        const BASE *data,
        igraph_long_t length);

/*-----------------------*/
/* Copying vectors       */
/*-----------------------*/

DECLDIR void FUNCTION(igraph_vector, copy_to)(const TYPE(igraph_vector) *v, BASE* to);
DECLDIR igraph_error_t FUNCTION(igraph_vector, update)(TYPE(igraph_vector) *to,
        const TYPE(igraph_vector) *from);
DECLDIR igraph_error_t FUNCTION(igraph_vector, append)(TYPE(igraph_vector) *to,
        const TYPE(igraph_vector) *from);
DECLDIR igraph_error_t FUNCTION(igraph_vector, swap)(TYPE(igraph_vector) *v1, TYPE(igraph_vector) *v2);

/*-----------------------*/
/* Exchanging elements   */
/*-----------------------*/

DECLDIR igraph_error_t FUNCTION(igraph_vector, swap_elements)(TYPE(igraph_vector) *v,
        igraph_long_t i, igraph_long_t j);
DECLDIR igraph_error_t FUNCTION(igraph_vector, reverse)(TYPE(igraph_vector) *v);
DECLDIR igraph_error_t FUNCTION(igraph_vector, shuffle)(TYPE(igraph_vector) *v);

/*-----------------------*/
/* Vector operations     */
/*-----------------------*/

DECLDIR void FUNCTION(igraph_vector, add_constant)(TYPE(igraph_vector) *v, BASE plus);
DECLDIR void FUNCTION(igraph_vector, scale)(TYPE(igraph_vector) *v, BASE by);
DECLDIR igraph_error_t FUNCTION(igraph_vector, add)(TYPE(igraph_vector) *v1,
        const TYPE(igraph_vector) *v2);
DECLDIR igraph_error_t FUNCTION(igraph_vector, sub)(TYPE(igraph_vector) *v1,
        const TYPE(igraph_vector) *v2);
DECLDIR igraph_error_t FUNCTION(igraph_vector, mul)(TYPE(igraph_vector) *v1,
        const TYPE(igraph_vector) *v2);
DECLDIR igraph_error_t FUNCTION(igraph_vector, div)(TYPE(igraph_vector) *v1,
        const TYPE(igraph_vector) *v2);
DECLDIR igraph_error_t FUNCTION(igraph_vector, cumsum)(TYPE(igraph_vector) *to,
        const TYPE(igraph_vector) *from);

#ifndef NOABS
    DECLDIR igraph_error_t FUNCTION(igraph_vector, abs)(TYPE(igraph_vector) *v);
#endif

/*------------------------------*/
/* Comparison                   */
/*------------------------------*/

DECLDIR igraph_bool_t FUNCTION(igraph_vector, all_e)(const TYPE(igraph_vector) *lhs,
        const TYPE(igraph_vector) *rhs);
DECLDIR igraph_bool_t FUNCTION(igraph_vector, all_l)(const TYPE(igraph_vector) *lhs,
        const TYPE(igraph_vector) *rhs);
DECLDIR igraph_bool_t FUNCTION(igraph_vector, all_g)(const TYPE(igraph_vector) *lhs,
        const TYPE(igraph_vector) *rhs);
DECLDIR igraph_bool_t FUNCTION(igraph_vector, all_le)(const TYPE(igraph_vector) *lhs,
        const TYPE(igraph_vector) *rhs);
DECLDIR igraph_bool_t FUNCTION(igraph_vector, all_ge)(const TYPE(igraph_vector) *lhs,
        const TYPE(igraph_vector) *rhs);

/*------------------------------*/
/* Finding minimum and maximum  */
/*------------------------------*/

DECLDIR BASE FUNCTION(igraph_vector, min)(const TYPE(igraph_vector)* v);
DECLDIR BASE FUNCTION(igraph_vector, max)(const TYPE(igraph_vector)* v);
DECLDIR igraph_error_t FUNCTION(igraph_vector, which_min)(const TYPE(igraph_vector)* v);
DECLDIR igraph_error_t FUNCTION(igraph_vector, which_max)(const TYPE(igraph_vector)* v);
DECLDIR igraph_error_t FUNCTION(igraph_vector, minmax)(const TYPE(igraph_vector) *v,
        BASE *min, BASE *max);
DECLDIR igraph_error_t FUNCTION(igraph_vector, which_minmax)(const TYPE(igraph_vector) *v,
        igraph_long_t *which_min, igraph_long_t *which_max);

/*-------------------*/
/* Vector properties */
/*-------------------*/

DECLDIR igraph_bool_t FUNCTION(igraph_vector, empty)     (const TYPE(igraph_vector)* v);
DECLDIR igraph_long_t FUNCTION(igraph_vector, size)      (const TYPE(igraph_vector)* v);
DECLDIR igraph_bool_t FUNCTION(igraph_vector, isnull)(const TYPE(igraph_vector) *v);
DECLDIR BASE FUNCTION(igraph_vector, sum)(const TYPE(igraph_vector) *v);
DECLDIR igraph_real_t FUNCTION(igraph_vector, sumsq)(const TYPE(igraph_vector) *v);
DECLDIR BASE FUNCTION(igraph_vector, prod)(const TYPE(igraph_vector) *v);
DECLDIR igraph_bool_t FUNCTION(igraph_vector, isininterval)(const TYPE(igraph_vector) *v,
        BASE low, BASE high);
DECLDIR igraph_bool_t FUNCTION(igraph_vector, any_smaller)(const TYPE(igraph_vector) *v,
        BASE limit);
DECLDIR igraph_bool_t FUNCTION(igraph_vector, is_equal)(const TYPE(igraph_vector) *lhs,
        const TYPE(igraph_vector) *rhs);
DECLDIR igraph_real_t FUNCTION(igraph_vector, maxdifference)(const TYPE(igraph_vector) *m1,
        const TYPE(igraph_vector) *m2);

/*------------------------*/
/* Searching for elements */
/*------------------------*/

DECLDIR igraph_bool_t FUNCTION(igraph_vector, contains)(const TYPE(igraph_vector) *v, BASE e);
DECLDIR igraph_bool_t FUNCTION(igraph_vector, search)(const TYPE(igraph_vector) *v,
        igraph_long_t from, BASE what,
        igraph_long_t *pos);
DECLDIR igraph_bool_t FUNCTION(igraph_vector, binsearch_slice)(const TYPE(igraph_vector) *v,
        BASE what, igraph_long_t *pos,
        igraph_long_t start, igraph_long_t end);
DECLDIR igraph_bool_t FUNCTION(igraph_vector, binsearch)(const TYPE(igraph_vector) *v,
        BASE what, igraph_long_t *pos);
DECLDIR igraph_bool_t FUNCTION(igraph_vector, binsearch2)(const TYPE(igraph_vector) *v,
        BASE what);

/*------------------------*/
/* Resizing operations    */
/*------------------------*/

DECLDIR void FUNCTION(igraph_vector, clear)(TYPE(igraph_vector)* v);
DECLDIR igraph_error_t FUNCTION(igraph_vector, resize)(TYPE(igraph_vector)* v, igraph_long_t newsize);
DECLDIR igraph_error_t FUNCTION(igraph_vector, resize_min)(TYPE(igraph_vector)*v);
DECLDIR igraph_error_t FUNCTION(igraph_vector, reserve)(TYPE(igraph_vector)* v, igraph_long_t size);
DECLDIR igraph_error_t FUNCTION(igraph_vector, push_back)(TYPE(igraph_vector)* v, BASE e);
DECLDIR BASE FUNCTION(igraph_vector, pop_back)(TYPE(igraph_vector)* v);
DECLDIR igraph_error_t FUNCTION(igraph_vector, insert)(TYPE(igraph_vector) *v, igraph_long_t pos, BASE value);
DECLDIR void FUNCTION(igraph_vector, remove)(TYPE(igraph_vector) *v, igraph_long_t elem);
DECLDIR void FUNCTION(igraph_vector, remove_section)(TYPE(igraph_vector) *v,
        igraph_long_t from, igraph_long_t to);

/*-----------*/
/* Sorting   */
/*-----------*/

DECLDIR void FUNCTION(igraph_vector, sort)(TYPE(igraph_vector) *v);
DECLDIR void FUNCTION(igraph_vector, reverse_sort)(TYPE(igraph_vector) *v);
DECLDIR igraph_error_t FUNCTION(igraph_vector, qsort_ind)(TYPE(igraph_vector) *v,
        igraph_vector_t *inds, igraph_bool_t descending);

/*-----------*/
/* Printing  */
/*-----------*/

DECLDIR igraph_error_t FUNCTION(igraph_vector, print)(const TYPE(igraph_vector) *v);
DECLDIR igraph_error_t FUNCTION(igraph_vector, printf)(const TYPE(igraph_vector) *v,
                                    const char *format);
DECLDIR igraph_error_t FUNCTION(igraph_vector, fprint)(const TYPE(igraph_vector) *v, FILE *file);

#ifdef BASE_COMPLEX

DECLDIR igraph_error_t igraph_vector_complex_real(const igraph_vector_complex_t *v,
                                       igraph_vector_t *real);
DECLDIR igraph_error_t igraph_vector_complex_imag(const igraph_vector_complex_t *v,
                                       igraph_vector_t *imag);
DECLDIR igraph_error_t igraph_vector_complex_realimag(const igraph_vector_complex_t *v,
        igraph_vector_t *real,
        igraph_vector_t *imag);
DECLDIR igraph_error_t igraph_vector_complex_create(igraph_vector_complex_t *v,
        const igraph_vector_t *real,
        const igraph_vector_t *imag);
DECLDIR igraph_error_t igraph_vector_complex_create_polar(igraph_vector_complex_t *v,
        const igraph_vector_t *r,
        const igraph_vector_t *theta);

#endif

DECLDIR igraph_error_t FUNCTION(igraph_vector, init_real)(TYPE(igraph_vector)*v, igraph_long_t no, ...);
DECLDIR igraph_error_t FUNCTION(igraph_vector, init_int)(TYPE(igraph_vector)*v, igraph_long_t no, ...);
DECLDIR igraph_error_t FUNCTION(igraph_vector, init_real_end)(TYPE(igraph_vector)*v, double endmark, ...);
DECLDIR igraph_error_t FUNCTION(igraph_vector, init_int_end)(TYPE(igraph_vector)*v, igraph_long_t endmark, ...);

DECLDIR igraph_error_t FUNCTION(igraph_vector, move_interval)(TYPE(igraph_vector) *v,
        igraph_long_t begin, igraph_long_t end, igraph_long_t to);
DECLDIR igraph_error_t FUNCTION(igraph_vector, move_interval2)(TYPE(igraph_vector) *v,
        igraph_long_t begin, igraph_long_t end, igraph_long_t to);
DECLDIR void FUNCTION(igraph_vector, permdelete)(TYPE(igraph_vector) *v,
        const igraph_vector_t *index,
        igraph_long_t nremove);
DECLDIR igraph_error_t FUNCTION(igraph_vector, filter_smaller)(TYPE(igraph_vector) *v, BASE elem);
DECLDIR igraph_error_t FUNCTION(igraph_vector, get_interval)(const TYPE(igraph_vector) *v,
        TYPE(igraph_vector) *res,
        igraph_long_t from, igraph_long_t to);
DECLDIR igraph_error_t FUNCTION(igraph_vector, difference_sorted)(const TYPE(igraph_vector) *v1,
        const TYPE(igraph_vector) *v2, TYPE(igraph_vector) *result);
DECLDIR igraph_error_t FUNCTION(igraph_vector, intersect_sorted)(const TYPE(igraph_vector) *v1,
        const TYPE(igraph_vector) *v2, TYPE(igraph_vector) *result);

DECLDIR igraph_error_t FUNCTION(igraph_vector, index)(const TYPE(igraph_vector) *v,
                                   TYPE(igraph_vector) *newv,
                                   const igraph_vector_t *idx);

DECLDIR igraph_error_t FUNCTION(igraph_vector, index_int)(TYPE(igraph_vector) *v,
                                       const igraph_vector_long_t *idx);
