/* This file was imported from a private scientific library
 * based on GSL coined Home Scientific Libray (HSL) by its author
 * Jerome Benoit; this very material is itself inspired from the
 * material written by G. Jungan and distributed by GSL.
 * Ultimately, some modifications were done in order to render the
 * imported material independent from the rest of GSL.
 */

/* `hsl/hsl_sf_zeta.h' C header file
// HSL - Home Scientific Library
// Copyright (C) 2005-2018  Jerome Benoit
//
// HSL is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

/* For futher details, see its source conterpart src/hzeta.c */

/* Author:  Jerome G. Benoit < jgmbenoit _at_ rezozer _dot_ net > */

#ifndef __HZETA_H__
#define __HZETA_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


/* Hurwitz Zeta Function
 * zeta(s,q) = Sum[ (k+q)^(-s), {k,0,Infinity} ]
 *
 * s > 1.0, q > 0.0
 */
double hsl_sf_hzeta(const double s, const double q);

/* First Derivative of Hurwitz Zeta Function
 * zeta'(s,q) = - Sum[ Ln(k+q)/(k+q)^(s), {k,0,Infinity} ]
 *
 * s > 1.0, q > 0.0
 */
double hsl_sf_hzeta_deriv(const double s, const double q);

/* Second Derivative of Hurwitz Zeta Function
 * zeta''(s,q) = + Sum[ Ln(k+q)^2/(k+q)^(s), {k,0,Infinity} ]
 *
 * s > 1.0, q > 0.0
 */
double hsl_sf_hzeta_deriv2(const double s, const double q);

/* Logarithm of Hurwitz Zeta Function
 * lnzeta(s,q) = ln(zeta(s,q))
 *
 * s > 1.0, q > 0.0 (and q >> 1)
 */
double hsl_sf_lnhzeta(const double s, const double q);

/* Logarithmic Derivative of Hurwitz Zeta Function
 * lnzeta'(s,q) = zeta'(s,q)/zeta(s,q)
 *
 * s > 1.0, q > 0.0 (and q >> 1)
 */
double hsl_sf_lnhzeta_deriv(const double s, const double q);

/* Logarithm and Logarithmic Derivative of Hurwitz Zeta Function:
 * nonredundant computation version:
 * - lnzeta(s,q) and lnzeta'(s,q) are stored in *deriv0 and *deriv1, respectively;
 * - the return value and the value stored in *deriv0 are the same;
 * - deriv0 and deriv1 must be effective pointers, that is, not the NULL pointer.
 *
 * s > 1.0, q > 0.0 (and q >> 1)
 */
double hsl_sf_lnhzeta_deriv_tuple(const double s, const double q, double * deriv0, double * deriv1);


__END_DECLS

#endif // __HZETA_H__
