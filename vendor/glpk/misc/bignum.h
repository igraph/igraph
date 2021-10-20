/* bignum.h (bignum arithmetic) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2006-2013 Free Software Foundation, Inc.
*  Written by Andrew Makhorin <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#ifndef BIGNUM_H
#define BIGNUM_H

#define bigmul _glp_bigmul
void bigmul(int n, int m, unsigned short x[], unsigned short y[]);
/* multiply unsigned integer numbers of arbitrary precision */

#define bigdiv _glp_bigdiv
void bigdiv(int n, int m, unsigned short x[], unsigned short y[]);
/* divide unsigned integer numbers of arbitrary precision */

#endif

/* eof */
