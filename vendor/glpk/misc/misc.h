/* misc.h (miscellaneous routines) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2000-2013 Free Software Foundation, Inc.
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

#ifndef MISC_H
#define MISC_H

#define str2int _glp_str2int
int str2int(const char *str, int *val);
/* convert character string to value of int type */

#define str2num _glp_str2num
int str2num(const char *str, double *val);
/* convert character string to value of double type */

#define strspx _glp_strspx
char *strspx(char *str);
/* remove all spaces from character string */

#define strtrim _glp_strtrim
char *strtrim(char *str);
/* remove trailing spaces from character string */

#define gcd _glp_gcd
int gcd(int x, int y);
/* find greatest common divisor of two integers */

#define gcdn _glp_gcdn
int gcdn(int n, int x[]);
/* find greatest common divisor of n integers */

#define round2n _glp_round2n
double round2n(double x);
/* round floating-point number to nearest power of two */

#define fp2rat _glp_fp2rat
int fp2rat(double x, double eps, double *p, double *q);
/* convert floating-point number to rational number */

#endif

/* eof */
