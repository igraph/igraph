/* glplib.h (miscellaneous library routines) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
*  2009, 2010 Andrew Makhorin, Department for Applied Informatics,
*  Moscow Aviation Institute, Moscow, Russia. All rights reserved.
*  E-mail: <mao@gnu.org>.
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

#ifndef GLPLIB_H
#define GLPLIB_H

#define bigmul _glp_lib_bigmul
void bigmul(int n, int m, unsigned short x[], unsigned short y[]);
/* multiply unsigned integer numbers of arbitrary precision */

#define bigdiv _glp_lib_bigdiv
void bigdiv(int n, int m, unsigned short x[], unsigned short y[]);
/* divide unsigned integer numbers of arbitrary precision */

#ifndef GLP_LONG_DEFINED
#define GLP_LONG_DEFINED
typedef struct { int lo, hi; } glp_long;
/* long integer data type */
#endif

typedef struct { glp_long quot, rem; } glp_ldiv;
/* result of long integer division */

#define xlset _glp_lib_xlset
glp_long xlset(int x);
/* expand integer to long integer */

#define xlneg _glp_lib_xlneg
glp_long xlneg(glp_long x);
/* negate long integer */

#define xladd _glp_lib_xladd
glp_long xladd(glp_long x, glp_long y);
/* add long integers */

#define xlsub _glp_lib_xlsub
glp_long xlsub(glp_long x, glp_long y);
/* subtract long integers */

#define xlcmp _glp_lib_xlcmp
int xlcmp(glp_long x, glp_long y);
/* compare long integers */

#define xlmul _glp_lib_xlmul
glp_long xlmul(glp_long x, glp_long y);
/* multiply long integers */

#define xldiv _glp_lib_xldiv
glp_ldiv xldiv(glp_long x, glp_long y);
/* divide long integers */

#define xltod _glp_lib_xltod
double xltod(glp_long x);
/* convert long integer to double */

#define xltoa _glp_lib_xltoa
char *xltoa(glp_long x, char *s);
/* convert long integer to character string */

#define str2int _glp_lib_str2int
int str2int(const char *str, int *val);
/* convert character string to value of int type */

#define str2num _glp_lib_str2num
int str2num(const char *str, double *val);
/* convert character string to value of double type */

#define strspx _glp_lib_strspx
char *strspx(char *str);
/* remove all spaces from character string */

#define strtrim _glp_lib_strtrim
char *strtrim(char *str);
/* remove trailing spaces from character string */

#define strrev _glp_lib_strrev
char *strrev(char *s);
/* reverse character string */

#define gcd _glp_lib_gcd
int gcd(int x, int y);
/* find greatest common divisor of two integers */

#define gcdn _glp_lib_gcdn
int gcdn(int n, int x[]);
/* find greatest common divisor of n integers */

#define lcm _glp_lib_lcm
int lcm(int x, int y);
/* find least common multiple of two integers */

#define lcmn _glp_lib_lcmn
int lcmn(int n, int x[]);
/* find least common multiple of n integers */

#define round2n _glp_lib_round2n
double round2n(double x);
/* round floating-point number to nearest power of two */

#define fp2rat _glp_lib_fp2rat
int fp2rat(double x, double eps, double *p, double *q);
/* convert floating-point number to rational number */

#define jday _glp_lib_jday
int jday(int d, int m, int y);
/* convert calendar date to Julian day number */

#define jdate _glp_lib_jdate
int jdate(int j, int *d, int *m, int *y);
/* convert Julian day number to calendar date */

#endif

/* eof */
