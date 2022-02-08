/* stdc.h (standard ANSI C headers) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2000-2017 Free Software Foundation, Inc.
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

#ifndef STDC_H
#define STDC_H

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef ENABLE_NON_SAFE /* 29/I-2017 */
/* disable using non-thread-safe functions directly */
#undef gmtime
#define gmtime ???
#undef strerror
#define strerror ???
#undef strtok
#define strtok ???
#endif

#if 1 /* 29/I-2017 */
/* provide replacements for these functions on a per-thread basis */
#define xgmtime _glp_xgmtime
struct tm *xgmtime(const time_t *);
#define xstrerr _glp_xstrerr
char *xstrerr(int);
#define xstrtok _glp_xstrtok
char *xstrtok(char *, const char *);
#endif

#if 1 /* 06/II-2018 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef __WOE__
#define CDECL
#else
#define CDECL __cdecl
#endif
#endif

#endif

/* eof */
