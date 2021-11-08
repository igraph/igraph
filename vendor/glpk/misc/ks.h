/* ks.h (0-1 knapsack problem) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2017-2018 Free Software Foundation, Inc.
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

#ifndef KS_H
#define KS_H

#define ks_enum _glp_ks_enum
int ks_enum(int n, const int a[/*1+n*/], int b, const int c[/*1+n*/],
      char x[/*1+n*/]);
/* solve 0-1 knapsack problem by complete enumeration */

#define ks_mt1 _glp_ks_mt1
int ks_mt1(int n, const int a[/*1+n*/], int b, const int c[/*1+n*/],
      char x[/*1+n*/]);
/* solve 0-1 knapsack problem with Martello & Toth algorithm */

#define ks_greedy _glp_ks_greedy
int ks_greedy(int n, const int a[/*1+n*/], int b, const int c[/*1+n*/],
      char x[/*1+n*/]);
/* solve 0-1 knapsack problem with greedy heuristic */

#endif

/* eof */
