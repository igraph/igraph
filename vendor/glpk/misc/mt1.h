/* mt1.h (0-1 knapsack problem; Martello & Toth algorithm) */

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

#ifndef MT1_H
#define MT1_H

#define mt1 _glp_mt1
int mt1(int n, int p[], int w[], int c, int x[], int jck, int xx[],
      int min[], int psign[], int wsign[], int zsign[]);
/* solve 0-1 single knapsack problem */

#endif

/* eof */
