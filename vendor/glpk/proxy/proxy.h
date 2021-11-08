/* proxy.h (proximity search heuristic algorithm) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2013 Free Software Foundation, Inc.
*  Written by Giorgio Sartor <0gioker0@gmail.com>.
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

#ifndef PROXY_H
#define PROXY_H

#define proxy _glp_proxy
int proxy(glp_prob *lp, double *zstar, double *xstar,
          const double initsol[], double rel_impr, int tlim,
          int verbose);

#endif

/* eof */
