/* simplex.h */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2015 Free Software Foundation, Inc.
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

#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "prob.h"

#define spx_primal _glp_spx_primal
int spx_primal(glp_prob *P, const glp_smcp *parm);
/* driver to the primal simplex method */

#define spy_dual _glp_spy_dual
int spy_dual(glp_prob *P, const glp_smcp *parm);
/* driver to the dual simplex method */

#endif

/* eof */
