/* spxprob.h */

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

#ifndef SPXPROB_H
#define SPXPROB_H

#include "prob.h"
#include "spxlp.h"

#define spx_init_lp _glp_spx_init_lp
void spx_init_lp(SPXLP *lp, glp_prob *P, int excl);
/* initialize working LP object */

#define spx_alloc_lp _glp_spx_alloc_lp
void spx_alloc_lp(SPXLP *lp);
/* allocate working LP arrays */

#define spx_build_lp _glp_spx_build_lp
void spx_build_lp(SPXLP *lp, glp_prob *P, int excl, int shift,
      int map[/*1+P->m+P->n*/]);
/* convert original LP to working LP */

#define spx_build_basis _glp_spx_build_basis
void spx_build_basis(SPXLP *lp, glp_prob *P, const int map[]);
/* convert original LP basis to working LP basis */

#define spx_store_basis _glp_spx_store_basis
void spx_store_basis(SPXLP *lp, glp_prob *P, const int map[],
      int daeh[/*1+n*/]);
/* convert working LP basis to original LP basis */

#define spx_store_sol _glp_spx_store_sol
void spx_store_sol(SPXLP *lp, glp_prob *P, int shift,
      const int map[], const int daeh[], const double beta[],
      const double pi[], const double d[]);
/* convert working LP solution to original LP solution */

#define spx_free_lp _glp_spx_free_lp
void spx_free_lp(SPXLP *lp);
/* deallocate working LP arrays */

#endif

/* eof */
