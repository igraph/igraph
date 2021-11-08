/* qmd.h (quotient minimum degree algorithm) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2001 Free Software Foundation, Inc.
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

#ifndef QMD_H
#define QMD_H

#define genqmd _glp_genqmd
void genqmd(int *neqns, int xadj[], int adjncy[], int perm[],
      int invp[], int deg[], int marker[], int rchset[], int nbrhd[],
      int qsize[], int qlink[], int *nofsub);
/* GENeral Quotient Minimum Degree algorithm */

#define qmdrch _glp_qmdrch
void qmdrch(int *root, int xadj[], int adjncy[], int deg[],
      int marker[], int *rchsze, int rchset[], int *nhdsze,
      int nbrhd[]);
/* Quotient MD ReaCHable set */

#define qmdqt _glp_qmdqt
void qmdqt(int *root, int xadj[], int adjncy[], int marker[],
      int *rchsze, int rchset[], int nbrhd[]);
/* Quotient MD Quotient graph Transformation */

#define qmdupd _glp_qmdupd
void qmdupd(int xadj[], int adjncy[], int *nlist, int list[],
      int deg[], int qsize[], int qlink[], int marker[], int rchset[],
      int nbrhd[]);
/* Quotient MD UPDate */

#define qmdmrg _glp_qmdmrg
void qmdmrg(int xadj[], int adjncy[], int deg[], int qsize[],
      int qlink[], int marker[], int *deg0, int *nhdsze, int nbrhd[],
      int rchset[], int ovrlp[]);
/* Quotient MD MeRGe */

#endif

/* eof */
