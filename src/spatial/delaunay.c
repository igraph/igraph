/*
    IGraph library.
    Copyright (C) 2025  The igraph development team <igraph@igraph.org>
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/



#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_matrix.h"

#include "qhull/libqhull_r/libqhull_r.h"

igraph_error_t igraph_delaunay_triangulation(igraph_t *graph, igraph_matrix_t *points_) {
  int curlong, totlong; /* used !qh_NOmem */
  int exitcode, numpoints, dim;
  //coordT *points;
  boolT ismalloc;
  qhT qh_qh;
  qhT *qh= &qh_qh;
  double points[] = {
    0, 1,
    2, 3,
    0, 0
  };

  QHULL_LIB_CHECK; /* Check for compatible library */

  qh_init_A(qh, stdin, stdout, stderr, 0, NULL);  /* sets qh->qhull_command */
  exitcode= setjmp(qh->errexit); /* simple statement for CRAY J916 */
  if (!exitcode) {
    qh->NOerrexit = False;
    qh_option(qh, "delaunay  Qbbound-last", NULL, NULL);
    qh->DELAUNAY= True;     /* 'd'   */
    qh->SCALElast= True;    /* 'Qbb' */
    qh->KEEPcoplanar= True; /* 'Qc', to keep coplanars in 'p' */
    qh_checkflags(qh, qh->qhull_command, "");
    qh_initflags(qh, qh->qhull_command);
    //points = qh_readpoints(qh, &numpoints, &dim, &ismalloc);
    qh_init_B(qh, points, numpoints, dim, ismalloc);
    qh_qhull(qh);
    qh_check_output(qh);
    qh_produce_output(qh);
    if (qh->VERIFYoutput && !qh->FORCEoutput && !qh->STOPpoint && !qh->STOPcone)
      qh_check_points(qh);
    exitcode= qh_ERRnone;
  }
  qh->NOerrexit= True;  /* no more setjmp */
  qh_freeqhull(qh, !qh_ALL);
  qh_memfreeshort(qh, &curlong, &totlong);
  if (curlong || totlong)
    qh_fprintf_stderr(7079, "qhull internal warning (main): did not free %d bytes of long memory(%d pieces)\n",
       totlong, curlong);
  return IGRAPH_SUCCESS;
}
