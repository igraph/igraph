/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <igraph.h>
#include <math.h>
#include <stdlib.h>

#define sqr(x) ((x)*(x))

int main() {
  float dist[8][8] = {
	  {0.00, 4.69, 6.79, 3.50, 3.11, 4.46, 5.57, 3.00},
	  {4.69, 0.00, 2.10, 2.27, 2.65, 2.36, 1.99, 1.74},
	  {6.79, 2.10, 0.00, 3.78, 4.53, 2.83, 2.44, 3.79},
	  {3.50, 2.27, 3.78, 0.00, 1.98, 4.35, 2.07, 0.53},
	  {3.11, 2.65, 4.53, 1.98, 0.00, 3.80, 3.31, 1.47},
	  {4.46, 2.36, 2.83, 4.35, 3.80, 0.00, 4.35, 3.82},
	  {5.57, 1.99, 2.44, 2.07, 3.31, 4.35, 0.00, 2.57},
	  {3.00, 1.74, 3.79, 0.53, 1.47, 3.82, 2.57, 0.00},
  };
  igraph_t g;
  igraph_matrix_t coords, dist_mat;
  igraph_real_t vc;
  igraph_arpack_options_t options;
  int i, j;
  srand(time(0));

  igraph_arpack_options_init(&options);

  igraph_tree(&g, 10, 2, IGRAPH_TREE_UNDIRECTED);
  igraph_matrix_init(&coords, 0, 0);
  igraph_layout_mds(&g, &coords, 0, 2, &options);
  if (MATRIX(coords, 0, 0) > 0) {
    for (i = 0; i < igraph_matrix_nrow(&coords); i++)
      MATRIX(coords, i, 0) *= -1;
  }
  if (MATRIX(coords, 0, 1) < 0) {
    for (i = 0; i < igraph_matrix_nrow(&coords); i++)
      MATRIX(coords, i, 1) *= -1;
  }
  igraph_matrix_print(&coords);
  igraph_matrix_destroy(&coords);
  igraph_destroy(&g);

  igraph_full(&g, 8, IGRAPH_UNDIRECTED, 0);
  igraph_matrix_init(&coords, 8, 2);
  igraph_matrix_init(&dist_mat, 8, 8);
  for (i = 0; i < 8; i++)
    for (j = 0; j < 2; j++)
      MATRIX(coords, i, j) = rand() % 1000;
  for (i = 0; i < 8; i++)
    for (j = i+1; j < 8; j++) {
      double dist_sq = 0.0;
      dist_sq += sqr(MATRIX(coords, i, 0)-MATRIX(coords, j, 0));
      dist_sq += sqr(MATRIX(coords, i, 1)-MATRIX(coords, j, 1));
      MATRIX(dist_mat, i, j) = sqrt(dist_sq);
      MATRIX(dist_mat, j, i) = sqrt(dist_sq);
	}
  igraph_layout_mds(&g, &coords, &dist_mat, 2, &options);
  for (i = 0; i < 8; i++)
    for (j = i+1; j < 8; j++) {
      double dist_sq = 0.0;
      dist_sq += sqr(MATRIX(coords, i, 0)-MATRIX(coords, j, 0));
      dist_sq += sqr(MATRIX(coords, i, 1)-MATRIX(coords, j, 1));
      if (fabs(sqrt(dist_sq) - MATRIX(dist_mat, i, j)) > 1e-2) {
        printf("dist(%d,%d) should be %.4f, but it is %.4f\n",
				i, j, MATRIX(dist_mat, i, j), sqrt(dist_sq));
        return 1;
      }
    }
  igraph_matrix_destroy(&dist_mat);
  igraph_matrix_destroy(&coords);
  igraph_destroy(&g);

  return 0;
}
