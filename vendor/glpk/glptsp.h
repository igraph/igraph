/* glptsp.h (TSP format) */

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

#ifndef GLPTSP_H
#define GLPTSP_H

typedef struct TSP TSP;

struct TSP
{     /* TSP (or related problem) instance in the format described in
         the report [G.Reinelt, TSPLIB 95] */
      /*--------------------------------------------------------------*/
      /* the specification part */
      char *name;
      /* identifies the data file */
      int type;
      /* specifies the type of data: */
#define TSP_UNDEF             0  /* undefined */
#define TSP_TSP               1  /* symmetric TSP */
#define TSP_ATSP              2  /* asymmetric TSP */
#define TSP_TOUR              3  /* collection of tours */
      char *comment;
      /* additional comments (usually the name of the contributor or
         creator of the problem instance is given here) */
      int dimension;
      /* for a TSP or ATSP, the dimension is the number of its nodes
         for a TOUR it is the dimension of the corresponding problem */
      int edge_weight_type;
      /* specifies how the edge weights (or distances) are given: */
#define TSP_UNDEF             0  /* undefined */
#define TSP_EXPLICIT          1  /* listed explicitly */
#define TSP_EUC_2D            2  /* Eucl. distances in 2-D */
#define TSP_CEIL_2D           3  /* Eucl. distances in 2-D rounded up */
#define TSP_GEO               4  /* geographical distances */
#define TSP_ATT               5  /* special distance function */
      int edge_weight_format;
      /* describes the format of the edge weights if they are given
         explicitly: */
#define TSP_UNDEF             0  /* undefined */
#define TSP_FUNCTION          1  /* given by a function */
#define TSP_FULL_MATRIX       2  /* given by a full matrix */
#define TSP_UPPER_ROW         3  /* upper triangulat matrix (row-wise
                                    without diagonal entries) */
#define TSP_LOWER_DIAG_ROW    4  /* lower triangular matrix (row-wise
                                    including diagonal entries) */
      int display_data_type;
      /* specifies how a graphical display of the nodes can be
         obtained: */
#define TSP_UNDEF             0  /* undefined */
#define TSP_COORD_DISPLAY     1  /* display is generated from the node
                                    coordinates */
#define TSP_TWOD_DISPLAY      2  /* explicit coordinates in 2-D are
                                    given */
      /*--------------------------------------------------------------*/
      /* data part */
      /* NODE_COORD_SECTION: */
      double *node_x_coord; /* double node_x_coord[1+dimension]; */
      double *node_y_coord; /* double node_y_coord[1+dimension]; */
      /* DISPLAY_DATA_SECTION: */
      double *dply_x_coord; /* double dply_x_coord[1+dimension]; */
      double *dply_y_coord; /* double dply_y_coord[1+dimension]; */
      /* TOUR_SECTION: */
      int *tour; /* int tour[1+dimension]; */
      /* EDGE_WEIGHT_SECTION: */
      int *edge_weight; /* int edge_weight[1+dimension*dimension]; */
};

#define tsp_read_data         _glp_tsp_read_data
#define tsp_free_data         _glp_tsp_free_data
#define tsp_distance          _glp_tsp_distance

TSP *tsp_read_data(char *fname);
/* read TSP instance data */

void tsp_free_data(TSP *tsp);
/* free TSP instance data */

int tsp_distance(TSP *tsp, int i, int j);
/* compute distance between two nodes */

#endif

/* eof */
