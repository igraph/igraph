/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include "igraph.h"
#include "math.h"

/* internal function */

int igraph_2dgrid_which(igraph_2dgrid_t *grid, real_t xc, real_t yc,
			long int *x, long int *y) {

  if (xc <= grid->minx) { 
    *x=0; 
  } else if (xc >= grid->maxx) {
    *x=grid->stepsx-1;
  } else {
    *x=floor((xc-(grid->minx))/(grid->deltax));
  }

  if (yc <= grid->miny) {
    *y=0;
  } else if (yc >= grid->maxy) {
    *y=grid->stepsy-1;
  } else {
    *y=floor((yc-(grid->miny))/(grid->deltay));
  }

  if (*x < 0 || *y < 0 || *x >= grid->stepsx || *y >= grid->stepsy) {
    printf("HOHOHO x: %li y: %li xc: %f yc: %f minx: %f maxx: %f miny: %f maxy: %f deltax: %f deltay: %f\n", x, y, xc, yc, grid->minx, grid->maxx, grid->miny, grid->maxy, grid->deltax, grid->deltay); 
  }

  return 0;
}

int igraph_2dgrid_init(igraph_2dgrid_t *grid, igraph_matrix_t *coords, 
		       real_t minx, real_t maxx, real_t deltax,
		       real_t miny, real_t maxy, real_t deltay) {
  long int i;
  
  grid->coords=coords;
  grid->minx=minx;
  grid->maxx=maxx;
  grid->deltax=deltax;
  grid->miny=miny;
  grid->maxy=maxy;
  grid->deltay=deltay;
  
  grid->stepsx=ceil((maxx-minx)/deltax);
  grid->stepsy=ceil((maxy-miny)/deltay);

  IGRAPH_CHECK(igraph_matrix_init(&grid->startidx, 
				  grid->stepsx, grid->stepsy));
  IGRAPH_FINALLY(igraph_matrix_destroy, &grid->startidx);
  IGRAPH_VECTOR_INIT_FINALLY(&grid->next, igraph_matrix_nrow(coords));
  IGRAPH_VECTOR_INIT_FINALLY(&grid->prev, igraph_matrix_nrow(coords));

  for (i=0; i<igraph_vector_size(&grid->next); i++) {
    VECTOR(grid->next)[i]=-1;
  }

  grid->massx=0;
  grid->massy=0;
  grid->vertices=0;
  
  IGRAPH_FINALLY_CLEAN(3);
  return 0;
}

void igraph_2dgrid_destroy(igraph_2dgrid_t *grid) {
  igraph_matrix_destroy(&grid->startidx);
  igraph_vector_destroy(&grid->next);
  igraph_vector_destroy(&grid->prev);
}

void igraph_2dgrid_add(igraph_2dgrid_t *grid, long int elem,
		       real_t xc, real_t yc) {
  long int x, y;
  long int first;

  MATRIX(*grid->coords, elem, 0)=xc;
  MATRIX(*grid->coords, elem, 1)=yc;

  /* add to cell */
  igraph_2dgrid_which(grid, xc, yc, &x, &y);
  first=MATRIX(grid->startidx, x, y);
  VECTOR(grid->prev)[elem]=0;
  VECTOR(grid->next)[elem]=first;
  if (first != 0) {
    VECTOR(grid->prev)[first-1]=elem+1;
  }
  MATRIX(grid->startidx, x, y)=elem+1;

  grid->massx += xc;
  grid->massy += yc;
  grid->vertices += 1;
}

void igraph_2dgrid_move(igraph_2dgrid_t *grid, long int elem, 
			real_t xc, real_t yc) {
  long int oldx, oldy;
  long int newx, newy;
  real_t oldxc=MATRIX(*grid->coords, elem, 0);
  real_t oldyc=MATRIX(*grid->coords, elem, 1);
  long int first;

  xc=oldxc+xc; yc=oldyc+yc;

  igraph_2dgrid_which(grid, oldxc, oldyc, &oldx, &oldy);
  /* remove from this cell */
  if (VECTOR(grid->prev)[elem] != 0) {
    VECTOR(grid->next) [ (long int) VECTOR(grid->prev)[elem]-1 ] =
      VECTOR(grid->next)[elem];
  } else {
    MATRIX(grid->startidx, oldx, oldy)=VECTOR(grid->next)[elem];
  }
  if (VECTOR(grid->next)[elem] != 0) {
    VECTOR(grid->prev)[ (long int) VECTOR(grid->next)[elem]-1 ] =
      VECTOR(grid->prev)[elem];
  }

  igraph_2dgrid_which(grid, xc, yc, &newx, &newy);
  /* add to this cell */
  first=MATRIX(grid->startidx, newx, newy);
  VECTOR(grid->prev)[elem]=0;
  VECTOR(grid->next)[elem]=first;
  if (first != 0) {
    VECTOR(grid->prev)[first-1]=elem+1;
  }
  MATRIX(grid->startidx, newx, newy)=elem+1;
  
  grid->massx += -oldxc+xc;
  grid->massy += -oldyc+yc;

  MATRIX(*grid->coords, elem, 0)=xc;
  MATRIX(*grid->coords, elem, 1)=yc;
  
}

void igraph_2dgrid_getcenter(const igraph_2dgrid_t *grid, 
			     real_t *massx, real_t *massy) {
  *massx = (grid->massx)/(grid->vertices);
  *massy = (grid->massy)/(grid->vertices);
}

bool_t igraph_2dgrid_in(const igraph_2dgrid_t *grid, long int elem) {
  return VECTOR(grid->next)[elem] != -1;
}

real_t igraph_2dgrid_dist(const igraph_2dgrid_t *grid, 
			  long int e1, long int e2) {
  real_t x1=MATRIX(*grid->coords, e1, 0);
  real_t x2=MATRIX(*grid->coords, e2, 0);
  real_t y1=MATRIX(*grid->coords, e1, 1);
  real_t y2=MATRIX(*grid->coords, e2, 1);
  
  return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

int igraph_i_2dgrid_addvertices(igraph_2dgrid_t *grid, igraph_vector_t *eids,
				integer_t vid, real_t r, 
				long int x, long int y) {
  long int act;

  act=MATRIX(grid->startidx, x, y);
  while (act != 0) {
    if (igraph_2dgrid_dist(grid, vid, act-1) < r) {
      IGRAPH_CHECK(igraph_vector_push_back(eids, act-1));
    }
    act=VECTOR(grid->next)[act-1];
  }
  return 0;
}

int igraph_2dgrid_neighbors(igraph_2dgrid_t *grid, igraph_vector_t *eids, 
			    integer_t vid, real_t r) {
  real_t xc=MATRIX(*grid->coords, (long int)vid, 0);
  real_t yc=MATRIX(*grid->coords, (long int)vid, 1);
  long int x, y;
  
  igraph_vector_clear(eids);
  
  igraph_2dgrid_which(grid, xc, yc, &x, &y);

  /* this cell */
  igraph_i_2dgrid_addvertices(grid, eids, vid, r, x, y);
  /* left */
  if (x!=0) {
    igraph_i_2dgrid_addvertices(grid, eids, vid, r, x-1, y);
  }
  /* right */
  if (x!=grid->stepsx-1) {
    igraph_i_2dgrid_addvertices(grid, eids, vid, r, x+1, y);
  }
  /* up */
  if (y!=0) {
    igraph_i_2dgrid_addvertices(grid, eids, vid, r, x, y-1);
  }
  /* down */
  if (y!=grid->stepsy-1) {
    igraph_i_2dgrid_addvertices(grid, eids, vid, r, x, y+1);
  }
  /* up & left */
  if (x != 0 && y != 0) {
    igraph_i_2dgrid_addvertices(grid, eids, vid, r, x-1, y-1);
  }
  /* up & right */
  if (x != grid->stepsx-1 && y != 0) {
    igraph_i_2dgrid_addvertices(grid, eids, vid, r, x+1, y-1);
  }
  /* down & left */
  if (x != 0 && y != grid->stepsy-1) {
    igraph_i_2dgrid_addvertices(grid, eids, vid, r, x-1, y+1);
  }
  /* down & right */
  if (x != grid->stepsx-1 && y != grid->stepsy-1) {
    igraph_i_2dgrid_addvertices(grid, eids, vid, r, x-1, y+1);
  }
  
  return 0;
}
