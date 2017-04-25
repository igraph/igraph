/* -*- mode: C -*-  */
/* 
   IGraph library R interface.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA
   
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA 

*/

#include "igraph.h"
#include "igraph_error.h"

#include "RayTracer.h"
#include "Sphere.h"

#include "config.h"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

using namespace igraph;

extern "C" {

SEXP R_igraph_getsphere(SEXP pos, SEXP radius, SEXP color, SEXP bgcolor,
			SEXP lightpos, SEXP lightcolor, SEXP width,
			SEXP height) {

  /* All error checking is done at the R level */

  int i;
  double *spos=REAL(pos);
  double *scolor=REAL(color);
  double *svgcolor=REAL(bgcolor);
  int no_lights=GET_LENGTH(lightpos); 
  RayTracer* p_ray_tracer;
  Sphere * sphere;
  int swidth=INTEGER(width)[0];
  int sheight=INTEGER(height)[0];
  int nopixels=swidth * sheight;
  SEXP result, dim;
  Image image;
  
  p_ray_tracer = new RayTracer();
  p_ray_tracer->EyePoint(Point(0,0,0));  
  
  for (i=0; i<no_lights; i++) {
    double *lpos=REAL(VECTOR_ELT(lightpos, i));
    double *lcol=REAL(VECTOR_ELT(lightcolor, i));
    Light *light = new Light(Point(lpos[0], lpos[1], lpos[2]));
    light->Intensity(1);
    light->LightColor(Color(lcol[0], lcol[1], lcol[2]));
    p_ray_tracer->AddLight(light);
  }
  
  sphere = new Sphere(Point(spos[0], spos[1], spos[2]), REAL(radius)[0]);
  sphere->ShapeColor(Color(scolor[0], scolor[1], scolor[2]));
  p_ray_tracer->AddShape(sphere);

  PROTECT(result=NEW_NUMERIC(nopixels * 4));
  PROTECT(dim=NEW_INTEGER(3));
  INTEGER(dim)[0]=swidth; INTEGER(dim)[1]=sheight; INTEGER(dim)[2]=4;
  SET_DIM(result, dim);
  
  image.width=swidth;
  image.height=sheight;
  image.red=REAL(result);
  image.green=image.red + nopixels;
  image.blue=image.green + nopixels;
  image.trans=image.blue + nopixels;

  p_ray_tracer->RayTrace(image);
  delete p_ray_tracer;
  
  UNPROTECT(2);
  return result;
}

} // extern C
