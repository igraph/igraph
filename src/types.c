/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_types.h"

int igraph_real_printf(igraph_real_t val) {
  if (igraph_finite(val)) {
    return printf("%g", val);
  } else if (igraph_is_nan(val)) {
    return printf("NaN");
  } else if (igraph_is_inf(val)) {
    if (val < 0) {
      return printf("-Inf");
    } else {
      return printf("Inf");
    }
  } else {
    /* fallback */
    return printf("%g", val);
  }
}

int igraph_real_fprintf(FILE *file, igraph_real_t val) {
  if (igraph_finite(val)) {
    return fprintf(file, "%g", val);
  } else if (igraph_is_nan(val)) {
    return fprintf(file, "NaN");
  } else if (igraph_is_inf(val)) {
    if (val < 0) {
      return fprintf(file, "-Inf");
    } else {
      return fprintf(file, "Inf");
    }
  } else {
    /* fallback */
    return fprintf(file, "%g", val);
  }
}
