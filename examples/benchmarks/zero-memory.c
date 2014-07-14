/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
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

#include <string.h>
#include <stdlib.h>
#include "bench.h"

#define SIZE (10000000)

int main() {

  int *buffer = malloc(SIZE * sizeof(int));
  int i;
  int sum = 0;

  if (!buffer) { return 1; }
  
  BENCH("1. Zero newly allocated", 100,
	memset(buffer, 0, SIZE);
	);

  BENCH("2. Sum                 ", 100,
	for (i = 0; i < SIZE; i++) { sum += buffer[i]; }
	);
  
  /* Need to include this, otherwise the compiler optimizes
     away everything */

  if (sum != 0) { exit(1); }

  free(buffer);

  return 0;
}
