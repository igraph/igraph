/*
Copyright (C) 2003-2006 Tommi Junttila

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License version 2
 as published by the Free Software Foundation.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* FSF address fixed in the above notice on 1 Oct 2009 by Tamas Nepusz */

#include <assert.h>
#include "bliss_utils.hh"

namespace igraph {

void print_permutation(FILE *fp,
		       const unsigned int N,
		       const unsigned int *perm)
{
  assert(N > 0);
  assert(perm);
  for(unsigned int i = 0; i < N; i++) {
    unsigned int j = perm[i];
    if(j == i)
      continue;
    bool is_first = true;
    while(j != i) {
      if(j < i) {
        is_first = false;
        break;
      }
      j = perm[j];
    }
    if(!is_first)
      continue;
    fprintf(fp, "(%u,", i);
    j = perm[i];
    while(j != i) {
      fprintf(fp, "%u", j);
      j = perm[j];
      if(j != i)
        fprintf(fp, ",");
    }
    fprintf(fp, ")");
  }
}

}
