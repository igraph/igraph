#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "defs.hh"
#include "heap.hh"

/* use 'and' instead of '&&' */
#if _MSC_VER
#include <ciso646>
#endif

/*
  Copyright (c) 2003-2015 Tommi Junttila
  Released under the GNU Lesser General Public License version 3.
  
  This file is part of bliss.
  
  bliss is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, version 3 of the License.

  bliss is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with bliss.  If not, see <http://www.gnu.org/licenses/>.
*/

namespace bliss {

Heap::~Heap()
{
  if(array)
    {
      free(array);
      array = 0;
      n = 0;
      N = 0;
    }
}

void Heap::upheap(unsigned int index)
{
  const unsigned int v = array[index];
  array[0] = 0;
  while(array[index/2] > v)
    {
      array[index] = array[index/2];
      index = index/2;
    }
  array[index] = v;
}

void Heap::downheap(unsigned int index)
{
  const unsigned int v = array[index];
  const unsigned int lim = n/2;
  while(index <= lim)
    {
      unsigned int new_index = index + index;
      if((new_index < n) and (array[new_index] > array[new_index+1]))
	new_index++;
      if(v <= array[new_index])
	break;
      array[index] = array[new_index];
      index = new_index;
    }
  array[index] = v;
}

void Heap::init(const unsigned int size)
{
  if(size > N)
    {
      if(array)
	free(array);
      array = (unsigned int*)malloc((size + 1) * sizeof(unsigned int));
      N = size;
    }
  n = 0;
}

void Heap::insert(const unsigned int v)
{
  array[++n] = v;
  upheap(n);
}

unsigned int Heap::remove()
{
  const unsigned int v = array[1];
  array[1] = array[n--];
  downheap(1);
  return v;
}

} // namespace bliss
