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

#include <cstdlib>
#include <cstdio>
#include <climits>
#include "bliss_defs.hh"
#include "bliss_heap.hh"

using namespace std;

namespace igraph {

Heap::~Heap()
{
  if(array)
    {
      free(array);
      array = 0;
    }
}

void Heap::upheap(unsigned int index)
{
  assert(n >= 1);
  assert(index >= 1 && index <= n);
  const unsigned int v = array[index];
  array[0] = UINT_MAX;
  while(array[index/2] <= v)
    {
      array[index] = array[index/2];
      index = index/2;
    }
  array[index] = v;
}

void Heap::downheap(unsigned int index)
{
  const unsigned int v = array[index];
  while(index <= n/2)
    {
      unsigned int new_index = index + index;
      if(new_index < n && array[new_index] < array[new_index+1])
	new_index++;
      if(v >= array[new_index])
	break;
      array[index] = array[new_index];
      index = new_index;
    }
  array[index] = v;
}

void Heap::init(unsigned int size)
{
  array = (unsigned int*)malloc((size + 1) * sizeof(unsigned int));
  n = 0;
#if defined(CONSISTENCY_CHECKS)
  assert(size > 0);
  N = size;
#endif
}

void Heap::insert(unsigned int v)
{
  DEBUG_ASSERT(n < N);
  array[++n] = v;
  upheap(n);
}

unsigned int Heap::remove()
{
  DEBUG_ASSERT(n >= 1 && n <= N);
  const unsigned int v = array[1];
  array[1] = array[n--];
  downheap(1);
  return v;
}

}
