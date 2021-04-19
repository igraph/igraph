#include "heap.hh"

#include <new>
#include <cassert>

/* Allow using 'and' instead of '&&' with MSVC */
#if _MSC_VER
#include <ciso646>
#endif

/*
  Copyright (c) 2003-2021 Tommi Junttila
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

Heap::Heap() {
  array = nullptr;
  n = 0;
  N = 0;
}

Heap::~Heap()
{
  delete[] array;
  array = nullptr;
  n = 0;
  N = 0;
}

void Heap::upheap(unsigned int index)
{
  assert(n >= 1);
  assert(index >= 1 and index <= n);
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
  assert(size > 0);
  if(size > N)
    {
      delete[] array;
      array = new unsigned int[size + 1];
      N = size;
    }
  n = 0;
}

void Heap::insert(const unsigned int v)
{
  assert(n < N);
  array[++n] = v;
  upheap(n);
}

unsigned int Heap::smallest() const
{
  assert(n >= 1 and n <= N);
  return array[1];
}

unsigned int Heap::remove()
{
  assert(n >= 1 and n <= N);
  const unsigned int v = array[1];
  array[1] = array[n--];
  downheap(1);
  return v;
}

} // namespace bliss
