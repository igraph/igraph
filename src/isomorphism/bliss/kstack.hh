#ifndef BLISS_KSTACK_HH
#define BLISS_KSTACK_HH

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

#include <new>
#include <cassert>

namespace bliss {

/**
 * \brief A simple implementation of a stack with fixed maximum capacity.
 */
template <class Type>
class KStack {
public:
  /**
   * Create a new stack with zero capacity.
   * The function init() should be called next.
   */
  KStack();

  /**
   * Create a new stack with the capacity to hold at most \a N elements.
   */
  KStack(int N);

  ~KStack();

  /**
   * Initialize the stack to have the capacity to hold at most \a N elements.
   */
  void init(int N);

  /**
   * Is the stack empty?
   */
  bool is_empty() const {return cursor == entries; }

  /**
   * Return (but don't remove) the top element of the stack.
   */
  Type top() const {assert(cursor > entries); return *cursor; }

  /**
   * Pop (remove) the top element of the stack.
   */
  Type pop()
  {
    assert(cursor > entries);
    return *cursor--;
  }

  /**
   * Push the element \a e in the stack.
   */
  void push(Type e)
  {
    assert(cursor < entries + kapacity);
    *(++cursor) = e;
  }

  /** Remove all the elements in the stack. */
  void clean() {cursor = entries; }

  /**
   * Get the number of elements in the stack.
   */
  unsigned int size() const {return cursor - entries; }

  /**
   * Return the i:th element in the stack, where \a i is in the range
   * 0,...,this.size()-1; the 0:th element is the bottom element
   * in the stack.
   */
  Type element_at(unsigned int i)
  {
    assert(i < size());
    return entries[i+1];
  }

  /** Return the capacity (NOT the number of elements) of the stack. */
  int capacity() const {return kapacity; }
private:
  int kapacity;
  Type *entries;
  Type *cursor;
};

template <class Type>
KStack<Type>::KStack()
{
  kapacity = 0;
  entries = nullptr;
  cursor = nullptr;
}

template <class Type>
KStack<Type>::KStack(int k)
{
  assert(k > 0);
  kapacity = k;
  entries = new Type[k+1];
  cursor = entries;
}

template <class Type>
void KStack<Type>::init(int k)
{
  assert(k > 0);
  delete[] entries;
  kapacity = k;
  entries = new Type[k+1];
  cursor = entries;
}

template <class Type>
KStack<Type>::~KStack()
{
  delete[] entries;
  kapacity = 0;
  entries = nullptr;
  cursor = nullptr;
}

} // namespace bliss

#endif // BLISS_KSTACK_HH
