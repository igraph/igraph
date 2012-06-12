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

#ifndef BLISS_KSTACK_H
#define BLISS_KSTACK_H

#include "bliss_defs.hh"
#include <cstdlib>		// malloc

namespace igraph {

/* 
 * A stack with fixed capacity
 */
template <class Type>
class KStack {
public:
  KStack();
  KStack(int k);
  ~KStack();
  void init(int k);

  bool is_empty() const {return(cursor == entries); }
  Type top() const {DEBUG_ASSERT(cursor > entries); return *cursor; }
  Type pop() {
    DEBUG_ASSERT(cursor > entries);
    Type obj = *cursor;
    cursor--;
    return obj;
  }
  void push(Type obj) {
    DEBUG_ASSERT(cursor < entries + kapacity);
    cursor++;
    *cursor = obj;
  }
  void clean() {cursor = entries; }
  unsigned int size() const {return(cursor - entries);
  }
  Type element_at(unsigned int i) {
    assert(i < size());
    return entries[i+1];
  }
  int capacity() {return kapacity; }
private:
  int kapacity;
  Type *entries;
  Type *cursor;
};

template <class Type>
KStack<Type>::KStack() {
  kapacity = 0;
  entries = 0;
  cursor = 0;
}

template <class Type>
KStack<Type>::KStack(int k) {
  assert(k > 0);
  kapacity = k;
  entries = (Type*)malloc((k+1) * sizeof(Type));
  cursor = entries;
}

template <class Type>
void KStack<Type>::init(int k) {
  assert(k > 0);
  if(entries)
    free(entries);
  kapacity = k;
  entries = (Type*)malloc((k+1) * sizeof(Type));
  cursor = entries;
}

template <class Type>
KStack<Type>::~KStack() {
  free(entries);
}

}

#endif
