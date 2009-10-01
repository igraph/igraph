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

#ifndef BLISS_KQUEUE_HH
#define BLISS_KQUEUE_HH

#include "bliss_defs.hh"
#include <cstdlib>		// malloc

namespace igraph {

/* 
 * A queue with fixed capacity
 */
template <class Type>
class KQueue
{
public:
  KQueue();
  ~KQueue();
  void init(const unsigned int k);
  
  bool is_empty() const;
  unsigned int size() const;
  unsigned int capacity() const;
  void clear();

  Type front() const;
  Type pop_front();
  void push_front(Type e);
  Type pop_back();
  void push_back(Type e);
private:
  Type *entries, *end;
  Type *head, *tail;
};

template <class Type>
KQueue<Type>::KQueue()
{
  entries = 0;
  end = 0;
  head = 0;
  tail = 0;
}

template <class Type>
KQueue<Type>::~KQueue()
{
  if(entries)
    free(entries);
}

template <class Type>
void KQueue<Type>::init(const unsigned int k)
{
  assert(k > 0);
  if(entries)
    free(entries);
  entries = (Type*)malloc((k + 1) * sizeof(Type));
  end = entries + k + 1;
  head = entries;
  tail = head;
}

template <class Type>
void KQueue<Type>::clear()
{
  head = entries;
  tail = head;
}

template <class Type>
bool KQueue<Type>::is_empty() const
{
  return(head == tail);
}

template <class Type>
unsigned int KQueue<Type>::size() const
{
  if(tail >= head)
    return(tail - head);
  return((end - head) + (tail - entries));
}

template <class Type>
Type KQueue<Type>::front() const
{
  DEBUG_ASSERT(head != tail);
  return *head;
}

template <class Type>
Type KQueue<Type>::pop_front()
{
  DEBUG_ASSERT(head != tail);
  Type *old_head = head;
  head++;
  if(head == end)
    head = entries;
  return *old_head;
}

template <class Type>
void KQueue<Type>::push_front(Type e)
{
  if(head == entries)
    head = end - 1;
  else
    head--;
  DEBUG_ASSERT(head != tail);
  *head = e;
}

template <class Type>
void KQueue<Type>::push_back(Type e)
{
  *tail = e;
  tail++;
  if(tail == end)
    tail = entries;
  DEBUG_ASSERT(head != tail);
}

}

#endif
