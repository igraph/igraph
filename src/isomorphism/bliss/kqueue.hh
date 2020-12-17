#ifndef BLISS_KQUEUE_HH
#define BLISS_KQUEUE_HH

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

#include "defs.hh"

namespace bliss {

/** \internal
 * \brief A very simple implementation of queues with fixed capacity.
 */

template <class Type>
class KQueue
{
public:
  /**
   * Create a new queue with capacity zero.
   * The function init() should be called next.
   */
  KQueue();

  ~KQueue();

  /**
   * Initialize the queue to have the capacity to hold at most \a N elements.
   */
  void init(const unsigned int N);
  
  /** Is the queue empty? */
  bool is_empty() const;

  /** Return the number of elements in the queue. */
  unsigned int size() const;

  /** Remove all the elements in the queue. */
  void clear();

  /** Return (but don't remove) the first element in the queue. */
  Type front() const;

  /** Remove and return the first element of the queue. */
  Type pop_front();

  /** Push the element \a e in the front of the queue. */
  void push_front(Type e);

  /** Remove and return the last element of the queue. */
  Type pop_back();

  /** Push the element \a e in the back of the queue. */
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
  return *head;
}

template <class Type>
Type KQueue<Type>::pop_front()
{
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
  *head = e;
}

template <class Type>
void KQueue<Type>::push_back(Type e)
{
  *tail = e;
  tail++;
  if(tail == end)
    tail = entries;
}

} // namespace bliss

#endif
