/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef VECTOR_MAP_H_
#define VECTOR_MAP_H_

#include <vector>
#include <limits>
#include <map>

namespace infomap {

template <typename T>
class VectorMap {
public:
  VectorMap(unsigned int capacity = 0)
      : m_capacity(capacity),
        m_values(capacity),
        m_redirect(capacity, 0),
        m_maxOffset(std::numeric_limits<unsigned int>::max() - 1 - capacity) { }

  void startRound()
  {
    if (m_size > 0) {
      m_offset += m_capacity;
      m_size = 0;
    }
    if (m_offset > m_maxOffset) {
      m_redirect.assign(m_capacity, 0);
      m_offset = 1;
    }
  }

  void add(unsigned int index, T value)
  {
    if (isSet(index)) {
      m_values[m_redirect[index] - m_offset] += value;
    } else {
      m_redirect[index] = m_offset + m_size;
      m_values[m_size] = value;
      ++m_size;
    }
  }

  bool isSet(unsigned int index)
  {
    return m_redirect[index] >= m_offset;
  }

  unsigned int size()
  {
    return m_size;
  }

  T& operator[](unsigned int index)
  {
    return m_values[m_redirect[index] - m_offset];
  }

  std::vector<T>& values()
  {
    return m_values;
  }

private:
  unsigned int m_capacity = 0;
  std::vector<T> m_values;
  std::vector<unsigned int> m_redirect;
  unsigned int m_maxOffset = std::numeric_limits<unsigned int>::max() - 1;
  unsigned int m_offset = 1;
  unsigned int m_size = 0;
};

} // namespace infomap

#endif // VECTOR_MAP_H_
