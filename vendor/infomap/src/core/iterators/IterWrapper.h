/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef ITER_WRAPPER_H_
#define ITER_WRAPPER_H_

#include <type_traits>
#include <utility>

namespace infomap {

template <typename Iter>
class IterWrapper {
  Iter m_begin, m_end;

public:
  IterWrapper(const Iter& begin, const Iter& end) : m_begin(begin), m_end(end) {}

  // template <typename Container>
  template <typename Container,
            typename = decltype(std::declval<Container>().begin()),
            typename = decltype(std::declval<Container>().end())>
  IterWrapper(Container& container) : m_begin(container.begin()), m_end(container.end()) {}

  template <typename Container,
            typename std::enable_if<std::is_convertible<Container&, IterWrapper<Iter>&>::value, int>::type = 0>
  IterWrapper(Container& container)
      : m_begin(static_cast<IterWrapper<Iter>&>(container).begin()),
        m_end(static_cast<IterWrapper<Iter>&>(container).end()) {}

  Iter begin() noexcept { return m_begin; };

  Iter end() noexcept { return m_end; };
};

} // namespace infomap

#endif // ITER_WRAPPER_H_
