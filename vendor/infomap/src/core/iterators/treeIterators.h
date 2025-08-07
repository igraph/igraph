/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef TREE_ITERATORS_H_
#define TREE_ITERATORS_H_

#include <cstddef>
#include <iterator>
#include <deque>

namespace infomap {

#ifndef ASSERT
#include <cassert>
#define ASSERT(x) assert(x)
#endif

using std::iterator_traits;

/**
 * Child iterator.
 */
template <typename NodePointerType> // pointer or const pointer
class ChildIterator {
  using iterator_category = std::bidirectional_iterator_tag;
  using value_type = typename iterator_traits<NodePointerType>::value_type;
  using difference_type = typename iterator_traits<NodePointerType>::difference_type;
  using reference = typename iterator_traits<NodePointerType>::reference;
  using pointer = typename iterator_traits<NodePointerType>::pointer;

protected:
  NodePointerType m_root = nullptr;
  NodePointerType m_current = nullptr;

public:
  ChildIterator() = default;

  ChildIterator(const NodePointerType& nodePointer)
      : m_root(nodePointer), m_current(nodePointer == nullptr ? nullptr : nodePointer->firstChild) { }

  ChildIterator(const ChildIterator& other)
      : m_root(other.m_root), m_current(other.m_current) { }

  ChildIterator& operator=(const ChildIterator& other)
  {
    m_root = other.m_root;
    m_current = other.m_current;
    return *this;
  }

  pointer current() const { return m_current; }

  reference operator*() const { return *m_current; }

  pointer operator->() const { return m_current; }

  bool operator==(const ChildIterator& rhs) const { return m_current == rhs.m_current; }

  bool operator!=(const ChildIterator& rhs) const { return !(m_current == rhs.m_current); }

  bool isEnd() const { return m_current == nullptr; }

  ChildIterator& operator++()
  {
    m_current = m_current->next;
    if (m_current != nullptr && m_current->parent != m_root) {
      m_current = nullptr;
    }
    return *this;
  }

  ChildIterator operator++(int)
  {
    ChildIterator copy(*this);
    ++(*this);
    return copy;
  }

  ChildIterator& operator--()
  {
    m_current = m_current->previous;
    if (m_current != nullptr && m_current->parent != m_root) {
      m_current = nullptr;
    }
    return *this;
  }

  ChildIterator operator--(int)
  {
    ChildIterator copy(*this);
    --(*this);
    return copy;
  }
};

/**
 * Tree iterator.
 */
template <typename NodePointerType> // pointer or const pointer
class TreeIterator {
  using iterator_category = std::forward_iterator_tag;
  using value_type = typename iterator_traits<NodePointerType>::value_type;
  using difference_type = typename iterator_traits<NodePointerType>::difference_type;
  using reference = typename iterator_traits<NodePointerType>::reference;
  using pointer = typename iterator_traits<NodePointerType>::pointer;

protected:
  NodePointerType m_root = nullptr;
  NodePointerType m_current = nullptr;
  int m_moduleIndexLevel = -1;
  unsigned int m_moduleIndex = 0;
  std::deque<unsigned int> m_path; // The child index path to current node
  unsigned int m_depth = 0;

public:
  TreeIterator() = default;

  TreeIterator(NodePointerType nodePointer, int moduleIndexLevel = -1)
      : m_root(nodePointer),
        m_current(nodePointer),
        m_moduleIndexLevel(moduleIndexLevel) { }

  TreeIterator(const TreeIterator& other)
      : m_root(other.m_root),
        m_current(other.m_current),
        m_moduleIndexLevel(other.m_moduleIndexLevel),
        m_moduleIndex(other.m_moduleIndex),
        m_path(other.m_path),
        m_depth(other.m_depth) { }

  virtual ~TreeIterator() = default;

  TreeIterator& operator=(const TreeIterator& other)
  {
    m_root = other.m_root;
    m_current = other.m_current;
    m_moduleIndexLevel = other.m_moduleIndexLevel;
    m_moduleIndex = other.m_moduleIndex;
    m_path = other.m_path;
    m_depth = other.m_depth;
    return *this;
  }

  pointer current() const { return m_current; }

  reference operator*() const { return *m_current; }

  pointer operator->() const { return m_current; }

  bool operator==(const TreeIterator& rhs) const { return m_current == rhs.m_current; }

  bool operator!=(const TreeIterator& rhs) const { return !(m_current == rhs.m_current); }

  const std::deque<unsigned int>& path() const { return m_path; }

  unsigned int moduleIndex() const { return m_moduleIndex; }

  unsigned int depth() const { return m_depth; }

  bool isEnd() const { return m_current == nullptr; }

  TreeIterator& operator++()
  {
    NodePointerType curr = m_current;
    NodePointerType infomapRoot = curr->getInfomapRoot();
    if (infomapRoot != nullptr) {
      curr = infomapRoot;
    }

    if (curr->firstChild != nullptr) {
      curr = curr->firstChild;
      ++m_depth;
      m_path.push_back(0);
    } else {
    // Current node is a leaf
    // Presupposes that the next pointer can't reach out from the current parent.
    tryNext:
      while (curr->next == nullptr) {
        if (curr->parent != nullptr) {
          curr = curr->parent;
          --m_depth;
          m_path.pop_back();
          if (curr == m_root) // Check if back to beginning
          {
            m_current = nullptr;
            return *this;
          }
          if (m_moduleIndexLevel < 0) {
            if (curr->isLeafModule()) // TODO: Generalize to -2 for second level to bottom
              ++m_moduleIndex;
          } else if (static_cast<unsigned int>(m_moduleIndexLevel) == m_depth)
            ++m_moduleIndex;
        } else {
          NodePointerType infomapOwner = curr->owner;
          if (infomapOwner != nullptr) {
            curr = infomapOwner;
            if (curr == m_root) // Check if back to beginning
            {
              m_current = nullptr;
              return *this;
            }
            goto tryNext;
          } else // null also if no children in first place
          {
            m_current = nullptr;
            return *this;
          }
        }
      }
      curr = curr->next;
      ++m_path.back();
    }
    m_current = curr;
    return *this;
  }

  TreeIterator operator++(int)
  {
    TreeIterator copy(*this);
    ++(*this);
    return copy;
  }

  TreeIterator& stepForward()
  {
    ++(*this);
    return *this;
  }
};

/**
 * Base node iterator.
 */
template <typename NodePointerType, typename iterator_tag = std::bidirectional_iterator_tag>
struct node_iterator_base {
  using iterator_category = iterator_tag;
  using value_type = typename iterator_traits<NodePointerType>::value_type;
  using difference_type = typename iterator_traits<NodePointerType>::difference_type;
  using reference = typename iterator_traits<NodePointerType>::reference;
  using pointer = typename iterator_traits<NodePointerType>::pointer;

  node_iterator_base() : m_current(nullptr) { }

  node_iterator_base(const NodePointerType& nodePointer) : m_current(nodePointer) { }

  node_iterator_base(const node_iterator_base& other) : m_current(other.m_current) { }

  node_iterator_base& operator=(const node_iterator_base& other)
  {
    m_current = other.m_current;
    return *this;
  }

  virtual ~node_iterator_base() = default;

  pointer base() const { return m_current; }

  reference operator*() const { return *m_current; }

  pointer operator->() const { return m_current; }

  bool operator==(const node_iterator_base& rhs) const { return m_current == rhs.m_current; }

  bool operator!=(const node_iterator_base& rhs) const { return !(m_current == rhs.m_current); }

  bool isEnd() const { return m_current == nullptr; }

protected:
  NodePointerType m_current;
};

template <typename NodePointerType>
class DepthFirstIteratorBase : public node_iterator_base<NodePointerType> {
  using Base = node_iterator_base<NodePointerType>;

public:
  DepthFirstIteratorBase() : Base(), m_root(nullptr), m_depth(0) { }

  DepthFirstIteratorBase(const NodePointerType& nodePointer) : Base(nodePointer), m_root(nodePointer), m_depth(0) { }

  DepthFirstIteratorBase(const DepthFirstIteratorBase& other) : Base(other), m_root(other.m_root), m_depth(other.m_depth) { }

  DepthFirstIteratorBase& operator=(const DepthFirstIteratorBase& other)
  {
    Base::operator=(other);
    m_root = other.m_root;
    m_depth = other.m_depth;
    return *this;
  }

  unsigned int depth() const { return m_depth; }

protected:
  NodePointerType m_root;
  unsigned int m_depth;
  using Base::m_current;
};

/**
 * Pre processing depth first iterator
 * Note:
 * This iterator presupposes that the next pointer of a node can't reach a node with a different parent.
 */
template <typename NodePointerType, bool pre_t = true>
class DepthFirstIterator : public DepthFirstIteratorBase<NodePointerType> {
  using Base = DepthFirstIteratorBase<NodePointerType>;

public:
  DepthFirstIterator() : Base() { }

  DepthFirstIterator(const NodePointerType& nodePointer) : Base(nodePointer) { }

  DepthFirstIterator(const DepthFirstIterator& other) : Base(other) { }

  DepthFirstIterator& operator=(const DepthFirstIterator& other)
  {
    Base::operator=(other);
    return *this;
  }

  DepthFirstIterator& operator++()
  {
    NodePointerType curr = Base::m_current;
    if (curr->firstChild != nullptr) {
      curr = curr->firstChild;
      ++Base::m_depth;
    } else {
      // Presupposes that the next pointer can't reach out from the current parent.
      while (curr->next == nullptr) {
        curr = curr->parent;
        --Base::m_depth;
        if (curr == Base::m_root || curr == nullptr) // 0 if no children in first place
        {
          Base::m_current = nullptr;
          return *this;
        }
      }
      curr = curr->next;
    }
    Base::m_current = curr;
    return *this;
  }

  DepthFirstIterator operator++(int)
  {
    auto copy(*this);
    ++(*this);
    return copy;
  }

  DepthFirstIterator next()
  {
    auto copy(*this);
    return ++copy;
  }
};

/**
 * Post processing depth first iterator
 * Note:
 * This iterator presupposes that the next pointer of a node can't reach a node with a different parent.
 */
template <typename NodePointerType>
class DepthFirstIterator<NodePointerType, false> : public DepthFirstIteratorBase<NodePointerType> {
  using Base = DepthFirstIteratorBase<NodePointerType>;

public:
  DepthFirstIterator() : Base() { }

  DepthFirstIterator(const NodePointerType& nodePointer) : Base(nodePointer) { init(); }

  DepthFirstIterator(const DepthFirstIterator& other) : Base(other) { }

  DepthFirstIterator& operator=(const DepthFirstIterator& other)
  {
    Base::operator=(other);
    return *this;
  }

  void init()
  {
    if (Base::m_current != nullptr) {
      while (Base::m_current->firstChild != nullptr) {
        Base::m_current = Base::m_current->firstChild;
        ++Base::m_depth;
      }
    }
  }

  DepthFirstIterator& operator++()
  {
    // The root should be the last node
    if (Base::m_current == Base::m_root) {
      Base::m_current = nullptr;
      return *this;
    }

    NodePointerType curr = Base::m_current;
    if (curr->next != nullptr) {
      curr = curr->next;
      while (curr->firstChild != nullptr) {
        curr = curr->firstChild;
        ++Base::m_depth;
      }
    } else {
      curr = curr->parent;
      --Base::m_depth;
    }

    Base::m_current = curr;

    return *this;
  }

  DepthFirstIterator operator++(int)
  {
    DepthFirstIterator copy(*this);
    ++(*this);
    return copy;
  }

  DepthFirstIterator
  next()
  {
    DepthFirstIterator copy(*this);
    return ++copy;
  }
};

/**
 * Leaf node iterator
 */
template <typename NodePointerType>
class LeafNodeIterator : public DepthFirstIteratorBase<NodePointerType> {
  using Base = DepthFirstIteratorBase<NodePointerType>;

public:
  LeafNodeIterator() : Base() { }

  LeafNodeIterator(const NodePointerType& nodePointer) : Base(nodePointer) { init(); }

  LeafNodeIterator(const LeafNodeIterator& other) : Base(other) { }

  LeafNodeIterator& operator=(const LeafNodeIterator& other)
  {
    Base::operator=(other);
    return *this;
  }

  void init()
  {
    if (Base::m_current != nullptr) {
      while (Base::m_current->firstChild != nullptr) {
        Base::m_current = Base::m_current->firstChild;
        ++Base::m_depth;
      }
    }
  }

  LeafNodeIterator& operator++()
  {
    ASSERT(Base::m_current != nullptr);
    while (Base::m_current->next == nullptr || Base::m_current->next->parent != Base::m_current->parent) {
      Base::m_current = Base::m_current->parent;
      --Base::m_depth;
      if (Base::m_current == nullptr)
        return *this;
    }

    Base::m_current = Base::m_current->next;

    if (Base::m_current != nullptr) {
      while (Base::m_current->firstChild != nullptr) {
        Base::m_current = Base::m_current->firstChild;
        ++Base::m_depth;
      }
    }
    return *this;
  }

  LeafNodeIterator operator++(int)
  {
    LeafNodeIterator copy(*this);
    ++(*this);
    return copy;
  }

  LeafNodeIterator next()
  {
    LeafNodeIterator copy(*this);
    return ++copy;
  }
};

/**
 * Leaf module iterator
 */
template <typename NodePointerType>
class LeafModuleIterator : public DepthFirstIteratorBase<NodePointerType> {
  using Base = DepthFirstIteratorBase<NodePointerType>;

public:
  LeafModuleIterator() : Base() { }

  LeafModuleIterator(const NodePointerType& nodePointer) : Base(nodePointer) { init(); }

  LeafModuleIterator(const LeafModuleIterator& other) : Base(other) { }

  LeafModuleIterator& operator=(const LeafModuleIterator& other)
  {
    Base::operator=(other);
    init();
    return *this;
  }

  void init()
  {
    if (Base::m_current != nullptr) {
      if (Base::m_current->firstChild == nullptr) {
        Base::m_current = nullptr; // End directly if no module
      } else {
        while (Base::m_current->firstChild->firstChild != nullptr) {
          Base::m_current = Base::m_current->firstChild;
          ++Base::m_depth;
        }
      }
    }
  }

  LeafModuleIterator& operator++()
  {
    ASSERT(Base::m_current != nullptr);
    while (Base::m_current->next == nullptr || Base::m_current->next->parent != Base::m_current->parent) {
      Base::m_current = Base::m_current->parent;
      --Base::m_depth;
      if (Base::m_current == nullptr)
        return *this;
    }

    Base::m_current = Base::m_current->next;

    if (Base::m_current != nullptr) {
      if (Base::m_current->firstChild == nullptr) {
        Base::m_current = Base::m_current->parent;
      } else {
        while (Base::m_current->firstChild->firstChild != nullptr) {
          Base::m_current = Base::m_current->firstChild;
          ++Base::m_depth;
        }
      }
    }
    return *this;
  }

  LeafModuleIterator operator++(int)
  {
    LeafModuleIterator copy(*this);
    ++(*this);
    return copy;
  }

  LeafModuleIterator next()
  {
    LeafModuleIterator copy(*this);
    return ++copy;
  }
};

/**
 * Sibling iterator.
 */
template <typename NodePointerType> // pointer or const pointer
class SiblingIterator : public node_iterator_base<NodePointerType> {
  using Base = node_iterator_base<NodePointerType>;

public:
  using self_type = SiblingIterator<NodePointerType>;

  SiblingIterator() : Base() { }

  SiblingIterator(const NodePointerType& nodePointer) : Base(nodePointer) { }

  SiblingIterator(const SiblingIterator& other) : Base(other) { }

  SiblingIterator& operator=(const SiblingIterator& other)
  {
    Base::operator=(other);
    return *this;
  }

  SiblingIterator& operator++()
  {
    ASSERT(Base::m_current != nullptr);
    Base::m_current = Base::m_current->next;
    return *this;
  }

  SiblingIterator operator++(int)
  {
    SiblingIterator copy(*this);
    ++(*this);
    return copy;
  }

  SiblingIterator& operator--()
  {
    ASSERT(Base::m_current != nullptr);
    Base::m_current = Base::m_current->previous;
    return *this;
  }

  SiblingIterator operator--(int)
  {
    SiblingIterator copy(*this);
    --(*this);
    return copy;
  }
};

} // namespace infomap

#endif // TREE_ITERATORS_H_
