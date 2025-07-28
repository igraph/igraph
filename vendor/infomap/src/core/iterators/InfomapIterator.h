/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFOMAP_ITERATOR_H_
#define INFOMAP_ITERATOR_H_

#include <deque>
#include <map>
#include <cmath>

namespace infomap {

class InfoNode;

/**
 * Pre processing depth first iterator that explores sub-Infomap instances
 * Note:
 * This iterator presupposes that the next pointer of a node can't reach a node with a different parent.
 */
struct InfomapIterator {
protected:
  InfoNode* m_root = nullptr;
  InfoNode* m_current = nullptr;
  int m_moduleIndexLevel = -1;
  unsigned int m_moduleIndex = 0;
  std::deque<unsigned int> m_path; // The tree path to current node (indexing starting from one!)
  unsigned int m_depth = 0;

public:
  InfomapIterator() = default;

  InfomapIterator(InfoNode* nodePointer, int moduleIndexLevel = -1)
      : m_root(nodePointer), m_current(nodePointer), m_moduleIndexLevel(moduleIndexLevel) { }

  virtual ~InfomapIterator() = default;
  InfomapIterator(const InfomapIterator&) = default;
  InfomapIterator& operator=(const InfomapIterator&) = default;
  InfomapIterator(InfomapIterator&&) noexcept = default;
  InfomapIterator& operator=(InfomapIterator&&) noexcept = default;

  InfoNode* current() noexcept { return m_current; }

  const InfoNode* current() const noexcept { return m_current; }

  InfoNode& operator*() noexcept { return *m_current; }

  const InfoNode& operator*() const noexcept { return *m_current; }

  InfoNode* operator->() noexcept { return m_current; }

  const InfoNode* operator->() const noexcept { return m_current; }

  bool operator==(const InfomapIterator& other) const noexcept { return m_current == other.m_current; }

  bool operator!=(const InfomapIterator& other) const noexcept { return m_current != other.m_current; }

  virtual InfomapIterator& operator++() noexcept;

  virtual InfomapIterator operator++(int) noexcept
  {
    InfomapIterator copy(*this);
    ++(*this);
    return copy;
  }

  virtual InfomapIterator& stepForward() noexcept
  {
    ++(*this);
    return *this;
  }

  const std::deque<unsigned int>& path() const noexcept { return m_path; }

  unsigned int moduleIndex() const noexcept { return m_moduleIndex; }

  unsigned int moduleId() const noexcept { return m_moduleIndex + 1; }

  unsigned int childIndex() const noexcept { return m_path.empty() ? 0 : m_path.back() - 1; }

  unsigned int depth() const noexcept { return m_depth; }

  double modularCentrality() const noexcept;

  bool isEnd() const noexcept { return m_current == nullptr; }
};

struct InfomapModuleIterator : public InfomapIterator {
public:
  InfomapModuleIterator() : InfomapIterator() { }

  InfomapModuleIterator(InfoNode* nodePointer, int moduleIndexLevel = -1) : InfomapIterator(nodePointer, moduleIndexLevel) { }

  ~InfomapModuleIterator() override = default;
  InfomapModuleIterator(const InfomapModuleIterator&) = default;
  InfomapModuleIterator& operator=(const InfomapModuleIterator&) = default;
  InfomapModuleIterator(InfomapModuleIterator&&) = default;
  InfomapModuleIterator& operator=(InfomapModuleIterator&&) = default;

  InfomapIterator& operator++() noexcept override;

  InfomapIterator operator++(int) noexcept override
  {
    InfomapModuleIterator copy(*this);
    ++(*this);
    return std::move(copy);
  }

  using InfomapIterator::childIndex;
  using InfomapIterator::current;
  using InfomapIterator::depth;
  using InfomapIterator::modularCentrality;
  using InfomapIterator::path;
};

struct InfomapLeafModuleIterator : public InfomapIterator {
public:
  InfomapLeafModuleIterator() : InfomapIterator() { }

  InfomapLeafModuleIterator(InfoNode* nodePointer, int moduleIndexLevel = -1)
      : InfomapIterator(nodePointer, moduleIndexLevel) { init(); }

  ~InfomapLeafModuleIterator() override = default;
  InfomapLeafModuleIterator(const InfomapLeafModuleIterator& other) : InfomapIterator(other) { init(); }
  InfomapLeafModuleIterator& operator=(const InfomapLeafModuleIterator&) = default;
  InfomapLeafModuleIterator(InfomapLeafModuleIterator&&) = default;
  InfomapLeafModuleIterator& operator=(InfomapLeafModuleIterator&&) = default;

  /**
   * Iterate to first leaf module
   */
  void init() noexcept;

  InfomapIterator& operator++() noexcept override;

  InfomapIterator operator++(int) noexcept override
  {
    InfomapLeafModuleIterator copy(*this);
    ++(*this);
    return std::move(copy);
  }

  using InfomapIterator::childIndex;
  using InfomapIterator::current;
  using InfomapIterator::depth;
  using InfomapIterator::modularCentrality;
  using InfomapIterator::path;
};

struct InfomapLeafIterator : public InfomapIterator {
public:
  InfomapLeafIterator() : InfomapIterator() { }

  InfomapLeafIterator(InfoNode* nodePointer, int moduleIndexLevel = -1)
      : InfomapIterator(nodePointer, moduleIndexLevel) { init(); }

  ~InfomapLeafIterator() override = default;
  InfomapLeafIterator(const InfomapLeafIterator& other) : InfomapIterator(other) { init(); }
  InfomapLeafIterator& operator=(const InfomapLeafIterator&) = default;
  InfomapLeafIterator(InfomapLeafIterator&&) = default;
  InfomapLeafIterator& operator=(InfomapLeafIterator&&) = default;

  /**
   * Iterate to first leaf node
   */
  void init() noexcept;

  InfomapIterator& operator++() noexcept override;

  InfomapIterator operator++(int) noexcept override
  {
    InfomapLeafIterator copy(*this);
    ++(*this);
    return std::move(copy);
  }

  using InfomapIterator::childIndex;
  using InfomapIterator::current;
  using InfomapIterator::depth;
  using InfomapIterator::modularCentrality;
  using InfomapIterator::path;
};

/**
 * Iterate over the whole tree, collecting physical nodes within same leaf modules
 * Note: The physical nodes are created when entering the parent module and removed
 * when leaving the module. The tree will not be modified.
 */
struct InfomapIteratorPhysical : public InfomapIterator {
protected:
  std::map<unsigned int, InfoNode> m_physNodes;
  std::map<unsigned int, InfoNode>::iterator m_physIter;
  InfomapIterator m_oldIter;

public:
  InfomapIteratorPhysical() : InfomapIterator() { }

  InfomapIteratorPhysical(InfoNode* nodePointer, int moduleIndexLevel = -1)
      : InfomapIterator(nodePointer, moduleIndexLevel) { }

  ~InfomapIteratorPhysical() override = default;
  InfomapIteratorPhysical(const InfomapIteratorPhysical&) = default;
  InfomapIteratorPhysical(const InfomapIterator& other) : InfomapIterator(other) { }

  InfomapIteratorPhysical(InfomapIteratorPhysical&&) = default;
  InfomapIteratorPhysical& operator=(const InfomapIteratorPhysical&) = default;

  // Don't allow moving from this iterator as we use the old iterator in operator++
  InfomapIteratorPhysical& operator=(InfomapIteratorPhysical&&) = delete;

  InfomapIteratorPhysical& operator=(const InfomapIterator& other)
  {
    InfomapIterator::operator=(other);
    return *this;
  }

  InfomapIterator& operator++() noexcept override;

  InfomapIterator operator++(int) noexcept override
  {
    InfomapIteratorPhysical copy(*this);
    ++(*this);
    return std::move(copy);
  }

  using InfomapIterator::childIndex;
  using InfomapIterator::current;
  using InfomapIterator::depth;
  using InfomapIterator::modularCentrality;
  using InfomapIterator::path;
};

/**
 * Iterate over all physical leaf nodes, joining physical nodes within same leaf modules
 * Note: The physical nodes are created when entering the parent module and removed
 * when leaving the module. The tree will not be modified.
 */
struct InfomapLeafIteratorPhysical : public InfomapIteratorPhysical {
public:
  InfomapLeafIteratorPhysical() : InfomapIteratorPhysical() { }

  InfomapLeafIteratorPhysical(InfoNode* nodePointer, int moduleIndexLevel = -1)
      : InfomapIteratorPhysical(nodePointer, moduleIndexLevel) { init(); }

  InfomapLeafIteratorPhysical(const InfomapLeafIteratorPhysical& other)
      : InfomapIteratorPhysical(other) { init(); }

  ~InfomapLeafIteratorPhysical() override = default;
  InfomapLeafIteratorPhysical(InfomapLeafIteratorPhysical&&) = default;
  InfomapLeafIteratorPhysical& operator=(const InfomapLeafIteratorPhysical&) = default;

  // Don't allow moving from this iterator as we use the old iterator in operator++
  InfomapLeafIteratorPhysical& operator=(InfomapLeafIteratorPhysical&&) = delete;

  /**
   * Iterate to first leaf node
   */
  void init() noexcept;

  InfomapIterator& operator++() noexcept override;

  InfomapIterator operator++(int) noexcept override
  {
    InfomapLeafIteratorPhysical copy(*this);
    ++(*this);
    return std::move(copy);
  }

  using InfomapIteratorPhysical::childIndex;
  using InfomapIteratorPhysical::current;
  using InfomapIteratorPhysical::depth;
  using InfomapIteratorPhysical::modularCentrality;
  using InfomapIteratorPhysical::path;
};

/**
 * Iterate parent by parent until it is nullptr,
 * moving up through possible sub infomap instances
 * on the way
 */
struct InfomapParentIterator {
protected:
  InfoNode* m_current = nullptr;

public:
  InfomapParentIterator() = default;

  InfomapParentIterator(InfoNode* nodePointer) : m_current(nodePointer) { }

  ~InfomapParentIterator() = default;

  InfomapParentIterator(const InfomapParentIterator&) = default;

  InfomapParentIterator& operator=(const InfomapParentIterator&) = default;

  InfomapParentIterator(InfomapParentIterator&&) = default;

  InfomapParentIterator& operator=(InfomapParentIterator&&) = default;

  InfoNode* current() noexcept { return m_current; }

  const InfoNode* current() const noexcept { return m_current; }

  InfoNode& operator*() noexcept { return *m_current; }

  const InfoNode& operator*() const noexcept { return *m_current; }

  InfoNode* operator->() noexcept { return m_current; }

  const InfoNode* operator->() const noexcept { return m_current; }

  bool operator==(const InfomapParentIterator& other) const noexcept { return m_current == other.m_current; }

  bool operator!=(const InfomapParentIterator& other) const noexcept { return m_current != other.m_current; }

  InfomapParentIterator& operator++() noexcept;

  InfomapParentIterator operator++(int) noexcept
  {
    InfomapParentIterator copy(*this);
    ++(*this);
    return copy;
  }

  InfomapParentIterator& stepForward() noexcept
  {
    ++(*this);
    return *this;
  }

  bool isEnd() const noexcept { return m_current == nullptr; }
};

} // namespace infomap

#endif // INFOMAP_ITERATOR_H_
