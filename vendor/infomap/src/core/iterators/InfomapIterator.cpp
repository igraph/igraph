/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "InfomapIterator.h"
#include "../InfoNode.h"
#include <iostream>
#include <utility> // std::pair

namespace infomap {

InfomapIterator& InfomapIterator::operator++() noexcept
{
  const auto root = m_current->getInfomapRoot();
  auto current = root ? root : m_current;

  if (current->firstChild) {
    current = current->firstChild;
    ++m_depth;
    m_path.push_back(1);
  } else {
    // Current node is a leaf
    // Presupposes that the next pointer can't reach out from the current parent.

    bool tryNext = true;

    while (tryNext) {
      tryNext = false;

      while (!current->next) {
        if (current->owner) {
          current = current->owner;

          if (current == m_root) { // Check if back to beginning
            m_current = nullptr;
            return *this;
          }

          tryNext = true;
          break;
        }
        if (current->parent) {
          current = current->parent;
          --m_depth;
          m_path.pop_back();

          if (current == m_root) { // Check if back to beginning
            m_current = nullptr;
            return *this;
          }
        } else { // null also if no children in first place
          m_current = nullptr;
          return *this;
        }
      }
    }

    current = current->next;

    if (!current->isLeaf() && (static_cast<unsigned int>(m_moduleIndexLevel) >= m_depth)) {
      ++m_moduleIndex;
    }

    ++m_path.back();
  }

  m_current = current;
  return *this;
}

double InfomapIterator::modularCentrality() const noexcept
{
  if (m_current->parent == nullptr) {
    // The root node has no modular centrality
    return 0.0;
  }

  const auto p_m = m_current->parent->data.flow;
  const auto p_u = m_current->data.flow;
  const auto p_diff = p_m - p_u;

  if (p_diff > 0.0) {
    return -p_diff * std::log2(p_diff / p_m);
  }

  return 0.0;
}

// -------------------------------------
// InfomapModuleIterator
// -------------------------------------

InfomapIterator& InfomapModuleIterator::operator++() noexcept
{
  InfomapIterator::operator++();
  while (!isEnd() && m_current->isLeaf()) {
    InfomapIterator::operator++();
  }
  return *this;
}

// -------------------------------------
// InfomapLeafModuleIterator
// -------------------------------------

void InfomapLeafModuleIterator::init() noexcept
{
  while (!isEnd() && !m_current->isLeafModule()) {
    InfomapIterator::operator++();
  }
}

InfomapIterator& InfomapLeafModuleIterator::operator++() noexcept
{
  InfomapIterator::operator++();
  while (!isEnd() && !m_current->isLeafModule()) {
    InfomapIterator::operator++();
  }
  return *this;
}

// -------------------------------------
// InfomapLeafIterator
// -------------------------------------

void InfomapLeafIterator::init() noexcept
{
  while (!isEnd() && !m_current->isLeaf()) {
    InfomapIterator::operator++();
  }
}

InfomapIterator& InfomapLeafIterator::operator++() noexcept
{
  InfomapIterator::operator++();
  while (!isEnd() && !m_current->isLeaf()) {
    InfomapIterator::operator++();
  }
  return *this;
}

// -------------------------------------
// InfomapIteratorPhysical
// -------------------------------------

InfomapIterator& InfomapIteratorPhysical::operator++() noexcept
{
  if (m_physNodes.empty()) {
    // Iterate modules
    InfomapIterator::operator++();
    if (isEnd()) {
      return *this;
    }
    if (m_current->isLeaf()) {
      // Copy current iterator to restart after iterating through the leaf nodes
      auto firstLeafIt = *this;
      // If on a leaf node, loop through and aggregate to physical nodes
      while (!isEnd() && m_current->isLeaf()) {
        auto ret = m_physNodes.insert(std::make_pair(m_current->physicalId, InfoNode(*m_current)));
        auto& physNode = ret.first->second;
        if (ret.second) {
          // New physical node, use same parent as the state leaf node
          physNode.parent = m_current->parent;
        } else {
          // Not inserted, add flow to existing physical node
          // TODO: If exitFlow should be correct, flow between memory nodes within same physical node should be subtracted.
          physNode.data += m_current->data;
        }
        physNode.stateNodes.push_back(m_current->stateId);
        InfomapIterator::operator++();
      }
      // Store current iterator to continue with after iterating physical leaf nodes
      m_oldIter = *this;
      // Reset path/depth/moduleIndex to values for first leaf node
      m_path = firstLeafIt.m_path;
      m_depth = firstLeafIt.m_depth;
      m_moduleIndex = firstLeafIt.m_moduleIndex;
      // Set current node to the first physical node
      m_physIter = m_physNodes.begin();
      m_current = &m_physIter->second;
    }
  } else {
    // Iterate physical nodes instead of leaf state nodes
    ++m_physIter;
    ++m_path.back();
    if (m_physIter == m_physNodes.end()) {
      // End of leaf nodes
      m_physNodes.clear();
      m_path.pop_back();
      // reset iterator to the one after the leaf nodes
      *this = m_oldIter;
    } else {
      // Set iterator node to the currently iterated physical node
      m_current = &m_physIter->second;
    }
  }
  return *this;
}

// -------------------------------------
// InfomapLeafIteratorPhysical
// -------------------------------------

void InfomapLeafIteratorPhysical::init() noexcept
{
  while (!isEnd() && !m_current->isLeaf()) {
    InfomapIteratorPhysical::operator++();
  }
}

InfomapIterator& InfomapLeafIteratorPhysical::operator++() noexcept
{
  InfomapIteratorPhysical::operator++();
  while (!isEnd() && !m_current->isLeaf()) {
    InfomapIteratorPhysical::operator++();
  }
  return *this;
}

// -------------------------------------
// InfomapParentIterator
// -------------------------------------

InfomapParentIterator& InfomapParentIterator::operator++() noexcept
{
  m_current = m_current->parent;
  if (m_current != nullptr && m_current->owner != nullptr) {
    m_current = m_current->owner;
  }
  return *this;
}

} // namespace infomap
