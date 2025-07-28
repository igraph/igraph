/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef META_COLLECTION_H_
#define META_COLLECTION_H_

#include "infomath.h"
#include <ostream>
#include <map>

namespace infomap {

struct FlowCount {
  FlowCount() = default;

  explicit FlowCount(double flow)
      : flow(flow), count(1) { }

  FlowCount& operator+=(const FlowCount& o)
  {
    flow += o.flow;
    count += o.count;
    return *this;
  }

  FlowCount& operator+=(double f)
  {
    flow += f;
    ++count;
    return *this;
  }

  FlowCount& operator-=(const FlowCount& o)
  {
    flow -= o.flow;
    count -= o.count;
    return *this;
  }

  FlowCount& operator-=(double f)
  {
    flow -= f;
    --count;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& out, const FlowCount& o)
  {
    return out << o.flow << "/" << o.count;
  }

  void reset()
  {
    flow = 0.0;
    count = 0;
  }

  bool empty() const { return count == 0; }
  double flow = 0.0;
  unsigned int count = 0;
};

using MetaToFlowCount = std::map<unsigned int, FlowCount>; // metaId -> (flow,count)

class MetaCollection {
protected:
  FlowCount m_total;
  MetaToFlowCount m_metaToFlowCount;

public:
  unsigned int size() const { return m_metaToFlowCount.size(); }

  bool empty() const { return m_metaToFlowCount.empty(); }

  MetaToFlowCount::iterator begin() { return m_metaToFlowCount.begin(); }
  MetaToFlowCount::iterator end() { return m_metaToFlowCount.end(); }
  MetaToFlowCount::const_iterator begin() const { return m_metaToFlowCount.begin(); }
  MetaToFlowCount::const_iterator end() const { return m_metaToFlowCount.end(); }

  void add(unsigned int meta, double flow = 1.0)
  {
    m_total += flow;
    m_metaToFlowCount[meta] += flow;
  }

  void add(unsigned int meta, const FlowCount& flow)
  {
    m_total += flow;
    m_metaToFlowCount[meta] += flow;
  }

  void add(const MetaCollection& other)
  {
    for (auto& it : other) {
      auto metaId = it.first;
      auto& flowCount = it.second;
      add(metaId, flowCount);
    }
  }

  void remove(unsigned int meta, const FlowCount& flow)
  {
    m_total -= flow;
    auto& metaFlowCount = m_metaToFlowCount[meta];
    metaFlowCount -= flow;
    if (metaFlowCount.empty())
      m_metaToFlowCount.erase(meta);
  }

  void remove(const MetaCollection& other)
  {
    for (auto& it : other) {
      auto metaId = it.first;
      auto& flowCount = it.second;
      remove(metaId, flowCount);
    }
  }

  double calculateEntropy()
  {
    double metaCodelength = 0.0;
    for (auto& it : m_metaToFlowCount) {
      metaCodelength -= infomath::plogp(it.second.flow / m_total.flow);
    }
    return m_total.flow * metaCodelength;
  }

  void clear()
  {
    m_total.reset();
    m_metaToFlowCount.clear();
  }

  friend std::ostream& operator<<(std::ostream& out, const MetaCollection& m)
  {
    out << "<( " << m.m_total << ": ";
    for (auto& it : m.m_metaToFlowCount) {
      out << "(" << it.first << "," << it.second << ") ";
    }
    out << ")>";
    return out;
  }
};

} // namespace infomap

#endif // META_COLLECTION_H_
