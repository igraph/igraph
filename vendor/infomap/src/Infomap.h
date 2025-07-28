/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFOMAP_H_
#define INFOMAP_H_

#include "core/InfomapBase.h"
#include "io/Config.h"

#include <string>
#include <utility>
#include <map>

namespace infomap {

// Wrapper class for the Python API
struct InfomapWrapper : public InfomapBase {
public:
  InfomapWrapper() : InfomapBase() { }
  InfomapWrapper(const std::string& flags) : InfomapBase(flags) { }
  InfomapWrapper(const Config& conf) : InfomapBase(conf) { }
  virtual ~InfomapWrapper() = default;

  // ===================================================
  // Wrapper methods
  // ===================================================

  void readInputData(std::string filename = "", bool accumulate = true) { m_network.readInputData(std::move(filename), accumulate); }

  void addNode(unsigned int id) { m_network.addNode(id); }
  void addNode(unsigned int id, std::string name) { m_network.addNode(id, std::move(name)); }
  void addNode(unsigned int id, double weight) { m_network.addNode(id, weight); }
  void addNode(unsigned int id, std::string name, double weight) { m_network.addNode(id, std::move(name), weight); }

  void addName(unsigned int id, const std::string& name) { m_network.addName(id, name); }
  std::string getName(unsigned int id) const
  {
    auto& names = m_network.names();
    auto it = names.find(id);
    return it != names.end() ? it->second : "";
  }

  const std::map<unsigned int, std::string>& getNames() const { return m_network.names(); }

  void addPhysicalNode(unsigned int id, const std::string& name = "") { m_network.addPhysicalNode(id, name); }
  void addStateNode(unsigned int id, unsigned int physId) { m_network.addStateNode(id, physId); }

  void addLink(unsigned int sourceId, unsigned int targetId, double weight = 1.0) { m_network.addLink(sourceId, targetId, weight); }
  void addLink(unsigned int sourceId, unsigned int targetId, unsigned long weight) { m_network.addLink(sourceId, targetId, weight); }
  void addMultilayerLink(unsigned int layer1, unsigned int n1, unsigned int layer2, unsigned int n2, double weight = 1.0) { m_network.addMultilayerLink(layer1, n1, layer2, n2, weight); }
  void addMultilayerIntraLink(unsigned int layer, unsigned int n1, unsigned int n2, double weight) { m_network.addMultilayerIntraLink(layer, n1, n2, weight); }
  void addMultilayerInterLink(unsigned int layer1, unsigned int n, unsigned int layer2, double interWeight) { m_network.addMultilayerInterLink(layer1, n, layer2, interWeight); }

  void setBipartiteStartId(unsigned int startId) { m_network.setBipartiteStartId(startId); }

  std::map<std::pair<unsigned int, unsigned int>, double> getLinks(bool flow) const
  {
    std::map<std::pair<unsigned int, unsigned int>, double> links;

    for (const auto& node : m_network.nodeLinkMap()) {
      const auto sourceId = node.first.id;

      for (const auto& link : node.second) {
        const auto targetId = link.first.id;
        links[{ sourceId, targetId }] = flow ? link.second.flow : link.second.weight;
      }
    }

    return links;
  }

  std::map<unsigned int, unsigned int> getModules(int level = 1, bool states = false)
  {
    if (haveMemory() && !states) {
      throw std::runtime_error("Cannot get modules on higher-order network without states.");
    }
    std::map<unsigned int, unsigned int> modules;
    for (auto it = iterTree(level); !it.isEnd(); ++it) {
      auto& node = *it;
      if (node.isLeaf()) {
        modules[states ? node.stateId : node.physicalId] = it.moduleId();
      }
    }
    return modules;
  }

  using InfomapBase::codelength;
  using InfomapBase::getEntropyRate;
  using InfomapBase::getMultilevelModules;
  using InfomapBase::iterLeafNodes;
  using InfomapBase::iterTree;
  using InfomapBase::run;
};

} // namespace infomap

#endif // INFOMAP_H_
