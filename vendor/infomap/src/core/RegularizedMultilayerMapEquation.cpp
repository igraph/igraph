/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "RegularizedMultilayerMapEquation.h"
#include "FlowData.h"
#include "InfoNode.h"
#include "../utils/format.h"

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace infomap {

// ===================================================
// IO
// ===================================================

std::ostream& RegularizedMultilayerMapEquation::print(std::ostream& out) const
{
  return out << indexCodelength << " + " << moduleCodelength << " = " << io::toPrecision(codelength);
}

std::ostream& operator<<(std::ostream& out, const RegularizedMultilayerMapEquation& mapEq)
{
  return mapEq.print(out);
}

// ===================================================
// Init
// ===================================================

void RegularizedMultilayerMapEquation::init(const Config& config)
{
  Log(3) << "RegularizedMultilayerMapEquation::init()...\n";
  Base::init(config);
}

void RegularizedMultilayerMapEquation::initNetwork(InfoNode& root)
{
  initPhysicalNodes(root);
}

void RegularizedMultilayerMapEquation::initSuperNetwork(InfoNode& /*root*/)
{
  // TODO: How use enterFlow instead of flow
}

void RegularizedMultilayerMapEquation::initSubNetwork(InfoNode& /*root*/)
{
}

void RegularizedMultilayerMapEquation::initPartition(std::vector<InfoNode*>& nodes)
{
  initPartitionOfPhysicalNodes(nodes);
  initPartitionLayerTeleFlowData(nodes);

  calculateCodelength(nodes);
}

void RegularizedMultilayerMapEquation::initPhysicalNodes(InfoNode& root)
{
  bool notInitiatedOnRoot = root.physicalNodes.empty();
  if (notInitiatedOnRoot) {
    // Assume leaf nodes directly under the root node
    std::unordered_map<unsigned int, double> physicalNodes;
    unsigned int maxPhysicalId = 0;
    unsigned int minPhysicalId = std::numeric_limits<unsigned int>::max();
    for (auto it(root.begin_leaf_nodes()); !it.isEnd(); ++it) {
      InfoNode& node = *it;
      physicalNodes[node.physicalId] += node.data.flow;
      minPhysicalId = std::min(minPhysicalId, node.physicalId);
      maxPhysicalId = std::max(maxPhysicalId, node.physicalId);
    }

    // Re-index physical nodes if necessary
    std::map<unsigned int, unsigned int> toZeroBasedIndex;
    if (maxPhysicalId - minPhysicalId + 1 > m_numPhysicalNodes) {
      unsigned int zeroBasedPhysicalId = 0;
      for (const auto& physNode : physicalNodes) {
        toZeroBasedIndex.insert(std::make_pair(physNode.first, zeroBasedPhysicalId++));
      }
    }

    for (const auto& physNode : physicalNodes) {
      unsigned int zeroBasedIndex = !toZeroBasedIndex.empty() ? toZeroBasedIndex[physNode.first] : (physNode.first - minPhysicalId);
      root.physicalNodes.emplace_back(zeroBasedIndex, physNode.second);
    }
  }
  auto firstLeafIt = root.begin_leaf_nodes();
  auto depth = firstLeafIt.depth();
  bool notInitiatedOnLeafNodes = firstLeafIt->physicalNodes.empty();
  if (notInitiatedOnLeafNodes) {
    Log(3) << "RegularizedMultilayerMapEquation::initPhysicalNodesOnOriginalNetwork()...\n";
    std::set<unsigned int> setOfPhysicalNodes;
    unsigned int maxPhysicalId = 0;
    unsigned int minPhysicalId = std::numeric_limits<unsigned int>::max();
    for (auto it(root.begin_leaf_nodes()); !it.isEnd(); ++it) {
      InfoNode& node = *it;
      setOfPhysicalNodes.insert(node.physicalId);
      maxPhysicalId = std::max(maxPhysicalId, node.physicalId);
      minPhysicalId = std::min(minPhysicalId, node.physicalId);
    }

    m_numPhysicalNodes = setOfPhysicalNodes.size();

    // Re-index physical nodes if necessary
    std::map<unsigned int, unsigned int> toZeroBasedIndex;
    if (maxPhysicalId - minPhysicalId + 1 > m_numPhysicalNodes) {
      unsigned int zeroBasedPhysicalId = 0;
      for (unsigned int physIndex : setOfPhysicalNodes) {
        toZeroBasedIndex.insert(std::make_pair(physIndex, zeroBasedPhysicalId++));
      }
    }

    for (auto it(root.begin_leaf_nodes()); !it.isEnd(); ++it) {
      InfoNode& node = *it;
      unsigned int zeroBasedIndex = !toZeroBasedIndex.empty() ? toZeroBasedIndex[node.physicalId] : (node.physicalId - minPhysicalId);
      node.physicalNodes.emplace_back(zeroBasedIndex, node.data.flow);
    }

    // If leaf nodes was not directly under root, make sure leaf modules have
    // physical nodes defined also
    if (depth > 1) {
      for (auto it(root.begin_leaf_modules()); !it.isEnd(); ++it) {
        InfoNode& module = *it;
        std::map<unsigned int, double> physToFlow;
        for (auto& node : module) {
          for (PhysData& physData : node.physicalNodes) {
            physToFlow[physData.physNodeIndex] += physData.sumFlowFromM2Node;
          }
        }
        for (auto& physFlow : physToFlow) {
          module.physicalNodes.emplace_back(physFlow.first, physFlow.second);
        }
      }
    }
  } else {
    // Either a sub-network (without modules) or the whole network with reconstructed tree
    if (depth == 1) {
      // new sub-network
      Log(3) << "RegularizedMultilayerMapEquation::initPhysicalNodesOnSubNetwork()...\n";
      std::set<unsigned int> setOfPhysicalNodes;
      unsigned int maxPhysNodeIndex = 0;
      unsigned int minPhysNodeIndex = std::numeric_limits<unsigned int>::max();

      // Collect all physical nodes in this sub network
      for (InfoNode& node : root) {
        for (PhysData& physData : node.physicalNodes) {
          setOfPhysicalNodes.insert(physData.physNodeIndex);
          maxPhysNodeIndex = std::max(maxPhysNodeIndex, physData.physNodeIndex);
          minPhysNodeIndex = std::min(minPhysNodeIndex, physData.physNodeIndex);
        }
      }

      m_numPhysicalNodes = setOfPhysicalNodes.size();

      // Re-index physical nodes if needed (not required when reconstructing tree)
      if (maxPhysNodeIndex >= m_numPhysicalNodes) {
        std::map<unsigned int, unsigned int> toZeroBasedIndex;
        if (maxPhysNodeIndex - minPhysNodeIndex + 1 > m_numPhysicalNodes) {
          unsigned int zeroBasedPhysicalId = 0;
          for (unsigned int physIndex : setOfPhysicalNodes) {
            toZeroBasedIndex.insert(std::make_pair(physIndex, zeroBasedPhysicalId++));
          }
        }

        for (InfoNode& node : root) {
          for (PhysData& physData : node.physicalNodes) {
            unsigned int zeroBasedIndex = !toZeroBasedIndex.empty() ? toZeroBasedIndex[physData.physNodeIndex] : (physData.physNodeIndex - minPhysNodeIndex);
            physData.physNodeIndex = zeroBasedIndex;
          }
        }
      }
    } else {
      // whole network with reconstructed tree
      for (auto it(root.begin_leaf_modules()); !it.isEnd(); ++it) {
        InfoNode& module = *it;
        std::map<unsigned int, double> physToFlow;
        for (auto& node : module) {
          for (PhysData& physData : node.physicalNodes) {
            physToFlow[physData.physNodeIndex] += physData.sumFlowFromM2Node;
          }
        }
        for (auto& physFlow : physToFlow) {
          module.physicalNodes.emplace_back(physFlow.first, physFlow.second);
        }
      }
    }
  }
}

void RegularizedMultilayerMapEquation::initPartitionOfPhysicalNodes(std::vector<InfoNode*>& nodes)
{
  m_physToModuleToMemNodes.clear();
  m_physToModuleToMemNodes.resize(m_numPhysicalNodes);

  for (auto& n : nodes) {
    InfoNode& node = *n;
    unsigned int moduleIndex = node.index; // Assume unique module index for all nodes in this initiation phase

    for (PhysData& physData : node.physicalNodes) {
      m_physToModuleToMemNodes[physData.physNodeIndex].emplace_back(moduleIndex, 1, physData.sumFlowFromM2Node);
    }
  }

  // Keep each physical node's entries sorted by module id; the move and
  // codelength code relies on this for lookups and summation order.
  for (auto& moduleToMemNodes : m_physToModuleToMemNodes) {
    std::sort(moduleToMemNodes.begin(), moduleToMemNodes.end(), [](const ModuleMemNodes& a, const ModuleMemNodes& b) { return a.module < b.module; });
  }

  m_memoryContributionsAdded = false;
}

void RegularizedMultilayerMapEquation::initPartitionLayerTeleFlowData(std::vector<InfoNode*>& nodes)
{
  m_moduleLayerTeleFlowData.clear();
  m_moduleLayerTeleFlowData.resize(nodes.size());

  for (auto& n : nodes) {
    InfoNode& node = *n;
    addLayerTeleFlow(node.index, node.layerTeleFlowData());
  }
}

// ===================================================
// Codelength
// ===================================================

void RegularizedMultilayerMapEquation::calculateCodelength(std::vector<InfoNode*>& nodes)
{
  calculateCodelengthTerms(nodes);

  calculateNodeFlow_log_nodeFlow();

  calculateCodelengthFromCodelengthTerms();
}

void RegularizedMultilayerMapEquation::calculateNodeFlow_log_nodeFlow()
{
  nodeFlow_log_nodeFlow = 0.0;
  for (unsigned int i = 0; i < m_numPhysicalNodes; ++i) {
    const ModuleToMemNodes& moduleToMemNodes = m_physToModuleToMemNodes[i];
    for (const auto& memNodes : moduleToMemNodes)
      nodeFlow_log_nodeFlow += infomath::plogp(memNodes.sumFlow);
  }
}

double RegularizedMultilayerMapEquation::calcCodelength(const InfoNode& parent) const
{
  if (parent.isLeafModule()) {
    return calcCodelengthOnModuleOfLeafNodes(parent);
  }
  // Use first-order model on index codebook
  return Base::calcCodelengthOnModuleOfModules(parent);
}

double RegularizedMultilayerMapEquation::calcCodelengthOnModuleOfLeafNodes(const InfoNode& parent) const
{
  if (parent.numPhysicalNodes() == 0) {
    return Base::calcCodelength(parent); // Infomap root node
  }

  // TODO: For ordinary networks, flow should be used instead of enter flow
  //  for leaf nodes, what about memory networks? sumFlowFromM2Node vs sumEnterFlowFromM2Node?
  double parentFlow = parent.data.flow;
  double parentExit = parent.data.exitFlow;
  double totalParentFlow = parentFlow + parentExit;
  if (totalParentFlow < 1e-16)
    return 0.0;

  double indexLength = 0.0;

  for (const PhysData& physData : parent.physicalNodes) {
    indexLength -= infomath::plogp(physData.sumFlowFromM2Node / totalParentFlow);
  }
  indexLength -= infomath::plogp(parentExit / totalParentFlow);

  indexLength *= totalParentFlow;

  return indexLength;
}

void RegularizedMultilayerMapEquation::addMemoryContributions(InfoNode& current,
                                                              DeltaFlowDataType& oldModuleDelta,
                                                              VectorMap<DeltaFlowDataType>& moduleDeltaFlow)
{
  // Overlapping modules
  /*
   * delta = old.first + new.first + old.second - new.second.
   * Two cases: (p(x) = plogp(x))
   * Moving to a module that already have that physical node: (old: p1, p2, new p3, moving p2 -> old:p1, new p2,p3)
   * Then old.second = new.second = plogp(physicalNodeSize) -> cancelation -> delta = p(p1) - p(p1+p2) + p(p2+p3) - p(p3)
   * Moving to a module that not have that physical node: (old: p1, p2, new -, moving p2 -> old: p1, new: p2)
   * Then new.first = new.second = 0 -> delta = p(p1) - p(p1+p2) + p(p2).
   */
  auto& physicalNodes = current.physicalNodes;
  unsigned int numPhysicalNodes = physicalNodes.size();
  for (unsigned int i = 0; i < numPhysicalNodes; ++i) {
    PhysData& physData = physicalNodes[i];
    ModuleToMemNodes& moduleToMemNodes = m_physToModuleToMemNodes[physData.physNodeIndex];
    for (const auto& memNodes : moduleToMemNodes) {
      unsigned int moduleIndex = memNodes.module;
      if (moduleIndex == current.index) // From where the multiple assigned node is moved
      {
        double oldPhysFlow = memNodes.sumFlow;
        double newPhysFlow = memNodes.sumFlow - physData.sumFlowFromM2Node;
        oldModuleDelta.sumDeltaPlogpPhysFlow += infomath::plogp(newPhysFlow) - infomath::plogp(oldPhysFlow);
        oldModuleDelta.sumPlogpPhysFlow += infomath::plogp(physData.sumFlowFromM2Node);
      } else // To where the multiple assigned node is moved
      {
        double oldPhysFlow = memNodes.sumFlow;
        double newPhysFlow = memNodes.sumFlow + physData.sumFlowFromM2Node;

        double sumDeltaPlogpPhysFlow = infomath::plogp(newPhysFlow) - infomath::plogp(oldPhysFlow);
        double sumPlogpPhysFlow = infomath::plogp(physData.sumFlowFromM2Node);
        moduleDeltaFlow.add(moduleIndex, DeltaFlowDataType(moduleIndex, 0.0, 0.0, sumDeltaPlogpPhysFlow, sumPlogpPhysFlow));
      }
    }
  }
  m_memoryContributionsAdded = true;
}

/**
 * Add teleportation flow for predefined move
 */
void RegularizedMultilayerMapEquation::addTeleportationFlow(InfoNode& current, const std::vector<FlowDataType>& moduleFlowData, DeltaFlowDataType& oldModuleDelta, DeltaFlowDataType& newModuleDelta)
{
  // TODO: Optimize: check if teleportFlow is zero and skip this (true for multilayer regularization now)
  Base::addTeleportationFlow(current, moduleFlowData, oldModuleDelta, newModuleDelta);

  auto& oldModuleLayerFlowData = m_moduleLayerTeleFlowData[oldModuleDelta.module];

  for (const auto& nodeLayerFlow : current.layerTeleFlowData()) {
    auto itModuleLayer = oldModuleLayerFlowData.find(nodeLayerFlow.layerId);
    if (itModuleLayer != oldModuleLayerFlowData.end()) {
      const auto& moduleLayerFlow = itModuleLayer->second;
      oldModuleDelta.deltaEnter += (moduleLayerFlow.teleportFlow - nodeLayerFlow.teleportFlow) * nodeLayerFlow.teleportWeight;
      oldModuleDelta.deltaExit += nodeLayerFlow.teleportFlow * (moduleLayerFlow.teleportWeight - nodeLayerFlow.teleportWeight);
    }
  }

  auto& newModuleLayerFlowData = m_moduleLayerTeleFlowData[newModuleDelta.module];

  for (const auto& nodeLayerFlow : current.layerTeleFlowData()) {
    auto itModuleLayer = newModuleLayerFlowData.find(nodeLayerFlow.layerId);
    if (itModuleLayer != newModuleLayerFlowData.end()) {
      const auto& moduleLayerFlow = itModuleLayer->second;
      newModuleDelta.deltaEnter += moduleLayerFlow.teleportFlow * nodeLayerFlow.teleportWeight;
      newModuleDelta.deltaExit += nodeLayerFlow.teleportFlow * moduleLayerFlow.teleportWeight;
    }
  }
}

void RegularizedMultilayerMapEquation::addTeleportationFlow(InfoNode& current, const std::vector<FlowDataType>& moduleFlowData, VectorMap<DeltaFlowDataType>& moduleDeltaFlow)
{
  Base::addTeleportationFlow(current, moduleFlowData, moduleDeltaFlow);

  const auto& moduleDeltaEnterExit = moduleDeltaFlow.values();

  for (unsigned int j = 0; j < moduleDeltaFlow.size(); ++j) {
    const auto& deltaEnterExit = moduleDeltaEnterExit[j];
    auto moduleIndex = deltaEnterExit.module;
    double deltaEnter = 0;
    double deltaExit = 0;

    if (moduleIndex == current.index) {
      auto& oldModuleLayerFlowData = m_moduleLayerTeleFlowData[moduleIndex];

      // A single node here can be multiple nodes moving (coarse-tune), so need to check all layers
      for (const auto& nodeLayerFlow : current.layerTeleFlowData()) {
        auto itModuleLayer = oldModuleLayerFlowData.find(nodeLayerFlow.layerId);
        if (itModuleLayer != oldModuleLayerFlowData.end()) {
          const auto& moduleLayerFlow = itModuleLayer->second;
          deltaEnter += (moduleLayerFlow.teleportFlow - nodeLayerFlow.teleportFlow) * nodeLayerFlow.teleportWeight;
          deltaExit += nodeLayerFlow.teleportFlow * (moduleLayerFlow.teleportWeight - nodeLayerFlow.teleportWeight);
        }
      }
    } else {
      auto& newModuleLayerFlowData = m_moduleLayerTeleFlowData[moduleIndex];

      for (const auto& nodeLayerFlow : current.layerTeleFlowData()) {
        auto itModuleLayer = newModuleLayerFlowData.find(nodeLayerFlow.layerId);
        if (itModuleLayer != newModuleLayerFlowData.end()) {
          const auto& moduleLayerFlow = itModuleLayer->second;
          deltaEnter += moduleLayerFlow.teleportFlow * nodeLayerFlow.teleportWeight;
          deltaExit += nodeLayerFlow.teleportFlow * moduleLayerFlow.teleportWeight;
        }
      }
    }
    moduleDeltaFlow.add(moduleIndex, DeltaFlowDataType(moduleIndex, deltaExit, deltaEnter));
  }
}

double RegularizedMultilayerMapEquation::getDeltaCodelengthOnMovingNode(InfoNode& current,
                                                                        DeltaFlowDataType& oldModuleDelta,
                                                                        DeltaFlowDataType& newModuleDelta,
                                                                        std::vector<FlowDataType>& moduleFlowData,
                                                                        std::vector<unsigned int>& moduleMembers)
{
  double deltaL = Base::getDeltaCodelengthOnMovingNode(current, oldModuleDelta, newModuleDelta, moduleFlowData, moduleMembers);

  double delta_nodeFlow_log_nodeFlow = oldModuleDelta.sumDeltaPlogpPhysFlow + newModuleDelta.sumDeltaPlogpPhysFlow + oldModuleDelta.sumPlogpPhysFlow - newModuleDelta.sumPlogpPhysFlow;

  return deltaL - delta_nodeFlow_log_nodeFlow;
}

double RegularizedMultilayerMapEquation::getDeltaCodelengthOnMovingNodeHoisted(InfoNode& current,
                                                                               DeltaFlowDataType& oldModuleDelta,
                                                                               const OldSideTerms& oldSide,
                                                                               DeltaFlowDataType& newModuleDelta,
                                                                               std::vector<FlowDataType>& moduleFlowData,
                                                                               std::vector<unsigned int>& moduleMembers)
{
  double deltaL = Base::getDeltaCodelengthOnMovingNodeHoisted(current, oldModuleDelta, oldSide, newModuleDelta, moduleFlowData, moduleMembers);

  double delta_nodeFlow_log_nodeFlow = oldModuleDelta.sumDeltaPlogpPhysFlow + newModuleDelta.sumDeltaPlogpPhysFlow + oldModuleDelta.sumPlogpPhysFlow - newModuleDelta.sumPlogpPhysFlow;

  return deltaL - delta_nodeFlow_log_nodeFlow;
}

// ===================================================
// Consolidation
// ===================================================

void RegularizedMultilayerMapEquation::updateCodelengthOnMovingNode(InfoNode& current,
                                                                    DeltaFlowDataType& oldModuleDelta,
                                                                    DeltaFlowDataType& newModuleDelta,
                                                                    std::vector<FlowDataType>& moduleFlowData,
                                                                    std::vector<unsigned int>& moduleMembers)
{
  Base::updateCodelengthOnMovingNode(current, oldModuleDelta, newModuleDelta, moduleFlowData, moduleMembers);
  if (m_memoryContributionsAdded)
    updatePhysicalNodes(current, oldModuleDelta.module, newModuleDelta.module);
  else
    addMemoryContributionsAndUpdatePhysicalNodes(current, oldModuleDelta, newModuleDelta);

  removeLayerTeleFlow(oldModuleDelta.module, current.layerTeleFlowData());
  addLayerTeleFlow(newModuleDelta.module, current.layerTeleFlowData());

  double delta_nodeFlow_log_nodeFlow = oldModuleDelta.sumDeltaPlogpPhysFlow + newModuleDelta.sumDeltaPlogpPhysFlow + oldModuleDelta.sumPlogpPhysFlow - newModuleDelta.sumPlogpPhysFlow;

  nodeFlow_log_nodeFlow += delta_nodeFlow_log_nodeFlow;
  moduleCodelength -= delta_nodeFlow_log_nodeFlow;
  codelength -= delta_nodeFlow_log_nodeFlow;

  // The flag marks that this move's memory contributions were already added to
  // the deltas by addMemoryContributions(). Consumed by this one move; reset so
  // a following standalone moveNodeToPredefinedModule() (fresh zero deltas)
  // recomputes them instead of leaving nodeFlow_log_nodeFlow stale. See the
  // matching note in MemMapEquation::updateCodelengthOnMovingNode.
  m_memoryContributionsAdded = false;
}

void RegularizedMultilayerMapEquation::addLayerTeleFlow(unsigned int moduleIndex, const std::vector<LayerTeleFlowData>& layerTeleFlowData)
{
  auto& moduleLayerTeleFlowData = m_moduleLayerTeleFlowData[moduleIndex];
  for (const auto& layerData : layerTeleFlowData) {
    auto& moduleLayerData = moduleLayerTeleFlowData[layerData.layerId];
    moduleLayerData.layerId = layerData.layerId;
    moduleLayerData += layerData;
  }
}

void RegularizedMultilayerMapEquation::removeLayerTeleFlow(unsigned int moduleIndex, const std::vector<LayerTeleFlowData>& layerTeleFlowData)
{
  auto& moduleLayerTeleFlowData = m_moduleLayerTeleFlowData[moduleIndex];
  for (const auto& layerData : layerTeleFlowData) {
    auto it = moduleLayerTeleFlowData.find(layerData.layerId);
    if (it == moduleLayerTeleFlowData.end())
      throw std::length_error("Couldn't find old module layer teleport data.");

    it->second -= layerData;
    if (it->second.isEmpty())
      moduleLayerTeleFlowData.erase(it);
  }
}

void RegularizedMultilayerMapEquation::updatePhysicalNodes(InfoNode& current, unsigned int oldModuleIndex, unsigned int bestModuleIndex)
{
  // For all multiple assigned nodes
  for (const auto& physData : current.physicalNodes) {
    ModuleToMemNodes& moduleToMemNodes = m_physToModuleToMemNodes[physData.physNodeIndex];

    // Remove contribution to old module
    auto overlapIt = findModuleMemNodes(moduleToMemNodes, oldModuleIndex);
    if (overlapIt == moduleToMemNodes.end() || overlapIt->module != oldModuleIndex)
      throw std::length_error(fmt::format(FMT_STRING("Couldn't find old module {} in physical node {}"), oldModuleIndex, physData.physNodeIndex));

    overlapIt->sumFlow -= physData.sumFlowFromM2Node;
    if (--overlapIt->numMemNodes == 0)
      moduleToMemNodes.erase(overlapIt);

    // Add contribution to new module
    overlapIt = findModuleMemNodes(moduleToMemNodes, bestModuleIndex);
    if (overlapIt == moduleToMemNodes.end() || overlapIt->module != bestModuleIndex) {
      moduleToMemNodes.insert(overlapIt, ModuleMemNodes(bestModuleIndex, 1, physData.sumFlowFromM2Node));
    } else {
      ++overlapIt->numMemNodes;
      overlapIt->sumFlow += physData.sumFlowFromM2Node;
    }
  }
}

void RegularizedMultilayerMapEquation::addMemoryContributionsAndUpdatePhysicalNodes(InfoNode& current, DeltaFlowDataType& oldModuleDelta, DeltaFlowDataType& newModuleDelta)
{
  unsigned int oldModuleIndex = oldModuleDelta.module;
  unsigned int bestModuleIndex = newModuleDelta.module;

  // For all multiple assigned nodes
  for (const auto& physData : current.physicalNodes) {
    ModuleToMemNodes& moduleToMemNodes = m_physToModuleToMemNodes[physData.physNodeIndex];

    // Remove contribution to old module
    auto overlapIt = findModuleMemNodes(moduleToMemNodes, oldModuleIndex);
    if (overlapIt == moduleToMemNodes.end() || overlapIt->module != oldModuleIndex)
      throw std::length_error("Couldn't find old module among physical node assignments.");

    double oldPhysFlow = overlapIt->sumFlow;
    double newPhysFlow = overlapIt->sumFlow - physData.sumFlowFromM2Node;
    oldModuleDelta.sumDeltaPlogpPhysFlow += infomath::plogp(newPhysFlow) - infomath::plogp(oldPhysFlow);
    oldModuleDelta.sumPlogpPhysFlow += infomath::plogp(physData.sumFlowFromM2Node);
    overlapIt->sumFlow -= physData.sumFlowFromM2Node;
    if (--overlapIt->numMemNodes == 0)
      moduleToMemNodes.erase(overlapIt);

    // Add contribution to new module
    overlapIt = findModuleMemNodes(moduleToMemNodes, bestModuleIndex);
    if (overlapIt == moduleToMemNodes.end() || overlapIt->module != bestModuleIndex) {
      moduleToMemNodes.insert(overlapIt, ModuleMemNodes(bestModuleIndex, 1, physData.sumFlowFromM2Node));
      oldPhysFlow = 0.0;
      newPhysFlow = physData.sumFlowFromM2Node;
      newModuleDelta.sumDeltaPlogpPhysFlow += infomath::plogp(newPhysFlow) - infomath::plogp(oldPhysFlow);
      newModuleDelta.sumPlogpPhysFlow += infomath::plogp(physData.sumFlowFromM2Node);
    } else {
      oldPhysFlow = overlapIt->sumFlow;
      newPhysFlow = overlapIt->sumFlow + physData.sumFlowFromM2Node;
      newModuleDelta.sumDeltaPlogpPhysFlow += infomath::plogp(newPhysFlow) - infomath::plogp(oldPhysFlow);
      newModuleDelta.sumPlogpPhysFlow += infomath::plogp(physData.sumFlowFromM2Node);
      ++overlapIt->numMemNodes;
      overlapIt->sumFlow += physData.sumFlowFromM2Node;
    }
  }
}

void RegularizedMultilayerMapEquation::consolidateModules(std::vector<InfoNode*>& modules)
{
  for (unsigned int i = 0; i < m_numPhysicalNodes; ++i) {
    ModuleToMemNodes& modToMemNodes = m_physToModuleToMemNodes[i];
    for (unsigned int j = 0; j < modToMemNodes.size(); ++j) {
      // Sorted strictly increasing module ids guarantee no (module, phys) duplicates
      if (j > 0 && modToMemNodes[j].module <= modToMemNodes[j - 1].module)
        throw std::domain_error("[InfomapGreedy::consolidateModules] Error updating physical nodes: duplication error");

      modules[modToMemNodes[j].module]->physicalNodes.emplace_back(i, modToMemNodes[j].sumFlow);
    }
  }

  for (unsigned int moduleIndex = 0; moduleIndex < modules.size(); ++moduleIndex) {
    if (modules[moduleIndex] == nullptr)
      continue;
    auto& moduleLayerTele = modules[moduleIndex]->ensureExtras().layerTeleFlowData;
    moduleLayerTele.clear();
    for (const auto& layerIt : m_moduleLayerTeleFlowData[moduleIndex]) {
      moduleLayerTele.push_back(layerIt.second);
    }
  }
}

// ===================================================
// Debug
// ===================================================

void RegularizedMultilayerMapEquation::printDebug() const
{
  Log().print("RegularizedMultilayerMapEquation::m_numPhysicalNodes: {}\n", m_numPhysicalNodes);
  Base::printDebug();
}

} // namespace infomap
