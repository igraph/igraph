/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "MemMapEquation.h"
#include "FlowData.h"
#include "InfoNode.h"
#include "../utils/Log.h"
#include "../utils/format.h"

#include <vector>
#include <set>
#include <map>
#include <limits>
#include <utility>
#include <algorithm>
#include <unordered_map>

namespace infomap {

// ===================================================
// IO
// ===================================================

std::ostream& MemMapEquation::print(std::ostream& out) const
{
  return out << indexCodelength << " + " << moduleCodelength << " = " << io::toPrecision(codelength);
}

std::ostream& operator<<(std::ostream& out, const MemMapEquation& mapEq)
{
  return mapEq.print(out);
}

// ===================================================
// Init
// ===================================================

void MemMapEquation::init(const Config& config)
{
  Log(3) << "MemMapEquation::init()...\n";
  Base::init(config);
}

void MemMapEquation::initNetwork(InfoNode& root)
{
  initPhysicalNodes(root);
}

void MemMapEquation::initSuperNetwork(InfoNode& /*root*/)
{
  // TODO: How use enterFlow instead of flow
}

void MemMapEquation::initSubNetwork(InfoNode& /*root*/)
{
}

void MemMapEquation::initPartition(std::vector<InfoNode*>& nodes)
{
  initPartitionOfPhysicalNodes(nodes);

  calculateCodelength(nodes);
}

void MemMapEquation::initPhysicalNodes(InfoNode& root)
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
    Log(3) << "MemMapEquation::initPhysicalNodesOnOriginalNetwork()...\n";
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
      Log(3) << "MemMapEquation::initPhysicalNodesOnSubNetwork()...\n";
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

void MemMapEquation::initPartitionOfPhysicalNodes(std::vector<InfoNode*>& nodes)
{
  Log(4) << "MemMapEquation::initPartitionOfPhysicalNodes()...\n";
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

// ===================================================
// Codelength
// ===================================================

void MemMapEquation::calculateCodelength(std::vector<InfoNode*>& nodes)
{
  calculateCodelengthTerms(nodes);

  calculateNodeFlow_log_nodeFlow();

  calculateCodelengthFromCodelengthTerms();
}

void MemMapEquation::calculateNodeFlow_log_nodeFlow()
{
  nodeFlow_log_nodeFlow = 0.0;
  for (unsigned int i = 0; i < m_numPhysicalNodes; ++i) {
    const ModuleToMemNodes& moduleToMemNodes = m_physToModuleToMemNodes[i];
    for (const auto& memNodes : moduleToMemNodes)
      nodeFlow_log_nodeFlow += infomath::plogp(memNodes.sumFlow);
  }
}

double MemMapEquation::calcCodelength(const InfoNode& parent) const
{
  if (parent.isLeafModule()) {
    return calcCodelengthOnModuleOfLeafNodes(parent);
  }
  // Use first-order model on index codebook
  return Base::calcCodelengthOnModuleOfModules(parent);
}

double MemMapEquation::calcCodelengthOnModuleOfLeafNodes(const InfoNode& parent) const
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

void MemMapEquation::addMemoryContributions(InfoNode& current,
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
    // plogp(sumFlowFromM2Node) is constant across every module holding this
    // physical node, so compute it once per physical node instead of per
    // (physical node, module). Bit-exact; the per-module plogp of the changed
    // flows (old/new) still varies and stays in the loop.
    const double physFlow = physData.sumFlowFromM2Node;
    const double plogpPhysFlow = infomath::plogp(physFlow);
    ModuleToMemNodes& moduleToMemNodes = m_physToModuleToMemNodes[physData.physNodeIndex];
    for (const auto& memNodes : moduleToMemNodes) {
      unsigned int moduleIndex = memNodes.module;
      if (moduleIndex == current.index) // From where the multiple assigned node is moved
      {
        double oldPhysFlow = memNodes.sumFlow;
        double newPhysFlow = memNodes.sumFlow - physFlow;
        oldModuleDelta.sumDeltaPlogpPhysFlow += infomath::plogp(newPhysFlow) - infomath::plogp(oldPhysFlow);
        oldModuleDelta.sumPlogpPhysFlow += plogpPhysFlow;
      } else // To where the multiple assigned node is moved
      {
        double oldPhysFlow = memNodes.sumFlow;
        double newPhysFlow = memNodes.sumFlow + physFlow;

        double sumDeltaPlogpPhysFlow = infomath::plogp(newPhysFlow) - infomath::plogp(oldPhysFlow);
        moduleDeltaFlow.add(moduleIndex, DeltaFlowDataType(moduleIndex, 0.0, 0.0, sumDeltaPlogpPhysFlow, plogpPhysFlow));
      }
    }
  }
  m_memoryContributionsAdded = true;
}

INFOMAP_HOT double MemMapEquation::getDeltaCodelengthOnMovingNode(InfoNode& current,
                                                                  DeltaFlowDataType& oldModuleDelta,
                                                                  DeltaFlowDataType& newModuleDelta,
                                                                  std::vector<FlowDataType>& moduleFlowData,
                                                                  std::vector<unsigned int>& moduleMembers)
{
  double deltaL = Base::getDeltaCodelengthOnMovingNode(current, oldModuleDelta, newModuleDelta, moduleFlowData, moduleMembers);

  double delta_nodeFlow_log_nodeFlow = oldModuleDelta.sumDeltaPlogpPhysFlow + newModuleDelta.sumDeltaPlogpPhysFlow + oldModuleDelta.sumPlogpPhysFlow - newModuleDelta.sumPlogpPhysFlow;

  return deltaL - delta_nodeFlow_log_nodeFlow;
}

INFOMAP_HOT double MemMapEquation::getDeltaCodelengthOnMovingNodeHoisted(InfoNode& current,
                                                                         DeltaFlowDataType& oldModuleDelta,
                                                                         const OldSideTerms& oldSide,
                                                                         DeltaFlowDataType& newModuleDelta,
                                                                         std::vector<FlowDataType>& moduleFlowData,
                                                                         std::vector<unsigned int>& moduleMembers)
{
  double deltaL = Base::getDeltaCodelengthOnMovingNodeHoisted(current, oldModuleDelta, oldSide, newModuleDelta, moduleFlowData, moduleMembers);

  // The physical-codebook old-side terms are already per-node constants:
  // addMemoryContributions() computed them once into oldModuleDelta.
  double delta_nodeFlow_log_nodeFlow = oldModuleDelta.sumDeltaPlogpPhysFlow + newModuleDelta.sumDeltaPlogpPhysFlow + oldModuleDelta.sumPlogpPhysFlow - newModuleDelta.sumPlogpPhysFlow;

  return deltaL - delta_nodeFlow_log_nodeFlow;
}

// ===================================================
// Consolidation
// ===================================================

void MemMapEquation::updateCodelengthOnMovingNode(InfoNode& current,
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

  double delta_nodeFlow_log_nodeFlow = oldModuleDelta.sumDeltaPlogpPhysFlow + newModuleDelta.sumDeltaPlogpPhysFlow + oldModuleDelta.sumPlogpPhysFlow - newModuleDelta.sumPlogpPhysFlow;

  nodeFlow_log_nodeFlow += delta_nodeFlow_log_nodeFlow;
  moduleCodelength -= delta_nodeFlow_log_nodeFlow;
  codelength -= delta_nodeFlow_log_nodeFlow;

  // The flag marks that the memory contributions for THIS move were already
  // accumulated onto the deltas by addMemoryContributions() (the optimize
  // candidate-evaluation path). It is consumed by this one move; reset it so a
  // following standalone moveNodeToPredefinedModule() — e.g. the single-node
  // "drag along" in tryMoveEachNodeIntoBestModule(), which passes fresh
  // zero-initialised deltas — recomputes its own contributions instead of
  // updating the physical-node map while leaving nodeFlow_log_nodeFlow stale.
  m_memoryContributionsAdded = false;
}

void MemMapEquation::updatePhysicalNodes(InfoNode& current, unsigned int oldModuleIndex, unsigned int bestModuleIndex)
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

void MemMapEquation::addMemoryContributionsAndUpdatePhysicalNodes(InfoNode& current, DeltaFlowDataType& oldModuleDelta, DeltaFlowDataType& newModuleDelta)
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

void MemMapEquation::consolidateModules(std::vector<InfoNode*>& modules)
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
}

// ===================================================
// Debug
// ===================================================

void MemMapEquation::printDebug() const
{
  Log().print("MemMapEquation::m_numPhysicalNodes: {}\n", m_numPhysicalNodes);
  Base::printDebug();
}

} // namespace infomap
