/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef MAPEQUATION_H_
#define MAPEQUATION_H_

#include "../utils/infomath.h"
#include "../utils/convert.h"
#include "../io/Config.h"
#include "../utils/Log.h"
#include "../utils/VectorMap.h"
#include "InfoNode.h"
#include "FlowData.h"
#include <vector>
#include <map>
#include <ostream>
#include <algorithm>

namespace infomap {

class InfoNode;

// One module entry of a physical node's flat module map, kept sorted by
// module id. The sort order is load-bearing: iteration order determines
// floating-point summation order in the codelength terms.
struct ModuleMemNodes {
  ModuleMemNodes(unsigned int module, unsigned int numMemNodes, double sumFlow)
      : module(module), numMemNodes(numMemNodes), sumFlow(sumFlow) {}
  unsigned int module;
  unsigned int numMemNodes; // use counter to check for zero to avoid round-off errors in sumFlow
  double sumFlow;
};

// Old-module plogp terms of getDeltaCodelengthOnMovingNode. For a fixed node
// visit the old module and the node's own flow are constant, so 6 of the 13
// plogp values the delta needs are identical for every candidate module the
// move loop probes for that node. hoistOldSide() computes them once per node;
// getDeltaCodelengthOnMovingNodeHoisted() consumes them per candidate, keeping
// the exact same expression structure as the unhoisted form. Plain doubles and
// namespace scope so the (privately-inheriting) objectives can pass it around
// by `auto` without accessibility juggling. Not bit-identical to the unhoisted
// path at the last ulp: splitting the single 13-wide plogp_batch into a 6-wide
// (once/node) and 7-wide (per candidate) batch regroups the SIMD lanes.
struct OldSideTerms {
  double deltaEnterExitOld; // oldModuleDelta.deltaEnter + deltaExit (with teleportation folded in)
  double curEnter, curExit, curFlow; // current.data, cached to keep the hoisted delta self-contained
  double plogpOldEnter, plogpOldEnterAfter;
  double plogpOldExit, plogpOldExitAfter;
  double plogpOldExitFlow, plogpOldExitFlowAfter;
};

/**
 * Base implementation of the map equation, shared by the concrete objectives
 * (BiasedMapEquation, MemMapEquation, MetaMapEquation, RegularizedMultilayerMapEquation).
 *
 * Dispatch is static, not dynamic: each objective is used as the concrete
 * `Objective` type parameter of InfomapOptimizer<Objective> and stored by value,
 * so no call here goes through a base pointer/reference. Subclasses inherit
 * privately and re-expose or redefine members by name; the methods below are
 * therefore deliberately non-virtual (a redefinition in a subclass hides the
 * base version). Internal cross-calls are qualified with `ME::` so they always
 * bind to this base implementation regardless of the most-derived type.
 */

inline std::vector<ModuleMemNodes>::iterator findModuleMemNodes(std::vector<ModuleMemNodes>& moduleToMemNodes, unsigned int module)
{
  return std::lower_bound(moduleToMemNodes.begin(), moduleToMemNodes.end(), module, [](const ModuleMemNodes& memNodes, unsigned int target) { return memNodes.module < target; });
}

template <typename FlowDataType = FlowData, typename DeltaFlowDataType = DeltaFlow>
class MapEquation {
  using ME = MapEquation<FlowDataType, DeltaFlowDataType>;

public:
  MapEquation() = default;

  MapEquation(const MapEquation& other) = default;

  MapEquation& operator=(const MapEquation& other) = default;

  MapEquation(MapEquation&& other) noexcept = default;

  MapEquation& operator=(MapEquation&& other) noexcept = default;

  ~MapEquation() = default;

  // ===================================================
  // Getters
  // ===================================================

  double getIndexCodelength() const { return indexCodelength; }

  double getModuleCodelength() const { return moduleCodelength; }

  double getCodelength() const { return codelength; }

  // ===================================================
  // IO
  // ===================================================

  std::ostream& print(std::ostream& out) const
  {
    return out << indexCodelength << " + " << moduleCodelength << " = " << io::toPrecision(codelength);
  }

  // ===================================================
  // Init
  // ===================================================

  void init(const Config& config)
  {
    Log(3) << "MapEquation::init()...\n";
    m_config = config;
  }

  void initNetwork(InfoNode& root)
  {
    Log(3) << "MapEquation::initNetwork()...\n";

    nodeFlow_log_nodeFlow = 0.0;
    for (InfoNode& node : root) {
      nodeFlow_log_nodeFlow += infomath::plogp(node.data.flow);
    }
    ME::initSubNetwork(root);
  }

  void initSuperNetwork(InfoNode& root)
  {
    Log(3) << "MapEquation::initSuperNetwork()...\n";

    nodeFlow_log_nodeFlow = 0.0;
    for (InfoNode& node : root) {
      nodeFlow_log_nodeFlow += infomath::plogp(node.data.enterFlow);
    }
  }

  void initSubNetwork(InfoNode& root)
  {
    exitNetworkFlow = root.data.exitFlow;
    exitNetworkFlow_log_exitNetworkFlow = infomath::plogp(exitNetworkFlow);
  }

  void initPartition(std::vector<InfoNode*>& nodes) { ME::calculateCodelength(nodes); }

  // ===================================================
  // Codelength
  // ===================================================

  double calcCodelength(const InfoNode& parent) const
  {
    return parent.isLeafModule() ? ME::calcCodelengthOnModuleOfLeafNodes(parent) : ME::calcCodelengthOnModuleOfModules(parent);
  }

  void addMemoryContributions(InfoNode& /*current*/, DeltaFlowDataType& /*oldModuleDelta*/, DeltaFlowDataType& /*newModuleDelta*/) {}

  void addMemoryContributions(InfoNode& /*current*/, DeltaFlowDataType& /*oldModuleDelta*/, VectorMap<DeltaFlowDataType>& /*moduleDeltaFlow*/) {}

  void addTeleportationFlow(InfoNode& current, const std::vector<FlowDataType>& moduleFlowData, DeltaFlowDataType& oldModuleDelta, DeltaFlowDataType& newModuleDelta);

  void addTeleportationFlow(InfoNode& current, const std::vector<FlowDataType>& moduleFlowData, VectorMap<DeltaFlowDataType>& moduleDeltaFlow);

  double getDeltaCodelengthOnMovingNode(InfoNode& current,
                                        DeltaFlowDataType& oldModuleDelta,
                                        DeltaFlowDataType& newModuleDelta,
                                        std::vector<FlowDataType>& moduleFlowData,
                                        std::vector<unsigned int>& /*moduleMembers*/);

  // Precompute the old-module (per-node-constant) plogp terms of the move delta.
  // Called once per node visit by the optimize move loop; the result feeds every
  // candidate through getDeltaCodelengthOnMovingNodeHoisted().
  OldSideTerms hoistOldSide(InfoNode& current,
                            DeltaFlowDataType& oldModuleDelta,
                            std::vector<FlowDataType>& moduleFlowData);

  double getDeltaCodelengthOnMovingNodeHoisted(InfoNode& current,
                                               DeltaFlowDataType& oldModuleDelta,
                                               const OldSideTerms& oldSide,
                                               DeltaFlowDataType& newModuleDelta,
                                               std::vector<FlowDataType>& moduleFlowData,
                                               std::vector<unsigned int>& /*moduleMembers*/);

  // ===================================================
  // Consolidation
  // ===================================================

  void updateCodelengthOnMovingNode(InfoNode& current,
                                    DeltaFlowDataType& oldModuleDelta,
                                    DeltaFlowDataType& newModuleDelta,
                                    std::vector<FlowDataType>& moduleFlowData,
                                    std::vector<unsigned int>& /*moduleMembers*/);

  // ===================================================
  // Debug
  // ===================================================

  void printDebug() const
  {
    Log().print("(enterFlow_log_enterFlow: {:g}, enter_log_enter: {:g}, exitNetworkFlow_log_exitNetworkFlow: {:g}) ",
                enterFlow_log_enterFlow,
                enter_log_enter,
                exitNetworkFlow_log_exitNetworkFlow);
  }

protected:
  // ===================================================
  // Protected member functions
  // ===================================================

  double calcCodelengthOnModuleOfLeafNodes(const InfoNode& parent) const;

  double calcCodelengthOnModuleOfModules(const InfoNode& parent) const;

  void calculateCodelength(std::vector<InfoNode*>& nodes)
  {
    ME::calculateCodelengthTerms(nodes);
    ME::calculateCodelengthFromCodelengthTerms();
  }

  void calculateCodelengthTerms(std::vector<InfoNode*>& nodes);

  void calculateCodelengthFromCodelengthTerms()
  {
    indexCodelength = enterFlow_log_enterFlow - enter_log_enter - exitNetworkFlow_log_exitNetworkFlow;
    moduleCodelength = -exit_log_exit + flow_log_flow - nodeFlow_log_nodeFlow;
    codelength = indexCodelength + moduleCodelength;
  }

public:
  // ===================================================
  // Public member variables
  // ===================================================

  double codelength = 0.0;
  double indexCodelength = 0.0;
  double moduleCodelength = 0.0;

protected:
  // ===================================================
  // Protected member variables
  // ===================================================

  Config m_config;

  double nodeFlow_log_nodeFlow = 0.0; // constant while the leaf network is the same
  double flow_log_flow = 0.0; // node.(flow + exitFlow)
  double exit_log_exit = 0.0;
  double enter_log_enter = 0.0;
  double enterFlow = 0.0;
  double enterFlow_log_enterFlow = 0.0;

  // For hierarchical
  double exitNetworkFlow = 0.0;
  double exitNetworkFlow_log_exitNetworkFlow = 0.0;
};

/**
 * Add teleportation flow for predefined move
 */
template <typename FlowDataType, typename DeltaFlowDataType>
void MapEquation<FlowDataType, DeltaFlowDataType>::addTeleportationFlow(InfoNode& current, const std::vector<FlowDataType>& moduleFlowData, DeltaFlowDataType& oldModuleDelta, DeltaFlowDataType& newModuleDelta)
{
  if (!m_config.recordedTeleportation)
    return;

  auto& oldModuleFlowData = moduleFlowData[oldModuleDelta.module];
  double deltaEnterOld = (oldModuleFlowData.teleportFlow - current.data.teleportFlow) * current.data.teleportWeight;
  double deltaExitOld = current.data.teleportFlow * (oldModuleFlowData.teleportWeight - current.data.teleportWeight);
  oldModuleDelta.deltaEnter += deltaEnterOld;
  oldModuleDelta.deltaExit += deltaExitOld;

  auto& newModuleFlowData = moduleFlowData[newModuleDelta.module];
  double deltaEnterNew = current.data.teleportFlow * newModuleFlowData.teleportWeight;
  double deltaExitNew = newModuleFlowData.teleportFlow * current.data.teleportWeight;
  newModuleDelta.deltaEnter += deltaEnterNew;
  newModuleDelta.deltaExit += deltaExitNew;
}

template <typename FlowDataType, typename DeltaFlowDataType>
void MapEquation<FlowDataType, DeltaFlowDataType>::addTeleportationFlow(InfoNode& current, const std::vector<FlowDataType>& moduleFlowData, VectorMap<DeltaFlowDataType>& moduleDeltaFlow)
{
  if (!m_config.recordedTeleportation)
    return;

  auto& moduleDeltaEnterExit = moduleDeltaFlow.values();

  for (unsigned int j = 0; j < moduleDeltaFlow.size(); ++j) {
    auto& deltaEnterExit = moduleDeltaEnterExit[j];
    auto moduleIndex = deltaEnterExit.module;

    if (moduleIndex == current.index) {
      auto& oldModuleFlowData = moduleFlowData[moduleIndex];
      double deltaEnterOld = (oldModuleFlowData.teleportFlow - current.data.teleportFlow) * current.data.teleportWeight;
      double deltaExitOld = current.data.teleportFlow * (oldModuleFlowData.teleportWeight - current.data.teleportWeight);
      moduleDeltaFlow.add(moduleIndex, DeltaFlowDataType(moduleIndex, deltaExitOld, deltaEnterOld));
    } else {
      auto& newModuleFlowData = moduleFlowData[moduleIndex];
      double deltaEnterNew = newModuleFlowData.teleportFlow * current.data.teleportWeight;
      double deltaExitNew = current.data.teleportFlow * newModuleFlowData.teleportWeight;
      moduleDeltaFlow.add(moduleIndex, DeltaFlowDataType(moduleIndex, deltaExitNew, deltaEnterNew));
    }
  }
}

template <typename FlowDataType, typename DeltaFlowDataType>
INFOMAP_HOT double MapEquation<FlowDataType, DeltaFlowDataType>::getDeltaCodelengthOnMovingNode(InfoNode& current, DeltaFlowDataType& oldModuleDelta, DeltaFlowDataType& newModuleDelta, std::vector<FlowDataType>& moduleFlowData, std::vector<unsigned int>&)
{
  using infomath::plogp_batch;
  unsigned int oldModule = oldModuleDelta.module;
  unsigned int newModule = newModuleDelta.module;
  double deltaEnterExitOldModule = oldModuleDelta.deltaEnter + oldModuleDelta.deltaExit;
  double deltaEnterExitNewModule = newModuleDelta.deltaEnter + newModuleDelta.deltaExit;

  FlowDataType& oldMfd = moduleFlowData[oldModule];
  FlowDataType& newMfd = moduleFlowData[newModule];
  double oldEnter = oldMfd.enterFlow;
  double newEnter = newMfd.enterFlow;
  double oldExit = oldMfd.exitFlow;
  double newExit = newMfd.exitFlow;
  double oldExitPlusFlow = oldMfd.exitFlow + oldMfd.flow;
  double newExitPlusFlow = newMfd.exitFlow + newMfd.flow;
  double curEnter = current.data.enterFlow;
  double curExit = current.data.exitFlow;
  double curFlow = current.data.flow;

  constexpr int kPlogpBatchN = 13;
  double args[kPlogpBatchN] = {
    enterFlow + deltaEnterExitOldModule - deltaEnterExitNewModule,
    oldEnter,
    newEnter,
    oldEnter - curEnter + deltaEnterExitOldModule,
    newEnter + curEnter - deltaEnterExitNewModule,
    oldExit,
    newExit,
    oldExit - curExit + deltaEnterExitOldModule,
    newExit + curExit - deltaEnterExitNewModule,
    oldExitPlusFlow,
    newExitPlusFlow,
    oldExitPlusFlow - curExit - curFlow + deltaEnterExitOldModule,
    newExitPlusFlow + curExit + curFlow - deltaEnterExitNewModule,
  };
  double pl[kPlogpBatchN];
  plogp_batch(args, pl, kPlogpBatchN);

  double delta_enter = pl[0] - enterFlow_log_enterFlow;
  double delta_enter_log_enter = -pl[1] - pl[2] + pl[3] + pl[4];
  double delta_exit_log_exit = -pl[5] - pl[6] + pl[7] + pl[8];
  double delta_flow_log_flow = -pl[9] - pl[10] + pl[11] + pl[12];

  double deltaL = delta_enter - delta_enter_log_enter - delta_exit_log_exit + delta_flow_log_flow;
  return deltaL;
}

template <typename FlowDataType, typename DeltaFlowDataType>
INFOMAP_HOT OldSideTerms MapEquation<FlowDataType, DeltaFlowDataType>::hoistOldSide(InfoNode& current, DeltaFlowDataType& oldModuleDelta, std::vector<FlowDataType>& moduleFlowData)
{
  using infomath::plogp_batch;
  unsigned int oldModule = oldModuleDelta.module;
  double deltaEnterExitOldModule = oldModuleDelta.deltaEnter + oldModuleDelta.deltaExit;

  FlowDataType& oldMfd = moduleFlowData[oldModule];
  double oldEnter = oldMfd.enterFlow;
  double oldExit = oldMfd.exitFlow;
  double oldExitPlusFlow = oldMfd.exitFlow + oldMfd.flow;
  double curEnter = current.data.enterFlow;
  double curExit = current.data.exitFlow;
  double curFlow = current.data.flow;

  // The 6 args are exactly the 6 old-module entries of the 13-wide batch in
  // getDeltaCodelengthOnMovingNode (indices 1,3,5,7,9,11 there).
  double args[6] = {
    oldEnter,
    oldEnter - curEnter + deltaEnterExitOldModule,
    oldExit,
    oldExit - curExit + deltaEnterExitOldModule,
    oldExitPlusFlow,
    oldExitPlusFlow - curExit - curFlow + deltaEnterExitOldModule,
  };
  double pl[6];
  plogp_batch(args, pl, 6);

  return OldSideTerms { deltaEnterExitOldModule, curEnter, curExit, curFlow, pl[0], pl[1], pl[2], pl[3], pl[4], pl[5] };
}

template <typename FlowDataType, typename DeltaFlowDataType>
INFOMAP_HOT double MapEquation<FlowDataType, DeltaFlowDataType>::getDeltaCodelengthOnMovingNodeHoisted(InfoNode& /*current*/, DeltaFlowDataType& /*oldModuleDelta*/, const OldSideTerms& oldSide, DeltaFlowDataType& newModuleDelta, std::vector<FlowDataType>& moduleFlowData, std::vector<unsigned int>&)
{
  using infomath::plogp_batch;
  unsigned int newModule = newModuleDelta.module;
  double deltaEnterExitNewModule = newModuleDelta.deltaEnter + newModuleDelta.deltaExit;
  double deltaEnterExitOldModule = oldSide.deltaEnterExitOld;

  FlowDataType& newMfd = moduleFlowData[newModule];
  double newEnter = newMfd.enterFlow;
  double newExit = newMfd.exitFlow;
  double newExitPlusFlow = newMfd.exitFlow + newMfd.flow;
  double curEnter = oldSide.curEnter;
  double curExit = oldSide.curExit;
  double curFlow = oldSide.curFlow;

  // The 7 candidate-dependent entries of the 13-wide batch (indices 0,2,4,6,8,10,12).
  constexpr int kPlogpBatchN = 7;
  double args[kPlogpBatchN] = {
    enterFlow + deltaEnterExitOldModule - deltaEnterExitNewModule,
    newEnter,
    newEnter + curEnter - deltaEnterExitNewModule,
    newExit,
    newExit + curExit - deltaEnterExitNewModule,
    newExitPlusFlow,
    newExitPlusFlow + curExit + curFlow - deltaEnterExitNewModule,
  };
  double pl[kPlogpBatchN];
  plogp_batch(args, pl, kPlogpBatchN);

  double delta_enter = pl[0] - enterFlow_log_enterFlow;
  double delta_enter_log_enter = -oldSide.plogpOldEnter - pl[1] + oldSide.plogpOldEnterAfter + pl[2];
  double delta_exit_log_exit = -oldSide.plogpOldExit - pl[3] + oldSide.plogpOldExitAfter + pl[4];
  double delta_flow_log_flow = -oldSide.plogpOldExitFlow - pl[5] + oldSide.plogpOldExitFlowAfter + pl[6];

  double deltaL = delta_enter - delta_enter_log_enter - delta_exit_log_exit + delta_flow_log_flow;
  return deltaL;
}

template <typename FlowDataType, typename DeltaFlowDataType>
void MapEquation<FlowDataType, DeltaFlowDataType>::updateCodelengthOnMovingNode(InfoNode& current, DeltaFlowDataType& oldModuleDelta, DeltaFlowDataType& newModuleDelta, std::vector<FlowDataType>& moduleFlowData, std::vector<unsigned int>&)
{
  using infomath::plogp;
  unsigned int oldModule = oldModuleDelta.module;
  unsigned int newModule = newModuleDelta.module;
  double deltaEnterExitOldModule = oldModuleDelta.deltaEnter + oldModuleDelta.deltaExit;
  double deltaEnterExitNewModule = newModuleDelta.deltaEnter + newModuleDelta.deltaExit;

  enterFlow -= moduleFlowData[oldModule].enterFlow + moduleFlowData[newModule].enterFlow;
  enter_log_enter -= plogp(moduleFlowData[oldModule].enterFlow) + plogp(moduleFlowData[newModule].enterFlow);
  exit_log_exit -= plogp(moduleFlowData[oldModule].exitFlow) + plogp(moduleFlowData[newModule].exitFlow);
  flow_log_flow -= plogp(moduleFlowData[oldModule].exitFlow + moduleFlowData[oldModule].flow) + plogp(moduleFlowData[newModule].exitFlow + moduleFlowData[newModule].flow);

  moduleFlowData[oldModule] -= current.data;
  moduleFlowData[newModule] += current.data;

  moduleFlowData[oldModule].enterFlow += deltaEnterExitOldModule;
  moduleFlowData[oldModule].exitFlow += deltaEnterExitOldModule;
  moduleFlowData[newModule].enterFlow -= deltaEnterExitNewModule;
  moduleFlowData[newModule].exitFlow -= deltaEnterExitNewModule;

  enterFlow += moduleFlowData[oldModule].enterFlow + moduleFlowData[newModule].enterFlow;
  enter_log_enter += plogp(moduleFlowData[oldModule].enterFlow) + plogp(moduleFlowData[newModule].enterFlow);
  exit_log_exit += plogp(moduleFlowData[oldModule].exitFlow) + plogp(moduleFlowData[newModule].exitFlow);
  flow_log_flow += plogp(moduleFlowData[oldModule].exitFlow + moduleFlowData[oldModule].flow) + plogp(moduleFlowData[newModule].exitFlow + moduleFlowData[newModule].flow);

  enterFlow_log_enterFlow = plogp(enterFlow);

  indexCodelength = enterFlow_log_enterFlow - enter_log_enter - exitNetworkFlow_log_exitNetworkFlow;
  moduleCodelength = -exit_log_exit + flow_log_flow - nodeFlow_log_nodeFlow;
  codelength = indexCodelength + moduleCodelength;
}

template <typename FlowDataType, typename DeltaFlowDataType>
double MapEquation<FlowDataType, DeltaFlowDataType>::calcCodelengthOnModuleOfLeafNodes(const InfoNode& parent) const
{
  double parentFlow = parent.data.flow;
  double parentExit = parent.data.exitFlow;
  double totalParentFlow = parentFlow + parentExit;
  if (totalParentFlow < 1e-16)
    return 0.0;

  double indexLength = 0.0;
  for (const auto& node : parent) {
    indexLength -= infomath::plogp(node.data.flow / totalParentFlow);
  }
  indexLength -= infomath::plogp(parentExit / totalParentFlow);

  indexLength *= totalParentFlow;

  return indexLength;
}

template <typename FlowDataType, typename DeltaFlowDataType>
double MapEquation<FlowDataType, DeltaFlowDataType>::calcCodelengthOnModuleOfModules(const InfoNode& parent) const
{
  double parentFlow = parent.data.flow;
  double parentExit = parent.data.exitFlow;
  if (parentFlow < 1e-16)
    return 0.0;

  // H(x) = -xlog(x), T = q + SUM(p), q = exitFlow, p = enterFlow
  // Normal format
  // L = q * -log(q/T) + SUM(p * -log(p/T))
  // Compact format
  // L = T * ( H(q/T) + SUM( H(p/T) ) )
  // Expanded format
  // L = q * -log(q) - q * -log(T) + SUM( p * -log(p) - p * -log(T) )
  // = T * log(T) - q*log(q) - SUM( p*log(p) )
  // = -H(T) + H(q) + SUM(H(p))
  // As T is not known, use expanded format to avoid two loops
  double sumEnter = 0.0;
  double sumEnterLogEnter = 0.0;
  for (const auto& node : parent) {
    sumEnter += node.data.enterFlow; // rate of enter to finer level
    sumEnterLogEnter += infomath::plogp(node.data.enterFlow);
  }
  // The possibilities from this module: Either exit to coarser level or enter one of its children
  double totalCodewordUse = parentExit + sumEnter;

  return infomath::plogp(totalCodewordUse) - sumEnterLogEnter - infomath::plogp(parentExit);
}

template <typename FlowDataType, typename DeltaFlowDataType>
void MapEquation<FlowDataType, DeltaFlowDataType>::calculateCodelengthTerms(std::vector<InfoNode*>& nodes)
{
  enter_log_enter = 0.0;
  flow_log_flow = 0.0;
  exit_log_exit = 0.0;
  enterFlow = 0.0;

  // For each module
  for (InfoNode* n : nodes) {
    InfoNode& node = *n;
    // own node/module codebook
    flow_log_flow += infomath::plogp(node.data.flow + node.data.exitFlow);

    // use of index codebook
    enter_log_enter += infomath::plogp(node.data.enterFlow);
    exit_log_exit += infomath::plogp(node.data.exitFlow);
    enterFlow += node.data.enterFlow;
  }
  enterFlow += exitNetworkFlow;
  enterFlow_log_enterFlow = infomath::plogp(enterFlow);
}

} // namespace infomap

#endif // MAPEQUATION_H_
