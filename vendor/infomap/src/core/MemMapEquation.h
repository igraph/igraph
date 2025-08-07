/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef MEM_MAPEQUATION_H_
#define MEM_MAPEQUATION_H_

#include "MapEquation.h"
#include "FlowData.h"
#include "../utils/Log.h"
#include <vector>
#include <set>
#include <map>
#include <utility>

namespace infomap {

class InfoNode;
struct MemNodeSet;

class MemMapEquation : private MapEquation<FlowData, MemDeltaFlow> {
  using Base = MapEquation<FlowData, MemDeltaFlow>;

public:
  using FlowDataType = FlowData;
  using DeltaFlowDataType = MemDeltaFlow;

  // ===================================================
  // Getters
  // ===================================================

  using Base::getIndexCodelength;

  using Base::getModuleCodelength;

  using Base::getCodelength;

  // ===================================================
  // IO
  // ===================================================

  std::ostream& print(std::ostream& out) const override;
  friend std::ostream& operator<<(std::ostream&, const MemMapEquation&);

  // ===================================================
  // Init
  // ===================================================

  void init(const Config& config) override;

  void initTree(InfoNode& /*root*/) override { }

  void initNetwork(InfoNode& root) override;

  void initSuperNetwork(InfoNode& root) override;

  void initSubNetwork(InfoNode& root) override;

  void initPartition(std::vector<InfoNode*>& nodes) override;

  // ===================================================
  // Codelength
  // ===================================================

  double calcCodelength(const InfoNode& parent) const override;

  void addMemoryContributions(InfoNode& current, MemDeltaFlow& oldModuleDelta, VectorMap<MemDeltaFlow>& moduleDeltaFlow) override;

  double getDeltaCodelengthOnMovingNode(InfoNode& current,
                                        MemDeltaFlow& oldModuleDelta,
                                        MemDeltaFlow& newModuleDelta,
                                        std::vector<FlowData>& moduleFlowData,
                                        std::vector<unsigned int>& moduleMembers) override;

  // ===================================================
  // Consolidation
  // ===================================================

  void updateCodelengthOnMovingNode(InfoNode& current,
                                    MemDeltaFlow& oldModuleDelta,
                                    MemDeltaFlow& newModuleDelta,
                                    std::vector<FlowData>& moduleFlowData,
                                    std::vector<unsigned int>& moduleMembers) override;

  void consolidateModules(std::vector<InfoNode*>& modules) override;

  // ===================================================
  // Debug
  // ===================================================

  void printDebug() const override;

private:
  // ===================================================
  // Private member functions
  // ===================================================
  double calcCodelengthOnModuleOfLeafNodes(const InfoNode& parent) const override;

  // ===================================================
  // Init
  // ===================================================

  void initPhysicalNodes(InfoNode& root);

  void initPartitionOfPhysicalNodes(std::vector<InfoNode*>& nodes);

  // ===================================================
  // Codelength
  // ===================================================

  void calculateCodelength(std::vector<InfoNode*>& nodes) override;

  using Base::calculateCodelengthTerms;

  using Base::calculateCodelengthFromCodelengthTerms;

  void calculateNodeFlow_log_nodeFlow();

  // ===================================================
  // Consolidation
  // ===================================================

  void updatePhysicalNodes(InfoNode& current, unsigned int oldModuleIndex, unsigned int bestModuleIndex);

  void addMemoryContributionsAndUpdatePhysicalNodes(InfoNode& current, MemDeltaFlow& oldModuleDelta, MemDeltaFlow& newModuleDelta);

public:
  // ===================================================
  // Public member variables
  // ===================================================

  using Base::codelength;
  using Base::indexCodelength;
  using Base::moduleCodelength;

private:
  // ===================================================
  // Private member variables
  // ===================================================

  using Base::enter_log_enter;
  using Base::enterFlow;
  using Base::enterFlow_log_enterFlow;
  using Base::exit_log_exit;
  using Base::flow_log_flow; // node.(flow + exitFlow)
  using Base::nodeFlow_log_nodeFlow; // constant while the leaf network is the same

  // For hierarchical
  using Base::exitNetworkFlow;
  using Base::exitNetworkFlow_log_exitNetworkFlow;

  using ModuleToMemNodes = std::map<unsigned int, MemNodeSet>;

  std::vector<ModuleToMemNodes> m_physToModuleToMemNodes; // vector[physicalNodeID] map<moduleID, {#memNodes, sumFlow}>
  unsigned int m_numPhysicalNodes = 0;
  bool m_memoryContributionsAdded = false;
};

struct MemNodeSet {
  MemNodeSet(unsigned int numMemNodes, double sumFlow) : numMemNodes(numMemNodes), sumFlow(sumFlow) { }
  unsigned int numMemNodes; // use counter to check for zero to avoid round-off errors in sumFlow
  double sumFlow;
};

} // namespace infomap

#endif // MEM_MAPEQUATION_H_
