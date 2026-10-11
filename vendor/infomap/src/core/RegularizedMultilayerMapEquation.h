/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef REGULARIZED_MULTILAYER_MAPEQUATION_H_
#define REGULARIZED_MULTILAYER_MAPEQUATION_H_

#include "MapEquation.h"
#include "FlowData.h"
#include "../utils/Log.h"
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace infomap {

class InfoNode;

class RegularizedMultilayerMapEquation final : private MapEquation<FlowData, MemDeltaFlow> {
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

  std::ostream& print(std::ostream& out) const;
  friend std::ostream& operator<<(std::ostream&, const RegularizedMultilayerMapEquation&);

  // ===================================================
  // Init
  // ===================================================

  void init(const Config& config);

  void initTree(InfoNode& /*root*/) {}

  void initNetwork(InfoNode& root);

  void initSuperNetwork(InfoNode& root);

  void initSubNetwork(InfoNode& root);

  void initPartition(std::vector<InfoNode*>& nodes);

  // ===================================================
  // Codelength
  // ===================================================

  double calcCodelength(const InfoNode& parent) const;

  void addMemoryContributions(InfoNode& current, DeltaFlowDataType& oldModuleDelta, VectorMap<DeltaFlowDataType>& moduleDeltaFlow);

  void addTeleportationFlow(InfoNode& current, const std::vector<FlowDataType>& moduleFlowData, DeltaFlowDataType& oldModuleDelta, DeltaFlowDataType& newModuleDelta);
  void addTeleportationFlow(InfoNode& current, const std::vector<FlowDataType>& moduleFlowData, VectorMap<DeltaFlowDataType>& moduleDeltaFlow);

  double getDeltaCodelengthOnMovingNode(InfoNode& current,
                                        DeltaFlowDataType& oldModuleDelta,
                                        DeltaFlowDataType& newModuleDelta,
                                        std::vector<FlowDataType>& moduleFlowData,
                                        std::vector<unsigned int>& moduleMembers);

  using Base::hoistOldSide;

  double getDeltaCodelengthOnMovingNodeHoisted(InfoNode& current,
                                               DeltaFlowDataType& oldModuleDelta,
                                               const OldSideTerms& oldSide,
                                               DeltaFlowDataType& newModuleDelta,
                                               std::vector<FlowDataType>& moduleFlowData,
                                               std::vector<unsigned int>& moduleMembers);

  // ===================================================
  // Consolidation
  // ===================================================

  void updateCodelengthOnMovingNode(InfoNode& current,
                                    DeltaFlowDataType& oldModuleDelta,
                                    DeltaFlowDataType& newModuleDelta,
                                    std::vector<FlowDataType>& moduleFlowData,
                                    std::vector<unsigned int>& moduleMembers);

  void consolidateModules(std::vector<InfoNode*>& modules);

  // ===================================================
  // Debug
  // ===================================================

  void printDebug() const;

private:
  // ===================================================
  // Private member functions
  // ===================================================
  double calcCodelengthOnModuleOfLeafNodes(const InfoNode& parent) const;

  // ===================================================
  // Init
  // ===================================================

  void initPhysicalNodes(InfoNode& root);

  void initPartitionOfPhysicalNodes(std::vector<InfoNode*>& nodes);
  void initPartitionLayerTeleFlowData(std::vector<InfoNode*>& nodes);

  // ===================================================
  // Codelength
  // ===================================================

  void calculateCodelength(std::vector<InfoNode*>& nodes);

  using Base::calculateCodelengthTerms;

  using Base::calculateCodelengthFromCodelengthTerms;

  void calculateNodeFlow_log_nodeFlow();

  // ===================================================
  // Consolidation
  // ===================================================

  void updatePhysicalNodes(InfoNode& current, unsigned int oldModuleIndex, unsigned int bestModuleIndex);

  void addMemoryContributionsAndUpdatePhysicalNodes(InfoNode& current, DeltaFlowDataType& oldModuleDelta, DeltaFlowDataType& newModuleDelta);

  void addLayerTeleFlow(unsigned int moduleIndex, const std::vector<LayerTeleFlowData>& layerTeleFlowData);
  void removeLayerTeleFlow(unsigned int moduleIndex, const std::vector<LayerTeleFlowData>& layerTeleFlowData);

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

  using ModuleToMemNodes = std::vector<ModuleMemNodes>; // sorted by module id
  using LayerTeleFlowMap = std::map<unsigned int, LayerTeleFlowData>;

  std::vector<ModuleToMemNodes> m_physToModuleToMemNodes; // vector[physicalNodeID] sorted vector of {moduleID, #memNodes, sumFlow}
  std::vector<LayerTeleFlowMap> m_moduleLayerTeleFlowData; // vector[moduleID] map<layerID, layer teleport flow>
  unsigned int m_numPhysicalNodes = 0;
  bool m_memoryContributionsAdded = false;
};

} // namespace infomap

#endif // REGULARIZED_MULTILAYER_MAPEQUATION_H_
