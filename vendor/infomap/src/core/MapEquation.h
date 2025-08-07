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
#include <iostream>

namespace infomap {

class InfoNode;

template <typename FlowDataType = FlowData, typename DeltaFlowDataType = DeltaFlow>
class MapEquation {
  using ME = MapEquation<FlowDataType, DeltaFlowDataType>;

public:
  MapEquation() = default;

  MapEquation(const MapEquation& other) = default;

  MapEquation& operator=(const MapEquation& other) = default;

  MapEquation(MapEquation&& other) noexcept = default;

  MapEquation& operator=(MapEquation&& other) noexcept = default;

  virtual ~MapEquation() = default;

  // ===================================================
  // Getters
  // ===================================================

  virtual double getIndexCodelength() const { return indexCodelength; }

  virtual double getModuleCodelength() const { return moduleCodelength; }

  virtual double getCodelength() const { return codelength; }

  // ===================================================
  // IO
  // ===================================================

  virtual std::ostream& print(std::ostream& out) const
  {
    return out << indexCodelength << " + " << moduleCodelength << " = " << io::toPrecision(codelength);
  }

  // ===================================================
  // Init
  // ===================================================

  virtual void init(const Config&)
  {
    Log(3) << "MapEquation::init()...\n";
  }

  virtual void initTree(InfoNode& /*root*/) = 0;

  virtual void initNetwork(InfoNode& root)
  {
    Log(3) << "MapEquation::initNetwork()...\n";

    nodeFlow_log_nodeFlow = 0.0;
    for (InfoNode& node : root) {
      nodeFlow_log_nodeFlow += infomath::plogp(node.data.flow);
    }
    ME::initSubNetwork(root);
  }

  virtual void initSuperNetwork(InfoNode& root)
  {
    Log(3) << "MapEquation::initSuperNetwork()...\n";

    nodeFlow_log_nodeFlow = 0.0;
    for (InfoNode& node : root) {
      nodeFlow_log_nodeFlow += infomath::plogp(node.data.enterFlow);
    }
  }

  virtual void initSubNetwork(InfoNode& root)
  {
    exitNetworkFlow = root.data.exitFlow;
    exitNetworkFlow_log_exitNetworkFlow = infomath::plogp(exitNetworkFlow);
  }

  virtual void initPartition(std::vector<InfoNode*>& nodes) { ME::calculateCodelength(nodes); }

  // ===================================================
  // Codelength
  // ===================================================

  virtual double calcCodelength(const InfoNode& parent) const
  {
    return parent.isLeafModule() ? ME::calcCodelengthOnModuleOfLeafNodes(parent) : ME::calcCodelengthOnModuleOfModules(parent);
  }

  virtual void addMemoryContributions(InfoNode& /*current*/, DeltaFlowDataType& /*oldModuleDelta*/, DeltaFlowDataType& /*newModuleDelta*/) { }

  virtual void addMemoryContributions(InfoNode& /*current*/, DeltaFlowDataType& /*oldModuleDelta*/, VectorMap<DeltaFlowDataType>& /*moduleDeltaFlow*/) { }

  virtual double getDeltaCodelengthOnMovingNode(InfoNode& current,
                                                DeltaFlowDataType& oldModuleDelta,
                                                DeltaFlowDataType& newModuleDelta,
                                                std::vector<FlowDataType>& moduleFlowData,
                                                std::vector<unsigned int>& /*moduleMembers*/);

  // ===================================================
  // Consolidation
  // ===================================================

  virtual void updateCodelengthOnMovingNode(InfoNode& current,
                                            DeltaFlowDataType& oldModuleDelta,
                                            DeltaFlowDataType& newModuleDelta,
                                            std::vector<FlowDataType>& moduleFlowData,
                                            std::vector<unsigned int>& /*moduleMembers*/);

  virtual void consolidateModules(std::vector<InfoNode*>& /*modules*/) = 0;

  // ===================================================
  // Debug
  // ===================================================

  virtual void printDebug() const
  {
    std::cout << "(enterFlow_log_enterFlow: " << enterFlow_log_enterFlow << ", "
              << "enter_log_enter: " << enter_log_enter << ", "
              << "exitNetworkFlow_log_exitNetworkFlow: " << exitNetworkFlow_log_exitNetworkFlow << ") ";
  }

protected:
  // ===================================================
  // Protected member functions
  // ===================================================

  virtual double calcCodelengthOnModuleOfLeafNodes(const InfoNode& parent) const;

  virtual double calcCodelengthOnModuleOfModules(const InfoNode& parent) const;

  virtual void calculateCodelength(std::vector<InfoNode*>& nodes)
  {
    ME::calculateCodelengthTerms(nodes);
    ME::calculateCodelengthFromCodelengthTerms();
  }

  virtual void calculateCodelengthTerms(std::vector<InfoNode*>& nodes);

  virtual void calculateCodelengthFromCodelengthTerms()
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

template <typename FlowDataType, typename DeltaFlowDataType>
double MapEquation<FlowDataType, DeltaFlowDataType>::getDeltaCodelengthOnMovingNode(InfoNode& current, DeltaFlowDataType& oldModuleDelta, DeltaFlowDataType& newModuleDelta, std::vector<FlowDataType>& moduleFlowData, std::vector<unsigned int>&)
{
  using infomath::plogp;
  unsigned int oldModule = oldModuleDelta.module;
  unsigned int newModule = newModuleDelta.module;
  double deltaEnterExitOldModule = oldModuleDelta.deltaEnter + oldModuleDelta.deltaExit;
  double deltaEnterExitNewModule = newModuleDelta.deltaEnter + newModuleDelta.deltaExit;

  double delta_enter = plogp(enterFlow + deltaEnterExitOldModule - deltaEnterExitNewModule) - enterFlow_log_enterFlow;

  double delta_enter_log_enter = -plogp(moduleFlowData[oldModule].enterFlow)
      - plogp(moduleFlowData[newModule].enterFlow)
      + plogp(moduleFlowData[oldModule].enterFlow - current.data.enterFlow + deltaEnterExitOldModule)
      + plogp(moduleFlowData[newModule].enterFlow + current.data.enterFlow - deltaEnterExitNewModule);

  double delta_exit_log_exit = -plogp(moduleFlowData[oldModule].exitFlow)
      - plogp(moduleFlowData[newModule].exitFlow)
      + plogp(moduleFlowData[oldModule].exitFlow - current.data.exitFlow + deltaEnterExitOldModule)
      + plogp(moduleFlowData[newModule].exitFlow + current.data.exitFlow - deltaEnterExitNewModule);

  double delta_flow_log_flow = -plogp(moduleFlowData[oldModule].exitFlow + moduleFlowData[oldModule].flow)
      - plogp(moduleFlowData[newModule].exitFlow + moduleFlowData[newModule].flow)
      + plogp(moduleFlowData[oldModule].exitFlow + moduleFlowData[oldModule].flow
              - current.data.exitFlow - current.data.flow + deltaEnterExitOldModule)
      + plogp(moduleFlowData[newModule].exitFlow + moduleFlowData[newModule].flow
              + current.data.exitFlow + current.data.flow - deltaEnterExitNewModule);

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
