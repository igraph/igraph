/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "BiasedMapEquation.h"
#include "FlowData.h"
#include "InfoNode.h"

#include <vector>
#include <utility>
#include <cstdlib>
#include "StateNetwork.h"

namespace infomap {

double BiasedMapEquation::s_totalDegree = 1;
unsigned int BiasedMapEquation::s_numNodes = 0;

void BiasedMapEquation::setNetworkProperties(const StateNetwork& network)
{
  s_totalDegree = network.sumWeightedDegree();
  // Negative entropy bias is based on discrete counts, if average weight is below 1, use unweighted total degree
  if (s_totalDegree < network.sumDegree()) {
    s_totalDegree = network.sumDegree();
  }
  s_numNodes = network.numNodes();
}

double BiasedMapEquation::getIndexCodelength() const
{
  return indexCodelength + indexEntropyBiasCorrection;
}

double BiasedMapEquation::getModuleCodelength() const
{
  return moduleCodelength + biasedCost + moduleEntropyBiasCorrection;
}

double BiasedMapEquation::getCodelength() const
{
  return codelength + biasedCost + getEntropyBiasCorrection();
}

double BiasedMapEquation::getEntropyBiasCorrection() const
{
  return indexEntropyBiasCorrection + moduleEntropyBiasCorrection;
}

// ===================================================
// IO
// ===================================================

std::ostream& BiasedMapEquation::print(std::ostream& out) const
{
  out << indexCodelength << " + " << moduleCodelength;
  if (preferredNumModules != 0) {
    out << " + " << biasedCost;
  }
  if (useEntropyBiasCorrection) {
    out << " + " << getEntropyBiasCorrection();
  }
  out << " = " << io::toPrecision(getCodelength());
  return out;
}

std::ostream& operator<<(std::ostream& out, const BiasedMapEquation& mapEq)
{
  return mapEq.print(out);
}

// ===================================================
// Init
// ===================================================

void BiasedMapEquation::init(const Config& config)
{
  Log(3) << "BiasedMapEquation::init()...\n";
  preferredNumModules = config.preferredNumberOfModules;
  useEntropyBiasCorrection = config.entropyBiasCorrection;
  entropyBiasCorrectionMultiplier = config.entropyBiasCorrectionMultiplier;
}

void BiasedMapEquation::initNetwork(InfoNode& root)
{
  Log(3) << "BiasedMapEquation::initNetwork()...\n";
  Base::initNetwork(root);
}

void BiasedMapEquation::initPartition(std::vector<InfoNode*>& nodes)
{
  calculateCodelength(nodes);
}

// ===================================================
// Codelength
// ===================================================

double BiasedMapEquation::calcNumModuleCost(unsigned int numModules) const
{
  if (preferredNumModules == 0) return 0;
  int deltaNumModules = numModules - preferredNumModules;
  return 1 * std::abs(deltaNumModules);
}

double BiasedMapEquation::calcIndexEntropyBiasCorrection(unsigned int numModules) const
{
  return useEntropyBiasCorrection ? entropyBiasCorrectionMultiplier * (numModules - 1) / (2 * s_totalDegree) : 0;
}

double BiasedMapEquation::calcModuleEntropyBiasCorrection() const
{
  return useEntropyBiasCorrection ? entropyBiasCorrectionMultiplier * s_numNodes / (2 * s_totalDegree) : 0;
}

double BiasedMapEquation::calcEntropyBiasCorrection(unsigned int numModules) const
{
  return useEntropyBiasCorrection ? entropyBiasCorrectionMultiplier * (numModules - 1 + s_numNodes) / (2 * s_totalDegree) : 0;
}

void BiasedMapEquation::calculateCodelength(std::vector<InfoNode*>& nodes)
{
  calculateCodelengthTerms(nodes);

  calculateCodelengthFromCodelengthTerms();

  currentNumModules = nodes.size();

  biasedCost = calcNumModuleCost(currentNumModules);

  indexEntropyBiasCorrection = calcIndexEntropyBiasCorrection(currentNumModules);
  moduleEntropyBiasCorrection = calcModuleEntropyBiasCorrection();
}

double BiasedMapEquation::calcCodelength(const InfoNode& parent) const
{
  return parent.isLeafModule()
      ? calcCodelengthOnModuleOfLeafNodes(parent)
      : calcCodelengthOnModuleOfModules(parent);
}

double BiasedMapEquation::calcCodelengthOnModuleOfModules(const InfoNode& parent) const
{
  double L = Base::calcCodelengthOnModuleOfModules(parent);
  if (!useEntropyBiasCorrection)
    return L;

  return L + entropyBiasCorrectionMultiplier * parent.childDegree() / (2 * s_totalDegree);
}

double BiasedMapEquation::calcCodelengthOnModuleOfLeafNodes(const InfoNode& parent) const
{
  double L = Base::calcCodelength(parent);
  if (!useEntropyBiasCorrection)
    return L;

  return L + entropyBiasCorrectionMultiplier * parent.childDegree() / (2 * s_totalDegree);
}

int BiasedMapEquation::getDeltaNumModulesIfMoving(unsigned int oldModule,
                                                  unsigned int newModule,
                                                  std::vector<unsigned int>& moduleMembers)
{
  bool removeOld = moduleMembers[oldModule] == 1;
  bool createNew = moduleMembers[newModule] == 0;
  int deltaNumModules = removeOld && !createNew ? -1 : (!removeOld && createNew ? 1 : 0);
  return deltaNumModules;
}

double BiasedMapEquation::getDeltaCodelengthOnMovingNode(InfoNode& current,
                                                         DeltaFlow& oldModuleDelta,
                                                         DeltaFlow& newModuleDelta,
                                                         std::vector<FlowData>& moduleFlowData,
                                                         std::vector<unsigned int>& moduleMembers)
{
  double deltaL = Base::getDeltaCodelengthOnMovingNode(current, oldModuleDelta, newModuleDelta, moduleFlowData, moduleMembers);

  if (preferredNumModules == 0)
    return deltaL;

  int deltaNumModules = getDeltaNumModulesIfMoving(oldModuleDelta.module, newModuleDelta.module, moduleMembers);

  double deltaBiasedCost = calcNumModuleCost(currentNumModules + deltaNumModules) - biasedCost;

  double deltaEntropyBiasCorrection = calcEntropyBiasCorrection(currentNumModules + deltaNumModules) - getEntropyBiasCorrection();

  return deltaL + deltaBiasedCost + deltaEntropyBiasCorrection;
}

// ===================================================
// Consolidation
// ===================================================

void BiasedMapEquation::updateCodelengthOnMovingNode(InfoNode& current,
                                                     DeltaFlow& oldModuleDelta,
                                                     DeltaFlow& newModuleDelta,
                                                     std::vector<FlowData>& moduleFlowData,
                                                     std::vector<unsigned int>& moduleMembers)
{
  Base::updateCodelengthOnMovingNode(current, oldModuleDelta, newModuleDelta, moduleFlowData, moduleMembers);

  if (preferredNumModules == 0)
    return;

  int deltaNumModules = getDeltaNumModulesIfMoving(oldModuleDelta.module, newModuleDelta.module, moduleMembers);

  currentNumModules += deltaNumModules;
  biasedCost = calcNumModuleCost(currentNumModules);
  indexEntropyBiasCorrection = calcIndexEntropyBiasCorrection(currentNumModules);
  moduleEntropyBiasCorrection = calcModuleEntropyBiasCorrection();
}

void BiasedMapEquation::consolidateModules(std::vector<InfoNode*>& modules)
{
  unsigned int numModules = 0;
  for (auto& module : modules) {
    if (module == nullptr)
      continue;
    ++numModules;
  }
  currentNumModules = numModules;
}

// ===================================================
// Debug
// ===================================================

void BiasedMapEquation::printDebug() const
{
  std::cout << "BiasedMapEquation\n";
  Base::printDebug();
}

} // namespace infomap
