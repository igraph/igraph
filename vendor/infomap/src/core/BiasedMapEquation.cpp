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
#include "../utils/Log.h"

#include <vector>
#include <utility>
#include <cstdlib>
#include "StateNetwork.h"

namespace infomap {

void BiasedMapEquation::setNetworkProperties(const StateNetwork& network)
{
  auto totalDegree = network.sumWeightedDegree();
  // Negative entropy bias is based on discrete counts, if average weight is below 1, use unweighted total degree
  if (totalDegree < network.sumDegree()) {
    totalDegree = network.sumDegree();
  }
  m_totalDegree = totalDegree;
  m_numNodes = network.numNodes();
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
  Base::init(config);
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
  return useEntropyBiasCorrection ? entropyBiasCorrectionMultiplier * (numModules - 1) / (2 * m_totalDegree) : 0;
}

double BiasedMapEquation::calcModuleEntropyBiasCorrection() const
{
  return useEntropyBiasCorrection ? entropyBiasCorrectionMultiplier * m_numNodes / (2 * m_totalDegree) : 0;
}

double BiasedMapEquation::calcEntropyBiasCorrection(unsigned int numModules) const
{
  return useEntropyBiasCorrection ? entropyBiasCorrectionMultiplier * (numModules - 1 + m_numNodes) / (2 * m_totalDegree) : 0;
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

  // NOTE (#830): this recompute and the tracked value disagree, by a now fully
  // measured amount. The tracked side uses (numModules - 1 + numNodes) / (2D),
  // i.e. one term for the root index codebook and one for the leaf codewords.
  // This recompute instead sums childDegree over every internal node, which
  // differs in exactly two places:
  //   * the root is charged childDegree here but childDegree - 1 there, and the
  //     latter is right: the root index codebook has no exit codeword;
  //   * every non-root module-of-modules contributes childDegree here and
  //     nothing at all there.
  // Measured with the root term corrected: the residual drift is exactly
  // (number of non-root module-of-modules) / (2 * totalDegree) -- 0 on two-level
  // trees (lossy_benchmark, ninetriangles: exact match) and 1/514 on
  // unbalanced_hierarchy, which has one. Both sides are left alone here so the
  // root term and the missing intermediate term land together; fixing only the
  // root would make two-level runs exact while leaving hierarchical runs off by
  // a tree-shape-dependent amount.
  return L + entropyBiasCorrectionMultiplier * parent.childDegree() / (2 * m_totalDegree);
}

double BiasedMapEquation::calcCodelengthOnModuleOfLeafNodes(const InfoNode& parent) const
{
  double L = Base::calcCodelength(parent);
  if (!useEntropyBiasCorrection)
    return L;

  return L + entropyBiasCorrectionMultiplier * parent.childDegree() / (2 * m_totalDegree);
}

// The post-move module count. deltaNumModules is -1, 0 or +1, and a removal
// requires the node's old module to be its last member -- which needs a second,
// non-empty module to move into -- so the count cannot drop below one. Computed
// through int to keep the unsigned-wrap question out of the reader's way.
unsigned int BiasedMapEquation::numModulesAfterMove(int deltaNumModules) const
{
  return static_cast<unsigned int>(static_cast<int>(currentNumModules) + deltaNumModules);
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

INFOMAP_HOT double BiasedMapEquation::getDeltaCodelengthOnMovingNode(InfoNode& current,
                                                                     DeltaFlow& oldModuleDelta,
                                                                     DeltaFlow& newModuleDelta,
                                                                     std::vector<FlowData>& moduleFlowData,
                                                                     std::vector<unsigned int>& moduleMembers)
{
  double deltaL = Base::getDeltaCodelengthOnMovingNode(current, oldModuleDelta, newModuleDelta, moduleFlowData, moduleMembers);

  // Both extra terms below depend on the module count, so they only matter when
  // the move changes it -- but each has its own switch. The guard used to test
  // preferredNumModules alone, which silently excluded the entropy-bias
  // correction from every move delta under the default
  // --preferred-number-of-modules 0: the search then optimized the uncorrected
  // map equation while reporting the corrected codelength (#830, #904).
  if (preferredNumModules == 0 && !useEntropyBiasCorrection)
    return deltaL;

  int deltaNumModules = getDeltaNumModulesIfMoving(oldModuleDelta.module, newModuleDelta.module, moduleMembers);

  const unsigned int numModules = numModulesAfterMove(deltaNumModules);

  double deltaBiasedCost = preferredNumModules == 0
      ? 0.0
      : calcNumModuleCost(numModules) - biasedCost;

  double deltaEntropyBiasCorrection = calcEntropyBiasCorrection(numModules) - getEntropyBiasCorrection();

  return deltaL + deltaBiasedCost + deltaEntropyBiasCorrection;
}

INFOMAP_HOT double BiasedMapEquation::getDeltaCodelengthOnMovingNodeHoisted(InfoNode& current,
                                                                            DeltaFlow& oldModuleDelta,
                                                                            const OldSideTerms& oldSide,
                                                                            DeltaFlow& newModuleDelta,
                                                                            std::vector<FlowData>& moduleFlowData,
                                                                            std::vector<unsigned int>& moduleMembers)
{
  double deltaL = Base::getDeltaCodelengthOnMovingNodeHoisted(current, oldModuleDelta, oldSide, newModuleDelta, moduleFlowData, moduleMembers);

  // See getDeltaCodelengthOnMovingNode: the two terms have independent switches.
  if (preferredNumModules == 0 && !useEntropyBiasCorrection)
    return deltaL;

  int deltaNumModules = getDeltaNumModulesIfMoving(oldModuleDelta.module, newModuleDelta.module, moduleMembers);

  const unsigned int numModules = numModulesAfterMove(deltaNumModules);

  double deltaBiasedCost = preferredNumModules == 0
      ? 0.0
      : calcNumModuleCost(numModules) - biasedCost;

  double deltaEntropyBiasCorrection = calcEntropyBiasCorrection(numModules) - getEntropyBiasCorrection();

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

  // Must match getDeltaCodelengthOnMovingNode's guard exactly: if the delta
  // accounts for a term, the accepted move has to update the state that term is
  // computed from. Widening the delta guard while leaving this one on
  // preferredNumModules would leave currentNumModules and
  // indexEntropyBiasCorrection frozen at their initPartition values, so every
  // later delta would be measured against a stale module count.
  if (preferredNumModules == 0 && !useEntropyBiasCorrection)
    return;

  int deltaNumModules = getDeltaNumModulesIfMoving(oldModuleDelta.module, newModuleDelta.module, moduleMembers);

  currentNumModules = numModulesAfterMove(deltaNumModules);
  if (preferredNumModules != 0)
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
  Log() << "BiasedMapEquation\n";
  Base::printDebug();
}

} // namespace infomap
