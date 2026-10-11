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
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "StateNetwork.h"

namespace infomap {

namespace {

  // Miller-Madow charges (K - 1) / (2n) *nats* for a codebook with K codewords
  // estimated from n observations, and the map equation is in bits, hence 1/ln2.
  const double LOG_2 = std::log(2.0);

  // Can the walker ever leave this module, i.e. does its codebook hold an exit
  // codeword? Leaving a module means emitting its exit codeword and then entering
  // a sibling somewhere further up, so a module that is its parent's only child
  // can only be left if that parent can be left. The root can never be left, so a
  // chain of only-children below the root has no exit codeword anywhere along it.
  bool hasExitCodeword(const InfoNode& module)
  {
    for (const InfoNode* node = &module; node->parent != nullptr; node = node->parent) {
      if (node->parent->childDegree() > 1)
        return true;
    }
    return false;
  }

} // namespace

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

// The network is the sample: an undirected link of weight w is w walk steps in
// each direction, so the sample size is the total degree D, and every codebook is
// a multinomial count vector over those steps. Codebook c is used u_c * D times
// and the map equation weights it by that same u_c, so the usage cancels and each
// codebook costs the same per free parameter:
//
//   sum_c u_c * (K_c - 1) / (2 u_c D ln2) = sum_c (K_c - 1) / (2 D ln2)
//
// Equivalently: L is the entropy rate of the Markov chain over codebooks, whose
// ML estimate is biased by -(free parameters)/(2D) nats.
double BiasedMapEquation::bitsPerFreeParameter() const
{
  return entropyBiasCorrectionMultiplier / (2 * m_totalDegree * LOG_2);
}

double BiasedMapEquation::calcIndexEntropyBiasCorrection(unsigned int numModules) const
{
  // One codeword per module, so numModules - 1 free parameters -- and none at all
  // for a single module, whose index codebook is never used.
  if (!useEntropyBiasCorrection)
    return 0;
  const double numCodewords = static_cast<double>(numModules);
  return bitsPerFreeParameter() * std::max(numCodewords - 1, 0.0);
}

double BiasedMapEquation::calcModuleEntropyBiasCorrection(unsigned int numModules) const
{
  // Each module codebook holds a codeword per node plus an exit codeword, i.e.
  // (n_i + 1) - 1 = n_i free parameters, summing to numNodes over the partition.
  // The alphabet is the one the partition declares, not the one this sample
  // happens to exercise: a node or a boundary link that no walk step touched is
  // unobserved, not impossible, which is the premise the correction exists to
  // serve. The single exception is impossible by construction rather than
  // unobserved -- one module holding the whole network can never be exited, so
  // that codeword does not exist and the partition has one parameter less.
  if (!useEntropyBiasCorrection)
    return 0;
  const double numParameters = numModules == 1 ? std::max(m_numNodes - 1.0, 0.0) : m_numNodes;
  return bitsPerFreeParameter() * numParameters;
}

double BiasedMapEquation::calcEntropyBiasCorrection(unsigned int numModules) const
{
  return calcIndexEntropyBiasCorrection(numModules) + calcModuleEntropyBiasCorrection(numModules);
}

void BiasedMapEquation::calculateCodelength(std::vector<InfoNode*>& nodes)
{
  calculateCodelengthTerms(nodes);

  calculateCodelengthFromCodelengthTerms();

  currentNumModules = nodes.size();

  biasedCost = calcNumModuleCost(currentNumModules);

  indexEntropyBiasCorrection = calcIndexEntropyBiasCorrection(currentNumModules);
  moduleEntropyBiasCorrection = calcModuleEntropyBiasCorrection(currentNumModules);
}

double BiasedMapEquation::calcCodelength(const InfoNode& parent) const
{
  return parent.isLeafModule()
      ? calcCodelengthOnModuleOfLeafNodes(parent)
      : calcCodelengthOnModuleOfModules(parent);
}

// Charged once for the whole tree, not per node. biasedCost is a single scalar in
// |K - K_pref|, so there is no share of it for calcCodelength to return, which is
// why the scored value used to lose the term entirely while the search kept it --
// up to 7 bits on a six-node fixture, unbounded in K_pref (#1021). Takes the count
// from the caller rather than from currentNumModules so a partition that was
// materialized rather than searched is charged for the modules it actually has, and
// so the one-level reference can price a module count its bare root does not show.
double BiasedMapEquation::calcTreeCodelengthCost(unsigned int numTopModules) const
{
  return calcNumModuleCost(numTopModules);
}

// Free parameters of the codebook this tree node owns: one codeword per child
// plus an exit codeword, minus one for the codebook's own normalisation. The root
// owns the index codebook, which has no exit codeword, and neither does any module
// that cannot be left (see hasExitCodeword).
double BiasedMapEquation::calcCodebookFreeParameters(const InfoNode& parent) const
{
  const double numCodewords = static_cast<double>(parent.childDegree()) + (hasExitCodeword(parent) ? 1.0 : 0.0);
  return std::max(numCodewords - 1, 0.0);
}

double BiasedMapEquation::calcCodelengthOnModuleOfModules(const InfoNode& parent) const
{
  double L = Base::calcCodelengthOnModuleOfModules(parent);
  if (!useEntropyBiasCorrection)
    return L;

  // NOTE (#830): this recompute and the tracked value now agree on two-level
  // trees -- the root term above was the whole of the disagreement there -- but
  // not on deeper ones. The tracked side charges (numModules - 1) + numNodes,
  // which is every codebook a two-level partition has; this recompute walks the
  // whole tree, so it also charges the intermediate module-of-modules codebooks,
  // one free parameter per non-root module-of-modules more. This side is the
  // correct one. The tracked side is left alone here because making it tree-aware
  // changes what the search minimises, not what it reports.
  return L + bitsPerFreeParameter() * calcCodebookFreeParameters(parent);
}

double BiasedMapEquation::calcCodelengthOnModuleOfLeafNodes(const InfoNode& parent) const
{
  double L = Base::calcCodelength(parent);
  if (!useEntropyBiasCorrection)
    return L;

  return L + bitsPerFreeParameter() * calcCodebookFreeParameters(parent);
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
  moduleEntropyBiasCorrection = calcModuleEntropyBiasCorrection(currentNumModules);
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
