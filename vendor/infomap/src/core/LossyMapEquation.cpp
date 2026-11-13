/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION

#include "LossyMapEquation.h"
#include "InfoEdge.h"
#include "InfoNode.h"
#include "../utils/Log.h"

#include <algorithm>
#include <unordered_map>
#include <vector>

namespace infomap {

// ===================================================
// IO
// ===================================================

std::ostream& LossyMapEquation::print(std::ostream& out) const
{
  return out << indexCodelength << " + " << (moduleCodelength - m_sumCorrection)
             << " = " << io::toPrecision(getCodelength());
}

std::ostream& operator<<(std::ostream& out, const LossyMapEquation& mapEq)
{
  return mapEq.print(out);
}

// ===================================================
// Init
// ===================================================

void LossyMapEquation::init(const Config& config)
{
  Log(3) << "LossyMapEquation::init()...\n";
  Base::init(config);
  m_lambda = config.lossyLambda;
}

void LossyMapEquation::initNetwork(InfoNode& root)
{
  Base::initNetwork(root);

  // Sub-Infomap networks are cloned from one module of the full network and lack the
  // links leaving the module, so recomputing h_alpha there would understate the
  // distortion; the cloned nodes already carry the full-network values.
  if (root.owner != nullptr)
    return;

  // Compute each leaf node's share of the Markov entropy rate, p_alpha * h_alpha,
  // and its plogp(p_alpha), from edge weights. h_alpha is the entropy of the walk's
  // transition over the node's distinct neighbours, so aggregate edge weights per
  // neighbour first: an undirected edge may reach the node as both an out- and an
  // in-edge (when the input lists the link in both directions) and genuine parallel
  // edges may repeat a neighbour; summing per neighbour makes h_alpha independent of
  // how the same undirected adjacency was written. Accumulate the totals on the root
  // so they survive tree rewrites that bypass consolidateModules.
  root.lossyEntropy = 0.0;
  root.lossyFlowLogFlow = 0.0;
  std::unordered_map<const InfoNode*, double> neighborWeight;
  for (InfoNode& node : root) {
    neighborWeight.clear();
    for (InfoEdge* e : node.outEdges())
      neighborWeight[e->target] += e->data.weight;
    for (InfoEdge* e : node.inEdges())
      neighborWeight[e->source] += e->data.weight;
    double sumWLogW = 0.0;
    double strength = 0.0;
    for (const auto& nw : neighborWeight) {
      sumWLogW += infomath::plogp(nw.second);
      strength += nw.second;
    }
    // h = -sum (w/s) log2 (w/s) = (plogp(s) - sum plogp(w)) / s
    double h = strength > 1e-16 ? (infomath::plogp(strength) - sumWLogW) / strength : 0.0;
    node.lossyEntropy = node.data.flow * h;
    node.lossyFlowLogFlow = infomath::plogp(node.data.flow);
    root.lossyEntropy += node.lossyEntropy;
    root.lossyFlowLogFlow += node.lossyFlowLogFlow;
  }
}

void LossyMapEquation::initPartition(std::vector<InfoNode*>& nodes)
{
  calculateCodelength(nodes);
}

// ===================================================
// Codelength
// ===================================================

double LossyMapEquation::calcLoss(double moduleFlow, double flowLogFlow) const
{
  return infomath::plogp(moduleFlow) - flowLogFlow;
}

double LossyMapEquation::calcCorrection(double moduleFlow, double flowLogFlow, double entropy) const
{
  return std::max(0.0, calcLoss(moduleFlow, flowLogFlow) - m_lambda * entropy);
}

void LossyMapEquation::calculateCodelength(std::vector<InfoNode*>& nodes)
{
  Base::calculateCodelength(nodes);

  // Rebuild per-module aggregates; at this point each active node is its own
  // module and node->index is its module id.
  auto numNodes = nodes.size();
  m_moduleEntropy.assign(numNodes, 0.0);
  m_moduleFlowLogFlow.assign(numNodes, 0.0);
  m_sumCorrection = 0.0;
  m_sumNoiseLoss = 0.0;
  m_sumNoiseEntropy = 0.0;

  for (InfoNode* n : nodes) {
    auto i = n->index;
    m_moduleEntropy[i] = n->lossyEntropy;
    m_moduleFlowLogFlow[i] = n->lossyFlowLogFlow;
    double c = calcCorrection(n->data.flow, n->lossyFlowLogFlow, n->lossyEntropy);
    if (c > 0.0) {
      m_sumCorrection += c;
      m_sumNoiseLoss += calcLoss(n->data.flow, n->lossyFlowLogFlow);
      m_sumNoiseEntropy += n->lossyEntropy;
    }
  }
}

double LossyMapEquation::calcCodelength(const InfoNode& parent) const
{
  double L = Base::calcCodelength(parent);

  // The one-level baseline (the whole network as a single module) is the lossless
  // map-equation code length L_1 = H(p_alpha): a fixed, lambda-independent upper
  // reference. Capture it before any noise correction so reporting and relative
  // savings use L_1 rather than the lambda-dependent corrected one-module objective
  // (which still drives the optimizer via m_oneLevelCodelength). Only the root whose
  // children are the leaf nodes equals L_1; once the root is a module-of-modules,
  // Base::calcCodelength(root) is an index codebook, so guard on isLeafModule() to
  // avoid overwriting L_1 with a partition-dependent value.
  if (parent.isRoot() && parent.isLeafModule())
    m_oneLevelLossless = L;

  // Corrections apply to modules that own a visit codebook, i.e. leaf modules in
  // the two-level partition. The root (and any module of modules) is an index
  // codebook and must not be treated as an anonymizable module.
  if (!parent.isLeafModule())
    return L;

  // Subtract this module's correction, computed from leaf sums carried by the children.
  double flowLogFlow = 0.0;
  double entropy = 0.0;
  for (const auto& node : parent) {
    flowLogFlow += node.lossyFlowLogFlow;
    entropy += node.lossyEntropy;
  }
  return L - calcCorrection(parent.data.flow, flowLogFlow, entropy);
}

INFOMAP_HOT double LossyMapEquation::getDeltaCodelengthOnMovingNode(InfoNode& current,
                                                                    DeltaFlow& oldModuleDelta,
                                                                    DeltaFlow& newModuleDelta,
                                                                    std::vector<FlowData>& moduleFlowData,
                                                                    std::vector<unsigned int>& moduleMembers)
{
  double deltaL = Base::getDeltaCodelengthOnMovingNode(current, oldModuleDelta, newModuleDelta, moduleFlowData, moduleMembers);

  unsigned int oldM = oldModuleDelta.module;
  unsigned int newM = newModuleDelta.module;
  if (oldM == newM)
    return deltaL;

  double curFlow = current.data.flow;
  double curFlf = current.lossyFlowLogFlow;
  double curEnt = current.lossyEntropy;

  double corrBefore = calcCorrection(moduleFlowData[oldM].flow, m_moduleFlowLogFlow[oldM], m_moduleEntropy[oldM])
      + calcCorrection(moduleFlowData[newM].flow, m_moduleFlowLogFlow[newM], m_moduleEntropy[newM]);
  double corrAfter = calcCorrection(moduleFlowData[oldM].flow - curFlow, m_moduleFlowLogFlow[oldM] - curFlf, m_moduleEntropy[oldM] - curEnt)
      + calcCorrection(moduleFlowData[newM].flow + curFlow, m_moduleFlowLogFlow[newM] + curFlf, m_moduleEntropy[newM] + curEnt);

  // J = L_full - sum of corrections, so a correction increase lowers the objective.
  return deltaL - (corrAfter - corrBefore);
}

INFOMAP_HOT double LossyMapEquation::getDeltaCodelengthOnMovingNodeHoisted(InfoNode& current,
                                                                           DeltaFlow& oldModuleDelta,
                                                                           const OldSideTerms& oldSide,
                                                                           DeltaFlow& newModuleDelta,
                                                                           std::vector<FlowData>& moduleFlowData,
                                                                           std::vector<unsigned int>& moduleMembers)
{
  double deltaL = Base::getDeltaCodelengthOnMovingNodeHoisted(current, oldModuleDelta, oldSide, newModuleDelta, moduleFlowData, moduleMembers);

  unsigned int oldM = oldModuleDelta.module;
  unsigned int newM = newModuleDelta.module;
  if (oldM == newM)
    return deltaL;

  double curFlow = current.data.flow;
  double curFlf = current.lossyFlowLogFlow;
  double curEnt = current.lossyEntropy;

  double corrBefore = calcCorrection(moduleFlowData[oldM].flow, m_moduleFlowLogFlow[oldM], m_moduleEntropy[oldM])
      + calcCorrection(moduleFlowData[newM].flow, m_moduleFlowLogFlow[newM], m_moduleEntropy[newM]);
  double corrAfter = calcCorrection(moduleFlowData[oldM].flow - curFlow, m_moduleFlowLogFlow[oldM] - curFlf, m_moduleEntropy[oldM] - curEnt)
      + calcCorrection(moduleFlowData[newM].flow + curFlow, m_moduleFlowLogFlow[newM] + curFlf, m_moduleEntropy[newM] + curEnt);

  // J = L_full - sum of corrections, so a correction increase lowers the objective.
  return deltaL - (corrAfter - corrBefore);
}

// ===================================================
// Consolidation
// ===================================================

void LossyMapEquation::updateCodelengthOnMovingNode(InfoNode& current,
                                                    DeltaFlow& oldModuleDelta,
                                                    DeltaFlow& newModuleDelta,
                                                    std::vector<FlowData>& moduleFlowData,
                                                    std::vector<unsigned int>& moduleMembers)
{
  unsigned int oldM = oldModuleDelta.module;
  unsigned int newM = newModuleDelta.module;

  auto removeModuleTerms = [&](unsigned int m) {
    double c = calcCorrection(moduleFlowData[m].flow, m_moduleFlowLogFlow[m], m_moduleEntropy[m]);
    if (c > 0.0) {
      m_sumCorrection -= c;
      m_sumNoiseLoss -= calcLoss(moduleFlowData[m].flow, m_moduleFlowLogFlow[m]);
      m_sumNoiseEntropy -= m_moduleEntropy[m];
    }
  };
  auto addModuleTerms = [&](unsigned int m) {
    double c = calcCorrection(moduleFlowData[m].flow, m_moduleFlowLogFlow[m], m_moduleEntropy[m]);
    if (c > 0.0) {
      m_sumCorrection += c;
      m_sumNoiseLoss += calcLoss(moduleFlowData[m].flow, m_moduleFlowLogFlow[m]);
      m_sumNoiseEntropy += m_moduleEntropy[m];
    }
  };

  removeModuleTerms(oldM);
  if (newM != oldM)
    removeModuleTerms(newM);

  Base::updateCodelengthOnMovingNode(current, oldModuleDelta, newModuleDelta, moduleFlowData, moduleMembers);

  m_moduleEntropy[oldM] -= current.lossyEntropy;
  m_moduleEntropy[newM] += current.lossyEntropy;
  m_moduleFlowLogFlow[oldM] -= current.lossyFlowLogFlow;
  m_moduleFlowLogFlow[newM] += current.lossyFlowLogFlow;

  addModuleTerms(oldM);
  if (newM != oldM)
    addModuleTerms(newM);
}

void LossyMapEquation::consolidateModules(std::vector<InfoNode*>& modules)
{
  // Carry the leaf sums onto the new module nodes so the next aggregation level
  // and calcCodelength see them.
  for (unsigned int i = 0; i < modules.size(); ++i) {
    if (modules[i] == nullptr)
      continue;
    modules[i]->lossyEntropy = m_moduleEntropy[i];
    modules[i]->lossyFlowLogFlow = m_moduleFlowLogFlow[i];
  }
}

// ===================================================
// Debug
// ===================================================

void LossyMapEquation::printDebug() const
{
  Log() << "LossyMapEquation\n";
  Base::printDebug();
}

} // namespace infomap

#endif // INFOMAP_FEATURE_LOSSY_MAP_EQUATION
