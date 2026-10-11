/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef LOSSY_MAPEQUATION_H_
#define LOSSY_MAPEQUATION_H_

#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION

#include "MapEquation.h"
#include "FlowData.h"
#include <vector>

namespace infomap {

class InfoNode;

// Lossy map equation (rate-distortion): each module is standard or noise; a noise
// module shares one visit codeword across its members. With per-module loss l_i
// (the naming overhead the shared codeword removes) and entropy share H_i (the
// members' share of the walk's Markov entropy rate),
//   J_lambda = L_full - sum_i max(0, l_i - lambda * H_i),
// so the optimal type follows in closed form from the membership alone:
// noise iff l_i > lambda * H_i. The base class maintains L_full.
class LossyMapEquation final : private MapEquation<> {
  using Base = MapEquation<>;

public:
  using FlowDataType = FlowData;
  using DeltaFlowDataType = DeltaFlow;

  // ===================================================
  // Getters
  // ===================================================

  double getIndexCodelength() const { return indexCodelength; }

  double getModuleCodelength() const { return moduleCodelength - m_sumCorrection; }

  double getCodelength() const { return codelength - m_sumCorrection; } // J_lambda

  double getLossyRate() const { return codelength - m_sumNoiseLoss; } // L_lossy

  double getLossyDistortion() const { return m_sumNoiseEntropy; } // D

  // L_1 = H(p_alpha): the lossless one-level reference, captured in calcCodelength
  // of the root. Independent of lambda; used for reporting and relative savings.
  double getOneLevelLossless() const { return m_oneLevelLossless; }

  // ===================================================
  // IO
  // ===================================================

  std::ostream& print(std::ostream& out) const;

  friend std::ostream& operator<<(std::ostream&, const LossyMapEquation&);

  // ===================================================
  // Init
  // ===================================================

  void init(const Config& config);

  void initTree(InfoNode& /*root*/) {}

  void initNetwork(InfoNode& root);

  using Base::initSuperNetwork;

  using Base::initSubNetwork;

  void initPartition(std::vector<InfoNode*>& nodes);

  // ===================================================
  // Codelength
  // ===================================================

  double calcCodelength(const InfoNode& parent) const;

  using Base::addMemoryContributions;

  using Base::addTeleportationFlow;

  double getDeltaCodelengthOnMovingNode(InfoNode& current,
                                        DeltaFlow& oldModuleDelta,
                                        DeltaFlow& newModuleDelta,
                                        std::vector<FlowData>& moduleFlowData,
                                        std::vector<unsigned int>& moduleMembers);

  using Base::hoistOldSide;

  double getDeltaCodelengthOnMovingNodeHoisted(InfoNode& current,
                                               DeltaFlow& oldModuleDelta,
                                               const OldSideTerms& oldSide,
                                               DeltaFlow& newModuleDelta,
                                               std::vector<FlowData>& moduleFlowData,
                                               std::vector<unsigned int>& moduleMembers);

  // ===================================================
  // Consolidation
  // ===================================================

  void updateCodelengthOnMovingNode(InfoNode& current,
                                    DeltaFlow& oldModuleDelta,
                                    DeltaFlow& newModuleDelta,
                                    std::vector<FlowData>& moduleFlowData,
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

  void calculateCodelength(std::vector<InfoNode*>& nodes);

  // l_i = -sum_{alpha in i} p_alpha log2(p_alpha / P_i) = plogp(P_i) - sum plogp(p_alpha),
  // for a module with member flow P_i and member-leaf sum F_i of plogp(p_alpha).
  double calcLoss(double moduleFlow, double flowLogFlow) const;

  // Correction c_i = max(0, l_i - lambda * H_i); the module is noise iff c_i > 0.
  double calcCorrection(double moduleFlow, double flowLogFlow, double entropy) const;

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
  using Base::flow_log_flow;
  using Base::nodeFlow_log_nodeFlow;

  // For hierarchical
  using Base::exitNetworkFlow;
  using Base::exitNetworkFlow_log_exitNetworkFlow;

  double m_lambda = 1.0;
  // Per-module aggregates, indexed like the optimizer's m_moduleFlowData.
  std::vector<double> m_moduleEntropy; // H_i
  std::vector<double> m_moduleFlowLogFlow; // F_i
  // Totals over all modules / noise modules.
  double m_sumCorrection = 0.0; // sum_i c_i  (J = codelength - m_sumCorrection)
  double m_sumNoiseLoss = 0.0; // sum over noise modules of l_i
  double m_sumNoiseEntropy = 0.0; // sum over noise modules of H_i
  mutable double m_oneLevelLossless = 0.0; // L_1 = H(p_alpha), set when scoring the root
};

} // namespace infomap

#endif // INFOMAP_FEATURE_LOSSY_MAP_EQUATION

#endif // LOSSY_MAPEQUATION_H_
