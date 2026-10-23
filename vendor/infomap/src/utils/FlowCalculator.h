/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef FLOW_CALCULATOR_H_
#define FLOW_CALCULATOR_H_

#include <unordered_map>
#include <string>
#include <vector>

namespace infomap {

struct Config;
class StateNetwork;

namespace detail {
  struct FlowLink {
    unsigned int source;
    unsigned int target;
    double flow;
  };
} // namespace detail

/**
 * Calculate flow on network based on different flow models
 */
class FlowCalculator {
public:
  FlowCalculator(StateNetwork&, const Config&);

private:
  void calcUndirectedFlow() noexcept;
  void calcDirectedFlow(const StateNetwork&, const Config&) noexcept;
  void calcDirectedRelaxToSelfFlow(const StateNetwork&, const Config&) noexcept;
  void calcUndirectedRegularizedFlow(const StateNetwork&, const Config&) noexcept;
  void calcUndirectedRegularizedMultilayerFlow(const StateNetwork&, const Config&);
  void calcDirectedRegularizedFlow(const StateNetwork&, const Config&) noexcept;
  // Not noexcept: this one throws on a flow imbalance it cannot normalize away,
  // and a throw escaping a noexcept function calls std::terminate rather than
  // reaching Infomap's error handler (#1019). The same reason
  // calcUndirectedRegularizedMultilayerFlow and usePrecomputedFlow are unmarked;
  // every flow method that cannot throw keeps the marking.
  void calcDirectedRegularizedMultilayerFlow(const StateNetwork&, const Config&);
  void calcDirectedBipartiteFlow(const StateNetwork&, const Config&) noexcept;
  void calcDirdirFlow(const Config&) noexcept;
  void calcRawdirFlow() noexcept;
  void usePrecomputedFlow(const StateNetwork&, const Config&);

  void finalize(StateNetwork&, const Config&, bool) noexcept;

  void addFlowNote(const std::string& note);
  void recordPageRank(unsigned int iterations, double error, bool converged) noexcept;

  double accumulateDanglingRank() const noexcept;

  unsigned int numNodes;
  //! One past the last dangling node, for the directed model's dangling-first
  //! ordering. Zero when that ordering was not applied -- see danglingIndices.
  unsigned int nonDanglingStartIndex = 0;
  unsigned int bipartiteStartIndex = 0;
  unsigned int bipartiteLinkStartIndex = 0;

  double sumLinkWeight = 0;
  double sumWeightedDegree = 0;
  double sumTeleFlow = 0;

  std::unordered_map<unsigned int, unsigned int> nodeIndexMap;
  std::vector<double> nodeFlow;
  std::vector<double> nodeTeleportWeights;
  std::vector<double> nodeTeleportFlow;
  std::vector<double> enterFlow;
  std::vector<double> exitFlow;
  std::vector<double> sumLinkOutWeight;
  std::vector<unsigned int> nodeOutDegree;
  //! Positions of the dangling nodes when they are not a prefix of nodeFlow.
  //! Bipartite input keeps the network's node order, because the feature side is
  //! identified by an index range, so the dangling-first ordering is unavailable.
  std::vector<unsigned int> danglingIndices;
  using FlowLink = detail::FlowLink;
  std::vector<FlowLink> flowLinks;

  std::string m_flowMethod;
  std::string m_teleportation;
  std::vector<std::string> m_flowNotes;
  unsigned int m_pageRankIterations = 0;
  double m_pageRankError = 0.0;
  bool m_havePageRank = false;
  bool m_pageRankConverged = true;
};

inline void calculateFlow(StateNetwork& network, const Config& config)
{
  FlowCalculator(network, config);
}

} // namespace infomap

#endif // FLOW_CALCULATOR_H_
