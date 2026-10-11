/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "FlowCalculator.h"
#include "../utils/Log.h"
#include "../utils/Console.h"
#include "../utils/convert.h"
#include "../utils/format.h"
#include "../utils/infomath.h"
#include "../core/StateNetwork.h"
#include "../io/InfomapError.h"
#include <cmath>
#include <numeric>
#include <limits>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <unordered_map>

namespace infomap {

namespace {

#if !INFOMAP_FEATURE_REGULARIZED_MULTILAYER
  const char* regularizedMultilayerFeatureError()
  {
    return "Regularized multilayer flow requires building with FEATURES=regularized-multilayer.";
  }
#endif

  // Converged means the error is within tolerance, not merely that the loop stopped
  // before the iteration limit. The two differ whenever the last allowed iteration is
  // the one that reaches tolerance: that was reported as a failure, with a warning to
  // match, while the error sat at 0. Shared by every power iteration in this file,
  // because the predicate existed in three copies and fixing one left the regularized
  // model exporting the old answer (#899).
  bool errorWithinFlowTolerance(const Config& config, double err) noexcept
  {
    return err <= config.flowTolerance;
  }

} // namespace

template <typename T>
inline void normalize(std::vector<T>& v, const T sum) noexcept
{
  // A zero or non-finite total has no meaningful normalization. Leaving the
  // vector untouched keeps the degenerate state recognizable (an all-zero
  // distribution) instead of turning every entry into NaN, and
  // InfomapBase::RunSession::checkFlowPostCondition -- run after initNetwork on
  // the run path -- reports it with context.
  if (sum == T {} || !std::isfinite(sum))
    return;
  for (auto& numerator : v) {
    numerator /= sum;
  }
}

template <typename T>
inline void normalize(std::vector<T>& v) noexcept
{
  const auto sum = std::accumulate(cbegin(v), cend(v), T {});
  normalize(v, sum);
}

FlowCalculator::FlowCalculator(StateNetwork& network, const Config& config)
    : numNodes(network.numNodes())
{
  // Prepare data in sequence containers for fast access of individual elements
  // Map to zero-based dense indexing
  nodeIndexMap.reserve(numNodes);
  nodeFlow.assign(numNodes, 0.0);
  nodeTeleportWeights.assign(numNodes, 0.0); // Fraction of teleportation flow landing on node i

  nodeOutDegree.assign(numNodes, 0);
  sumLinkOutWeight.assign(numNodes, 0.0);

  unsigned int nodeIndex = 0;
  network.ensureFinalized(); // CSR is the consumed link store; build it if needed

  if (network.isBipartite()) {
    // Preserve node order
    for (const auto& node : network.nodes()) {
      const auto nodeId = node.second.id;
      nodeIndexMap[nodeId] = nodeIndex++;
    }

    auto bipartiteStartId = network.bipartiteStartId();
    bipartiteStartIndex = nodeIndexMap[bipartiteStartId];

  } else {
    if (config.flowModel != FlowModel::directed) {
      // Preserve node order
      for (const auto& node : network.nodes()) {
        const auto nodeId = node.second.id;
        nodeIndexMap[nodeId] = nodeIndex++;
      }
    } else {
      // Store dangling nodes out-of-order,
      // with dangling nodes first to optimize calculation of dangling rank

      for (const auto& node : network.nodes()) {
        const auto isDangling = network.isDangling(network.indexOfId(node.second.id));
        if (!isDangling) continue;

        const auto& nodeId = node.second.id;
        nodeIndexMap[nodeId] = nodeIndex++;
      }

      nonDanglingStartIndex = nodeIndex;

      for (const auto& node : network.nodes()) {
        const auto isDangling = network.isDangling(network.indexOfId(node.second.id));
        if (isDangling) continue;

        const auto& nodeId = node.second.id;
        nodeIndexMap[nodeId] = nodeIndex++;
      }
    }
  }

  flowLinks.resize(network.numLinks(), { 0, 0, 0.0 });
  sumLinkWeight = network.sumLinkWeight();
  sumWeightedDegree = network.sumWeightedDegree();

  if (network.isBipartite()) {
    const auto bipartiteStartId = network.bipartiteStartId();

    for (unsigned int s = 0; s < numNodes; ++s) {
      const auto sourceIsFeature = network.nodeId(s) >= bipartiteStartId;
      if (sourceIsFeature) continue;
      bipartiteLinkStartIndex += network.outDegree(s);
    }
  }

  unsigned int linkIndex = 0;
  unsigned int featureLinkIndex = bipartiteLinkStartIndex; // bipartite case
  double undirectedLinkNormalization = 2 * sumLinkWeight - network.sumSelfLinkWeight();

  // CSR iterates links in (source, target) id order, identical to the nested
  // map, so flowLinks indices match. Translate CSR index -> id -> FlowCalculator
  // internal index (the directed model reorders dangling nodes first).
  network.forEachLink([&](unsigned int srcIdx, unsigned int tgtIdx, double linkWeight, double&) {
    const auto sourceId = network.nodeId(srcIdx);
    const auto sourceIndex = nodeIndexMap[sourceId];
    const auto targetIndex = nodeIndexMap[network.nodeId(tgtIdx)];

    ++nodeOutDegree[sourceIndex];
    sumLinkOutWeight[sourceIndex] += linkWeight;
    nodeFlow[sourceIndex] += linkWeight / undirectedLinkNormalization;

    if (network.isBipartite() && sourceId >= network.bipartiteStartId()) {
      // Link from feature node to ordinary node
      flowLinks[featureLinkIndex].source = sourceIndex;
      flowLinks[featureLinkIndex].target = targetIndex;
      flowLinks[featureLinkIndex].flow = linkWeight;
      ++featureLinkIndex;
    } else {
      // Ordinary link, or unipartite
      flowLinks[linkIndex].source = sourceIndex;
      flowLinks[linkIndex].target = targetIndex;
      flowLinks[linkIndex].flow = linkWeight;
      ++linkIndex;
    }

    if (sourceIndex != targetIndex) {
      if (config.isUndirectedFlow()) {
        ++nodeOutDegree[targetIndex];
        sumLinkOutWeight[targetIndex] += linkWeight;
      }
      if (config.flowModel != FlowModel::outdirdir) {
        nodeFlow[targetIndex] += linkWeight / undirectedLinkNormalization;
      }
    }
  });

  // The directed model normally finds the dangling nodes as a prefix of nodeFlow,
  // having ordered them first above. Bipartite input keeps the network's node order,
  // since the feature side is identified by an index range, so nonDanglingStartIndex
  // stays 0 and their positions have to be collected here instead -- out-degree is
  // only known now that the links have been counted.
  if (network.isBipartite() && config.flowModel == FlowModel::directed) {
    for (unsigned int i = 0; i < numNodes; ++i) {
      if (nodeOutDegree[i] == 0) {
        danglingIndices.push_back(i);
      }
    }
  }

  bool normalizeNodeFlow = false;

  switch (config.flowModel) {
  case FlowModel::undirected:
    if (config.regularized) {
      if (config.isMultilayerNetwork()) {
#if INFOMAP_FEATURE_REGULARIZED_MULTILAYER
        calcUndirectedRegularizedMultilayerFlow(network, config);
#else
        throw std::runtime_error(regularizedMultilayerFeatureError());
#endif
      } else {
        calcUndirectedRegularizedFlow(network, config);
      }
    } else {
      calcUndirectedFlow();
    }
    break;
  case FlowModel::directed:
    if (network.isBipartite() && config.bipartiteTeleportation) {
      // The two-step walk needs a step back: the option's second half sends flow from
      // the feature side to the primary side along explicit feature -> primary links.
      // A canonical bipartite file has none -- every link runs primary -> feature -- so
      // that half is a no-op and what comes out is not this flow model. Depending on the
      // teleportation flags it was zero flow, NaN, or, with --to-nodes, the uniform
      // teleport distribution with a plausible codelength and no warning: flow that
      // ignores the network entirely (#899). Refused with the precondition named.
      if (bipartiteLinkStartIndex == flowLinks.size()) {
        throw InfomapError(ExitCode::InputError,
                           fmt::format(FMT_STRING("--bipartite-teleportation needs links from the feature side back to the "
                                                  "primary side, and all {} links in this network run from the primary side "
                                                  "to the feature side. Without them the second step of the two-step walk "
                                                  "has nothing to follow and the resulting flow does not describe the "
                                                  "network. Drop --bipartite-teleportation to use the two-step projection "
                                                  "instead, which handles this input."),
                                       flowLinks.size()));
      }
      calcDirectedBipartiteFlow(network, config);
    } else {
      if (config.regularized) {
        if (config.isMultilayerNetwork()) {
#if INFOMAP_FEATURE_REGULARIZED_MULTILAYER
          calcDirectedRegularizedMultilayerFlow(network, config);
#else
          throw std::runtime_error(regularizedMultilayerFeatureError());
#endif
        } else {
          calcDirectedRegularizedFlow(network, config);
        }
      } else if (config.isMultilayerNetwork() && config.multilayerRelaxToSelf) {
        calcDirectedRelaxToSelfFlow(network, config);
      } else {
        calcDirectedFlow(network, config);
      }
    }
    break;
  case FlowModel::undirdir:
  case FlowModel::outdirdir:
    calcDirdirFlow(config);
    normalizeNodeFlow = true;
    break;
  case FlowModel::rawdir:
    calcRawdirFlow();
    normalizeNodeFlow = true;
    break;
  case FlowModel::precomputed:
    usePrecomputedFlow(network, config);
    normalizeNodeFlow = true;
    break;
  }

  finalize(network, config, normalizeNodeFlow);
}

void FlowCalculator::addFlowNote(const std::string& note)
{
  m_flowNotes.push_back(note);
}

void FlowCalculator::recordPageRank(unsigned int iterations, double error, bool converged) noexcept
{
  m_havePageRank = true;
  m_pageRankIterations = iterations;
  m_pageRankError = error;
  m_pageRankConverged = converged;
}

void FlowCalculator::calcUndirectedFlow() noexcept
{
  m_flowMethod = "undirected links";

  // Flow is outgoing transition probability times source node flow
  // = w_ij / s_ij * s_ij / sum(s_ij) = w_ij / sum(s_ij)
  // Count twice for non-loops to cover flow in both directions
  // Assuming convention to treat self-links as directed

  for (auto& link : flowLinks) {
    link.flow /= sumWeightedDegree;
    if (link.source != link.target) {
      link.flow *= 2;
    }
  }
}

void FlowCalculator::calcDirdirFlow(const Config& config) noexcept
{
  m_flowMethod = config.flowModel == FlowModel::outdirdir ? "ingoing links only" : "undirected links, directed steady state";

  // Take one last power iteration
  const std::vector<double> nodeFlowSteadyState(nodeFlow);
  nodeFlow.assign(numNodes, 0.0);

  for (const auto& link : flowLinks) {
    nodeFlow[link.target] += nodeFlowSteadyState[link.source] * link.flow / sumLinkOutWeight[link.source];
  }

  double sumNodeFlow = std::accumulate(cbegin(nodeFlow), cend(nodeFlow), 0.0);

  // Update link data to represent flow instead of weight
  for (auto& link : flowLinks) {
    link.flow *= nodeFlowSteadyState[link.source] / sumLinkOutWeight[link.source] / sumNodeFlow;
  }
}

void FlowCalculator::calcRawdirFlow() noexcept
{
  m_flowMethod = "directed links with raw flow";
  addFlowNote(fmt::format(FMT_STRING("Total link weight {:g}"), sumLinkWeight));

  // Treat the link weights as flow (after global normalization) and
  // do one power iteration to set the node flow
  nodeFlow.assign(numNodes, 0.0);

  for (auto& link : flowLinks) {
    link.flow /= sumLinkWeight;
    nodeFlow[link.target] += link.flow;
  }
}

void FlowCalculator::usePrecomputedFlow(const StateNetwork& network, const Config&)
{
  m_flowMethod = "precomputed directed flow";
  addFlowNote(fmt::format(FMT_STRING("Total link flow {:g}"), sumLinkWeight));

  if (network.haveFileInput()) {
    if (network.haveMemoryInput() && !network.haveStateNodeWeights()) {
      throw std::runtime_error("Missing node flow in input data. Should be passed as a third field under a *States section.");
    }
    if (!network.haveMemoryInput() && !network.haveNodeWeights()) {
      throw std::runtime_error("Missing node flow in input data. Should be passed as a third field under a *Vertices section.");
    }
  }

  // Treat the link weights as flow
  nodeFlow.assign(numNodes, 0.0);
  double sumFlow = 0.0;

  for (const auto& nodeIt : network.nodes()) {
    auto& node = nodeIt.second;
    nodeFlow[nodeIndexMap[node.id]] = node.weight;
    sumFlow += node.weight;
  }
  addFlowNote(fmt::format(FMT_STRING("Total node flow {:g}"), sumFlow));

  if (infomath::isEqual(sumFlow, 0)) {
    throw std::runtime_error("Missing node flow. Set it on the node weight property.");
  }
  if (!infomath::isEqual(sumFlow, 1)) {
    if (infomath::isEqual(sumFlow, numNodes) && infomath::isEqual(nodeFlow[0], 1)) {
      Log() << "\n";
      Console::warn(0, "Node flow sums to the number of nodes; was node flow provided, or are default node weights being used? Normalizing.");
    } else {
      Log() << "\n";
      Console::warn(0, "Node flow sums to {:g}, normalizing.", sumFlow);
    }
    for (unsigned int i = 0; i < numNodes; ++i) {
      nodeFlow[i] /= sumFlow;
    }
  }
}

struct IterationResult {
  double alpha;
  double beta;
  unsigned int iterations;
  double error;
  bool converged;
};

template <typename Iteration>
IterationResult powerIterate(const Config& config, double alpha, Iteration&& iter)
{
  unsigned int iterations = 0;
  double beta = 1.0 - alpha;
  double err = 0.0;

  do {
    double oldErr = err;
    err = iter(iterations, alpha, beta);

    // Perturb the system if equilibrium
    if (std::abs(err - oldErr) < 1e-17) {
      alpha += 1.0e-12;
      beta = 1.0 - alpha;
    }

    ++iterations;
  } while (iterations < config.maxFlowIterations && (err > config.flowTolerance || iterations < config.minFlowIterations));

  const bool converged = errorWithinFlowTolerance(config, err);
  if (!converged) {
    Log() << "\n";
    Console::warn(0, "PageRank calculation did not converge after {} iterations with error {:g}.", iterations, err);
  }

  return { alpha, beta, iterations, err, converged };
}

double FlowCalculator::accumulateDanglingRank() const noexcept
{
  // Two ways of locating the dangling nodes, one fast and one general. The directed
  // model orders them first, making their flow a prefix; bipartite input cannot use
  // that ordering, so their positions were collected explicitly. An empty list means
  // the prefix form applies -- including the case of no dangling nodes at all, where
  // both forms sum to zero.
  if (danglingIndices.empty()) {
    return std::accumulate(cbegin(nodeFlow), cbegin(nodeFlow) + nonDanglingStartIndex, 0.0);
  }

  double sum = 0.0;
  for (const auto i : danglingIndices) {
    sum += nodeFlow[i];
  }
  return sum;
}

void FlowCalculator::calcDirectedFlow(const StateNetwork& network, const Config& config) noexcept
{
  m_flowMethod = "directed links";
  m_teleportation = fmt::format(FMT_STRING("{}, to {}"), config.recordedTeleportation ? "recorded" : "unrecorded", config.teleportToNodes ? "nodes" : "links");

  // Calculate the teleport rate distribution
  if (config.teleportToNodes) {
    double sumNodeWeights = 0.0;

    for (const auto& nodeIt : network.nodes()) {
      auto& node = nodeIt.second;
      nodeTeleportWeights[nodeIndexMap[node.id]] = node.weight;
      sumNodeWeights += node.weight;
    }

    normalize(nodeTeleportWeights, sumNodeWeights);
  } else {
    // Teleport to links

    // Teleport proportionally to out-degree, or in-degree if recorded teleportation.
    for (const auto& link : flowLinks) {
      auto toNode = config.recordedTeleportation ? link.target : link.source;
      nodeTeleportWeights[toNode] += link.flow / sumLinkWeight;
    }
  }

  // Normalize link weights with respect to its source nodes total out-link weight;
  for (auto& link : flowLinks) {
    if (sumLinkOutWeight[link.source] > 0) {
      link.flow /= sumLinkOutWeight[link.source];
    }
  }

  std::vector<double> nodeFlowTmp(numNodes, 0.0);
  double danglingRank;

  // Calculate PageRank
  const auto iteration = [&](const auto iter, const double alpha, const double beta) {
    danglingRank = accumulateDanglingRank();

    // Flow from teleportation
    const auto teleportationFlow = alpha + beta * danglingRank;
    for (unsigned int i = 0; i < numNodes; ++i) {
      nodeFlowTmp[i] = teleportationFlow * nodeTeleportWeights[i];
    }

    // Flow from links
    for (const auto& link : flowLinks) {
      nodeFlowTmp[link.target] += beta * link.flow * nodeFlow[link.source];
    }

    // Update node flow from the power iteration above and check if converged
    double nodeFlowDiff = -1.0; // Start with -1.0 so we don't have to subtract it later
    double error = 0.0;
    for (unsigned int i = 0; i < numNodes; ++i) {
      nodeFlowDiff += nodeFlowTmp[i];
      error += std::abs(nodeFlowTmp[i] - nodeFlow[i]);
    }

    nodeFlow = nodeFlowTmp;

    // Normalize if needed
    if (std::abs(nodeFlowDiff) > 1.0e-10) {
      Console::detail(1, "normalizing flow after {} power iterations with error {:g}", iter, nodeFlowDiff);
      normalize(nodeFlow, nodeFlowDiff + 1.0);
    }

    return error;
  };

  const auto result = powerIterate(config, config.teleportationProbability, iteration);
  recordPageRank(result.iterations, result.error, result.converged);

  double sumNodeRank = 1.0;
  double beta = result.beta;

  if (!config.recordedTeleportation) {
    // Take one last power iteration excluding the teleportation
    // and normalize node flow
    sumNodeRank = 1.0 - danglingRank;
    nodeFlow.assign(numNodes, 0.0);

    for (const auto& link : flowLinks) {
      nodeFlow[link.target] += link.flow * nodeFlowTmp[link.source] / sumNodeRank;
    }

    beta = 1.0;
  }

  // Update the links with their global flow from the PageRank values.
  // Note: beta is set to 1 if unrecorded teleportation
  for (auto& link : flowLinks) {
    link.flow *= beta * nodeFlowTmp[link.source] / sumNodeRank;
  }
}

void FlowCalculator::calcDirectedRelaxToSelfFlow(const StateNetwork& network, const Config& config) noexcept
{
  // Node flow for --multilayer-relax-to-self. This mirrors calcDirectedFlow but
  // uses a two-step transition: an intra-layer link is an ordinary one-hop, while
  // an inter-layer link (to the same physical node in another layer) is transient
  // -- its flow is relayed through the target's intra-layer out-links (the deferred
  // relax intra-step) within the same iteration, so the inter-layer target accrues
  // no visit. The fused inter+intra step equals the default spread model's
  // transition, so the stationary node flow is exactly spread's on the compact
  // O(L*k) network -- the same two-step trick as calcDirectedBipartiteFlow. The
  // matching link flow is produced in finalize().
  m_flowMethod = "directed multilayer relax-to-self (two-step)";
  m_teleportation = fmt::format(FMT_STRING("{}, to {}"), config.recordedTeleportation ? "recorded" : "unrecorded", config.teleportToNodes ? "nodes" : "links");

  // Calculate the teleport rate distribution
  if (config.teleportToNodes) {
    double sumNodeWeights = 0.0;

    for (const auto& nodeIt : network.nodes()) {
      auto& node = nodeIt.second;
      nodeTeleportWeights[nodeIndexMap[node.id]] = node.weight;
      sumNodeWeights += node.weight;
    }

    normalize(nodeTeleportWeights, sumNodeWeights);
  } else {
    // Teleport to links

    // Teleport proportionally to out-degree, or in-degree if recorded teleportation.
    for (const auto& link : flowLinks) {
      auto toNode = config.recordedTeleportation ? link.target : link.source;
      nodeTeleportWeights[toNode] += link.flow / sumLinkWeight;
    }
  }

  // Normalize link weights with respect to its source nodes total out-link weight;
  for (auto& link : flowLinks) {
    if (sumLinkOutWeight[link.source] > 0) {
      link.flow /= sumLinkOutWeight[link.source];
    }
  }

  // Classify links (inter-layer vs intra-layer) and accumulate per-node intra-layer out-mass.
  std::vector<unsigned int> physId(numNodes, 0);
  for (const auto& nodeIt : network.nodes()) {
    physId[nodeIndexMap[nodeIt.second.id]] = nodeIt.second.physicalId;
  }
  std::vector<char> isInterLayer(flowLinks.size(), 0);
  std::vector<double> intraOutSum(numNodes, 0.0);
  for (unsigned int k = 0; k < flowLinks.size(); ++k) {
    const auto& link = flowLinks[k];
    if (physId[link.source] == physId[link.target] && link.source != link.target) {
      isInterLayer[k] = 1;
    } else {
      intraOutSum[link.source] += link.flow;
    }
  }

  std::vector<double> nodeFlowTmp(numNodes, 0.0);
  std::vector<double> interArrived(numNodes, 0.0);
  double danglingRank;

  // One fused transition step. First pass: every intra-layer link makes the ordinary
  // one-hop dst += beta * P(src->dst) * src[src], while flow arriving on an inter-layer
  // link is held back in interArrived (its target is not visited). Second pass: push
  // that held-back flow on through the target's intra-layer out-links, split by their
  // transition probability -- the deferred relax intra-step.
  const auto twoStep = [&](const double beta, const std::vector<double>& src, std::vector<double>& dst) {
    std::fill(interArrived.begin(), interArrived.end(), 0.0);
    for (unsigned int k = 0; k < flowLinks.size(); ++k) {
      const auto& link = flowLinks[k];
      if (isInterLayer[k]) {
        interArrived[link.target] += beta * link.flow * src[link.source];
      } else {
        dst[link.target] += beta * link.flow * src[link.source];
      }
    }
    for (unsigned int k = 0; k < flowLinks.size(); ++k) {
      const auto& link = flowLinks[k];
      if (!isInterLayer[k] && intraOutSum[link.source] > 0.0) {
        dst[link.target] += interArrived[link.source] * link.flow / intraOutSum[link.source];
      }
    }
  };

  const auto iteration = [&](const auto iter, const double alpha, const double beta) {
    danglingRank = accumulateDanglingRank();
    const auto teleportationFlow = alpha + beta * danglingRank;
    for (unsigned int i = 0; i < numNodes; ++i) {
      nodeFlowTmp[i] = teleportationFlow * nodeTeleportWeights[i];
    }
    twoStep(beta, nodeFlow, nodeFlowTmp);

    double nodeFlowDiff = -1.0;
    double error = 0.0;
    for (unsigned int i = 0; i < numNodes; ++i) {
      nodeFlowDiff += nodeFlowTmp[i];
      error += std::abs(nodeFlowTmp[i] - nodeFlow[i]);
    }
    nodeFlow = nodeFlowTmp;
    if (std::abs(nodeFlowDiff) > 1.0e-10) {
      Console::detail(1, "normalizing flow after {} power iterations with error {:g}", iter, nodeFlowDiff);
      normalize(nodeFlow, nodeFlowDiff + 1.0);
    }
    return error;
  };

  const auto result = powerIterate(config, config.teleportationProbability, iteration);
  recordPageRank(result.iterations, result.error, result.converged);

  double sumNodeRank = 1.0;
  double beta = result.beta;

  if (!config.recordedTeleportation) {
    // Take one last (un-teleported) two-step iteration and normalize node flow.
    sumNodeRank = 1.0 - danglingRank;
    nodeFlow.assign(numNodes, 0.0);
    twoStep(1.0, nodeFlowTmp, nodeFlow);
    if (sumNodeRank > 0.0) {
      normalize(nodeFlow, sumNodeRank);
    }
    beta = 1.0;
  }

  // Update the links with their global flow from the PageRank values.
  // (Inter-layer links get their transit flow here; finalize() also relays it onto
  // the target's intra-layer links for the matching two-step link flow.)
  for (auto& link : flowLinks) {
    link.flow *= beta * nodeFlowTmp[link.source] / sumNodeRank;
  }
}

void FlowCalculator::calcDirectedRegularizedFlow(const StateNetwork& network, const Config& config) noexcept
{
  m_flowMethod = "directed regularized flow";
  m_teleportation = "recorded, Bayesian prior to nodes";

  // Calculate node weights w_i = s_i/k_i, where s_i is the node strength (weighted degree) and k_i the (unweighted) degree
  unsigned int N = network.numNodes();

  std::vector<unsigned int> k_out(N, 0);
  std::vector<unsigned int> k_in(N, 0);
  std::vector<double> s_out(N, 0);
  std::vector<double> s_in(N, 0);
  double sum_s = sumWeightedDegree;
  unsigned int sum_k = network.sumDegree();
  if (sum_k == 0) {
    const auto uniformFlow = 1.0 / N;
    std::fill(nodeFlow.begin(), nodeFlow.end(), uniformFlow);
    std::fill(nodeTeleportWeights.begin(), nodeTeleportWeights.end(), uniformFlow);
    nodeTeleportFlow.assign(numNodes, uniformFlow);
    return;
  }
  double average_weight = sum_s / sum_k;

  for (auto& link : flowLinks) {
    k_out[link.source] += 1;
    s_out[link.source] += link.flow;
    k_in[link.target] += 1;
    s_in[link.target] += link.flow;
  }

  double min_u_out = std::numeric_limits<double>::max();
  double min_u_in = std::numeric_limits<double>::max();
  for (unsigned int i = 0; i < N; ++i) {
    if (k_out[i] > 0) {
      min_u_out = std::min(min_u_out, s_out[i] / k_out[i]);
    }
    if (k_in[i] > 0) {
      min_u_in = std::min(min_u_in, s_in[i] / k_in[i]);
    }
  }

  auto u_out = [&s_out, &k_out, min_u_out](auto i) { return k_out[i] == 0 ? min_u_out : s_out[i] / k_out[i]; };
  auto u_in = [&s_in, &k_in, min_u_in](auto i) { return k_in[i] == 0 ? min_u_in : s_in[i] / k_in[i]; };

  unsigned int numNodesAsTeleportationTargets = config.noSelfLinks ? N - 1 : N;
  double lambda = config.regularizationStrength * std::log(N) / numNodesAsTeleportationTargets;
  if (network.numPhysicalNodes() < N) {
    // Adjust for higher-order network, lnN/N^2 recovers lnN/N on physical network, but divide by N again because physically unconstrained teleportation
    unsigned int Np = network.numPhysicalNodes();
    lambda = config.regularizationStrength * std::log(Np) / (Np * Np * Np);
  }
  double u_t = average_weight;

  double sum_u_in = 0.0;
  for (unsigned int i = 0; i < N; ++i) {
    sum_u_in += u_in(i);
  }

  for (unsigned int i = 0; i < N; ++i) {
    nodeTeleportWeights[i] = u_in(i) / sum_u_in;
  }

  std::function<double(unsigned int)> t_out_withoutSelfLinks = [lambda, u_t, u_out, u_in, sum_u_in](unsigned int i) { return lambda / u_t * u_out(i) * (sum_u_in - u_in(i)); };
  std::function<double(unsigned int)> t_out_withSelfLinks = [lambda, u_t, u_out, sum_u_in](unsigned int i) { return lambda / u_t * u_out(i) * sum_u_in; };
  auto t_out = config.noSelfLinks ? t_out_withoutSelfLinks : t_out_withSelfLinks;

  std::vector<double> alpha(N, 0);
  for (unsigned int i = 0; i < N; ++i) {
    auto t_i = t_out(i);
    alpha[i] = t_i / (s_out[i] + t_i); // = 1 for dangling nodes
    if (config.noSelfLinks) {
      // Inflate to adjust for no self-teleportation
      // TODO: Check possible side-effects
      alpha[i] /= 1 - nodeTeleportWeights[i];
    }
  }

  // Normalize link weights with respect to its source nodes total out-link weight;
  for (auto& link : flowLinks) {
    if (sumLinkOutWeight[link.source] > 0) {
      link.flow /= sumLinkOutWeight[link.source];
    }
  }

  std::vector<double> nodeFlowTmp(numNodes, 0.0);

  // Calculate PageRank
  const auto iteration = [&](const auto iter) {
    double teleTmp = 0.0;
    for (unsigned int i = 0; i < N; ++i) {
      teleTmp += alpha[i] * nodeFlow[i];
    }

    for (unsigned int i = 0; i < N; ++i) {
      nodeFlowTmp[i] = nodeTeleportWeights[i] * (config.noSelfLinks ? (teleTmp - alpha[i] * nodeFlow[i]) : teleTmp);
    }

    // Flow from links
    for (const auto& link : flowLinks) {
      double beta = 1 - alpha[link.source] * (config.noSelfLinks ? 1 - nodeTeleportWeights[link.source] : 1);
      nodeFlowTmp[link.target] += beta * link.flow * nodeFlow[link.source];
    }

    // Update node flow from the power iteration above and check if converged
    double nodeFlowDiff = -1.0; // Start with -1.0 so we don't have to subtract it later
    double error = 0.0;
    for (unsigned int i = 0; i < numNodes; ++i) {
      nodeFlowDiff += nodeFlowTmp[i];
      error += std::abs(nodeFlowTmp[i] - nodeFlow[i]);
    }

    nodeFlow = nodeFlowTmp;

    // Normalize if needed
    if (std::abs(nodeFlowDiff) > 1.0e-10) {
      Console::detail(1, "normalizing ranks after {} power iterations with error {:g}", iter, nodeFlowDiff);
      normalize(nodeFlow, nodeFlowDiff + 1.0);
    }

    return error;
  };

  unsigned int iterations = 0;
  double err = 0.0;

  do {
    err = iteration(iterations);

    ++iterations;
  } while (iterations < config.maxFlowIterations && (err > config.flowTolerance || iterations < config.minFlowIterations));

  const bool converged = errorWithinFlowTolerance(config, err);
  recordPageRank(iterations, err, converged);
  if (!converged) {
    Log() << "\n";
    Console::warn(0, "PageRank calculation did not converge after {} iterations with error {:g}.", iterations, err);
  }

  double sumNodeRank = 1.0;
  for (auto& link : flowLinks) {
    double beta = 1 - alpha[link.source] * (config.noSelfLinks ? 1 - nodeTeleportWeights[link.source] : 1);
    link.flow *= beta * nodeFlow[link.source] / sumNodeRank;
  }

  nodeTeleportFlow.assign(numNodes, 0.0);
  for (unsigned int i = 0; i < N; ++i) {
    nodeTeleportFlow[i] = nodeFlow[i] * alpha[i];
  }
}

#if INFOMAP_FEATURE_REGULARIZED_MULTILAYER
void FlowCalculator::calcDirectedRegularizedMultilayerFlow(const StateNetwork& network, const Config& config) noexcept
{
  // Calculate node weights w_i = s_i/k_i, where s_i is the node strength (weighted degree) and k_i the (unweighted) degree
  unsigned int N = network.numNodes();
  unsigned int N_phys = network.numPhysicalNodes();
  unsigned int L = network.numLayers();
  // unsigned int N_states = network.numNodes();
  // double nodeWeight = 1.0 / N;
  double intraOutWeight = config.regularizationStrength * std::log(N_phys) / L;
  // double interOutWeight = config.regularizationStrength * std::log(L);

  // Log(1) << "\n N: " << N_phys << ", N_states: " << N << ", L: " << L << "\n";
  // Log(1) << "ln(N)/(NL): " << std::log(N_phys) / (N_phys * L) << "\n";
  // Log(1) << "ln(N)/(L): " << std::log(N_phys) / (L) << "\n";

  std::vector<unsigned int> layerIds(N, 0);
  std::vector<unsigned int> physicalIds(N, 0);

  std::unordered_map<unsigned int, unsigned int> layerIdToIndex;
  unsigned int layerIndex = 0;
  for (unsigned int layerId : network.layers()) {
    layerIdToIndex[layerId] = layerIndex++;
    // Log(1) << "Layer " << layerId << " -> index " << layerIdToIndex[layerId] << "\n";
  }
  // Log(1) << "\n -> " << layerIdToIndex.size() << " layers...\n";

  std::vector<bool> isInterLink(flowLinks.size(), false);
  std::vector<unsigned int> layerIndices(N);

  for (const auto& node : network.nodes()) {
    const auto nodeIndex = nodeIndexMap[node.second.id];
    layerIds[nodeIndex] = node.second.layerId;
    physicalIds[nodeIndex] = node.second.physicalId;
    layerIndices[nodeIndex] = layerIdToIndex[node.second.layerId];
    // nodeTeleportWeights[nodeIndexMap[nodeId]] = node.weight;
    // if (layerIdToIndex.count(node.second.layerId) == 0) {
    //   layerIdToIndex[node.second.layerId] = layerIndex++;
    // }
    // Log(1) << "Node (physId: " << node.second.physicalId << ", layerId: " << node.second.layerId << ") -> index: " << nodeIndexMap[node.second.id] << "\n";
  }

  unsigned int linkIndex = 0;

  // Log(1) << "\nLinks:\n";

  for (const auto& link : flowLinks) {
    isInterLink[linkIndex] = physicalIds[link.source] == physicalIds[link.target];
    // Log(1) << linkIndex << ": (" << layerIds[link.source] << "," << physicalIds[link.source] << ") -> (" << layerIds[link.target] << "," << physicalIds[link.target] << ") is inter: " << isInterLink[linkIndex] << "\n";
    ++linkIndex;
  }

  std::vector<unsigned int> k_out(N, 0);
  std::vector<unsigned int> k_in(N, 0);
  std::vector<double> s_out(N, 0);
  std::vector<double> s_in(N, 0);
  std::vector<double> inter_out(N, 0);
  // double sum_s = sumWeightedDegree;
  // unsigned int sum_k = network.sumDegree();
  // double average_weight = sum_s / sum_k;

  linkIndex = 0;
  for (auto& link : flowLinks) {
    if (isInterLink[linkIndex++]) {
      inter_out[link.source] += link.flow;
    } else {
      k_out[link.source] += 1;
      s_out[link.source] += link.flow;
      k_in[link.target] += 1;
      s_in[link.target] += link.flow;
      // if (link.source == 0) {
      //   Log(1) << link.source << " -> " << link.target << " => k_out[0] -> " << k_out[link.source] << "\n";
      // }
    }
  }

  // auto gamma = [s_out, intraOutWeight, interOutWeight](auto i) { return 1 + interOutWeight / (s_out[i] + intraOutWeight); };

  // double min_u_out = std::numeric_limits<double>::max();
  // double min_u_in = std::numeric_limits<double>::max();
  // for (unsigned int i = 0; i < N; ++i) {
  //   if (k_out[i] > 0) {
  //     min_u_out = std::min(min_u_out, s_out[i] / k_out[i]);
  //   }
  //   if (k_in[i] > 0) {
  //     min_u_in = std::min(min_u_in, s_in[i] / k_in[i]);
  //   }
  // }

  // auto u_out = [s_out, k_out, min_u_out](auto i) { return k_out[i] == 0 ? min_u_out : s_out[i] / k_out[i]; };
  // auto u_in = [s_in, k_in, min_u_in](auto i) { return k_in[i] == 0 ? min_u_in : s_in[i] / k_in[i]; };

  // unsigned int numNodesAsTeleportationTargets = config.noSelfLinks ? N - 1 : N;
  // double lambda = config.regularizationStrength * std::log(N) / numNodesAsTeleportationTargets;
  // double u_t = average_weight;

  // double sum_u_in = 0.0;
  // for (unsigned int i = 0; i < N; ++i) {
  //   sum_u_in += u_in(i);
  // }

  // Log(1) << "\nNodes:\n";
  for (unsigned int i = 0; i < N; ++i) {
    // nodeTeleportWeights[i] = u_in(i) / sum_u_in;
    // nodeTeleportWeights[i] = 1 / N;
    nodeTeleportWeights[i] = 1.0 / N_phys;
    // nodeTeleportWeights[i] = 1.0 / (config.noSelfLinks ? N_phys - 1 : N_phys);
    // Log(1) << i << ") k_out: " << k_out[i] << ", s_out: " << s_out[i] << ", inter_out: " << inter_out[i] << "\n";
  }

  // std::function<double(unsigned int)> t_out_withoutSelfLinks = [lambda, u_t, u_out, u_in, sum_u_in](unsigned int i) { return lambda / u_t * u_out(i) * (sum_u_in - u_in(i)); };
  // std::function<double(unsigned int)> t_out_withSelfLinks = [lambda, u_t, u_out, sum_u_in](unsigned int i) { return lambda / u_t * u_out(i) * sum_u_in; };
  // auto t_out = config.noSelfLinks ? t_out_withoutSelfLinks : t_out_withSelfLinks;

  auto intraLayerTeleRate = [s_out, k_out, intraOutWeight](auto i) { return k_out[i] == 0 ? 1 : intraOutWeight / (intraOutWeight + s_out[i]); };

  std::vector<double> alpha(N, 0);
  std::vector<double> alphaInter(N, 0);
  // Log(1) << "\nTele probabilities:\n";
  for (unsigned int i = 0; i < N; ++i) {
    // auto t_i = t_out(i);
    alpha[i] = intraLayerTeleRate(i); // = 1 for dangling nodes
    alphaInter[i] = inter_out[i] / (inter_out[i] + s_out[i] + intraOutWeight);
    if (config.noSelfLinks) {
      // Inflate to adjust for no self-teleportation
      // TODO: Check possible side-effects
      alpha[i] /= 1 - nodeTeleportWeights[i];
      // alphaInter[i] /= 1 - nodeTeleportWeights[i];
    }
    // Log(1) << i << ": intra: " << alpha[i] << ", inter: " << alphaInter[i] << "\n";
    // Log(1) << i << ": intra: " << alpha[i] << ", inter: " << alphaInter[i] << " (inter_out: " << inter_out[i] << ", s_out: " << s_out[i] << ", intra_prior_out: " << intraOutWeight << ")\n";
  }

  // Log(1) << "\nLink probabilities:\n";
  // Normalize link weights to probabilities, separate for intra and inter links
  linkIndex = 0;
  for (auto& link : flowLinks) {
    // if (sumLinkOutWeight[link.source] > 0) {
    //   link.flow /= sumLinkOutWeight[link.source];
    // }
    if (isInterLink[linkIndex++]) {
      link.flow /= inter_out[link.source];
    } else {
      if (k_out[link.source] > 0) {
        link.flow /= s_out[link.source];
      }
    }
    // Log(1) << link.source << " -> " << link.target << ": " << link.flow << "\n";
  }

  std::vector<double> unrecordedInterFlow(N, 0);
  std::vector<double> nodeFlowTmp(numNodes, 0.0);
  std::vector<double> layerTeleFlow(L, 0.0);

  for (unsigned int i = 0; i < N; ++i) {
    nodeFlow[i] = 1.0 / N;
  }

  // Calculate two-step PageRank:
  const auto iteration = [&](const auto iter) {
    // Log(1) << "\nIter " << iter << ":\n";

    // 1. Unrecorded inter-layer step: push fraction of flow on inter-layer links to temporary location
    linkIndex = 0;
    unrecordedInterFlow.assign(N, 0.0);
    for (auto& link : flowLinks) {
      if (!isInterLink[linkIndex++]) {
        continue;
      }
      unrecordedInterFlow[link.target] += alphaInter[link.source] * nodeFlow[link.source] * link.flow;
      // unrecordedInterFlow[link.target] += alphaInter[link.source] * nodeFlow[link.source] * link.flow * (config.noSelfLinks ? 1 - nodeTeleportWeights[link.source] : 1);
      // Log(1) << "  " << link.source << " -> " << link.target << ": unrecorded[" << link.target << "] += " << nodeFlow[link.source] << " * " << alphaInter[link.source] << " * " << link.flow << "\n";
    }

    // double sumFlow = 0.0;
    // double sumUnrecordedFlow = 0.0;
    // for (unsigned int i = 0; i < N; ++i) {
    //   sumFlow += nodeFlow[i];
    //   sumUnrecordedFlow += unrecordedInterFlow[i];
    // }
    // Log(1) << "  Sum flow: " << sumFlow << "\n";
    // Log(1) << "  Sum unrecorded flow: " << sumUnrecordedFlow << "\n";

    // 2. Recorded intra-layer step: push rest of flow plus temporarily stored flow to intra-layer with intra-layer teleportation
    layerTeleFlow.assign(L, 0.0);
    for (unsigned int i = 0; i < N; ++i) {
      layerTeleFlow[layerIndices[i]] += alpha[i] * ((1 - alphaInter[i]) * nodeFlow[i] + unrecordedInterFlow[i]);
      // Log(1) << "Node " << i << " in layer " << layerIndices[i] << ": alpha: " << alpha[i] << ", non-inter flow: " << (1 - alphaInter[i]) * nodeFlow[i] << ", unrecorded flow: " << unrecordedInterFlow[i] << ", += " << alpha[i] * ((1 - alphaInter[i]) * nodeFlow[i] + unrecordedInterFlow[i]) << " -> " << layerTeleFlow[layerIndices[i]] << "\n";
    }

    // for (unsigned int i = 0; i < layerTeleFlow.size(); ++i) {
    //   Log(1) << "Layer " << i << " tele flow: " << layerTeleFlow[i] << "\n";
    // }

    for (unsigned int i = 0; i < N; ++i) {
      nodeFlowTmp[i] = nodeTeleportWeights[i] * (layerTeleFlow[layerIndices[i]] - (config.noSelfLinks ? (alpha[i] * nodeFlow[i]) : 0));
      // nodeFlowTmp[i] = nodeTeleportWeights[i] * layerTeleFlow[layerIndices[i]];
      // Log(1) << i << ": tele flow: " << nodeFlowTmp[i] << "\n";
    }

    // Flow from links
    linkIndex = 0;
    for (const auto& link : flowLinks) {
      if (isInterLink[linkIndex++]) {
        continue;
      }
      double beta = 1 - alpha[link.source] * (config.noSelfLinks ? 1 - nodeTeleportWeights[link.source] : 1);
      // double beta = 1 - alpha[link.source];
      nodeFlowTmp[link.target] += beta * link.flow * ((1 - alphaInter[link.source]) * nodeFlow[link.source] + unrecordedInterFlow[link.source]);
    }

    // Update node flow from the power iteration above and check if converged
    double nodeFlowDiff = -1.0; // Start with -1.0 so we don't have to subtract it later
    double error = 0.0;
    for (unsigned int i = 0; i < numNodes; ++i) {
      nodeFlowDiff += nodeFlowTmp[i];
      error += std::abs(nodeFlowTmp[i] - nodeFlow[i]);
    }

    nodeFlow = nodeFlowTmp;

    // Normalize if needed
    if (std::abs(nodeFlowDiff) > 1.0e-10) {
      Console::detail(1, "normalizing ranks after {} power iterations with error {:g}", iter, nodeFlowDiff);
      if (std::abs(nodeFlowDiff) > 1.0e-4) {
        throw std::runtime_error(fmt::format(FMT_STRING("Total flow differs from 1 by {} after {} iterations. Please report the issue.\n"), nodeFlowDiff, iter));
      }
      normalize(nodeFlow, nodeFlowDiff + 1.0);
    }

    return error;
  };

  unsigned int iterations = 0;
  double err = 0.0;

  do {
    err = iteration(iterations);

    ++iterations;
  } while (iterations < config.maxFlowIterations && (err > config.flowTolerance || iterations < config.minFlowIterations));

  // Recorded here too, so the exported outcome means the same thing for every flow
  // model: this loop warned but reported nothing, leaving the run report without a flow
  // object at all for regularized multilayer input.
  recordPageRank(iterations, err, errorWithinFlowTolerance(config, err));
  if (!errorWithinFlowTolerance(config, err)) {
    Console::warn(0, "PageRank calculation stopped after the maximum of {} iterations with diff {:g}.", iterations, err);
  }

  for (unsigned int i = 0; i < layerTeleFlow.size(); ++i) {
    sumTeleFlow += layerTeleFlow[i];
  }

  linkIndex = 0;
  enterFlow.assign(numNodes, 0.0);
  exitFlow.assign(numNodes, 0.0);
  for (auto& link : flowLinks) {
    if (isInterLink[linkIndex++]) {
      link.flow = alphaInter[link.source] * nodeFlow[link.source] * link.flow;
      // Need to add enter/exit flow to eventually collapse
      exitFlow[link.source] += link.flow;
      enterFlow[link.target] += link.flow;
    } else {
      double beta = 1 - alpha[link.source];
      link.flow = beta * link.flow * ((1 - alphaInter[link.source]) * nodeFlow[link.source] + unrecordedInterFlow[link.source]);
      exitFlow[link.source] += link.flow;
      enterFlow[link.target] += link.flow;
    }
  }

  nodeTeleportFlow.assign(numNodes, 0.0);
  for (unsigned int i = 0; i < N; ++i) {
    nodeTeleportFlow[i] = alpha[i] * ((1 - alphaInter[i]) * nodeFlow[i] + unrecordedInterFlow[i]);

    exitFlow[i] += nodeTeleportFlow[i] * (1 - nodeTeleportWeights[i]); // + node.intraLayerTeleFlow * (1 - node.intraLayerTeleWeight);
    enterFlow[i] += (layerTeleFlow[layerIndices[i]] - nodeTeleportFlow[i]) * nodeTeleportWeights[i];
  }
}
#endif // INFOMAP_FEATURE_REGULARIZED_MULTILAYER

void FlowCalculator::calcUndirectedRegularizedFlow(const StateNetwork& network, const Config& config) noexcept
{
  m_flowMethod = "undirected regularized flow";
  m_teleportation = "recorded, Bayesian prior to nodes";

  // Calculate node weights w_i = s_i/k_i, where s_i is the node strength (weighted degree) and k_i the (unweighted) degree
  unsigned int N = network.numNodes();
  std::vector<unsigned int> k(N, 0);
  std::vector<double> s(N, 0);
  double sum_s = sumWeightedDegree;
  unsigned int sum_k = network.sumDegree();
  if (sum_k == 0) {
    const auto uniformFlow = 1.0 / N;
    std::fill(nodeFlow.begin(), nodeFlow.end(), uniformFlow);
    std::fill(nodeTeleportWeights.begin(), nodeTeleportWeights.end(), uniformFlow);
    nodeTeleportFlow.assign(numNodes, uniformFlow);
    return;
  }
  double average_weight = sum_s / sum_k;

  for (auto& link : flowLinks) {
    k[link.source] += 1;
    s[link.source] += link.flow;
    if (link.source != link.target) {
      k[link.target] += 1;
      s[link.target] += link.flow;
    }
  }

  double min_u = std::numeric_limits<double>::max();
  for (unsigned int i = 0; i < N; ++i) {
    if (k[i] > 0) {
      min_u = std::min(min_u, s[i] / k[i]);
    }
  }

  auto u = [&s, &k, min_u](auto i) { return k[i] == 0 ? min_u : s[i] / k[i]; };

  unsigned int numNodesAsTeleportationTargets = config.noSelfLinks ? N - 1 : N;
  double lambda = config.regularizationStrength * std::log(N) / numNodesAsTeleportationTargets;
  if (network.numPhysicalNodes() < N) {
    // Adjust for higher-order network, lnN/N^2 recovers lnN/N on physical network, but divide by N again because physically unconstrained teleportation
    unsigned int Np = network.numPhysicalNodes();
    lambda = config.regularizationStrength * std::log(Np) / (Np * Np * Np);
  }
  double u_t = average_weight;

  double sum_u = 0.0;
  for (unsigned int i = 0; i < N; ++i) {
    sum_u += u(i);
  }

  // nodeTeleportWeights is the fraction of teleportation flow landing on each node. This is proportional to u_in
  for (unsigned int i = 0; i < N; ++i) {
    nodeTeleportWeights[i] = u(i) / sum_u;
  }

  std::function<double(unsigned int)> t_withoutSelfLinks = [lambda, u_t, u, sum_u](unsigned int i) { return lambda / u_t * u(i) * (sum_u - u(i)); };
  std::function<double(unsigned int)> t_withSelfLinks = [lambda, u_t, u, sum_u](unsigned int i) { return lambda / u_t * u(i) * sum_u; };
  auto t = config.noSelfLinks ? t_withoutSelfLinks : t_withSelfLinks;

  std::vector<double> alpha(N, 0);
  double sum_t = 0.0;
  for (unsigned int i = 0; i < N; ++i) {
    auto t_i = t(i);
    alpha[i] = t_i / (s[i] + t_i);
    if (config.noSelfLinks) {
      // Inflate to adjust for no self-teleportation
      // TODO: No later side effects of cheating here? Need to normalize targets instead?
      alpha[i] /= 1 - nodeTeleportWeights[i];
    }
    sum_t += t_i;
  }

  for (auto& link : flowLinks) {
    if (sumLinkOutWeight[link.source] > 0) {
      link.flow /= sumLinkOutWeight[link.source];
    }
  }

  for (unsigned int i = 0; i < N; ++i) {
    nodeFlow[i] = (s[i] + t(i)) / (sum_s + sum_t);
  }

  nodeTeleportFlow.assign(numNodes, 0.0);
  for (unsigned int i = 0; i < N; ++i) {
    nodeTeleportFlow[i] = nodeFlow[i] * alpha[i];
  }

  for (auto& link : flowLinks) {
    // TODO: Side effect from inflating alpha, need real alpha here.
    double beta = 1 - alpha[link.source] * (config.noSelfLinks ? 1 - nodeTeleportWeights[link.source] : 1);
    link.flow *= beta * nodeFlow[link.source] * 2;
  }
}

#if INFOMAP_FEATURE_REGULARIZED_MULTILAYER
void FlowCalculator::calcUndirectedRegularizedMultilayerFlow(const StateNetwork& network, const Config& config)
{
  (void)network;
  (void)config;
  throw std::runtime_error("Undirected regularized multilayer flow is not implemented.");
}
#endif // INFOMAP_FEATURE_REGULARIZED_MULTILAYER

void FlowCalculator::calcDirectedBipartiteFlow(const StateNetwork& network, const Config& config) noexcept
{
  m_flowMethod = "directed bipartite links";
  m_teleportation = fmt::format(FMT_STRING("{}, to {}"), config.recordedTeleportation ? "recorded" : "unrecorded", config.teleportToNodes ? "nodes" : "links");

  const auto bipartiteStartId = network.bipartiteStartId();

  if (config.teleportToNodes) {
    for (const auto& nodeIt : network.nodes()) {
      auto& node = nodeIt.second;
      if (node.id < bipartiteStartId) {
        nodeTeleportWeights[nodeIndexMap[node.id]] = node.weight;
      }
    }
  } else {
    // Teleport proportionally to out-degree, or in-degree if recorded teleportation.
    // Two-step degree: sum of products between incoming and outgoing links from bipartite nodes

    if (config.recordedTeleportation) {
      for (auto link = begin(flowLinks) + bipartiteLinkStartIndex; link != end(flowLinks); ++link) {
        // target is an ordinary node
        nodeTeleportWeights[link->target] += link->flow;
      }
    } else {
      // Unrecorded teleportation

      for (auto link = begin(flowLinks); link != begin(flowLinks) + bipartiteLinkStartIndex; ++link) {
        // source is an ordinary node
        nodeTeleportWeights[link->source] += link->flow;
      }
    }
  }

  normalize(nodeTeleportWeights);

  nodeFlow = nodeTeleportWeights;

  // Normalize link weights with respect to its source nodes total out-link weight;
  for (auto& link : flowLinks) {
    if (sumLinkOutWeight[link.source] > 0) {
      link.flow /= sumLinkOutWeight[link.source];
    }
  }

  std::vector<double> nodeFlowTmp(numNodes, 0.0);
  double danglingRank;

  // Calculate two-step PageRank
  const auto iteration = [&](const auto iter, const double alpha, const double beta) {
    danglingRank = accumulateDanglingRank();

    // Flow from teleportation
    const auto teleportationFlow = alpha + beta * danglingRank;
    for (unsigned int i = 0; i < bipartiteStartIndex; ++i) {
      nodeFlowTmp[i] = teleportationFlow * nodeTeleportWeights[i];
    }

    for (unsigned int i = bipartiteStartIndex; i < numNodes; ++i) {
      nodeFlowTmp[i] = 0.0;
    }

    // Flow from links
    // First step
    for (auto link = begin(flowLinks); link != begin(flowLinks) + bipartiteLinkStartIndex; ++link) {
      nodeFlow[link->target] += beta * link->flow * nodeFlow[link->source];
    }

    // Second step back to primary nodes
    for (auto link = begin(flowLinks) + bipartiteLinkStartIndex; link != end(flowLinks); ++link) {
      nodeFlowTmp[link->target] += link->flow * nodeFlow[link->source];
    }

    // Update node flow from the power iteration above and check if converged
    double nodeFlowDiff = -1.0;
    double error = 0.0;
    for (unsigned int i = 0; i < bipartiteStartIndex; ++i) {
      nodeFlowDiff += nodeFlowTmp[i];
      error += std::abs(nodeFlowTmp[i] - nodeFlow[i]);
    }

    nodeFlow = nodeFlowTmp;

    // Normalize if needed
    if (std::abs(nodeFlowDiff) > 1.0e-10) {
      Console::detail(1, "normalizing ranks after {} power iterations with error {:g}", iter, nodeFlowDiff);
      normalize(nodeFlow, nodeFlowDiff + 1.0);
    }

    return error;
  };

  const auto result = powerIterate(config, config.teleportationProbability, iteration);
  recordPageRank(result.iterations, result.error, result.converged);

  double sumNodeRank = 1.0;
  double beta = result.beta;

  if (!config.recordedTeleportation) {
    // Take one last power iteration excluding the teleportation (and normalize node flow to sum 1.0)
    sumNodeRank = 1.0 - danglingRank;
    nodeFlow.assign(numNodes, 0.0);

    for (auto link = begin(flowLinks); link != begin(flowLinks) + bipartiteLinkStartIndex; ++link) {
      nodeFlowTmp[link->target] += link->flow * nodeFlowTmp[link->source];
    }
    // Second step back to primary nodes
    for (auto link = begin(flowLinks) + bipartiteLinkStartIndex; link != end(flowLinks); ++link) {
      nodeFlow[link->target] += link->flow * nodeFlowTmp[link->source];
    }

    beta = 1.0;
  }

  // Update the links with their global flow from the PageRank values.
  // Note: beta is set to 1 if unrecorded teleportation
  for (auto& link : flowLinks) {
    link.flow *= beta * nodeFlowTmp[link.source] / sumNodeRank;
  }
}

void FlowCalculator::finalize(StateNetwork& network, const Config& config, bool normalizeNodeFlow) noexcept
{
  // Hand the power iteration's outcome to the network so it outlives this calculator.
  network.m_haveFlowConvergence = m_havePageRank;
  network.m_flowConverged = m_pageRankConverged;
  network.m_flowIterations = m_pageRankIterations;
  network.m_flowError = m_pageRankError;

  // TODO: Skip bipartite flow adjustment for directed / rawdir / .. ?
  if (network.isBipartite()) {
    addFlowNote("Using bipartite links");

    if (!config.skipAdjustBipartiteFlow && !config.bipartiteTeleportation) {
      // A node that is not coded has no visit rate, and therefore no teleportation into it
      // either. Clearing the node's flow while leaving its teleport flow and weight behind
      // broke the invariant the enter/exit flow below rests on -- that a node's
      // self-teleportation teleFlow * weight is part of its flow -- so the subtraction ran
      // past zero: with --regularized on examples/networks/bipartite.net the two feature
      // nodes came out at flow 0 and enter/exit -0.02901104907 = -(0.1232969585 *
      // 0.2352941176), and every module containing one inherited a negative enter flow.
      // The two-level index codelength then evaluated to 1.1e-16, charging nothing at all
      // for entering either module (#957).
      const auto uncode = [this](unsigned int nodeIndex) {
        nodeFlow[nodeIndex] = 0.0;
        if (!nodeTeleportFlow.empty())
          nodeTeleportFlow[nodeIndex] = 0.0;
        if (!nodeTeleportWeights.empty())
          nodeTeleportWeights[nodeIndex] = 0.0;
      };

      // Only links between ordinary nodes and feature nodes in bipartite network
      // Don't code feature nodes -> distribute all flow from those to ordinary nodes
      for (auto& link : flowLinks) {
        auto sourceIsFeature = link.source >= bipartiteStartIndex;

        if (sourceIsFeature) {
          nodeFlow[link.target] += link.flow;
          uncode(link.source); // Doesn't matter if done multiple times on each node.
        } else {
          nodeFlow[link.source] += link.flow;
          uncode(link.target); // Doesn't matter if done multiple times on each node.
        }
        // TODO: Should flow double before moving to nodes, does it cancel out in normalization?

        // Markov time 2 on the full network will correspond to markov time 1 between the real nodes.
        link.flow *= 2;
      }
      // TODO: Should flow double before moving to nodes, does it cancel out in normalization?

      normalizeNodeFlow = true;

    } else if (config.bipartiteTeleportation) {
      for (auto& link : flowLinks) {
        // Markov time 2 on the full network will correspond to markov time 1 between the real nodes.
        link.flow *= 2;
      }
    }
  }

  // Link flow for --multilayer-relax-to-self. calcDirectedRelaxToSelfFlow gave the
  // nodes spread's flow with a two-step walk; the links get the matching two-step
  // (Markov-time-2) treatment so module-exit flows -- and thus the codelength -- match
  // spread. The relax flow from (i,n) crosses to layer j via the inter-layer link
  // (i,n)->(j,n) and then continues along (j,n)'s intra-layer out-links, so keep BOTH:
  // the inter-layer link keeps its own flow (the layer switch, charged to (i,n)'s
  // module when a node's copies are split across modules), and that flow is also
  // relayed onto (j,n)'s intra out-links (the onward step, charged within the target
  // layer). The codelength depends only on node flows and module-exit flows, so the
  // extra within-module relay volume is harmless. This is exact when (j,n) is
  // co-modular with its target-layer neighbours -- always for coherent partitions, and
  // for simple overlaps (e.g. one node shared between two layer-local communities). It
  // is approximate only when those neighbours are cross-community (the relay then
  // charges a second crossing the one-hop spread model does not); reproducing that
  // exactly would need spread's O(L^2*k) links. (Same Markov-time-2 idea as the
  // bipartite handling above.)
  if (config.multilayerRelaxToSelf && !network.isBipartite()) {
    std::vector<unsigned int> physId(numNodes, 0);
    for (const auto& nodeIt : network.nodes()) {
      physId[nodeIndexMap[nodeIt.second.id]] = nodeIt.second.physicalId;
    }
    std::vector<std::vector<unsigned int>> outLinks(numNodes);
    for (unsigned int k = 0; k < flowLinks.size(); ++k) {
      outLinks[flowLinks[k].source].push_back(k);
    }
    std::vector<double> delta(flowLinks.size(), 0.0);
    for (unsigned int k = 0; k < flowLinks.size(); ++k) {
      const auto& link = flowLinks[k];
      const bool interLayer = physId[link.source] == physId[link.target] && link.source != link.target;
      if (!interLayer) {
        continue;
      }
      const unsigned int t = link.target;
      double sumIntra = 0.0;
      for (const auto l : outLinks[t]) {
        if (physId[flowLinks[l].target] != physId[t]) {
          sumIntra += flowLinks[l].flow;
        }
      }
      if (sumIntra <= 0.0) {
        continue; // dangling target: leave the inter-layer link as is
      }
      const double f = link.flow;
      for (const auto l : outLinks[t]) {
        if (physId[flowLinks[l].target] != physId[t]) {
          delta[l] += f * flowLinks[l].flow / sumIntra;
        }
      }
      // The inter-layer link keeps its own flow (the layer switch); it is not dropped.
    }
    for (unsigned int k = 0; k < flowLinks.size(); ++k) {
      flowLinks[k].flow += delta[k];
    }
  }

  if (config.useNodeWeightsAsFlow) {
    addFlowNote("Using node weights as flow");

    for (auto& nodeIt : network.nodes()) {
      auto& node = nodeIt.second;
      nodeFlow[nodeIndexMap[node.id]] = node.weight;
    }

    normalizeNodeFlow = true;
  }

  if (normalizeNodeFlow) {
    normalize(nodeFlow);
  }

  // Write back flow to network
  double sumNodeFlow = 0.0;
  double sumLinkFlow = 0.0;
  unsigned int linkIndex = 0;
  auto featureLinkIndex = bipartiteLinkStartIndex;

  network.forEachLink([&](unsigned int srcIdx, unsigned int, double, double& flow) {
    if (network.isBipartite() && network.nodeId(srcIdx) >= network.bipartiteStartId()) {
      flow = flowLinks[featureLinkIndex++].flow;
    } else {
      flow = flowLinks[linkIndex++].flow;
    }
    sumLinkFlow += flow;
  });

  double fractionIntraFlow = config.isMultilayerNetwork() && config.regularized ? 1 : 0;

  sumTeleFlow = 0.0;

  for (auto& nodeIt : network.m_nodes) {
    auto& node = nodeIt.second;
    const auto nodeIndex = nodeIndexMap[node.id];
    node.flow = nodeFlow[nodeIndex];
    node.weight = nodeTeleportWeights[nodeIndex];
    node.teleFlow = !nodeTeleportFlow.empty() ? nodeTeleportFlow[nodeIndex] : nodeFlow[nodeIndex] * (nodeOutDegree[nodeIndex] == 0 ? 1 : config.teleportationProbability);
    node.intraLayerTeleFlow = fractionIntraFlow * node.teleFlow;
    node.teleFlow *= 1 - fractionIntraFlow;
    node.intraLayerTeleWeight = nodeTeleportWeights[nodeIndex];
    node.enterFlow = node.flow;
    node.exitFlow = node.flow;

    if (!config.noSelfLinks) {
      // Remove self-teleportation flow
      node.enterFlow -= node.teleFlow * node.weight;
      node.exitFlow -= node.teleFlow * node.weight;

      // Remove self-link flow
      unsigned int norm = config.isUndirectedFlow() ? 2 : 1;
      const auto srcIdx = network.indexOfId(node.id);
      for (unsigned int e = network.m_linkOffsets[srcIdx]; e < network.m_linkOffsets[srcIdx + 1]; ++e) {
        if (network.m_linkTargets[e] == srcIdx) { // self-link: target index == source index
          node.enterFlow -= network.m_linkFlows[e] / norm;
          node.exitFlow -= network.m_linkFlows[e] / norm;
          break;
        }
      }
    }

    sumTeleFlow += node.teleFlow;
    sumNodeFlow += node.flow;
  }

  // Enter/exit flow
  if (enterFlow.empty()) {
    if (!config.isUndirectedClustering() && !config.regularized) {
      enterFlow.assign(numNodes, 0);
      exitFlow.assign(numNodes, 0);
      double alpha = config.teleportationProbability;
      double sumDanglingFlow = 0.0;
      for (unsigned int i = 0; i < numNodes; ++i) {
        if (nodeOutDegree[i] == 0) {
          sumDanglingFlow += nodeFlow[i];
        }
      }
      for (auto& nodeIt : network.m_nodes) {
        auto& node = nodeIt.second;
        const auto sourceIndex = nodeIndexMap[node.id];
        const auto srcIdx = network.indexOfId(node.id);
        double danglingFlow = network.isDangling(srcIdx) ? node.flow : 0.0;
        if (config.recordedTeleportation) {
          // Don't let self-teleportation add to the enter/exit flow (i.e. multiply with (1.0 - node.data.teleportWeight))
          exitFlow[sourceIndex] += alpha * node.flow * (1.0 - node.weight);
          enterFlow[sourceIndex] += (alpha * (1.0 - node.flow) + (1 - alpha) * (sumDanglingFlow - danglingFlow)) * node.weight;
        }
        for (unsigned int e = network.m_linkOffsets[srcIdx]; e < network.m_linkOffsets[srcIdx + 1]; ++e) {
          const auto targetIndex = nodeIndexMap[network.nodeId(network.m_linkTargets[e])];
          exitFlow[sourceIndex] += network.m_linkFlows[e];
          enterFlow[targetIndex] += network.m_linkFlows[e];
        }
      }
    }
  }

  // Save enter/exit flow on nodes
  if (!enterFlow.empty()) {
    for (auto& nodeIt : network.m_nodes) {
      auto& node = nodeIt.second;
      const auto nodeIndex = nodeIndexMap[node.id];
      node.enterFlow = enterFlow[nodeIndex];
      node.exitFlow = exitFlow[nodeIndex];
    }
  }

  Console console;
  console.section("Flow");
  console.metric("Model", io::stringify(config.flowModel));
  console.metric("Method", m_flowMethod.empty() ? "standard" : m_flowMethod);
  if (!m_teleportation.empty())
    console.metric("Teleportation", m_teleportation);
  if (m_havePageRank) {
    console.metric(m_pageRankConverged ? "PageRank" : "PageRank warning",
                   fmt::format(FMT_STRING("{} iterations, error {}"), m_pageRankIterations, io::toPrecision(m_pageRankError)));
  }
  for (const auto& note : m_flowNotes)
    console.status("Note", note);
  console.metric("Node flow sum", io::toPrecision(sumNodeFlow));
  console.metric("Link flow sum", io::toPrecision(sumLinkFlow));
}

} // namespace infomap
