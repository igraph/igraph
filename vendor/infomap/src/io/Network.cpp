/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "Network.h"
#include "InfomapError.h"
#include "NetworkBuilder.h"
#include "../utils/Log.h"
#include "../utils/Console.h"
#include "../utils/convert.h"
#include "../utils/format.h"

#include <cmath>
#include <algorithm>
#include <limits>
#include <set>
#include <stdexcept>

namespace infomap {

using std::make_pair;

namespace {

#if !INFOMAP_FEATURE_REGULARIZED_MULTILAYER
  const char* regularizedMultilayerFeatureError()
  {
    return "Regularized multilayer flow requires building with FEATURES=regularized-multilayer.";
  }
#endif

} // namespace

void Network::init()
{
  updateDerivedConfig();
}

void Network::updateDerivedConfig()
{
  m_multilayerStateIdBitShift = m_config.matchableMultilayerIds == 0
      ? 0
      : static_cast<unsigned int>(std::ceil(std::log2(m_config.matchableMultilayerIds)));
}

void Network::clear()
{
  StateNetwork::clear();
  m_networks.clear();
  m_interLinks.clear();
  m_layerNodeToStateId.clear();
  m_sumIntraOutWeight.clear();
  m_layers.clear();
  m_numInterLayerLinks = 0;
  m_numIntraLayerLinks = 0;

  // Bipartite
  m_bipartiteStartId = 0;

  // Meta data
  m_metaData.clear();
  m_numMetaDataColumns = 0;

  // Embedded initial partition
  m_initialPartitionPaths.clear();
}

void Network::readInputData(std::string filename, bool accumulate)
{
  if (!accumulate) {
    clear();
  }
  if (filename.empty())
    filename = m_config.networkFile;
  if (filename.empty()) {
    throw std::runtime_error("No input file to read network");
  }
  try {
    buildNetworkFromInput(*this, filename);
  } catch (const InfomapError&) {
    throw;
  } catch (const std::exception& e) {
    throw InfomapError(ExitCode::InputError, e.what());
  }
  printSummary();
}

void Network::postProcessInputData()
{
  if (!m_networks.empty()) {
    generateStateNetworkFromMultilayer();
  }

  if (!haveMemoryInput()) {
    // If no memory input, add physical nodes as state nodes to not miss unconnected nodes.
    // Only a weight given for the physical node replaces the state node's, so a
    // weight set through addNode(id, weight) survives (#1183).
    for (const auto& it : m_physNodes) {
      const auto& physNode = it.second;
      if (physNode.hasWeight) {
        addNode(physNode.physId, physNode.weight);
      } else {
        addNode(physNode.physId);
      }
    }
    // The physical nodes now mirror the state nodes 1:1, and nothing reads
    // m_physNodes' contents past this point (only numPhysicalNodes(), which is
    // backed by a counter). Release the map -- a per-source std::map -- to cut
    // the build peak. std::map::clear() deallocates its nodes.
    m_physNodes.clear();
  }
}

void Network::readMetaData(const std::string& filename)
{
  if (!m_metaData.empty()) {
    Console::note(0, "External meta-data '{}' overrides {} meta-data entries already on the network (e.g. embedded JSON 'meta').", filename, m_metaData.size());
  }
  try {
    buildMetaDataFromInput(*this, filename);
  } catch (const InfomapError&) {
    throw;
  } catch (const std::exception& e) {
    throw InfomapError(ExitCode::InputError, e.what());
  }
}

void Network::printSummary()
{
  Console console;
  console.section("Network");
  console.metric("Input", m_config.networkFile);
  console.metric("Direction", m_config.isUndirectedClustering() ? "undirected" : "directed");
  if (haveMemoryInput()) {
    console.metric("Type", isMultilayerNetwork() ? "higher-order multilayer" : "higher-order state");
    console.metric("State nodes", fmt::to_string(numNodes()));
    console.metric("Physical nodes", fmt::to_string(numPhysicalNodes()));
  } else if (m_bipartiteStartId > 0) {
    console.metric("Type", "bipartite first-order");
    console.metric("Bipartite start id", fmt::to_string(m_bipartiteStartId));
    console.metric("Nodes", fmt::to_string(numNodes()));
  } else {
    console.metric("Type", "first-order");
    console.metric("Nodes", fmt::to_string(numNodes()));
  }
  if (isMultilayerNetwork()) {
    console.metric("Layers", fmt::to_string(m_layers.size()));
    console.metric("Layer links", fmt::format(FMT_STRING("{} ({} intra, {} inter)"), m_numIntraLayerLinks + m_numInterLayerLinks, m_numIntraLayerLinks, m_numInterLayerLinks));
  }
  console.metric("Links", fmt::to_string(numLinks()));
  console.metric("Total weight", fmt::format(FMT_STRING("{:g}"), m_totalLinkWeightAdded));
  if (m_numLinksIgnoredByWeightThreshold > 0) {
    console.metric("Ignored by threshold", fmt::format(FMT_STRING("{} links, weight {:g} ({})"), m_numLinksIgnoredByWeightThreshold, m_totalLinkWeightIgnored, Console::percent(m_totalLinkWeightIgnored / (m_totalLinkWeightIgnored + m_totalLinkWeightAdded) * 100)));
  }
}

void Network::addMultilayerLink(unsigned int layer1, unsigned int n1, unsigned int layer2, unsigned int n2, double weight)
{
  m_higherOrderInputMethodCalled = true;
  if (weight < m_config.weightThreshold) {
    ++m_numLinksIgnoredByWeightThreshold;
    m_totalLinkWeightIgnored += weight;
    return;
  }
  unsigned int stateId1 = addMultilayerNode(layer1, n1);
  unsigned int stateId2 = addMultilayerNode(layer2, n2);
  addMultilayerLink(stateId1, layer1, n1, stateId2, layer2, n2, weight);
}

void Network::addMultilayerLink(unsigned int stateId1, unsigned int layer1, unsigned int n1, unsigned int stateId2, unsigned int layer2, unsigned int, double weight)
{
  // if (stateId1 == stateId2) {
  // TODO: Handle self-links?
  // }

  if (layer1 == layer2) {
    ++m_numIntraLayerLinks;
    m_sumIntraOutWeight[layer1][n1] += weight; // TODO: Not used? Add on target also if undirected (not inter/intra format)?
  } else {
    ++m_numInterLayerLinks;
  }

  // The main expanded network builds through the flat buffer (mode A), exactly
  // like ordinary/state networks: this is a write-only sink during expansion
  // (the only main-map read was already dead code) and during explicit
  // *Multilayer parsing, so aggregation can be deferred to finalizeLinks()
  // instead of paying a per-source std::map. Mode is left at its default (A);
  // only the transient per-layer networks stay mode B (see addMultilayerIntraLink).
  addLink(stateId1, stateId2, weight);
}

void Network::addMultilayerLinks(const std::vector<unsigned int>& sourceLayerIds,
                                 const std::vector<unsigned int>& sourceNodeIds,
                                 const std::vector<unsigned int>& targetLayerIds,
                                 const std::vector<unsigned int>& targetNodeIds,
                                 const std::vector<double>& weights)
{
  if (sourceLayerIds.size() != sourceNodeIds.size() || sourceLayerIds.size() != targetLayerIds.size() || sourceLayerIds.size() != targetNodeIds.size() || sourceLayerIds.size() != weights.size()) {
    throw std::invalid_argument("sourceLayerIds, sourceNodeIds, targetLayerIds, targetNodeIds, and weights must have the same length");
  }

  for (std::size_t i = 0; i < sourceLayerIds.size(); ++i) {
    addMultilayerLink(sourceLayerIds[i], sourceNodeIds[i], targetLayerIds[i], targetNodeIds[i], weights[i]);
  }
}

void Network::generateStateNetworkFromMultilayer()
{
#if !INFOMAP_FEATURE_REGULARIZED_MULTILAYER
  if (m_config.regularized) {
    throw std::runtime_error(regularizedMultilayerFeatureError());
  }
#endif

  // As inter-layer links is directed to neighbouring nodes in target layer,
  // the symmetry is broken so we need directed links for inter-layer flow
  if (!m_networks.empty() && m_config.isUndirectedFlow()) {
    m_haveDirectedInput = true;
    // undirdir is replaced by the directed model here and warns in configureNetworkMode (#1233).
    // outdirdir and rawdir are not undirected flow and keep the intra-layer links as given.
    // Expand each undirected intra-layer link to two opposite directed links
    Console::detail(1, "expanding undirected links to directed");
    for (auto& layerIt : m_networks) {
      auto& network = layerIt.second;
      network.undirectedToDirected();
    }
  }

  if (!m_interLinks.empty()) {
    generateStateNetworkFromMultilayerWithInterLinks();
  } else if (m_config.regularized) {
    generateStateNetworkFromMultilayerWithSimulatedInterLinksBasedOnNodeStrengthRegularized();
  } else {
    generateStateNetworkFromMultilayerWithSimulatedInterLinks();
  }

  // Update physId and layerId on link map
  // TODO: Fix this here?
  // for (const auto& node : m_nodeLinkMap) {
  //   for (const auto& link : node.second) {
  //     auto& s2 = m_nodes[link.first.id];
  //     link.first.physicalId = s2.physicalId;
  //     link.first.layerId = s2.layerId;
  //   }
  // }

  m_networks.clear();
  m_interLinks.clear();
}

void Network::generateStateNetworkFromMultilayerWithInterLinks()
{
  Console::detail(1, "generating state network from multilayer networks with inter-layer links");
  Log(1) << std::flush;
  // First add intra-layer links
  for (auto& layerIt : m_networks) {
    unsigned int layer1 = layerIt.first;
    auto& network = layerIt.second;
    for (auto& linkIt : network.nodeLinkMap()) {
      // Per-layer networks key links by physical node id
      auto& source = linkIt.first;
      const auto& subLinks = linkIt.second;
      for (auto& subIt : subLinks) {
        auto& target = subIt.first;
        double linkWeight = subIt.second.weight;
        addMultilayerLink(layer1, source, layer1, target, linkWeight);
      }
    }
  }

  // Spread drops an inter-layer link to a layer where the node has no intra-layer
  // out-links, so the to-self links below are dropped there too: the relax-to-self
  // flow relays the link's flow onward along those out-links, and without any the
  // flow would be lost (#1314).
  auto hasIntraOutLinks = [this](unsigned int layer, unsigned int physId) {
    auto networkIt = m_networks.find(layer);
    if (networkIt == m_networks.end()) {
      return false;
    }
    const auto& linkMap = networkIt->second.nodeLinkMap();
    auto linksIt = linkMap.find(physId);
    return linksIt != linkMap.end() && !linksIt->second.empty();
  };

  Console::detail(1, "connecting layers");
  // Connect layers with inter-layer links spread out in target layer
  for (auto& it : m_interLinks) {
    auto& layerNode = it.first;
    unsigned int layer1 = layerNode.layer;
    unsigned int physId = layerNode.node;
    unsigned int stateId1 = addMultilayerNode(layer1, physId);
    for (auto& it2 : it.second) {
      unsigned int layer2 = it2.first;
      double interWeight = it2.second;

      if (m_config.multilayerRelaxToSelf) {
        if (!hasIntraOutLinks(layer2, physId)) {
          continue;
        }
        unsigned int stateId2 = addMultilayerNode(layer2, physId);
        addLink(stateId1, stateId2, interWeight);
        ++m_numInterLayerLinks;
        continue;
      }

      auto& targetNetwork = m_networks[layer2];

      auto& targetLinks = targetNetwork.nodeLinkMap();
      auto& outlinks = targetLinks[physId];
      if (outlinks.empty()) {
        continue;
      }
      auto& targetOutWeights = targetNetwork.outWeights();
      double sumIntraOutWeightTargetLayer = targetOutWeights[physId];

      for (auto& outLink : outlinks) {
        auto& targetPhysId = outLink.first;
        auto& linkData = outLink.second;
        double intraWeight = linkData.weight;
        unsigned int stateId2i = addMultilayerNode(layer2, targetPhysId);

        double weight = sumIntraOutWeightTargetLayer == 0.0 ? 0.0 : interWeight * intraWeight / sumIntraOutWeightTargetLayer;

        addLink(stateId1, stateId2i, weight);
        ++m_numInterLayerLinks; // TODO: Count all as one?
      }
    }
  }
  if (m_config.isUndirectedFlow()) {
    // For undirected inter-layer links, expand and add in other direction too
    for (auto& it : m_interLinks) {
      auto& layerNode = it.first;
      unsigned int layer2 = layerNode.layer;
      unsigned int physId = layerNode.node;

      if (m_config.multilayerRelaxToSelf) {
        if (!hasIntraOutLinks(layer2, physId)) {
          continue;
        }
        unsigned int stateId2 = addMultilayerNode(layer2, physId);
        for (auto& it2 : it.second) {
          unsigned int layer1 = it2.first;
          double interWeight = it2.second;
          unsigned int stateId1 = addMultilayerNode(layer1, physId);
          addLink(stateId1, stateId2, interWeight);
          ++m_numInterLayerLinks;
        }
        continue;
      }

      auto& targetNetwork = m_networks[layer2];
      auto& targetLinks = targetNetwork.nodeLinkMap();
      auto& outlinks = targetLinks[physId];
      if (outlinks.empty()) {
        continue;
      }
      auto& targetOutWeights = targetNetwork.outWeights();
      double sumIntraOutWeightTargetLayer = targetOutWeights[physId];
      for (auto& it2 : it.second) {
        unsigned int layer1 = it2.first;
        double interWeight = it2.second;
        unsigned int stateId1 = addMultilayerNode(layer1, physId);

        for (auto& outLink : outlinks) {
          auto& targetPhysId = outLink.first;
          auto& linkData = outLink.second;
          double intraWeight = linkData.weight;
          unsigned int stateId2i = addMultilayerNode(layer2, targetPhysId);

          double weight = sumIntraOutWeightTargetLayer == 0.0 ? 0.0 : interWeight * intraWeight / sumIntraOutWeightTargetLayer;

          addLink(stateId1, stateId2i, weight);
          ++m_numInterLayerLinks; // TODO: Count all as one?
        }
      }
    }
  }
}

void Network::generateStateNetworkFromMultilayerWithSimulatedInterLinks()
{
  if (m_config.multilayerRelaxByJensenShannonDivergence) {
    return generateStateNetworkFromMultilayerWithSimulatedInterLinksBasedOnLayerSimilarity();
  }
  return generateStateNetworkFromMultilayerWithSimulatedInterLinksBasedOnNodeStrength();
}

void Network::generateStateNetworkFromMultilayerWithSimulatedInterLinksBasedOnLayerSimilarity()
{
  Console::detail(1, "generating state network from multilayer networks with simulated inter-layer links based on layer similarity");
  Log(1) << std::flush;
  double relaxRate = m_config.multilayerRelaxRate;

  // A negative limit means no limit, also across gaps in the layer ids
  const unsigned int noRelaxLimit = std::numeric_limits<unsigned int>::max();
  auto toRelaxLimit = [](int limit) { return limit < 0 ? std::numeric_limits<unsigned int>::max() : static_cast<unsigned int>(limit); };
  unsigned int relaxLimitSymmetric = toRelaxLimit(m_config.multilayerRelaxLimit);
  unsigned int relaxLimitDown = std::min(relaxLimitSymmetric, toRelaxLimit(m_config.multilayerRelaxLimitDown));
  unsigned int relaxLimitUp = std::min(relaxLimitSymmetric, toRelaxLimit(m_config.multilayerRelaxLimitUp));
  auto haveUpOrDownLimit = m_config.multilayerRelaxLimitDown >= 0 || m_config.multilayerRelaxLimitUp >= 0;

  Console::detail(1, "{} networks", m_networks.size());
  Console::detail(1, "relax rate: {:g}", relaxRate);
  if (haveUpOrDownLimit) {
    Console::detail(1, "relax limit up: {}", relaxLimitUp == noRelaxLimit ? "no limit" : std::to_string(relaxLimitUp));
    Console::detail(1, "relax limit down: {}", relaxLimitDown == noRelaxLimit ? "no limit" : std::to_string(relaxLimitDown));
  } else if (m_config.multilayerRelaxLimit >= 0) {
    Console::detail(1, "relax limit: {}", m_config.multilayerRelaxLimit);
  }

  // Whether a state in layer1 may relax to layer2. A layer is always within reach
  // of itself, so the intra-layer links are kept with limit 0.
  auto withinRelaxLimit = [relaxLimitDown, relaxLimitUp](unsigned int layer1, unsigned int layer2) {
    return layer1 >= layer2 ? layer1 - layer2 <= relaxLimitDown : layer2 - layer1 <= relaxLimitUp;
  };
  // Layers further below than this are out of reach in both directions
  unsigned int relaxReach = std::max(relaxLimitDown, relaxLimitUp);

  if (m_config.multilayerRelaxByJensenShannonDivergence) {
    Console::detail(1, "using Jensen-Shannon divergence");

    for (unsigned int nodeId = 0; nodeId <= m_maxNodeIdInIntraLayerNetworks; ++nodeId) {
      // Calculate Jensen-Shannon similarity between all layers such that layer1 >= layer2,
      // and then use its symmetry for layer2 > layer1. Keep only the directions within
      // the relax limits.
      std::map<unsigned int, std::map<unsigned int, double>> jsRelaxWeights;
      std::map<unsigned int, double> jsTotWeight;

      for (auto& layer1It : m_networks) {
        unsigned int layer1 = layer1It.first;
        auto& layer1OutLinks = layer1It.second.nodeLinkMap()[nodeId];
        // Skip dangling nodes, because they have no information to calculate similarity
        if (layer1OutLinks.empty())
          continue;

        double sumOutLinkWeightLayer1 = layer1It.second.outWeights()[nodeId];

        auto layer2It = m_networks.lower_bound(layer1 > relaxReach ? layer1 - relaxReach : 0);
        for (; layer2It != m_networks.end() && layer2It->first <= layer1; ++layer2It) {
          unsigned int layer2 = layer2It->first;
          bool down = withinRelaxLimit(layer1, layer2);
          bool up = withinRelaxLimit(layer2, layer1);
          if (!down && !up)
            continue;

          auto& layer2OutLinks = layer2It->second.nodeLinkMap()[nodeId];
          if (layer2OutLinks.empty())
            continue;

          double sumOutLinkWeightLayer2 = layer2It->second.outWeights()[nodeId];

          bool intersect;
          double div = calculateJensenShannonDivergence(intersect, layer1OutLinks, sumOutLinkWeightLayer1, layer2OutLinks, sumOutLinkWeightLayer2);
          double jsWeight = 1.0 - div;
          if (intersect && (jsWeight >= m_config.multilayerJSRelaxLimit)) {
            if (down) {
              jsTotWeight[layer1] += jsWeight;
              jsRelaxWeights[layer1][layer2] = jsWeight;
            }
            if (up && layer1 != layer2) {
              jsTotWeight[layer2] += jsWeight;
              jsRelaxWeights[layer2][layer1] = jsWeight;
            }
          }
        }
      }

      // Second loop over the pairs of layers within the relax limits
      for (auto& jsRelaxWeightsLayer1It : jsRelaxWeights) {
        unsigned int layer1 = jsRelaxWeightsLayer1It.first;
        double sumOutLinkWeightLayer1 = m_networks[layer1].outWeights()[nodeId];
        double jsTotWeightLayer1 = jsTotWeight[layer1];

        // Create inter-links to the intra-connected nodes in other layers
        for (auto& jsRelaxWeightsIt : jsRelaxWeightsLayer1It.second) {
          unsigned int layer2 = jsRelaxWeightsIt.first;
          bool isIntra = layer2 == layer1;

          // Create inter-links to the outgoing nodes in the target layer
          double linkWeightNormalizationFactor;
          if (isIntra) {
            linkWeightNormalizationFactor = 1;
          } else {
            linkWeightNormalizationFactor = jsRelaxWeightsIt.second * relaxRate / (1.0 - relaxRate) * sumOutLinkWeightLayer1 / jsTotWeightLayer1;
          }

          auto& targetLinks = m_networks[layer2].nodeLinkMap();
          auto& targetOutlinks = targetLinks[nodeId];
          if (targetOutlinks.empty()) {
            continue;
          }
          for (auto& outLink : targetOutlinks) {
            auto& n2 = outLink.first;
            auto& linkData = outLink.second;
            double intraWeight = linkData.weight;
            // Add intra link weight as teleport weight to source node
            // TODO: Why, and only first intra link that sets the teleport weight?
            unsigned int stateId1 = addMultilayerNode(layer1, nodeId, intraWeight);
            unsigned int stateId2i = addMultilayerNode(layer2, n2, 0.0);

            double weight = intraWeight == 0.0 ? 0.0 : linkWeightNormalizationFactor * intraWeight;
            addLink(stateId1, stateId2i, weight);
            ++m_numInterLayerLinks;
          }
        }
      }
    }
  }
}

void Network::generateStateNetworkFromMultilayerWithSimulatedInterLinksBasedOnNodeStrength()
{
  Console::detail(1, "generating state network from multilayer networks with simulated inter-layer links");
  Log(1) << std::flush;
  double relaxRate = m_config.multilayerRelaxRate;

  // A negative limit means no limit. Layer ids can have gaps, so no limit is
  // not the number of layers but the largest difference (#1322).
  const unsigned int maxRelaxLimit = std::numeric_limits<unsigned int>::max();
  auto toRelaxLimit = [](int limit) { return limit < 0 ? std::numeric_limits<unsigned int>::max() : static_cast<unsigned int>(limit); };
  unsigned int relaxLimitSymmetric = toRelaxLimit(m_config.multilayerRelaxLimit);
  unsigned int relaxLimitDown = std::min(relaxLimitSymmetric, toRelaxLimit(m_config.multilayerRelaxLimitDown));
  unsigned int relaxLimitUp = std::min(relaxLimitSymmetric, toRelaxLimit(m_config.multilayerRelaxLimitUp));
  auto haveUpOrDownLimit = m_config.multilayerRelaxLimitDown >= 0 || m_config.multilayerRelaxLimitUp >= 0;

  Console::detail(1, "{} networks", m_networks.size());
  Console::detail(1, "relax rate: {:g}", relaxRate);
  if (haveUpOrDownLimit) {
    Console::detail(1, "relax limit up: {}", relaxLimitUp == maxRelaxLimit ? "no limit" : std::to_string(relaxLimitUp));
    Console::detail(1, "relax limit down: {}", relaxLimitDown == maxRelaxLimit ? "no limit" : std::to_string(relaxLimitDown));
  } else if (m_config.multilayerRelaxLimit >= 0) {
    Console::detail(1, "relax limit: {}", m_config.multilayerRelaxLimit);
  }

  auto withinRelaxLimit = [relaxLimitDown, relaxLimitUp](unsigned int layer1, unsigned int layer2) {
    return layer1 >= layer2 ? layer1 - layer2 <= relaxLimitDown : layer2 - layer1 <= relaxLimitUp;
  };

  for (auto& it1 : m_networks) {
    auto layer1 = it1.first;
    auto& network1 = it1.second;

    for (auto& n1It : network1.nodes()) {
      auto& n1 = n1It.first;
      unsigned int stateId1 = addMultilayerNode(layer1, n1);

      double sumOutLinkWeightLayer1 = network1.outWeights()[n1];
      double sumOutWeightAllLayers = 0.0;

      for (auto& it2 : m_networks) {
        auto layer2 = it2.first;
        if (!withinRelaxLimit(layer1, layer2)) {
          continue;
        }
        auto& network2 = it2.second;
        sumOutWeightAllLayers += network2.outWeights()[n1];
      }

      if (sumOutWeightAllLayers <= 0) {
        continue;
      }

      for (auto& it2 : m_networks) {
        auto layer2 = it2.first;
        if (!withinRelaxLimit(layer1, layer2)) {
          continue;
        }
        auto& network2 = it2.second;
        bool isIntra = layer2 == layer1;

        if (m_config.multilayerRelaxToSelf && !isIntra) {
          // Couple only to the same physical node in the target layer, carrying
          // the same per-layer relaxation mass the spread model distributes over
          // the node's neighbours there (r * s_j(n) / S_n). The home layer is
          // left to the spread model below.
          double targetOutWeight = network2.outWeights()[n1];
          if (targetOutWeight <= 0) {
            continue; // n1 absent/dangling in the target layer
          }
          double weight = relaxRate * targetOutWeight / sumOutWeightAllLayers;
          if (weight < 1e-16) {
            continue;
          }
          unsigned int stateId2 = addMultilayerNode(layer2, n1);
          addLink(stateId1, stateId2, weight);
          ++m_numInterLayerLinks;
          continue;
        }

        auto& targetLinks = network2.nodeLinkMap();
        auto& targetOutlinks = targetLinks[n1];
        if (targetOutlinks.empty()) {
          continue;
        }
        // After the empty check: a node without intra-layer out-links has sumOutLinkWeightLayer1 == 0.
        double linkWeightNormalizationFactor = relaxRate / sumOutWeightAllLayers;
        if (isIntra) {
          linkWeightNormalizationFactor += (1.0 - relaxRate) / sumOutLinkWeightLayer1;
        }
        for (auto& outLink : targetOutlinks) {
          auto& n2 = outLink.first;
          auto& linkData = outLink.second;
          double intraWeight = linkData.weight;
          unsigned int stateId2i = addMultilayerNode(layer2, n2);

          double weight = linkWeightNormalizationFactor * intraWeight;
          if (weight < 1e-16) {
            continue;
          }
          addLink(stateId1, stateId2i, weight);
          ++m_numInterLayerLinks; // TODO: Count all as one?
        }
      }
    }
  }
}

void Network::generateStateNetworkFromMultilayerWithSimulatedInterLinksBasedOnNodeStrengthRegularized()
{
  Console::detail(1, "generating regularized state network from multilayer networks with simulated inter-layer links");
  Log(1) << std::flush;
  // double relaxRate = m_config.multilayerRelaxRate;

  unsigned int L = m_networks.size();
  // unsigned int N = numPhysicalNodes();
  // A negative limit means no limit, also across gaps in the layer ids (#1322).
  const unsigned int maxRelaxLimit = std::numeric_limits<unsigned int>::max();
  auto toRelaxLimit = [](int limit) { return limit < 0 ? std::numeric_limits<unsigned int>::max() : static_cast<unsigned int>(limit); };
  unsigned int relaxLimitSymmetric = toRelaxLimit(m_config.multilayerRelaxLimit);
  unsigned int relaxLimitDown = std::min(relaxLimitSymmetric, toRelaxLimit(m_config.multilayerRelaxLimitDown));
  unsigned int relaxLimitUp = std::min(relaxLimitSymmetric, toRelaxLimit(m_config.multilayerRelaxLimitUp));
  auto haveUpOrDownLimit = m_config.multilayerRelaxLimitDown >= 0 || m_config.multilayerRelaxLimitUp >= 0;

  // double interLinkStrength = std::log(L);
  // double interLinkWeight = std::log(L) / L;
  // The inter-layer prior spreads ln(L_i) over the L_i layers that hold physical node i,
  // so its per-link weight is only a constant while every node is in every layer.
  // double intraLinkStrength = std::log(N) / L;
  if (haveUpOrDownLimit) {
    // TODO: create a map with aggregated inter-flow for each layer with limits.
  }

  Console::detail(1, "{} networks and {} physical nodes", m_networks.size(), m_physNodes.size());
  // Log(1) << "-> Relax rate: " << relaxRate << "\n";
  if (haveUpOrDownLimit) {
    Console::detail(1, "relax limit up: {}", relaxLimitUp == maxRelaxLimit ? "no limit" : std::to_string(relaxLimitUp));
    Console::detail(1, "relax limit down: {}", relaxLimitDown == maxRelaxLimit ? "no limit" : std::to_string(relaxLimitDown));
    throw std::runtime_error("Relax limits not implemented for regularized flow");
  } else if (m_config.multilayerRelaxLimit >= 0) {
    Console::detail(1, "relax limit: {}", m_config.multilayerRelaxLimit);
    throw std::runtime_error("Relax limits not implemented for regularized flow");
  }

  auto withinRelaxLimit = [relaxLimitDown, relaxLimitUp](unsigned int layer1, unsigned int layer2) {
    return layer1 >= layer2 ? layer1 - layer2 <= relaxLimitDown : layer2 - layer1 <= relaxLimitUp;
  };

  // Which physical nodes each layer holds, and which layers hold each physical node.
  // With --multilayer-skip-absent-nodes a node belongs to a layer only if that layer's
  // links mention it, so the prior spans what the data shows; by default a node is
  // assumed to exist in every layer, unobserved where its links are missing. Both maps
  // stay proportional to the input: the rosters together hold one entry per (node,
  // layer) pair that the prior covers.
  std::map<unsigned int, std::vector<unsigned int>> nodesInLayer;
  std::map<unsigned int, std::vector<unsigned int>> layersOfNode;

  for (auto& it : m_networks) {
    auto layer = it.first;
    auto& roster = nodesInLayer[layer];

    if (m_config.multilayerSkipAbsentNodes) {
      std::set<unsigned int> present;
      for (const auto& outLinks : it.second.nodeLinkMap()) {
        present.insert(outLinks.first);
        for (const auto& outLink : outLinks.second) {
          present.insert(outLink.first);
        }
      }
      roster.assign(present.begin(), present.end());
    } else {
      roster.reserve(m_physNodes.size());
      for (const auto& physNodeIt : m_physNodes) {
        roster.push_back(physNodeIt.first);
      }
    }

    for (auto n : roster) {
      layersOfNode[n].push_back(layer);
    }
  }

  for (auto& it1 : m_networks) {
    auto layer1 = it1.first;
    auto& network1 = it1.second;

    // Loop over the physical nodes this layer holds, even if they lack intra links, to
    // add inter links
    for (auto n1 : nodesInLayer[layer1]) {
      unsigned int stateId1 = addMultilayerNode(layer1, n1);

      // double sumOutLinkWeightLayer1 = network1.outWeights()[n1];
      // Multiply observed link weight with this factor to account for the physical part of inter-link flow
      // double gamma = 1 + interLinkStrength / (sumOutLinkWeightLayer1 + intraLinkStrength);
      // Log(1) << "(" << layer1 << "," << n1 << "): gamma = 1 + " << interLinkStrength << " / (" << sumOutLinkWeightLayer1 << " + " << intraLinkStrength << ") = " << gamma << "\n";
      // Need original weight to calculate correct inter-layer relax rate.

      // Add intra links
      auto targetOutLinksIt = network1.nodeLinkMap().find(n1);
      if (targetOutLinksIt != network1.nodeLinkMap().end()) {
        for (auto& outLink : targetOutLinksIt->second) {
          auto& n2 = outLink.first;
          auto& linkData = outLink.second;
          // double weight = linkData.weight * gamma;
          double weight = linkData.weight;
          unsigned int stateId2i = addMultilayerNode(layer1, n2);

          if (weight < 1e-16) {
            continue;
          }
          addLink(stateId1, stateId2i, weight);
          ++m_numIntraLayerLinks;
        }
      }

      // Add inter links. The prior couples node n1 across the L_i layers that hold it
      // with total strength ln(L_i), which is zero for a node confined to one layer:
      // nothing to spread, and the per-link weight below would divide by zero.
      const auto& layers = layersOfNode[n1];
      if (layers.size() < 2) {
        continue;
      }
      const auto numInterTargets = m_config.noSelfLinks ? layers.size() - 1 : layers.size();
      const double interLinkWeight = m_config.regularizationStrength * m_config.interRegularizationStrength * std::log(static_cast<double>(layers.size())) / numInterTargets;

      for (auto layer2 : layers) {
        if (!withinRelaxLimit(layer1, layer2) || (m_config.noSelfLinks && layer1 == layer2)) {
          continue;
        }

        unsigned int stateId2 = addMultilayerNode(layer2, n1);
        addLink(stateId1, stateId2, interLinkWeight);
        ++m_numInterLayerLinks;
      }
    }
  }
}

double Network::calculateJensenShannonDivergence(bool& intersect, const OutLinkMap& layer1OutLinks, double sumOutLinkWeightLayer1, const OutLinkMap& layer2OutLinks, double sumOutLinkWeightLayer2)
{
  intersect = false;
  double h1 = 0.0; // The entropy rate of the node in the first layer
  double h2 = 0.0; // The entropy rate of the node in the second layer
  double h12 = 0.0; // The entropy rate of the lumped node
  // The out-link weights of the nodes
  double ow1 = sumOutLinkWeightLayer1;
  double ow2 = sumOutLinkWeightLayer2;
  // Normalized weights over node in layer 1 and 2
  double pi1 = ow1 / (ow1 + ow2);
  double pi2 = ow2 / (ow1 + ow2);

  auto layer1OutLinkIt = layer1OutLinks.begin();
  auto layer2OutLinkIt = layer2OutLinks.begin();
  auto layer1OutLinkItEnd = layer1OutLinks.end();
  auto layer2OutLinkItEnd = layer2OutLinks.end();
  while (layer1OutLinkIt != layer1OutLinkItEnd && layer2OutLinkIt != layer2OutLinkItEnd) {
    if (layer1OutLinkIt->first < layer2OutLinkIt->first) {
      // If the first state node has a link that the second has not
      double p1 = layer1OutLinkIt->second.weight / ow1;
      h1 -= p1 * log2(p1);
      double p12 = pi1 * layer1OutLinkIt->second.weight / ow1;
      h12 -= p12 * log2(p12);
      layer1OutLinkIt++;
    } else if (layer2OutLinkIt->first < layer1OutLinkIt->first) {
      // If the second state node has a link that the second has not
      double p2 = layer2OutLinkIt->second.weight / ow2;
      h2 -= p2 * log2(p2);
      double p12 = pi2 * layer2OutLinkIt->second.weight / ow2;
      h12 -= p12 * log2(p12);
      layer2OutLinkIt++;
    } else { // If both state nodes have the link
      intersect = true;
      double p1 = layer1OutLinkIt->second.weight / ow1;
      h1 -= p1 * log2(p1);
      double p2 = layer2OutLinkIt->second.weight / ow2;
      h2 -= p2 * log2(p2);
      double p12 = pi1 * layer1OutLinkIt->second.weight / ow1 + pi2 * layer2OutLinkIt->second.weight / ow2;
      h12 -= p12 * log2(p12);
      layer1OutLinkIt++;
      layer2OutLinkIt++;
    }
  }

  while (layer1OutLinkIt != layer1OutLinkItEnd) {
    // If the first state node has a link that the second has not
    double p1 = layer1OutLinkIt->second.weight / ow1;
    h1 -= p1 * log2(p1);
    double p12 = pi1 * layer1OutLinkIt->second.weight / ow1;
    h12 -= p12 * log2(p12);
    layer1OutLinkIt++;
  }

  while (layer2OutLinkIt != layer2OutLinkItEnd) {
    // If the second state node has a link that the second has not
    double p2 = layer2OutLinkIt->second.weight / ow2;
    h2 -= p2 * log2(p2);
    double p12 = pi2 * layer2OutLinkIt->second.weight / ow2;
    h12 -= p12 * log2(p12);
    layer2OutLinkIt++;
  }

  double div = (pi1 + pi2) * h12 - pi1 * h1 - pi2 * h2;

  // Fix precision problems
  if (div < 0.0)
    div = 0.0;
  else if (div > 1.0)
    div = 1.0;

  return div;
}

void Network::simulateInterLayerLinks()
{
}

void Network::addMultilayerIntraLink(unsigned int layer, unsigned int n1, unsigned int n2, double weight)
{
  m_higherOrderInputMethodCalled = true;
  // The main network gets its expansion links later and builds via the flat
  // buffer (mode A) like ordinary/state networks, so it is left at its default
  // mode. Only the transient per-layer network stays mode B: its nested map +
  // outWeights are read directly during expansion, so set before addLink so it
  // does not enter the buffer path.
  auto& layerNetwork = m_networks[layer];
  layerNetwork.m_useMapBuild = true;
  bool added = layerNetwork.addLink(n1, n2, weight);
  if (added) {
    ++m_numIntraLayerLinks;
    m_maxNodeIdInIntraLayerNetworks = std::max(m_maxNodeIdInIntraLayerNetworks, std::max(n1, n2));
  }
  addPhysicalNode(n1);
  addPhysicalNode(n2);
}

void Network::addMultilayerIntraLinks(const std::vector<unsigned int>& layerIds,
                                      const std::vector<unsigned int>& sourceNodeIds,
                                      const std::vector<unsigned int>& targetNodeIds,
                                      const std::vector<double>& weights)
{
  if (layerIds.size() != sourceNodeIds.size() || layerIds.size() != targetNodeIds.size() || layerIds.size() != weights.size()) {
    throw std::invalid_argument("layerIds, sourceNodeIds, targetNodeIds, and weights must have the same length");
  }

  for (std::size_t i = 0; i < layerIds.size(); ++i) {
    addMultilayerIntraLink(layerIds[i], sourceNodeIds[i], targetNodeIds[i], weights[i]);
  }
}

void Network::addMultilayerInterLink(unsigned int layer1, unsigned int n, unsigned int layer2, double interWeight)
{
  if (layer1 == layer2) {
    throw std::runtime_error(fmt::format(FMT_STRING("Inter-layer link (layer1, node, layer2): {}, {}, {} must have layer1 != layer2"), layer1, n, layer2));
  }
  m_higherOrderInputMethodCalled = true;
  // Inter-layer links are buffered in m_interLinks and only realized on the main
  // network during expansion, which builds via the flat buffer (mode A). The main
  // network is left at its default mode here. (File input expands before any
  // numLinks() read; an in-memory numLinks() read before expansion just finalizes
  // an empty buffer, which the first expansion addLink re-opens via definalize().)
  auto& interLinks = m_interLinks[LayerNode(layer1, n)];
  auto it = interLinks.find(layer2);

  if (it == interLinks.end()) {
    ++m_numInterLayerLinks;
  }
  interLinks[layer2] += interWeight;
}

void Network::addMultilayerInterLinks(const std::vector<unsigned int>& sourceLayerIds,
                                      const std::vector<unsigned int>& nodeIds,
                                      const std::vector<unsigned int>& targetLayerIds,
                                      const std::vector<double>& weights)
{
  if (sourceLayerIds.size() != nodeIds.size() || sourceLayerIds.size() != targetLayerIds.size() || sourceLayerIds.size() != weights.size()) {
    throw std::invalid_argument("sourceLayerIds, nodeIds, targetLayerIds, and weights must have the same length");
  }

  for (std::size_t i = 0; i < sourceLayerIds.size(); ++i) {
    addMultilayerInterLink(sourceLayerIds[i], nodeIds[i], targetLayerIds[i], weights[i]);
  }
}

unsigned int Network::addMultilayerNode(unsigned int layerId, unsigned int physicalId, double weight)
{
  m_higherOrderInputMethodCalled = true;

  // Create state node if not already exist, return state node id
  auto& layerIt = m_layerNodeToStateId[layerId];
  auto it = layerIt.find(physicalId);

  if (it != layerIt.end()) {
    return it->second;
  }

  bool matchableMultilayerIds = m_config.matchableMultilayerIds != 0;

  if (matchableMultilayerIds && layerId > m_config.matchableMultilayerIds) {
    throw std::runtime_error(fmt::format(FMT_STRING("Cannot add node with layer {} to network with matchable multilayer ids using largest layer id {}"), layerId, m_config.matchableMultilayerIds));
  }

  auto ret = matchableMultilayerIds
      ? addStateNodeWithDeterministicId(physicalId, layerId, m_multilayerStateIdBitShift)
      : addStateNodeWithAutogeneratedId(physicalId);
  auto& stateNode = ret.first->second;
  stateNode.layerId = layerId;
  stateNode.weight = weight;
  m_layerNodeToStateId[layerId][physicalId] = stateNode.id;
  m_layers.insert(layerId);
  return stateNode.id;
}

unsigned int Network::addMultilayerNode(unsigned int stateId, unsigned int layerId, unsigned int physicalId, double weight)
{
  return addMultilayerNode(stateId, layerId, physicalId, weight, std::string());
}

unsigned int Network::addMultilayerNode(unsigned int stateId, unsigned int layerId, unsigned int physicalId, double weight, std::string name)
{
  m_higherOrderInputMethodCalled = true;

  // Create state node if not already exist, return state node id
  auto& layerIt = m_layerNodeToStateId[layerId];
  auto it = layerIt.find(physicalId);

  if (it != layerIt.end()) {
    return it->second;
  }

  auto ret = addStateNode(stateId, physicalId, std::move(name));
  auto& stateNode = ret.first->second;
  stateNode.layerId = layerId;
  stateNode.weight = weight;
  m_layerNodeToStateId[layerId][physicalId] = stateNode.id;
  m_layers.insert(layerId);
  return stateNode.id;
}

void Network::addInitialPartitionPath(unsigned int nodeId, const Path& path, TreeLeafIdType idType)
{
  m_initialPartitionPaths.push_back(TreePath { nodeId, path, idType });
}

void Network::addMetaData(unsigned int nodeId, int meta)
{
  std::vector<int> metaData(1, meta);
  addMetaData(nodeId, metaData);
}

void Network::addMetaData(unsigned int nodeId, const std::vector<int>& metaData)
{
  m_metaData[nodeId] = metaData;
  if (m_numMetaDataColumns == 0) {
    m_numMetaDataColumns = metaData.size();
  } else if (metaData.size() != m_numMetaDataColumns) {
    throw std::runtime_error(fmt::format(FMT_STRING("Must have same number of dimensions in meta data, error trying to add meta data '{}' on node {}."), io::stringify(metaData, ","), nodeId));
  }
}

} // namespace infomap
