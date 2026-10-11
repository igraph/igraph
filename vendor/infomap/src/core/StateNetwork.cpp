/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "StateNetwork.h"
#include "../utils/FlowCalculator.h"
#include "../utils/Log.h"
#include "../utils/format.h"
#include "../io/SafeFile.h"
#include <algorithm>
#include <cmath>
#include <deque>
#include <stdexcept>
#include <utility>

namespace infomap {

std::pair<StateNetwork::NodeMap::iterator, bool> StateNetwork::addStateNode(const StateNode& node)
{
  auto ret = m_nodes.insert(StateNetwork::NodeMap::value_type(node.id, node));
  if (ret.second) {
    // If state node didn't exist, also create the associated physical node
    addPhysicalNode(node.physicalId);
    if (node.id != node.physicalId) {
      m_haveMemoryInput = true;
    }
  }
  return ret;
}

std::pair<StateNetwork::NodeMap::iterator, bool> StateNetwork::addStateNode(unsigned int id, unsigned int physId)
{
  m_higherOrderInputMethodCalled = true;
  return addStateNode(StateNode(id, physId));
}

std::pair<StateNetwork::NodeMap::iterator, bool> StateNetwork::addStateNode(unsigned int id, unsigned int physId, std::string name)
{
  m_higherOrderInputMethodCalled = true;
  return addStateNode(StateNode(id, physId, std::move(name)));
}

std::pair<StateNetwork::NodeMap::iterator, bool> StateNetwork::addNode(unsigned int id)
{
  return addStateNode(StateNode(id));
}

std::pair<StateNetwork::NodeMap::iterator, bool> StateNetwork::addNode(unsigned int id, double weight)
{
  auto ret = addNode(id);
  auto& node = ret.first->second;
  node.weight = weight;
  return ret;
}

std::pair<StateNetwork::NodeMap::iterator, bool> StateNetwork::addNode(unsigned int id, std::string name)
{
  m_names[id] = std::move(name);
  return addNode(id);
}

std::pair<StateNetwork::NodeMap::iterator, bool> StateNetwork::addNode(unsigned int id, std::string name, double weight)
{
  m_names[id] = std::move(name);
  return addNode(id, weight);
}

StateNetwork::PhysNode& StateNetwork::addPhysicalNode(unsigned int physId)
{
  auto result = m_physNodes.emplace(physId, PhysNode(physId));
  if (result.second) {
    ++m_numPhysicalNodesFound; // count uniques so numPhysicalNodes() survives a freed map
  }
  auto& physNode = result.first->second;
  physNode.physId = physId;
  m_sumNodeWeight += 1.0;
  return physNode;
}

StateNetwork::PhysNode& StateNetwork::addPhysicalNode(unsigned int physId, double weight)
{
  auto& physNode = addPhysicalNode(physId);
  physNode.weight = weight;
  m_sumNodeWeight += weight;
  return physNode;
}

StateNetwork::PhysNode& StateNetwork::addPhysicalNode(unsigned int physId, const std::string& name)
{
  auto& physNode = addPhysicalNode(physId);
  m_names[physId] = name;
  m_sumNodeWeight += 1.0;
  return physNode;
}

StateNetwork::PhysNode& StateNetwork::addPhysicalNode(unsigned int physId, double weight, const std::string& name)
{
  auto& physNode = addPhysicalNode(physId);
  physNode.weight = weight;
  m_sumNodeWeight += weight;
  m_names[physId] = name;
  return physNode;
}

std::pair<std::map<unsigned int, std::string>::iterator, bool> StateNetwork::addName(unsigned int id, const std::string& name)
{
  return m_names.insert(std::make_pair(id, name));
}

bool StateNetwork::addLink(unsigned int sourceId, unsigned int targetId, double weight)
{
  // Reject ill-defined weights at ingestion. Link weights feed a flow
  // distribution, so negative, NaN and infinite values are never valid for the
  // map equation -- without this guard the engine silently computes a
  // meaningless result. Zero is well-defined (no flow) and stays allowed; it is
  // filtered by the weight threshold below. This is the single funnel for the
  // in-memory API, file parsing and multilayer expansion, so the language
  // bindings (R, JS, native CLI) no longer need to re-implement the check.
  if (!std::isfinite(weight) || weight < 0) {
    throw std::invalid_argument(
        fmt::format(FMT_STRING("Link weight must be finite and non-negative, got {} for link ({}, {})"), weight, sourceId, targetId));
  }

  // weight is finite and >= 0 past the guard above, so this drops zero-weight
  // links (no flow) and anything below the configured threshold.
  if (weight < m_config.weightThreshold || weight == 0) {
    ++m_numLinksIgnoredByWeightThreshold;
    m_totalLinkWeightIgnored += weight;
    return false;
  }
  m_totalLinkWeightAdded += weight;

  if (sourceId == targetId) {
    ++m_numSelfLinksFound;
    if (m_config.noSelfLinks) {
      return false;
    }
    ++m_numSelfLinks;
    m_sumSelfLinkWeight += weight;
  }

  addNode(sourceId);
  addNode(targetId);

  m_sumLinkWeight += weight;

  if (m_useMapBuild) {
    return addLinkToMap(sourceId, targetId, weight); // mode B (multilayer)
  }

  // mode A (first-order): deferred storage -- append now, aggregate in
  // finalizeLinks(). New-vs-aggregated isn't known until then, so the return
  // value becomes "accepted" (passed the filters above) rather than "new unique".
  if (m_linksFinalized) {
    definalize(); // re-entered build (e.g. accumulate=true read, addLink after run)
  }
  ++m_rawLinkCount;
  m_linkBuffer.push_back({ sourceId, targetId, weight });
  return true;
}

bool StateNetwork::addLinkToMap(unsigned int sourceId, unsigned int targetId, double weight)
{
  // Mode B (multilayer): eager aggregation into the nested map, exactly as the
  // pre-CSR build did. Multilayer expansion reads this map and outWeights()
  // mid-build, so both must stay populated with the original arrival-order sums.
  if (m_linksFinalized) {
    // A link added after a CSR build: the map is still the source of truth in
    // mode B, so just invalidate the cached CSR for the next finalize.
    m_linksFinalized = false;
  }

  ++m_numLinks;

  bool addedNewLink = true;

  // Aggregate link weights if they are defined more than once
  // TODO: Sort source/target order if undirected?
  auto& outLinks = m_nodeLinkMap[sourceId];
  if (outLinks.empty()) {
    outLinks[targetId] = weight;
  } else {
    auto ret = outLinks.insert(std::make_pair(targetId, LinkData(weight)));
    auto& linkData = ret.first->second;
    if (!ret.second) {
      linkData.weight += weight;
      ++m_numAggregatedLinks;
      --m_numLinks;
      addedNewLink = false;
    } else {
      linkData.weight = weight;
    }
  }
  m_outWeights[sourceId] += weight;
  return addedNewLink;
}

bool StateNetwork::addLink(unsigned int sourceId, unsigned int targetId, unsigned long weight)
{
  return addLink(sourceId, targetId, static_cast<double>(weight));
}

void StateNetwork::addLinks(const std::vector<unsigned int>& sourceIds, const std::vector<unsigned int>& targetIds, const std::vector<double>& weights)
{
  if (sourceIds.size() != targetIds.size() || sourceIds.size() != weights.size()) {
    throw std::invalid_argument("sourceIds, targetIds, and weights must have the same length");
  }

  for (std::size_t i = 0; i < sourceIds.size(); ++i) {
    addLink(sourceIds[i], targetIds[i], weights[i]);
  }
}

bool StateNetwork::removeLink(unsigned int sourceId, unsigned int targetId)
{
  // removeLink operates on the nested map. For a mode-A (flat-buffer) network the
  // map is empty, so materialize it from CSR and switch to map build first;
  // otherwise this would silently no-op for ordinary/state networks.
  ensureMapBuild();
  auto itSource = m_nodeLinkMap.find(sourceId);
  if (itSource == m_nodeLinkMap.end()) {
    return false;
  }

  auto& targetMap = itSource->second;
  auto itTarget = targetMap.find(targetId);

  if (itTarget == targetMap.end()) {
    return false;
  }

  double weight = itTarget->second.weight;

  targetMap.erase(itTarget);
  if (targetMap.empty()) {
    m_nodeLinkMap.erase(itSource);
  } else if (m_numAggregatedLinks > 0) {
    // Guard the unsigned counter: it's 0 for networks with no aggregated
    // duplicates (per-edge occurrence counts aren't tracked, so this can't be
    // exact, but must not underflow now that numAggregatedLinks() is exposed).
    --m_numAggregatedLinks;
  }

  --m_numLinks;
  m_sumLinkWeight -= weight;
  m_totalLinkWeightAdded -= weight;

  if (sourceId == targetId) {
    --m_numSelfLinksFound;
    --m_numSelfLinks;
    m_sumSelfLinkWeight -= weight;
  }

  m_outWeights[sourceId] -= weight;

  return true;
}

bool StateNetwork::undirectedToDirected()
{
  // Collect links in separate data structure to not risk iterating newly added links
  if (!m_config.isUndirectedFlow()) {
    return false;
  }
  // Reads/extends the nested map; materialize it from CSR for mode-A networks
  // (no-op for the multilayer sub-networks that already build the map).
  ensureMapBuild();
  std::deque<StateLink> oppositeLinks;
  for (auto& linkIt : m_nodeLinkMap) {
    unsigned int sourceId = linkIt.first;
    const auto& subLinks = linkIt.second;
    for (auto& subIt : subLinks) {
      unsigned int targetId = subIt.first;
      if (targetId == sourceId) {
        continue; // Self-links are treated as directed on undirected networks
      }
      double weight = subIt.second.weight;
      oppositeLinks.emplace_back(targetId, sourceId, weight);
    }
  }
  for (auto& link : oppositeLinks) {
    addLink(link.source, link.target, link.weight);
  }
  return true;
}

void StateNetwork::clearLinks()
{
  // Release ALL link storage so none of it survives into the optimize phase.
  // releaseInputLinksIfCli() calls this before runTrials(); the pre-CSR build
  // freed its nested map here, so the CSR (and outWeights) must be freed too
  // (vector::clear keeps capacity, so swap with empty).
  NodeLinkMap().swap(m_nodeLinkMap);
  std::vector<LinkTriple>().swap(m_linkBuffer);
  std::vector<unsigned int>().swap(m_nodeIds);
  std::vector<unsigned int>().swap(m_linkOffsets);
  std::vector<unsigned int>().swap(m_linkTargets);
  std::vector<double>().swap(m_linkWeights);
  std::vector<double>().swap(m_linkFlows);
  std::map<unsigned int, double>().swap(m_outWeights);
  // Mark finalized (CSR is "built" -- it's just empty) and zero the raw count so a
  // lazy finalize triggered after clearLinks() can't rebuild a bogus CSR or compute
  // m_numAggregatedLinks = m_rawLinkCount - 0. m_numLinks/m_numAggregatedLinks keep
  // their pre-clear values so a post-clear numLinks() still reports the count (as
  // before); forEachLink() yields nothing by design (empty node list).
  m_rawLinkCount = 0;
  m_linksFinalized = true;
}

void StateNetwork::clear()
{
  m_nodes.clear();
  NodeLinkMap().swap(m_nodeLinkMap);
  std::vector<LinkTriple>().swap(m_linkBuffer);
  std::vector<unsigned int>().swap(m_nodeIds);
  std::vector<unsigned int>().swap(m_linkOffsets);
  std::vector<unsigned int>().swap(m_linkTargets);
  std::vector<double>().swap(m_linkWeights);
  std::vector<double>().swap(m_linkFlows);
  m_physNodes.clear();
  m_outWeights.clear();
  m_names.clear();

  m_linksFinalized = false;
  m_rawLinkCount = 0;
  m_useMapBuild = false; // reset build mode so a reused instance defaults to mode A
  m_haveDirectedInput = false;
  m_haveMemoryInput = false;
  m_numStateNodesFound = 0;
  m_numPhysicalNodesFound = 0;
  m_numLinks = 0;
  m_numSelfLinksFound = 0;
  m_sumLinkWeight = 0.0;
  m_numSelfLinks = 0;
  m_sumSelfLinkWeight = 0.0;
  m_numAggregatedLinks = 0;
  m_totalLinkWeightAdded = 0.0;
  m_numLinksIgnoredByWeightThreshold = 0;
  m_totalLinkWeightIgnored = 0.0;
}

void StateNetwork::writeStateNetwork(const std::string& filename) const
{
  if (filename.empty())
    throw std::runtime_error("writeStateNetwork called with empty filename");

  SafeOutFile outFile(filename, std::ios_base::out, m_config.overwriteOutput());

  outFile << "# v" << INFOMAP_VERSION << "\n"
          << "# ./Infomap " << m_config.parsedString << "\n";

  if (!m_names.empty()) {
    outFile << "*Vertices\n";
    for (auto& nameIt : m_names) {
      auto& physId = nameIt.first;
      auto& name = nameIt.second;
      outFile << physId << " \"" << name << "\"\n";
    }
  }

  outFile << "*States\n";
  outFile << "# stateId physicalId\n";
  for (const auto& nodeIt : nodes()) {
    const auto& node = nodeIt.second;
    outFile << node.id << " " << node.physicalId;
    // Optional name
    if (!node.name.empty())
      outFile << " \"" << node.name << "\"";
    outFile << "\n";
  }

  outFile << "*Links\n";
  forEachLink([&](unsigned int s, unsigned int t, double weight, double&) {
    outFile << nodeId(s) << " " << nodeId(t) << " " << weight << "\n";
  });
  outFile.commit();
}

void StateNetwork::writePajekNetwork(const std::string& filename, bool printFlow) const
{
  if (filename.empty())
    throw std::runtime_error("writePajekNetwork called with empty filename");

  SafeOutFile outFile(filename, std::ios_base::out, m_config.overwriteOutput());

  outFile << "# v" << INFOMAP_VERSION << "\n"
          << "# ./Infomap " << m_config.parsedString << "\n";
  if (haveMemoryInput())
    outFile << "# State network as physical network\n";

  outFile << "*Vertices\n";
  outFile << "#id name " << (printFlow ? "flow" : "weight") << "\n";
  for (const auto& nodeIt : nodes()) {
    const auto& node = nodeIt.second;
    outFile << node.id << " \"";
    // Name, default to id
    const auto& nameIt = haveMemoryInput() ? m_names.end() : m_names.find(node.id);
    if (nameIt != m_names.end())
      outFile << nameIt->second;
    else
      outFile << node.id;
    outFile << "\" " << (printFlow ? node.flow : node.weight) << "\n";
  }

  outFile << (m_config.printAsUndirected() ? "*Edges" : "*Arcs") << "\n";
  outFile << "#source target " << (printFlow ? "flow" : "weight") << "\n";
  forEachLink([&](unsigned int s, unsigned int t, double weight, double& flow) {
    outFile << nodeId(s) << " " << nodeId(t) << " " << (printFlow ? flow : weight) << "\n";
  });
  outFile.commit();
}

std::pair<StateNetwork::NodeMap::iterator, bool> StateNetwork::addStateNodeWithAutogeneratedId(unsigned int physId)
{
  // Keys sorted with std::less comparator, so last key is the largest
  unsigned int stateId = m_nodes.empty() ? 0 : m_nodes.crbegin()->first + 1;
  return addStateNode(stateId, physId);
}

std::pair<StateNetwork::NodeMap::iterator, bool> StateNetwork::addStateNodeWithDeterministicId(unsigned int physId, unsigned int layerId, unsigned int numLayersLog2)
{
  // The shift has to fit in 32 bits. Only layerId was checked (Network.cpp), so a
  // physical id at or above this bound wrapped and produced a state id already
  // taken by a smaller id in the same layer -- 2147483748 << 2 is 400, and so is
  // 100 << 2. Network::addMultilayerNode then aliased the second physical node
  // onto the first, so a node disappeared from the partition and every flow moved,
  // at exit 0 with no warning. With --matchable-multilayer-ids 2 the safe range is
  // only [0, 2^30), which ordinary ids can exceed.
  const unsigned int shift = numLayersLog2 + 1;
  // Its own error: with shift >= 32 the bound below is 0, so every id would be
  // reported as "must be below 0", which describes neither the cause nor a fix.
  // Config::adaptDefaults rejects a largest-layer-id this big, but a caller that
  // builds a Network directly never goes through it.
  if (shift >= 32) {
    const unsigned int maxLargestLayerId = Config::maxMatchableMultilayerIds;
    throw std::runtime_error(fmt::format(FMT_STRING("Matchable multilayer ids need to shift the physical id left by {} bits, which does not fit a 32-bit id. Reduce the largest layer id to at most {}."), shift, maxLargestLayerId));
  }

  const unsigned int maxPhysId = 1u << (32 - shift);
  if (physId >= maxPhysId) {
    throw std::runtime_error(fmt::format(FMT_STRING("Physical node id {} is too large for matchable multilayer ids: the id is shifted left by {} bits to make room for the layer, so ids must be below {}. Use smaller physical ids, or drop --matchable-multilayer-ids."), physId, shift, maxPhysId));
  }

  unsigned int stateId = physId << shift | layerId;
  return addStateNode(stateId, physId);
}

unsigned int StateNetwork::indexOfId(unsigned int id) const
{
  // Precondition: `id` must exist in m_nodeIds (i.e. be one of nodes()). Every
  // link endpoint is addNode'd before finalize and all callers query ids drawn
  // from nodes(), so this holds; the result is then a valid index in [0,numNodes).
  // For an absent id this returns the lower_bound insertion point (possibly
  // numNodes), which the CSR accessors would read out of bounds -- don't call it
  // with an unknown id.
  return static_cast<unsigned int>(
      std::lower_bound(m_nodeIds.begin(), m_nodeIds.end(), id) - m_nodeIds.begin());
}

void StateNetwork::finalizeLinks()
{
  // Idempotent: build the consumed CSR exactly once. A second call would rebuild
  // from a freed buffer (mode A -> empty CSR) or re-read flow-less map entries
  // (mode B -> wiped flows), so mutations after finalize go through definalize()
  // / addLinkToMap() which reset m_linksFinalized instead of rebuilding here.
  if (m_linksFinalized) {
    return;
  }
  if (m_useMapBuild) {
    buildCsrFromMap(); // mode B (multilayer): sorted+aggregated map -> CSR
  } else {
    buildCsrFromBuffer(); // mode A (first-order): sort + merge the flat buffer -> CSR
  }
  m_linksFinalized = true;
}

void StateNetwork::buildCsrFromBuffer()
{
  // Node index = sorted position in m_nodes (every link endpoint was addNode'd
  // during addLink), so index->id is direct and id->index is a binary search.
  m_nodeIds.clear();
  m_nodeIds.reserve(m_nodes.size());
  for (const auto& n : m_nodes) {
    m_nodeIds.push_back(n.first);
  }

  // Sort an index permutation rather than the 16-byte triples: std::sort is
  // in-place (no N-element temporary, unlike stable_sort), and the index array
  // is 4 B/link instead of 16, cutting the finalize transient peak. The
  // comparator orders by (source, target) and breaks ties on the arrival index,
  // a total order whose ties preserve arrival order -- so the left-to-right
  // weight merge below sums duplicates exactly as the old map's '+=' did,
  // bit-for-bit identical.
  const std::size_t numTriples = m_linkBuffer.size();
  std::vector<unsigned int> order(numTriples);
  for (std::size_t i = 0; i < numTriples; ++i) {
    order[i] = static_cast<unsigned int>(i);
  }
  const auto& buf = m_linkBuffer;
  std::sort(order.begin(), order.end(), [&buf](unsigned int a, unsigned int b) {
    if (buf[a].source != buf[b].source) return buf[a].source < buf[b].source;
    if (buf[a].target != buf[b].target) return buf[a].target < buf[b].target;
    return a < b; // arrival-order tiebreak keeps the weight sum bit-exact
  });

  const unsigned int numNodes = static_cast<unsigned int>(m_nodeIds.size());
  m_linkOffsets.assign(numNodes + 1, 0);
  m_linkTargets.clear();
  m_linkWeights.clear();
  m_linkTargets.reserve(numTriples);
  m_linkWeights.reserve(numTriples);

  std::size_t k = 0;
  unsigned int curSourceIndex = 0;
  while (k < numTriples) {
    const unsigned int src = buf[order[k]].source;
    const unsigned int tgt = buf[order[k]].target;
    double w = buf[order[k]].weight;
    std::size_t m = k + 1;
    while (m < numTriples && buf[order[m]].source == src && buf[order[m]].target == tgt) {
      w += buf[order[m]].weight; // arrival order via the index tiebreak
      ++m;
    }
    const unsigned int srcIdx = indexOfId(src);
    while (curSourceIndex < srcIdx) {
      m_linkOffsets[++curSourceIndex] = static_cast<unsigned int>(m_linkTargets.size());
    }
    m_linkTargets.push_back(indexOfId(tgt));
    m_linkWeights.push_back(w);
    k = m;
  }
  while (curSourceIndex < numNodes) {
    m_linkOffsets[++curSourceIndex] = static_cast<unsigned int>(m_linkTargets.size());
  }

  // Free the build buffer + permutation before sizing the flow array, to keep
  // the transient peak down.
  std::vector<LinkTriple>().swap(m_linkBuffer);
  std::vector<unsigned int>().swap(order);

  m_numLinks = static_cast<unsigned int>(m_linkTargets.size());
  m_numAggregatedLinks = m_rawLinkCount - m_numLinks;
  m_linkFlows.assign(m_linkTargets.size(), 0.0);

  // m_outWeights is deliberately NOT derived in mode A: first-order consumers
  // don't read it (only the multilayer mode-B expansion does, where it's kept
  // eager in addLinkToMap). Leaving it empty avoids a ~per-source std::map.
  m_outWeights.clear();
}

void StateNetwork::definalize()
{
  // Move CSR rows back into the flat buffer so building can continue. Merged
  // weights become single occurrences; a later finalizeLinks() re-sorts and
  // re-merges them (and any newly appended links).
  // NOTE: interleaving addLink with finalizing reads (numLinks()/outWeights()/
  // forEachLink()) re-runs the full sort+merge on every addLink-after-finalize.
  // The parse path appends all links first and finalizes once, so this is cheap
  // there; callers doing repeated query/append cycles pay O(L log L) per append.
  m_linkBuffer.clear();
  m_linkBuffer.reserve(m_linkTargets.size()); // final size is known -> single allocation
  m_rawLinkCount = 0;
  for (unsigned int s = 0; s < m_nodeIds.size(); ++s) {
    for (unsigned int e = m_linkOffsets[s]; e < m_linkOffsets[s + 1]; ++e) {
      m_linkBuffer.push_back({ m_nodeIds[s], m_nodeIds[m_linkTargets[e]], m_linkWeights[e] });
      ++m_rawLinkCount;
    }
  }
  m_linkOffsets.clear();
  m_linkTargets.clear();
  m_linkWeights.clear();
  m_linkFlows.clear();
  m_linksFinalized = false;
}

void StateNetwork::buildCsrFromMap()
{
  m_nodeIds.clear();
  m_nodeIds.reserve(m_nodes.size());
  for (const auto& n : m_nodes) {
    m_nodeIds.push_back(n.first);
  }

  const unsigned int numNodes = static_cast<unsigned int>(m_nodeIds.size());
  m_linkOffsets.assign(numNodes + 1, 0);
  m_linkTargets.clear();
  m_linkWeights.clear();
  m_linkFlows.clear();
  m_linkTargets.reserve(m_numLinks);
  m_linkWeights.reserve(m_numLinks);
  m_linkFlows.reserve(m_numLinks);

  unsigned int curSourceIndex = 0;
  unsigned int linkCount = 0;
  for (const auto& node : m_nodeLinkMap) {
    const unsigned int srcIdx = indexOfId(node.first);
    while (curSourceIndex < srcIdx) {
      m_linkOffsets[++curSourceIndex] = linkCount;
    }
    for (const auto& link : node.second) {
      m_linkTargets.push_back(indexOfId(link.first));
      m_linkWeights.push_back(link.second.weight);
      m_linkFlows.push_back(link.second.flow);
      ++linkCount;
    }
  }
  while (curSourceIndex < numNodes) {
    m_linkOffsets[++curSourceIndex] = linkCount;
  }
}

void StateNetwork::ensureMapBuild()
{
  // Materialize the nested map (mode B) from the consumed CSR so the map-based
  // mutation APIs (removeLink, undirectedToDirected) work on a mode-A network.
  // No-op when already building the map (e.g. multilayer sub-networks), so those
  // paths are unchanged. The network then stays in map-build mode; a later
  // finalizeLinks() rebuilds the CSR via buildCsrFromMap().
  if (m_useMapBuild) {
    return;
  }
  ensureFinalized();
  NodeLinkMap().swap(m_nodeLinkMap);
  m_outWeights.clear();
  for (unsigned int s = 0; s < m_nodeIds.size(); ++s) {
    if (m_linkOffsets[s] == m_linkOffsets[s + 1]) {
      continue; // dangling: the historic nested map only held sources with out-links
    }
    const unsigned int srcId = m_nodeIds[s];
    auto& outLinks = m_nodeLinkMap[srcId];
    for (unsigned int e = m_linkOffsets[s]; e < m_linkOffsets[s + 1]; ++e) {
      LinkData linkData(m_linkWeights[e]);
      linkData.flow = m_linkFlows[e];
      outLinks[m_nodeIds[m_linkTargets[e]]] = linkData;
      m_outWeights[srcId] += m_linkWeights[e];
    }
  }
  m_useMapBuild = true;
  m_linksFinalized = false; // map is now the build rep; CSR will be rebuilt on finalize
}

void StateNetwork::deriveOutWeightsIfNeeded()
{
  // Mode A doesn't keep m_outWeights during build (first-order consumers don't
  // read it). Derive it on demand from CSR when the public outWeights() getter is
  // called, so the bound API still returns per-source out-weights. Cleared on
  // every (re)finalize, so it stays consistent.
  if (m_useMapBuild || !m_outWeights.empty()) {
    return;
  }
  ensureFinalized(); // mode A: m_numLinks/CSR are only valid after finalize
  for (unsigned int s = 0; s < m_nodeIds.size(); ++s) {
    double w = 0.0;
    for (unsigned int e = m_linkOffsets[s]; e < m_linkOffsets[s + 1]; ++e) {
      w += m_linkWeights[e];
    }
    if (w != 0.0) {
      m_outWeights[m_nodeIds[s]] = w;
    }
  }
}

} // namespace infomap
