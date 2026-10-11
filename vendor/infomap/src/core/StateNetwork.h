/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef STATE_NETWORK_H_
#define STATE_NETWORK_H_

#include "../io/Config.h"
#include <string>
#include <map>
#include <utility>
#include <vector>
#include <set>
#include <utility>

namespace infomap {

class StateNetwork {
public:
  struct StateNode {
    unsigned int id = 0;
    unsigned int physicalId = 0;
    std::string name;
    unsigned int layerId = 0;
    double weight = 1.0;
    double flow = 0.0;
    double enterFlow = 0.0;
    double exitFlow = 0.0;
    double teleFlow = 0.0;
    double intraLayerTeleFlow = 0.0;
    double intraLayerTeleWeight = 0.0;

    StateNode(unsigned int id = 0) : id(id), physicalId(id) {}

    StateNode(unsigned int id, unsigned int physicalId) : id(id), physicalId(physicalId) {}

    StateNode(unsigned int id, unsigned int physicalId, std::string name) : id(id), physicalId(physicalId), name(std::move(name)) {}

    bool operator==(const StateNode& rhs) const { return id == rhs.id; }
    bool operator!=(const StateNode& rhs) const { return id != rhs.id; }
    bool operator<(const StateNode& rhs) const { return id < rhs.id; }
  };

  struct PhysNode {
    unsigned int physId = 0;
    double weight = 1.0;
    PhysNode(unsigned int physId) : physId(physId) {}
    PhysNode(unsigned int physId, double weight) : physId(physId), weight(weight) {}
    PhysNode(double weight = 1.0) : weight(weight) {}
  };

  struct LinkData {
    double weight = 1.0;
    double flow = 0.0;

    LinkData(double weight = 1.0) : weight(weight) {}

    LinkData& operator+=(double w)
    {
      weight += w;
      return *this;
    }
  };

  struct StateLink {
    StateLink(unsigned int sourceIndex = 0, unsigned int targetIndex = 0, double weight = 0.0)
        : source(sourceIndex),
          target(targetIndex),
          weight(weight),
          flow(weight) {}

    unsigned int source;
    unsigned int target;
    double weight;
    double flow;
  };

  // Unique state id to state node
  using NodeMap = std::map<unsigned int, StateNode>;
  // Links keyed by state id. StateNode keys were only ever compared by id
  // (and never updated past construction), so plain ids give the same
  // iteration order at a fraction of the per-link footprint.
  using OutLinkMap = std::map<unsigned int, LinkData>;
  using NodeLinkMap = std::map<unsigned int, OutLinkMap>;

  // Flat build buffer for first-order (mode A) link intake.
  struct LinkTriple {
    unsigned int source;
    unsigned int target;
    double weight;
  };

protected:
  friend class FlowCalculator;
  // Config
  Config m_config;
  // Network
  bool m_haveDirectedInput = false;
  bool m_haveMemoryInput = false;
  bool m_higherOrderInputMethodCalled = false;
  NodeMap m_nodes; // Nodes indexed by state id (equal physical id for first-order networks)
  // --- Link build representations (one active per instance) ---
  NodeLinkMap m_nodeLinkMap; // mode B (multilayer) build rep
  mutable std::vector<LinkTriple> m_linkBuffer; // mode A (first-order) build rep
  bool m_useMapBuild = false; // true => mode B
  mutable bool m_linksFinalized = false;
  mutable unsigned int m_rawLinkCount = 0; // pre-aggregation occurrences (mode A)
  // --- Consumed CSR representation (valid after finalizeLinks) ---
  mutable std::vector<unsigned int> m_nodeIds; // sorted unique ids; index->id
  mutable std::vector<unsigned int> m_linkOffsets; // size numNodes+1
  mutable std::vector<unsigned int> m_linkTargets; // dense target indices
  mutable std::vector<double> m_linkWeights;
  mutable std::vector<double> m_linkFlows;
  unsigned int m_numStateNodesFound = 0;
  double m_sumNodeWeight = 0.0;
  unsigned int m_numLinks = 0;
  double m_sumLinkWeight = 0.0;
  unsigned int m_numSelfLinksFound = 0;
  unsigned int m_numSelfLinks = 0;
  double m_sumSelfLinkWeight = 0.0;
  unsigned int m_numAggregatedLinks = 0;
  double m_totalLinkWeightAdded = 0.0;
  unsigned int m_numLinksIgnoredByWeightThreshold = 0;
  double m_totalLinkWeightIgnored = 0.0;
  std::map<unsigned int, double> m_outWeights;
  bool m_haveNodeWeights = false;
  bool m_haveStateNodeWeights = false;
  bool m_haveFileInput = false;
  // Attributes
  std::map<unsigned int, std::string> m_names;
  std::map<unsigned int, PhysNode> m_physNodes;
  // Unique physical-id count, tracked independently of m_physNodes so the map
  // can be released once redundant (postProcessInputData, first-order input)
  // while numPhysicalNodes() stays correct.
  unsigned int m_numPhysicalNodesFound = 0;

  // Flow calculation outcome. Written by FlowCalculator, which is a friend, so the
  // power iteration's result outlives the calculator and can reach the run reports --
  // it used to exist only as a console line, invisible to every automated consumer
  // (#899).
  bool m_haveFlowConvergence = false;
  bool m_flowConverged = true;
  unsigned int m_flowIterations = 0;
  double m_flowError = 0.0;

  // Bipartite
  unsigned int m_bipartiteStartId = 0;

  // Multilayer
  std::set<unsigned int> m_layers;

public:
  StateNetwork() {}
  StateNetwork(Config config) : m_config(std::move(config)) {}
  virtual ~StateNetwork() = default;

  StateNetwork(const StateNetwork&) = delete;
  StateNetwork& operator=(const StateNetwork&) = delete;
  StateNetwork(StateNetwork&&) = delete;
  StateNetwork& operator=(StateNetwork&&) = delete;

  // Config
  void setConfig(const Config& config) { m_config = config; }

  // Mutators
  std::pair<NodeMap::iterator, bool> addStateNode(const StateNode& node);
  std::pair<NodeMap::iterator, bool> addStateNode(unsigned int id, unsigned int physId);
  std::pair<NodeMap::iterator, bool> addStateNode(unsigned int id, unsigned int physId, std::string name);
  std::pair<NodeMap::iterator, bool> addNode(unsigned int id);
  std::pair<NodeMap::iterator, bool> addNode(unsigned int id, std::string name);
  std::pair<NodeMap::iterator, bool> addNode(unsigned int id, double weight);
  std::pair<NodeMap::iterator, bool> addNode(unsigned int id, std::string, double weight);
  PhysNode& addPhysicalNode(unsigned int physId);
  PhysNode& addPhysicalNode(unsigned int physId, double weight);
  PhysNode& addPhysicalNode(unsigned int physId, const std::string& name);
  PhysNode& addPhysicalNode(unsigned int physId, double weight, const std::string& name);
  std::pair<std::map<unsigned int, std::string>::iterator, bool> addName(unsigned int id, const std::string&);
  bool addLink(unsigned int sourceId, unsigned int targetId, double weight = 1.0);
  bool addLink(unsigned int sourceId, unsigned int targetId, unsigned long weight);
  void addLinks(const std::vector<unsigned int>& sourceIds, const std::vector<unsigned int>& targetIds, const std::vector<double>& weights);

  /**
   * Remove link
   * Note: It will not remove nodes if they become dangling
   */
  bool removeLink(unsigned int sourceId, unsigned int targetId);

  // Expand each undirected link to two opposite directed links
  bool undirectedToDirected();

  /**
   * Clear all network data and reset to default state.
   */
  virtual void clear();

  /**
   * Clear link data but keep node data.
   */
  virtual void clearLinks();

  // Getters
  const NodeMap& nodes() const { return m_nodes; }
  unsigned int numNodes() const { return m_nodes.size(); }
  unsigned int numPhysicalNodes() const { return m_numPhysicalNodesFound; }
  double sumNodeWeight() const { return m_sumNodeWeight; }
#ifndef SWIG
  // Mode B (multilayer) build representation, consumed only by the multilayer
  // expansion in Network. Hidden from the bindings: external callers read links
  // through getLinks()/getLinkResults(), which now serve the consumed CSR store.
  const NodeLinkMap& nodeLinkMap() const { return m_nodeLinkMap; }
  NodeLinkMap& nodeLinkMap() { return m_nodeLinkMap; }
#endif

#ifndef SWIG
  // --- CSR link storage access (internal; valid after finalize). Hidden from the
  // bindings: these expose CSR index mechanics, so external callers read links
  // through getLinks()/getLinkResults(), which serve the consumed CSR store. ---
  void finalizeLinks();
  void ensureFinalized() const
  {
    if (!m_linksFinalized) const_cast<StateNetwork*>(this)->finalizeLinks();
  }
  unsigned int nodeId(unsigned int index) const { return m_nodeIds[index]; }
  unsigned int indexOfId(unsigned int id) const;
  unsigned int outDegree(unsigned int index) const { return m_linkOffsets[index + 1] - m_linkOffsets[index]; }
  bool isDangling(unsigned int index) const { return outDegree(index) == 0; }
  // fn(srcIndex, targetIndex, weight, double& flow) in (src,tgt) order. flow is
  // writable (CSR arrays are mutable) so FlowCalculator can write flow back.
  template <typename Fn>
  void forEachLink(Fn&& fn) const
  {
    ensureFinalized();
    for (unsigned int s = 0; s < m_nodeIds.size(); ++s)
      for (unsigned int e = m_linkOffsets[s]; e < m_linkOffsets[s + 1]; ++e)
        fn(s, m_linkTargets[e], m_linkWeights[e], m_linkFlows[e]);
  }
#endif

  // Mode A (first-order and the main multilayer network) defers dedup to
  // finalizeLinks(), so these counts are only valid afterwards -- ensure it. An
  // early read (e.g. an in-memory numLinks() before the multilayer expansion
  // populates the main network) just finalizes early; the next addLink re-opens
  // the buffer via definalize(). Mode B (the
  // transient per-layer networks) keeps the counts eager and must NOT finalize:
  // those networks are read through nodeLinkMap()/outWeights() and discarded,
  // never consumed as CSR. The ensureFinalized() call is guarded so SWIG never
  // parses a reference to the hidden method (the compiled library the wrapper
  // calls still lazy-finalizes).
  unsigned int numAggregatedLinks() const
  {
#ifndef SWIG
    if (!m_useMapBuild) ensureFinalized();
#endif
    return m_numAggregatedLinks;
  }
  unsigned int numLinks() const
  {
#ifndef SWIG
    if (!m_useMapBuild) ensureFinalized();
#endif
    return m_numLinks;
  }
  double sumLinkWeight() const { return m_sumLinkWeight; }
  unsigned int numSelfLinks() const { return m_numSelfLinks; }
  double sumSelfLinkWeight() const { return m_sumSelfLinkWeight; }
  // Use convention of counting self-links only once, treating them as directed
  double sumWeightedDegree() const { return 2 * sumLinkWeight() - sumSelfLinkWeight(); }
  unsigned int sumDegree() const { return 2 * numLinks() - numSelfLinks(); }
  std::map<unsigned int, double>& outWeights()
  {
#ifndef SWIG
    deriveOutWeightsIfNeeded(); // mode A derives it lazily; mode B keeps it eager
#endif
    return m_outWeights;
  }
  std::map<unsigned int, std::string>& names() { return m_names; }
  const std::map<unsigned int, std::string>& names() const { return m_names; }
  //! Whether a power iteration ran and recorded its outcome.
  bool haveFlowConvergence() const { return m_haveFlowConvergence; }
  //! Whether the power iteration's final error is within the flow tolerance. True also
  //! when tolerance is reached on the last allowed iteration, and when no iteration ran.
  bool flowConverged() const { return m_flowConverged; }
  unsigned int flowIterations() const { return m_flowIterations; }
  double flowError() const { return m_flowError; }

  bool haveNodeWeights() const { return m_haveNodeWeights; }
  bool haveStateNodeWeights() const { return m_haveStateNodeWeights; }
  bool haveFileInput() const { return m_haveFileInput; }

  virtual const std::map<unsigned int, std::vector<int>>& metaData() const = 0;

  bool haveDirectedInput() const { return m_haveDirectedInput; }
  bool haveMemoryInput() const { return m_haveMemoryInput; }
  bool higherOrderInputMethodCalled() const { return m_higherOrderInputMethodCalled; }
  // Bipartite
  bool isBipartite() const { return m_bipartiteStartId > 0; }
  unsigned int bipartiteStartId() const { return m_bipartiteStartId; }
  void setBipartiteStartId(unsigned int value) { m_bipartiteStartId = value; }
  // Multilayer
#ifndef SWIG
  unsigned int numLayers() const { return m_layers.size(); }
  const std::set<unsigned int>& layers() const { return m_layers; }
#endif

  /**
   * Write state network to file.
   */
  void writeStateNetwork(const std::string& filename) const;

  /**
   * Write state network as first-order Pajek network, where
   * state nodes are treated as physical nodes.
   * For a non-memory input, the state nodes are equivalent to
   * physical nodes.
   */
  void writePajekNetwork(const std::string& filename, bool printFlow = false) const;

protected:
  std::pair<NodeMap::iterator, bool> addStateNodeWithAutogeneratedId(unsigned int physId);

  std::pair<NodeMap::iterator, bool> addStateNodeWithDeterministicId(unsigned int physId, unsigned int layerId, unsigned int numLayersLog2);

  // Build the consumed CSR arrays. Mode B (multilayer): from the already-sorted,
  // already-aggregated m_nodeLinkMap. Mode A (first-order): from m_linkBuffer.
  void buildCsrFromMap();
  void buildCsrFromBuffer();
  // Mode B eager aggregation into the nested map (the pre-CSR build path).
  bool addLinkToMap(unsigned int sourceId, unsigned int targetId, double weight);
  // Materialize the nested map from CSR so map-based mutation APIs work on a
  // mode-A network (no-op in map-build mode).
  void ensureMapBuild();
  // Lazily fill m_outWeights from CSR for mode-A networks when outWeights() is read.
  void deriveOutWeightsIfNeeded();
  // Move CSR rows back into the flat buffer so a finalized network can keep
  // building (mode A accumulate / addLink-after-run).
  void definalize();
};

} // namespace infomap

#endif // STATE_NETWORK_H_
