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

    StateNode(unsigned int id = 0) : id(id), physicalId(id) { }

    StateNode(unsigned int id, unsigned int physicalId) : id(id), physicalId(physicalId) { }

    StateNode(unsigned int id, unsigned int physicalId, std::string name) : id(id), physicalId(physicalId), name(std::move(name)) { }

    bool operator==(const StateNode& rhs) const { return id == rhs.id; }
    bool operator!=(const StateNode& rhs) const { return id != rhs.id; }
    bool operator<(const StateNode& rhs) const { return id < rhs.id; }
  };

  struct PhysNode {
    unsigned int physId = 0;
    double weight = 1.0;
    PhysNode(unsigned int physId) : physId(physId) { }
    PhysNode(unsigned int physId, double weight) : physId(physId), weight(weight) { }
    PhysNode(double weight = 1.0) : weight(weight) { }
  };

  struct LinkData {
    double weight = 1.0;
    double flow = 0.0;

    LinkData(double weight = 1.0) : weight(weight) { }

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
          flow(weight) { }

    unsigned int source;
    unsigned int target;
    double weight;
    double flow;
  };

  // Unique state id to state node
  using NodeMap = std::map<unsigned int, StateNode>;
  using OutLinkMap = std::map<StateNode, LinkData>;
  using NodeLinkMap = std::map<StateNode, OutLinkMap>;

protected:
  friend class FlowCalculator;
  // Config
  Config m_config;
  // Network
  bool m_haveDirectedInput = false;
  bool m_haveMemoryInput = false;
  bool m_higherOrderInputMethodCalled = false;
  NodeMap m_nodes; // Nodes indexed by state id (equal physical id for first-order networks)
  NodeLinkMap m_nodeLinkMap;
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

  // Bipartite
  unsigned int m_bipartiteStartId = 0;

public:
  StateNetwork() : m_config(Config()) { }
  StateNetwork(Config config) : m_config(std::move(config)) { }
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
  unsigned int numPhysicalNodes() const { return m_physNodes.size(); }
  double sumNodeWeight() const { return m_sumNodeWeight; }
  const NodeLinkMap& nodeLinkMap() const { return m_nodeLinkMap; }
  NodeLinkMap& nodeLinkMap() { return m_nodeLinkMap; }
  unsigned int numLinks() const { return m_numLinks; }
  double sumLinkWeight() const { return m_sumLinkWeight; }
  unsigned int numSelfLinks() const { return m_numSelfLinks; }
  double sumSelfLinkWeight() const { return m_sumSelfLinkWeight; }
  // Use convention of counting self-links only once, treating them as directed
  double sumWeightedDegree() const { return 2 * sumLinkWeight() - (m_config.isUndirectedFlow() ? sumSelfLinkWeight() : 0); }
  unsigned int sumDegree() const { return 2 * numLinks() - (m_config.isUndirectedFlow() ? numSelfLinks() : 0); }
  std::map<unsigned int, double>& outWeights() { return m_outWeights; }
  std::map<unsigned int, std::string>& names() { return m_names; }
  const std::map<unsigned int, std::string>& names() const { return m_names; }
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
};

} // namespace infomap

#endif // STATE_NETWORK_H_
