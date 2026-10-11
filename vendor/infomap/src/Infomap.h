/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFOMAP_H_
#define INFOMAP_H_

#include "core/InfomapBase.h"
#include "io/Config.h"

#include <string>
#include <utility>
#include <map>
#include <vector>
#include <stdexcept>

namespace infomap {

#ifndef SWIGPYTHON
struct LinkResult {
  LinkResult() = default;
  LinkResult(unsigned int source, unsigned int target, double weight, double flow)
      : source(source),
        target(target),
        weight(weight),
        flow(flow) {}

  unsigned int source = 0;
  unsigned int target = 0;
  double weight = 0.0;
  double flow = 0.0;
};
#endif

// Single-traversal bulk node data for the Python (and JS) bindings. Excluded
// from the R binding (#ifndef SWIGR): R wraps the struct as an undocumented S4
// class, which trips an R CMD check WARNING, and R has no use for the Python
// extractor. The C++ compiler and the Python/JS SWIG passes leave SWIGR
// undefined, so the struct stays available there. Parallel std::vector columns
// reuse the existing vector_uint/vector_double SWIG templates; the ragged path
// is stored CSR-style (flat values + lengths).
#ifndef SWIGR
struct NodeData {
  std::vector<unsigned int> node_id; // physicalId
  std::vector<unsigned int> state_id; // stateId
  std::vector<unsigned int> module_id; // moduleId at `level`
  std::vector<double> flow;
  std::vector<unsigned int> depth;
  std::vector<unsigned int> layer_id;
  std::vector<unsigned int> child_index;
  std::vector<double> modular_centrality;
  std::vector<unsigned int> path_flat; // CSR values: all path entries concatenated
  std::vector<unsigned int> path_len; // CSR lengths: per-node path length
};
#endif

// Wrapper class for the Python API
struct InfomapWrapper : public InfomapBase {
public:
  InfomapWrapper() {}
  InfomapWrapper(const std::string& flags) : InfomapBase(flags) {}
  InfomapWrapper(const Config& conf) : InfomapBase(conf) {}
  ~InfomapWrapper() override = default;

  // ===================================================
  // Wrapper methods
  // ===================================================

  void readInputData(std::string filename = "", bool accumulate = true)
  {
    // The bindings never set networkFile, so this is the only record of what a
    // run read (#1026). Recorded after the read, so a failed one leaves no
    // claim; cleared when accumulate is false, since the network is replaced.
    // Several accumulated files are a union no single path identifies, so it
    // gets no identity rather than the last file's.
    //
    // The identity is hashed here, from the bytes just parsed -- not at run
    // start. A file edited or deleted between read_file() and run() would
    // otherwise be hashed in the name of a partition computed from the old
    // content.
    const std::string path = filename;
    if (!accumulate) {
      // Before the call, not after: readInputData() clears the network first
      // (Network.cpp), so the previous sources stop describing it the moment
      // the read begins, whether or not the parse then succeeds.
      m_readInputPaths.clear();
      m_networkHasInMemoryContent = false;
      m_inputIdentityFromRead = InputIdentity();
    }
    try {
      m_network.readInputData(std::move(filename), accumulate);
    } catch (...) {
      // The parser writes straight into the network, so a malformed file can
      // leave part of itself behind before it throws -- and nothing tells us
      // how far it got. Node and link counts do not: a bad file whose ids
      // overlap what is already there moves neither. So a failed read makes the
      // network unidentifiable, and stays that way until accumulate=false
      // replaces it. That costs a false negative on a read that failed before
      // touching anything, such as a mistyped path, which is the safe
      // direction: an artifact naming no file beats one naming the wrong file.
      m_readInputPaths.clear();
      m_networkHasInMemoryContent = true;
      m_inputIdentityFromRead = InputIdentity();
      throw;
    }
    if (!path.empty())
      m_readInputPaths.push_back(path);
    m_inputIdentityFromRead = (m_readInputPaths.size() == 1 && !m_networkHasInMemoryContent)
        ? infomap::inputIdentity(m_readInputPaths.front())
        : InputIdentity();
  }

  // See provenanceJson() in io/Output.h: everything a caller needs to log or
  // reproduce this run, as one JSON string the bindings parse (#1026).
  std::string provenanceJson() const { return infomap::provenanceJson(*this); }

#ifndef SWIG
  // Every in-memory build call below goes through this first: a network the
  // caller extended after reading a file is a mixture that file alone does not
  // describe, so the file stops identifying the run (#1026).
  //
  // Guarded from SWIG like baseSeed() above: it is bookkeeping the mutators do
  // for themselves, not an operation to offer the bindings. The generated
  // wrapper still compiles against the real header, where it is visible.
  void noteInMemoryMutation()
  {
    m_readInputPaths.clear();
    m_inputIdentityFromRead = InputIdentity();
    // Sticky until accumulate=false replaces the network: a read that follows
    // an in-memory build with the default accumulate=true keeps both sources,
    // so the file does not describe the network then either.
    m_networkHasInMemoryContent = true;
  }
#endif

  // Reached from Python's remove_link() and set_meta_data(), which used to call
  // straight through to the network and so left a file identity standing over a
  // network the call had changed (#1026).
  bool removeLink(unsigned int sourceId, unsigned int targetId)
  {
    noteInMemoryMutation();
    return m_network.removeLink(sourceId, targetId);
  }

  void addMetaData(unsigned int nodeId, int meta)
  {
    noteInMemoryMutation();
    m_network.addMetaData(nodeId, meta);
  }

  void addNode(unsigned int id)
  {
    noteInMemoryMutation();
    m_network.addNode(id);
  }
  void addNode(unsigned int id, std::string name)
  {
    noteInMemoryMutation();
    m_network.addNode(id, std::move(name));
  }
  void addNode(unsigned int id, double weight)
  {
    noteInMemoryMutation();
    m_network.addNode(id, weight);
  }
  void addNode(unsigned int id, std::string name, double weight)
  {
    noteInMemoryMutation();
    m_network.addNode(id, std::move(name), weight);
  }

  void addName(unsigned int id, const std::string& name)
  {
    noteInMemoryMutation();
    m_network.addName(id, name);
  }
  std::string getName(unsigned int id) const
  {
    auto& names = m_network.names();
    auto it = names.find(id);
    return it != names.end() ? it->second : "";
  }

  const std::map<unsigned int, std::string>& getNames() const { return m_network.names(); }

  // State-id -> name for higher-order (state/memory) networks. Physical names
  // live in getNames(), keyed by physical id; state nodes carry their own
  // optional name parsed from the *States section (StateNode::name). Only
  // non-empty names are returned, so a first-order network (no per-state names)
  // yields an empty map and callers fall back to the physical name.
  std::map<unsigned int, std::string> getStateNames() const
  {
    std::map<unsigned int, std::string> stateNames;
    for (const auto& entry : m_network.nodes()) {
      if (!entry.second.name.empty())
        stateNames[entry.first] = entry.second.name;
    }
    return stateNames;
  }

  // Name of a single state node, or "" if the id is unknown or unnamed. The
  // scalar counterpart of getStateNames(), mirroring getName(); the R binding
  // walks state nodes calling this since SWIG-R exposes std::map opaquely.
  std::string getStateName(unsigned int stateId) const
  {
    const auto& nodes = m_network.nodes();
    auto it = nodes.find(stateId);
    return it != nodes.end() ? it->second.name : "";
  }

  void addPhysicalNode(unsigned int id, const std::string& name = "")
  {
    noteInMemoryMutation();
    m_network.addPhysicalNode(id, name);
  }
  void addStateNode(unsigned int id, unsigned int physId)
  {
    noteInMemoryMutation();
    m_network.addStateNode(id, physId);
  }
  void addStateNode(unsigned int id, unsigned int physId, std::string name)
  {
    noteInMemoryMutation();
    m_network.addStateNode(id, physId, std::move(name));
  }

  void addLink(unsigned int sourceId, unsigned int targetId, double weight = 1.0)
  {
    noteInMemoryMutation();
    m_network.addLink(sourceId, targetId, weight);
  }
  void addLink(unsigned int sourceId, unsigned int targetId, unsigned long weight)
  {
    noteInMemoryMutation();
    m_network.addLink(sourceId, targetId, weight);
  }
  void addLinks(const std::vector<unsigned int>& sourceIds, const std::vector<unsigned int>& targetIds, const std::vector<double>& weights)
  {
    noteInMemoryMutation();
    m_network.addLinks(sourceIds, targetIds, weights);
  }
  void addMultilayerLink(unsigned int layer1, unsigned int n1, unsigned int layer2, unsigned int n2, double weight = 1.0)
  {
    noteInMemoryMutation();
    m_network.addMultilayerLink(layer1, n1, layer2, n2, weight);
  }
  void addMultilayerLinks(const std::vector<unsigned int>& sourceLayerIds,
                          const std::vector<unsigned int>& sourceNodeIds,
                          const std::vector<unsigned int>& targetLayerIds,
                          const std::vector<unsigned int>& targetNodeIds,
                          const std::vector<double>& weights)
  {
    noteInMemoryMutation();
    m_network.addMultilayerLinks(sourceLayerIds, sourceNodeIds, targetLayerIds, targetNodeIds, weights);
  }
  void addMultilayerIntraLink(unsigned int layer, unsigned int n1, unsigned int n2, double weight)
  {
    noteInMemoryMutation();
    m_network.addMultilayerIntraLink(layer, n1, n2, weight);
  }
  void addMultilayerIntraLinks(const std::vector<unsigned int>& layerIds,
                               const std::vector<unsigned int>& sourceNodeIds,
                               const std::vector<unsigned int>& targetNodeIds,
                               const std::vector<double>& weights)
  {
    noteInMemoryMutation();
    m_network.addMultilayerIntraLinks(layerIds, sourceNodeIds, targetNodeIds, weights);
  }
  void addMultilayerInterLink(unsigned int layer1, unsigned int n, unsigned int layer2, double interWeight)
  {
    noteInMemoryMutation();
    m_network.addMultilayerInterLink(layer1, n, layer2, interWeight);
  }
  void addMultilayerInterLinks(const std::vector<unsigned int>& sourceLayerIds,
                               const std::vector<unsigned int>& nodeIds,
                               const std::vector<unsigned int>& targetLayerIds,
                               const std::vector<double>& weights)
  {
    noteInMemoryMutation();
    m_network.addMultilayerInterLinks(sourceLayerIds, nodeIds, targetLayerIds, weights);
  }

  // Resolve a physical (layer_id, node_id) to the internally generated state id.
  // Argument order is layer-first, consistent with MultilayerNode and the
  // add_multilayer_* / setMultilayerInitialPartition APIs. The network must be
  // built first (state nodes are created on the fly). Throws std::out_of_range
  // if the combination is not present (issue #616).
  unsigned int getMultilayerStateId(unsigned int layerId, unsigned int nodeId) const
  {
    const auto& map = m_network.layerNodeToStateId();
    auto layerIt = map.find(layerId);
    if (layerIt != map.end()) {
      auto nodeIt = layerIt->second.find(nodeId);
      if (nodeIt != layerIt->second.end())
        return nodeIt->second;
    }
    throw std::out_of_range("No state id for (layer " + std::to_string(layerId) + ", node " + std::to_string(nodeId) + "); build the multilayer network before requesting it");
  }

  void setBipartiteStartId(unsigned int startId) { m_network.setBipartiteStartId(startId); }

  std::map<std::pair<unsigned int, unsigned int>, double> getLinks(bool flow) const
  {
    std::map<std::pair<unsigned int, unsigned int>, double> links;

    m_network.forEachLink([&](unsigned int s, unsigned int t, double weight, double& linkFlow) {
      links[{ m_network.nodeId(s), m_network.nodeId(t) }] = flow ? linkFlow : weight;
    });

    return links;
  }

#ifndef SWIGPYTHON
  std::vector<LinkResult> getLinkResults() const
  {
    std::vector<LinkResult> links;
    links.reserve(m_network.numLinks());

    m_network.forEachLink([&](unsigned int s, unsigned int t, double weight, double& linkFlow) {
      links.emplace_back(m_network.nodeId(s), m_network.nodeId(t), weight, linkFlow);
    });

    return links;
  }
#endif

  std::map<unsigned int, unsigned int> getModules(int level = 1, bool states = false)
  {
    if (haveMemory() && !states) {
      throw std::runtime_error("Cannot get modules on higher-order network without states.");
    }
    std::map<unsigned int, unsigned int> modules;
    for (auto it = iterTree(level); !it.isEnd(); ++it) {
      auto& node = *it;
      if (node.isLeaf()) {
        modules[states ? node.stateId : node.physicalId] = it.moduleId();
      }
    }
    return modules;
  }

  // Single-traversal bulk node data; NodeData is defined at namespace scope
  // (above) so SWIG wraps its parallel-vector members for Python. Excluded from
  // the R binding together with NodeData (see the struct's #ifndef SWIGR note).
#ifndef SWIGR
  NodeData getNodeData(int level = 1, bool states = false)
  {
    // Mirror _results.get_nodes: physical iterator for higher-order networks
    // unless state-level data is requested.
    if (haveMemory() && !states) {
      return collectLeafData(iterTreePhysical(level));
    }
    return collectLeafData(iterLeafNodes(level));
  }
#endif

private:
#ifndef SWIGR
  template <typename Iter>
  NodeData collectLeafData(Iter it) const
  {
    NodeData d;
    for (; !it.isEnd(); ++it) {
      auto& node = *it;
      if (!node.isLeaf()) continue;
      d.node_id.push_back(node.physicalId);
      d.state_id.push_back(node.stateId);
      d.module_id.push_back(it.moduleId());
      d.flow.push_back(node.data.flow);
      d.depth.push_back(it.depth());
      d.layer_id.push_back(node.layerId);
      d.child_index.push_back(it.childIndex());
      d.modular_centrality.push_back(it.modularCentrality());
      const auto& p = it.path();
      d.path_len.push_back(static_cast<unsigned int>(p.size()));
      d.path_flat.insert(d.path_flat.end(), p.begin(), p.end());
    }
    return d;
  }
#endif

public:
  using InfomapBase::codelength;
  using InfomapBase::getEntropyRate;
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
  using InfomapBase::getLossyDistortion;
  using InfomapBase::getLossyRate;
#endif
  using InfomapBase::getMultilevelModules;
  using InfomapBase::iterLeafNodes;
  using InfomapBase::iterTree;
  using InfomapBase::run;
};

} // namespace infomap

#endif // INFOMAP_H_
