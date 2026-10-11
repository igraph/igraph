/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef NETWORK_H_
#define NETWORK_H_

#include "Config.h"
#ifndef SWIG
#include "ClusterMap.h"
#endif
#include "../core/StateNetwork.h"

#include <limits>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace infomap {

struct LayerNode;
class NetworkIntakeAdapter;

class Network : public StateNetwork {
  friend class NetworkIntakeAdapter;

private:
  // Multilayer
  std::map<unsigned int, Network> m_networks; // intra-layer links
  std::map<LayerNode, std::map<unsigned int, double>> m_interLinks;
  // { layer -> { physId -> stateId }}
  std::map<unsigned int, std::map<unsigned int, unsigned int>> m_layerNodeToStateId;
  std::map<unsigned int, std::map<unsigned int, double>> m_sumIntraOutWeight;
  unsigned int m_numInterLayerLinks = 0;
  unsigned int m_numIntraLayerLinks = 0;
  unsigned int m_maxNodeIdInIntraLayerNetworks = 0;

  unsigned int m_multilayerStateIdBitShift = 0;

  // Meta data
  std::map<unsigned int, std::vector<int>> m_metaData;
  unsigned int m_numMetaDataColumns = 0;

#ifndef SWIG
  // Embedded initial partition paths (from JSON nodes[].path / states[].path).
  // Hierarchical, reusing the same TreePaths model as a .tree cluster file.
  // Internal plumbing (JSON parser -> initTrialPartition); not a binding API.
  TreePaths m_initialPartitionPaths;
#endif

public:
  Network() { init(); }
  explicit Network(const Config& config) : StateNetwork(config) { init(); }
  explicit Network(const std::string& flags) : StateNetwork(Config(flags)) { init(); }
  ~Network() override = default;

  Network(const Network&) = delete;
  Network& operator=(const Network&) = delete;
  Network(Network&&) = delete;
  Network& operator=(Network&&) = delete;

  void setConfig(const Config& config)
  {
    StateNetwork::setConfig(config);
    updateDerivedConfig();
  }

  void clear() override;

  /**
   * Parse network data from file and generate network
   * @param filename input network
   * @param accumulate add to possibly existing network data (default), else clear before.
   */
  virtual void readInputData(std::string filename = "", bool accumulate = true);

  /**
   * Init categorical meta data on all nodes from a file with the following format:
   * # nodeId metaData
   * 1 1
   * 2 1
   * 3 2
   * 4 2
   * 5 3
   * @param filename input filename for metadata
   */
  virtual void readMetaData(const std::string& filename);

  unsigned int numMetaDataColumns() const { return m_numMetaDataColumns; }
  const std::map<unsigned int, std::vector<int>>& metaData() const override { return m_metaData; }

#ifndef SWIG
  void addInitialPartitionPath(unsigned int nodeId, const Path& path, TreeLeafIdType idType);
  const TreePaths& initialPartitionPaths() const { return m_initialPartitionPaths; }
#endif

  bool isMultilayerNetwork() const { return !m_layerNodeToStateId.empty(); }
  const std::map<unsigned int, std::map<unsigned int, unsigned int>>& layerNodeToStateId() const { return m_layerNodeToStateId; }

  void postProcessInputData();
  void generateStateNetworkFromMultilayer();
  void generateStateNetworkFromMultilayerWithInterLinks();
  void generateStateNetworkFromMultilayerWithSimulatedInterLinks();
  void simulateInterLayerLinks();

  /**
   * Create state node corresponding to this multilayer node if not already exist
   * @return state node id
   */
  unsigned int addMultilayerNode(unsigned int layerId, unsigned int physicalId, double weight = 1.0);

  unsigned int addMultilayerNode(unsigned int stateId, unsigned int layerId, unsigned int physicalId, double weight);

  void addMultilayerLink(unsigned int layer1, unsigned int n1, unsigned int layer2, unsigned int n2, double weight);
  void addMultilayerLink(unsigned int stateId1, unsigned int layer1, unsigned int n1, unsigned int stateId2, unsigned int layer2, unsigned int n2, double weight);
  void addMultilayerLinks(const std::vector<unsigned int>& sourceLayerIds,
                          const std::vector<unsigned int>& sourceNodeIds,
                          const std::vector<unsigned int>& targetLayerIds,
                          const std::vector<unsigned int>& targetNodeIds,
                          const std::vector<double>& weights);

private:
  void generateStateNetworkFromMultilayerWithSimulatedInterLinksBasedOnLayerSimilarity();
  void generateStateNetworkFromMultilayerWithSimulatedInterLinksBasedOnNodeStrength();
  void generateStateNetworkFromMultilayerWithSimulatedInterLinksBasedOnNodeStrengthRegularized();

public:
  /**
   * Create an intra-layer link
   */
  void addMultilayerIntraLink(unsigned int layer, unsigned int n1, unsigned int n2, double weight);
  void addMultilayerIntraLinks(const std::vector<unsigned int>& layerIds,
                               const std::vector<unsigned int>& sourceNodeIds,
                               const std::vector<unsigned int>& targetNodeIds,
                               const std::vector<double>& weights);

  /**
   * Create links between (layer1,n) and (layer2,m) for all m connected to n in layer 2.
   * The weight is distributed proportionally.
   * TODO: This is done later..
   */
  void addMultilayerInterLink(unsigned int layer1, unsigned int n, unsigned int layer2, double interWeight);
  void addMultilayerInterLinks(const std::vector<unsigned int>& sourceLayerIds,
                               const std::vector<unsigned int>& nodeIds,
                               const std::vector<unsigned int>& targetLayerIds,
                               const std::vector<double>& weights);

  void addMetaData(unsigned int nodeId, int meta);

  void addMetaData(unsigned int nodeId, const std::vector<int>& metaData);

private:
  void init();
  void updateDerivedConfig();

  static double calculateJensenShannonDivergence(bool& intersect, const OutLinkMap& layer1OutLinks, double sumOutLinkWeightLayer1, const OutLinkMap& layer2OutLinks, double sumOutLinkWeightLayer2);

  void printSummary();
};

struct LayerNode {
  unsigned int layer, node;
  explicit LayerNode(unsigned int layer = 0, unsigned int node = 0) : layer(layer), node(node) {}

  bool operator<(const LayerNode other) const
  {
    return layer == other.layer ? node < other.node : layer < other.layer;
  }
};

} // namespace infomap

#endif // NETWORK_H_
