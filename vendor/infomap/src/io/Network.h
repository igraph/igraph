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
#include "../core/StateNetwork.h"

#include <string>
#include <map>
#include <utility>
#include <vector>
#include <set>
#include <utility>
#include <limits>
#include <sstream>
#include <locale>

namespace infomap {

struct LayerNode;

class Network : public StateNetwork {
private:
  // Helpers
  std::istringstream m_extractor;

  // Multilayer
  std::map<unsigned int, Network> m_networks; // intra-layer links
  std::map<LayerNode, std::map<unsigned int, double>> m_interLinks;
  // { layer -> { physId -> stateId }}
  std::map<unsigned int, std::map<unsigned int, unsigned int>> m_layerNodeToStateId;
  std::map<unsigned int, std::map<unsigned int, double>> m_sumIntraOutWeight;
  std::set<unsigned int> m_layers;
  unsigned int m_numInterLayerLinks = 0;
  unsigned int m_numIntraLayerLinks = 0;
  unsigned int m_maxNodeIdInIntraLayerNetworks = 0;

  unsigned int m_multilayerStateIdBitShift = 0;

  // Meta data
  std::map<unsigned int, std::vector<int>> m_metaData;
  unsigned int m_numMetaDataColumns = 0;

  using InsensitiveStringSet = std::set<std::string, io::InsensitiveCompare>;

  std::map<std::string, InsensitiveStringSet> m_ignoreHeadings;
  std::map<std::string, InsensitiveStringSet> m_validHeadings; // {
  // 	{ "pajek", {"*Vertices", "*Edges", "*Arcs"} },
  // 	{ "link-list", {"*Links"} },
  // 	{ "bipartite", {"*Vertices", "*Bipartite"} },
  // 	{ "general", {"*Vertices", "*States", "*Edges", "*Arcs", "*Links", "*Context"} }
  // };

public:
  Network() : StateNetwork() { init(); }
  explicit Network(const Config& config) : StateNetwork(config) { init(); }
  explicit Network(const std::string& flags) : StateNetwork(Config(flags)) { init(); }
  ~Network() override = default;

  Network(const Network&) = delete;
  Network& operator=(const Network&) = delete;
  Network(Network&&) = delete;
  Network& operator=(Network&&) = delete;

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

  void addMultilayerLink(unsigned int layer1, unsigned int n1, unsigned int layer2, unsigned int n2, double weight);

  /**
   * Create an intra-layer link
   */
  void addMultilayerIntraLink(unsigned int layer, unsigned int n1, unsigned int n2, double weight);

  /**
   * Create links between (layer1,n) and (layer2,m) for all m connected to n in layer 2.
   * The weight is distributed proportionally.
   * TODO: This is done later..
   */
  void addMultilayerInterLink(unsigned int layer1, unsigned int n, unsigned int layer2, double interWeight);

  void addMetaData(unsigned int nodeId, int meta);

  void addMetaData(unsigned int nodeId, const std::vector<int>& metaData);

private:
  void init();
  void initValidHeadings();

  void parseNetwork(const std::string& filename);
  void parseNetwork(const std::string& filename, const InsensitiveStringSet& validHeadings, const InsensitiveStringSet& ignoreHeadings, const std::string& startHeading = "");

  // Helper methods

  /**
   * Parse vertices under the heading
   * @return The line after the vertices
   */
  std::string parseVertices(std::ifstream& file, const std::string& heading);
  std::string parseStateNodes(std::ifstream& file, const std::string& heading);

  std::string parseLinks(std::ifstream& file);

  /**
   * Parse multilayer links from a *multilayer section
   */
  std::string parseMultilayerLinks(std::ifstream& file);

  /**
   * Parse multilayer links from an *intra section
   */
  std::string parseMultilayerIntraLinks(std::ifstream& file);

  /**
   * Parse multilayer links from an *inter section
   */
  std::string parseMultilayerInterLinks(std::ifstream& file);

  std::string parseBipartiteLinks(std::ifstream& file, const std::string& heading);

  static std::string ignoreSection(std::ifstream& file, const std::string& heading);

  void parseStateNode(const std::string& line, StateNetwork::StateNode& stateNode);

  /**
   * Parse a string of link data.
   * If no weight data can be extracted, the default value 1.0 will be used.
   * @throws an error if not both node ids can be extracted.
   */
  void parseLink(const std::string& line, unsigned int& n1, unsigned int& n2, double& weight);

  /**
   * Parse a string of multilayer link data.
   * If no weight data can be extracted, the default value 1.0 will be used.
   * @throws an error if not both node and layer ids can be extracted.
   */
  void parseMultilayerLink(const std::string& line, unsigned int& layer1, unsigned int& n1, unsigned int& layer2, unsigned int& n2, double& weight);

  /**
   * Parse a string of intra-multilayer link data.
   * If no weight data can be extracted, the default value 1.0 will be used.
   * @throws an error if not both node and layer ids can be extracted.
   */
  void parseMultilayerIntraLink(const std::string& line, unsigned int& layer, unsigned int& n1, unsigned int& n2, double& weight);

  /**
   * Parse a string of inter-multilayer link data.
   * If no weight data can be extracted, the default value 1.0 will be used.
   * @throws an error if not both node and layer ids can be extracted.
   */
  void parseMultilayerInterLink(const std::string& line, unsigned int& layer1, unsigned int& n, unsigned int& layer2, double& weight);

  static double calculateJensenShannonDivergence(bool& intersect, const OutLinkMap& layer1OutLinks, double sumOutLinkWeightLayer1, const OutLinkMap& layer2OutLinks, double sumOutLinkWeightLayer2);

  void printSummary();
};

struct LayerNode {
  unsigned int layer, node;
  explicit LayerNode(unsigned int layer = 0, unsigned int node = 0) : layer(layer), node(node) { }

  bool operator<(const LayerNode other) const
  {
    return layer == other.layer ? node < other.node : layer < other.layer;
  }
};

} // namespace infomap

#endif // NETWORK_H_
