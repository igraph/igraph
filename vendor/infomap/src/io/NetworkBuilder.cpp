/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "NetworkBuilder.h"
#include "Network.h"
#include "NetworkInputParser.h"
#include "JsonNetworkInputParser.h"
#include "../utils/FileURI.h"
#include "../utils/format.h"

#include <stdexcept>
#include <vector>

namespace infomap {

namespace {

  void addParsedVertex(Network& network, const input::ParsedVertex& vertex)
  {
    if (vertex.hasWeight) {
      network.addPhysicalNode(vertex.id, vertex.weight, vertex.name);
    } else {
      network.addPhysicalNode(vertex.id, vertex.name);
    }
  }

} // namespace

class NetworkIntakeAdapter {
public:
  explicit NetworkIntakeAdapter(Network& network) : m_network(network) {}

  input::NetworkInputOptions inputOptions() const
  {
    return { m_network.m_config.isUndirectedFlow(), m_network.m_config.matchableMultilayerIds };
  }

  void onFileInput() { m_network.m_haveFileInput = true; }

  void onNetworkParsed()
  {
    m_network.postProcessInputData();
  }

  void onPhysicalNode(const input::ParsedVertex& vertex)
  {
    if (vertex.hasWeight) {
      m_network.m_haveNodeWeights = true;
    }
    addParsedVertex(m_network, vertex);
  }

  void onStateNode(const input::ParsedStateNode& parsed)
  {
    const auto& stateNode = parsed.node;
    m_network.m_higherOrderInputMethodCalled = true;
    if (parsed.hasWeight) {
      m_network.m_haveStateNodeWeights = true;
    }
    m_network.addStateNode(stateNode);
    m_network.addPhysicalNode(stateNode.physicalId);
    ++m_network.m_numStateNodesFound;
  }

  void onLink(const input::ParsedLink& link)
  {
    m_network.addLink(link.source, link.target, link.weight);
  }

  void onMultilayerLink(const input::ParsedMultilayerLink& link)
  {
    if (m_network.m_config.regularized) {
      throw std::runtime_error("Regularized multilayer flow requires *Intra/*Inter input; explicit *Multilayer input is not supported with --regularized.");
    }
    m_network.addMultilayerLink(link.sourceLayer, link.sourceNode, link.targetLayer, link.targetNode, link.weight);
  }

  void onMultilayerIntraLink(const input::ParsedMultilayerIntraLink& link)
  {
    m_network.addMultilayerIntraLink(link.layer, link.sourceNode, link.targetNode, link.weight);
  }

  void onMultilayerInterLink(const input::ParsedMultilayerInterLink& link)
  {
    m_network.addMultilayerInterLink(link.sourceLayer, link.node, link.targetLayer, link.weight);
  }

  void onBipartiteStart(unsigned int bipartiteStartId)
  {
    m_network.m_bipartiteStartId = bipartiteStartId;
    m_network.m_config.bipartite = true;
  }

  void onBipartiteLink(const std::string& line, const input::ParsedLink& link)
  {
    bool sourceIsFeature = link.source >= m_network.m_bipartiteStartId;
    bool targetIsFeature = link.target >= m_network.m_bipartiteStartId;
    if (sourceIsFeature == targetIsFeature) {
      throw std::runtime_error(fmt::format(FMT_STRING("Bipartite link '{}' must cross bipartite start id {}."), line, m_network.m_bipartiteStartId));
    }
    m_network.addLink(link.source, link.target, link.weight);
  }

  void onMetaData(unsigned int nodeId, const std::vector<int>& metaData)
  {
    m_network.addMetaData(nodeId, metaData);
  }

  void onInitialPartitionPath(unsigned int id, const std::vector<unsigned int>& path, bool stateLevel)
  {
    m_network.addInitialPartitionPath(id, path, stateLevel ? TreeLeafIdType::state : TreeLeafIdType::physical);
  }

private:
  Network& m_network;
};

void buildNetworkFromInput(Network& network, const std::string& filename)
{
  FileURI { filename, false };

  NetworkIntakeAdapter sink(network);
  sink.onFileInput();
  if (input::looksLikeJsonNetwork(filename)) {
    input::parseJsonNetworkInput(filename, sink, sink.inputOptions());
  } else {
    input::parseNetworkInput(filename, sink, sink.inputOptions());
  }
  sink.onNetworkParsed();
}

void buildMetaDataFromInput(Network& network, const std::string& filename)
{
  NetworkIntakeAdapter sink(network);
  input::parseMetaDataInput(filename, sink);
}

} // namespace infomap
