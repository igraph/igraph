/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef CLUSTER_MAP_H_
#define CLUSTER_MAP_H_

#include <string>
#include <map>
#include <vector>
#include <deque>

namespace infomap {

using Path = std::deque<unsigned int>; // 1-based indexing

using NodePath = std::pair<unsigned int, Path>;

using NodePaths = std::vector<NodePath>;

class ClusterMap {
public:
  void readClusterData(const std::string& filename, bool includeFlow = false, const std::map<unsigned int, std::map<unsigned int, unsigned int>>* layerNodeToStateId = nullptr);

  const std::map<unsigned int, unsigned int>& clusterIds() const noexcept
  {
    return m_clusterIds;
  }

  const NodePaths& nodePaths() const noexcept { return m_nodePaths; }

  const std::string& extension() const noexcept { return m_extension; }

private:
  void readTree(const std::string& filename, bool includeFlow, const std::map<unsigned int, std::map<unsigned int, unsigned int>>* layerNodeToStateId = nullptr);
  void readClu(const std::string& filename, bool includeFlow, const std::map<unsigned int, std::map<unsigned int, unsigned int>>* layerNodeToStateId = nullptr);

  std::map<unsigned int, unsigned int> m_clusterIds;
  std::map<unsigned int, double> m_flowData;
  NodePaths m_nodePaths;
  std::string m_extension;
  bool m_isHigherOrder = false;
};

} // namespace infomap

#endif // CLUSTER_MAP_H_
