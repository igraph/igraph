/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef CLUSTER_MAP_H_
#define CLUSTER_MAP_H_

#include <cstdint>
#include <string>
#include <map>
#include <vector>

namespace infomap {

using Path = std::vector<unsigned int>; // 1-based indexing

using NodePath = std::pair<unsigned int, Path>;

using NodePaths = std::vector<NodePath>;

enum class TreeLeafIdType : std::uint8_t {
  physical,
  state,
};

struct TreePath {
  unsigned int nodeId = 0; // Either a physical id or a state id, see idType
  Path path;
  TreeLeafIdType idType = TreeLeafIdType::physical;
};

using TreePaths = std::vector<TreePath>;

class ClusterMap {
public:
  void readClusterData(const std::string& filename, bool includeFlow = false, const std::map<unsigned int, std::map<unsigned int, unsigned int>>* layerNodeToStateId = nullptr);

  const std::map<unsigned int, unsigned int>& clusterIds() const noexcept
  {
    return m_clusterIds;
  }

  const TreePaths& treePaths() const noexcept { return m_treePaths; }

  const std::string& extension() const noexcept { return m_extension; }

  TreeLeafIdType treeLeafIdType() const noexcept { return m_treeLeafIdType; }

private:
  void readTree(const std::string& filename, bool includeFlow, const std::map<unsigned int, std::map<unsigned int, unsigned int>>* layerNodeToStateId = nullptr);
  void readClu(const std::string& filename, bool includeFlow, const std::map<unsigned int, std::map<unsigned int, unsigned int>>* layerNodeToStateId = nullptr);

  std::map<unsigned int, unsigned int> m_clusterIds;
  std::map<unsigned int, double> m_flowData;
  TreePaths m_treePaths;
  std::string m_extension;
  bool m_isHigherOrder = false;
  bool m_hasTreeLeafIdType = false;
  TreeLeafIdType m_treeLeafIdType = TreeLeafIdType::physical;
};

} // namespace infomap

#endif // CLUSTER_MAP_H_
