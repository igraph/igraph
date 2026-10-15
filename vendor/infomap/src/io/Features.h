#ifndef INFOMAP_FEATURES_H
#define INFOMAP_FEATURES_H

#include <string>
#include <vector>

namespace infomap {

inline std::vector<std::string> enabledFeatures()
{
  std::vector<std::string> features;
#ifdef INFOMAP_ENABLED_FEATURES
  const std::string encodedFeatures { INFOMAP_ENABLED_FEATURES };
  std::string::size_type begin = 0;
  while (begin <= encodedFeatures.size()) {
    const auto end = encodedFeatures.find(',', begin);
    const auto feature = encodedFeatures.substr(begin, end - begin);
    if (!feature.empty())
      features.push_back(feature);
    if (end == std::string::npos)
      break;
    begin = end + 1;
  }
#endif
  return features;
}

} // namespace infomap

#endif // INFOMAP_FEATURES_H
