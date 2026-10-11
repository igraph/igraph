/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef RUNMETADATA_H_
#define RUNMETADATA_H_

#include <string>

namespace infomap {

struct Config;

std::string canonicalConfigJson(const Config& config);
std::string configFingerprint(const Config& config);
std::string inputFingerprintJson(const std::string& path);
// Node-stable content hash of the input network (no size/mtime). Used by the
// distributed-trial merge guard to confirm all shards ran on the same network.
std::string networkFingerprint(const std::string& path);
std::string runManifestJson(const Config& config);

} // namespace infomap

#endif // RUNMETADATA_H_
