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

// Identity of one input file: path as given, size, mtime, and a hash of the
// whole content (FNV-1a 64, one streaming pass). Empty when no file is known
// (an in-memory network built through the bindings).
struct InputIdentity {
  std::string path;
  unsigned long long size = 0;
  long long mtime = 0;
  std::string hash;

  bool known() const { return !hash.empty(); }
};

// The identity of the file at ``path``. An empty path, or one that cannot be
// stat'ed, yields an unknown identity rather than an error: the reader that
// opens the file reports a missing file with its own classified error, and
// identity capture must not pre-empt that. A file that exists but cannot be
// read still throws.
InputIdentity inputIdentity(const std::string& path);
// ``{"path","size","mtime","hash"}`` as JSON text, or ``null`` when unknown.
std::string inputIdentityJson(const InputIdentity& identity);

// The identities a run captured when it started, for the artifacts it writes.
// An empty fingerprint means "compute it from the config".
struct RunIdentities {
  InputIdentity input;
  InputIdentity clusterData;
  InputIdentity metaData;
  std::string configFingerprint;
};

std::string canonicalConfigJson(const Config& config);
std::string configFingerprint(const Config& config);
std::string inputFingerprintJson(const std::string& path);
// Node-stable content hash of the input network (no size/mtime). Used by the
// distributed-trial merge guard to confirm all shards ran on the same network.
std::string networkFingerprint(const std::string& path);
// Writes the identities the run captured. The Config-only overload re-reads
// the paths in the config and exists for callers with no captured run.
std::string runManifestJson(const Config& config, const RunIdentities& identities);
std::string runManifestJson(const Config& config);

} // namespace infomap

#endif // RUNMETADATA_H_
