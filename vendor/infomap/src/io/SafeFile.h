/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef SAFEFILE_H_
#define SAFEFILE_H_

#include "InfomapError.h"
#include "../utils/format.h"
#include <atomic>
#include <cerrno>
#include <ctime>
#include <fstream>
#include <ios>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#include <process.h>
#else
#include <unistd.h>
#endif

namespace infomap {

using std::ifstream;
using std::ofstream;

/**
 * A wrapper for the C++ file stream class that automatically closes
 * the file stream when the destructor is called. Allocate it on the
 * stack to have it automatically closed when going out of scope.
 *
 * Note:
 * In C++, the only code that can be guaranteed to be executed after an
 * exception is thrown are the destructors of objects residing on the stack.
 *
 * You can exploit that fact to avoid resource leaks by tying all resources
 * to the lifespan of an object allocated on the stack. This technique is
 * called Resource Acquisition Is Initialization (RAII).
 *
 */
class SafeInFile : public ifstream {
public:
  SafeInFile(const std::string& filename, ios_base::openmode mode = ios_base::in)
      : ifstream(filename, mode)
  {
    if (fail())
      throw std::runtime_error(fmt::format(FMT_STRING("Cannot open file '{}'. Check that the path points to a file and that you have read permissions."), filename));
  }

  ~SafeInFile() override
  {
    if (is_open())
      close();
  }
};

inline bool pathExists(const std::string& path)
{
  struct stat info;
  return stat(path.c_str(), &info) == 0;
}

// The process working directory, or an empty string if it cannot be read.
// The buffer doubles on ERANGE up to 64 KiB, so a deep CI path is not a silent
// failure -- callers that anchor relative paths to this would otherwise quietly
// lose the anchoring. The cap is deliberate rather than a limit to raise: it sits
// past every platform's own path limit, Windows long paths included at 32767, so
// reaching it means the call is failing for a reason growing cannot fix, and
// giving up beats doubling until the allocation does.
inline std::string currentWorkingDirectory()
{
  for (std::size_t size = 512; size <= 65536; size *= 2) {
    std::vector<char> buffer(size);
#ifdef _WIN32
    const char* result = _getcwd(buffer.data(), static_cast<int>(buffer.size()));
#else
    const char* result = getcwd(buffer.data(), buffer.size());
#endif
    if (result != nullptr)
      return std::string(buffer.data());
    if (errno != ERANGE)
      break;
  }
  return std::string();
}

inline bool isDirectory(const std::string& path)
{
  struct stat info;
  if (stat(path.c_str(), &info) != 0)
    return false;
  return (info.st_mode & S_IFDIR) != 0;
}

inline std::string stripTrailingPathSeparators(std::string path)
{
  while (path.size() > 1 && (path.back() == '/' || path.back() == '\\'))
    path.pop_back();
  return path;
}

#ifdef _WIN32
inline bool isWindowsPathRoot(const std::string& path)
{
  // Drive root, e.g. "C:".
  if (path.size() == 2 && path[1] == ':')
    return true;
  // UNC root, e.g. "\\server" or "\\server\share": leading "\\" with at most one
  // path separator after it. Recursing/creating above this point is invalid.
  if (path.size() >= 2 && (path[0] == '\\' || path[0] == '/') && (path[1] == '\\' || path[1] == '/')) {
    const auto rest = path.substr(2);
    const auto firstSep = rest.find_first_of("/\\");
    if (firstSep == std::string::npos)
      return true;
    return rest.find_first_of("/\\", firstSep + 1) == std::string::npos;
  }
  return false;
}
#endif

inline void ensureDirectoryExists(const std::string& directory)
{
  if (directory.empty() || isDirectory(directory))
    return;

  const auto normalized = stripTrailingPathSeparators(directory);
#ifdef _WIN32
  // Never try to create a drive or UNC share root; assume it exists. If it truly
  // does not, creating a child below will fail with a clear error instead.
  if (isWindowsPathRoot(normalized))
    return;
#endif
  const auto separator = normalized.find_last_of("/\\");
  if (separator != std::string::npos) {
    const auto parent = normalized.substr(0, separator);
    if (!parent.empty())
      ensureDirectoryExists(parent);
  }

#ifdef _WIN32
  const auto created = _mkdir(normalized.c_str());
#else
  const auto created = mkdir(normalized.c_str(), 0777);
#endif
  if (created != 0 && !isDirectory(normalized)) {
    throw InfomapError(ExitCode::OutputError, fmt::format(FMT_STRING("Can't create output directory '{}'."), directory));
  }
}

inline std::string makeTemporaryOutputFilename(const std::string& filename)
{
  static std::atomic<unsigned int> counter { 0 };
#ifdef _WIN32
  const auto processId = static_cast<unsigned int>(_getpid());
#else
  const auto processId = static_cast<unsigned int>(getpid());
#endif
  return fmt::format(FMT_STRING("{}.tmp.{}.{}.{}"), filename, processId, static_cast<unsigned long long>(std::time(nullptr)), counter.fetch_add(1));
}

class SafeOutFile : public ofstream {
public:
  SafeOutFile(const std::string& filename, ios_base::openmode mode = ios_base::out, bool overwrite = true)
      : m_filename(filename),
        m_temporaryFilename(makeTemporaryOutputFilename(filename)),
        m_overwrite(overwrite)
  {
    if (!m_overwrite && pathExists(m_filename)) {
      throw InfomapError(ExitCode::OutputError, fmt::format(FMT_STRING("Output file already exists: '{}'"), m_filename));
    }

    open(m_temporaryFilename, mode);
    if (fail())
      throw InfomapError(ExitCode::OutputError, fmt::format(FMT_STRING("Cannot open file '{}'. Check that the directory you are writing to exists and that you have write permissions."), filename));
  }

  ~SafeOutFile() override
  {
    if (is_open())
      close();
    if (!m_committed)
      std::remove(m_temporaryFilename.c_str());
  }

  void commit()
  {
    if (m_committed)
      return;

    flush();
    if (fail()) {
      close();
      std::remove(m_temporaryFilename.c_str());
      throw InfomapError(ExitCode::OutputError, fmt::format(FMT_STRING("Cannot write file '{}'."), m_filename));
    }

    close();
    if (fail()) {
      std::remove(m_temporaryFilename.c_str());
      throw InfomapError(ExitCode::OutputError, fmt::format(FMT_STRING("Cannot close file '{}'."), m_filename));
    }

    if (!m_overwrite && pathExists(m_filename)) {
      std::remove(m_temporaryFilename.c_str());
      throw InfomapError(ExitCode::OutputError, fmt::format(FMT_STRING("Output file already exists: '{}'"), m_filename));
    }

#ifdef _WIN32
    if (m_overwrite && pathExists(m_filename)) {
      std::remove(m_filename.c_str());
    }
#endif

    if (std::rename(m_temporaryFilename.c_str(), m_filename.c_str()) != 0) {
      std::remove(m_temporaryFilename.c_str());
      throw InfomapError(ExitCode::OutputError, fmt::format(FMT_STRING("Cannot replace file '{}'."), m_filename));
    }

    m_committed = true;
  }

private:
  std::string m_filename;
  std::string m_temporaryFilename;
  bool m_overwrite = true;
  bool m_committed = false;
};

inline bool isDirectoryWritable(const std::string& dir)
{
  std::string path = fmt::format(FMT_STRING("{}_1nf0m4p_.tmp"), dir);
  bool ok = true;
  try {
    SafeOutFile out(path);
  } catch (const std::runtime_error&) {
    ok = false;
  }
  if (ok)
    std::remove(path.c_str());
  return ok;
}

} // namespace infomap

#endif // SAFEFILE_H_
