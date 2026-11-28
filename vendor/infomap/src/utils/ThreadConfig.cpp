/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "ThreadConfig.h"

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <string>
#include <thread>

#if defined(__linux__)
#include <sched.h>
#endif

namespace infomap {

ThreadBudget resolveThreadBudget(const ThreadSources& s)
{
  const auto pick = [](unsigned int value, ThreadSource source, ThreadBudget& out) {
    if (value > 0) {
      out.threads = value;
      out.source = source;
      return true;
    }
    return false;
  };

  ThreadBudget budget;
  if (pick(s.explicitThreads, ThreadSource::Explicit, budget)
      || pick(s.infomapEnv, ThreadSource::InfomapEnv, budget)
      || pick(s.slurmCpusPerTask, ThreadSource::Slurm, budget)
      || pick(s.ompEnv, ThreadSource::Omp, budget)
      || pick(s.cpusetCount, ThreadSource::Cpuset, budget)) {
    return budget;
  }
  budget.threads = std::max(1u, s.hardwareConcurrency);
  budget.source = ThreadSource::Hardware;
  return budget;
}

const char* threadSourceName(ThreadSource source)
{
  switch (source) {
  case ThreadSource::Explicit:
    return "--num-threads";
  case ThreadSource::InfomapEnv:
    return "INFOMAP_NUM_THREADS";
  case ThreadSource::Slurm:
    return "SLURM_CPUS_PER_TASK";
  case ThreadSource::Omp:
    return "OMP_NUM_THREADS";
  case ThreadSource::Cpuset:
    return "cpuset";
  case ThreadSource::Hardware:
    return "hardware_concurrency";
  }
  // unreachable — all enumerators handled above
  return "unknown";
}

namespace {
  std::string withoutSurroundingSpace(const std::string& text)
  {
    const auto first = text.find_first_not_of(" \t\n\v\f\r");
    if (first == std::string::npos) {
      return std::string();
    }
    const auto last = text.find_last_not_of(" \t\n\v\f\r");
    return text.substr(first, last - first + 1);
  }

  // A single positive count, or 0 when the text is not one. Trailing garbage
  // ("4x") is rejected rather than parsed loosely; surrounding space is not
  // garbage, since a value assembled by a shell script easily carries it.
  unsigned int parseCount(const std::string& raw)
  {
    const std::string text = withoutSurroundingSpace(raw);
    if (text.empty()) {
      return 0;
    }
    try {
      std::size_t consumed = 0;
      const long value = std::stol(text, &consumed);
      if (consumed == text.size() && value > 0 && value <= static_cast<long>(std::numeric_limits<unsigned int>::max())) {
        return static_cast<unsigned int>(value);
      }
      return 0;
    } catch (...) {
      return 0;
    }
  }

  unsigned int parseEnvCount(const char* name)
  {
    const char* raw = std::getenv(name);
    return raw == nullptr ? 0 : parseCount(std::string(raw));
  }

  // OMP_NUM_THREADS is specified as a comma-separated list of positive integers,
  // one per nesting level, so "2,1" is a valid way to ask for two threads in the
  // outermost parallel region. Reading it as a single integer rejected the whole
  // value and left the budget to fall through to the cpuset or the core count --
  // which then got pushed into the runtime, raising the count above the request
  // rather than merely ignoring it. Take the first element, which governs the
  // region we care about. Every element still has to parse, so a malformed list
  // is treated as unset instead of half-read.
  unsigned int parseOmpEnvCount(const char* name)
  {
    const char* raw = std::getenv(name);
    if (raw == nullptr) {
      return 0;
    }

    const std::string text(raw);
    unsigned int firstCount = 0;
    std::size_t start = 0;
    while (true) {
      const auto separator = text.find(',', start);
      const std::string element = text.substr(start, separator == std::string::npos ? std::string::npos : separator - start);
      const unsigned int count = parseCount(element);
      if (count == 0) {
        return 0;
      }
      if (firstCount == 0) {
        firstCount = count;
      }
      if (separator == std::string::npos) {
        break;
      }
      start = separator + 1;
    }
    return firstCount;
  }

  unsigned int readCpusetCount()
  {
#if defined(__linux__)
    // Fast path: the fixed-size set covers machines up to CPU_SETSIZE (typically 1024) CPUs.
    {
      cpu_set_t set;
      CPU_ZERO(&set);
      if (sched_getaffinity(0, sizeof(set), &set) == 0) {
        return static_cast<unsigned int>(CPU_COUNT(&set));
      }
      if (errno != EINVAL) {
        return 0; // a non-size error won't be fixed by a larger mask
      }
    }
    // Large machines (> CPU_SETSIZE CPUs) make the fixed-size call fail with EINVAL;
    // grow a dynamically-allocated mask until sched_getaffinity succeeds.
    for (unsigned int numCpus = CPU_SETSIZE * 2; numCpus <= CPU_SETSIZE * 1024u; numCpus *= 2) {
      cpu_set_t* set = CPU_ALLOC(numCpus);
      if (set == nullptr) {
        break;
      }
      const std::size_t size = CPU_ALLOC_SIZE(numCpus);
      CPU_ZERO_S(size, set);
      const int rc = sched_getaffinity(0, size, set);
      if (rc == 0) {
        const unsigned int count = static_cast<unsigned int>(CPU_COUNT_S(size, set));
        CPU_FREE(set);
        return count;
      }
      const int err = errno;
      CPU_FREE(set);
      if (err != EINVAL) {
        break;
      }
    }
#endif
    return 0;
  }
} // namespace

ThreadSources readThreadSourcesFromEnv()
{
  ThreadSources s;
  s.infomapEnv = parseEnvCount("INFOMAP_NUM_THREADS");
  s.slurmCpusPerTask = parseEnvCount("SLURM_CPUS_PER_TASK");
  s.ompEnv = parseOmpEnvCount("OMP_NUM_THREADS");
  s.cpusetCount = readCpusetCount();
  s.hardwareConcurrency = std::max(1u, std::thread::hardware_concurrency());
  return s;
}

} // namespace infomap
