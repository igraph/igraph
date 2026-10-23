/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef THREAD_CONFIG_H_
#define THREAD_CONFIG_H_

#include <cstdint>

namespace infomap {

// A thread count of 0 means "absent / not specified" (C++14 — no std::optional).
struct ThreadSources {
  unsigned int explicitThreads = 0; // --num-threads N (0 = auto)
  unsigned int infomapEnv = 0; // INFOMAP_NUM_THREADS
  unsigned int slurmCpusPerTask = 0; // SLURM_CPUS_PER_TASK
  unsigned int ompEnv = 0; // OMP_NUM_THREADS
  unsigned int cpusetCount = 0; // sched_getaffinity
  unsigned int hardwareConcurrency = 1; // always >= 1, final fallback
};

enum class ThreadSource : std::uint8_t { Explicit,
                                         InfomapEnv,
                                         Slurm,
                                         Omp,
                                         Cpuset,
                                         Hardware };

struct ThreadBudget {
  unsigned int threads = 1;
  ThreadSource source = ThreadSource::Hardware;
};

// Pure: first non-zero source in precedence order wins, clamped to >= 1.
ThreadBudget resolveThreadBudget(const ThreadSources& sources);

// Human-readable source label for logs / JSON (e.g. "SLURM_CPUS_PER_TASK").
const char* threadSourceName(ThreadSource source);

// Production wiring: reads env vars + sched_getaffinity + hardware_concurrency.
ThreadSources readThreadSourcesFromEnv();

} // namespace infomap

#endif // THREAD_CONFIG_H_
