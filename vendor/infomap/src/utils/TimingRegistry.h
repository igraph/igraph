/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef TIMINGREGISTRY_H_
#define TIMINGREGISTRY_H_

#include "Stopwatch.h"

#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace infomap {

struct TrialTimingRecord {
  unsigned int trial = 0;
  int thread = 0;
  unsigned long seed = 0;
  double timeSec = 0.0;
  double codelength = 0.0;
  unsigned int topModules = 0;
  unsigned int numLevels = 0;
  bool valid = false;
};

class TimingRegistry {
public:
  class Scope {
  public:
    Scope(TimingRegistry& registry, std::string phase)
        : m_registry(&registry),
          m_phase(std::move(phase)),
          m_timer(true),
          m_active(true)
    {
    }

    Scope(Scope&& other) noexcept
        : m_registry(other.m_registry),
          m_phase(std::move(other.m_phase)),
          m_timer(other.m_timer),
          m_active(other.m_active)
    {
      other.m_active = false;
    }

    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;
    Scope& operator=(Scope&&) = delete;

    ~Scope()
    {
      if (m_active && m_registry != nullptr) {
        m_registry->addPhase(m_phase, m_timer.getElapsedTimeInSec());
      }
    }

  private:
    TimingRegistry* m_registry;
    std::string m_phase;
    Stopwatch m_timer;
    bool m_active;
  };

  explicit TimingRegistry(unsigned int numTrials = 0) : m_trials(numTrials) {}

  Scope scope(const std::string& phase) { return Scope(*this, phase); }

  void addPhase(const std::string& phase, double seconds)
  {
    std::lock_guard<std::mutex> lock(m_mutex);
    for (auto& entry : m_phases) {
      if (entry.first == phase) {
        entry.second += seconds;
        return;
      }
    }
    m_phases.emplace_back(phase, seconds);
  }

  void setPhase(const std::string& phase, double seconds)
  {
    std::lock_guard<std::mutex> lock(m_mutex);
    for (auto& entry : m_phases) {
      if (entry.first == phase) {
        entry.second = seconds;
        return;
      }
    }
    m_phases.emplace_back(phase, seconds);
  }

  void resetTrials(unsigned int numTrials)
  {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_trials.assign(numTrials, TrialTimingRecord());
  }

  void recordTrial(unsigned int trialIndex, int thread, unsigned long seed, double timeSec, double codelength, unsigned int topModules, unsigned int numLevels)
  {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (trialIndex >= m_trials.size()) {
      m_trials.resize(trialIndex + 1);
    }
    TrialTimingRecord record;
    record.trial = trialIndex + 1;
    record.thread = thread;
    record.seed = seed;
    record.timeSec = timeSec;
    record.codelength = codelength;
    record.topModules = topModules;
    record.numLevels = numLevels;
    record.valid = true;
    m_trials[trialIndex] = record;
  }

  std::vector<std::pair<std::string, double>> phases() const
  {
    std::lock_guard<std::mutex> lock(m_mutex);
    return m_phases;
  }

  std::vector<TrialTimingRecord> trials() const
  {
    std::lock_guard<std::mutex> lock(m_mutex);
    return m_trials;
  }

private:
  mutable std::mutex m_mutex;
  std::vector<std::pair<std::string, double>> m_phases;
  std::vector<TrialTimingRecord> m_trials;
};

} // namespace infomap

#endif // TIMINGREGISTRY_H_
