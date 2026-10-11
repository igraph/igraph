/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "InfomapBase.h"
#include "InfomapConfig.h"
#include "InfoNode.h"
#include "FlowData.h"
#include "BiasedMapEquation.h"
#include "MemMapEquation.h"
#include "MetaMapEquation.h"
#if INFOMAP_FEATURE_REGULARIZED_MULTILAYER
#include "RegularizedMultilayerMapEquation.h"
#endif
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
#include "LossyMapEquation.h"
#endif
#include "InfomapOptimizer.h"
#include "../io/InfomapError.h"
#include "../io/RunReport.h"
#include "../io/RunMetadata.h"
#include "../io/SafeFile.h"
#include "../io/OutputPlan.h"
#include "../io/TrialResults.h"
#include "../utils/FileURI.h"
#include "../utils/FlowCalculator.h"
#include "../utils/MemoryUsage.h"
#include "../utils/Console.h"
#include "../utils/ThreadConfig.h"
#include "../utils/TimingRegistry.h"
#include "../utils/convert.h"
#include "../utils/format.h"
#include "../utils/formatters.h"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <exception>
#include <mutex>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace infomap {

namespace {

  // Depth-prior bias (bits) for a hierarchy-stage accept test: the log-odds of a
  // Bernoulli depth prior, added to the codelength comparison. > 0 below the target
  // depth (toward deeper), < 0 above. Disabled for target 0 or strength <= 0 (a
  // negative would invert it; the public Config field is binding-settable).
  // Derivation: issue #308.
  double prefLevelsBias(unsigned int target, double strength, unsigned int currentDepth)
  {
    if (target == 0 || strength <= 0.0)
      return 0.0;
    constexpr double ln2 = 0.6931471805599453;
    const double delta = static_cast<double>(static_cast<int>(target) - static_cast<int>(currentDepth));
    return strength * (delta + 0.5) / ln2;
  }

  std::string compressionSummary(const std::vector<std::string>& compression)
  {
    if (compression.empty())
      return "0%";

    return io::stringify(compression, ", ");
  }

  void printPrettyStart(const Config& config, const Date& startDate, unsigned int numInitialModuleIds)
  {
    Console console;
    Log().print("{}Infomap v{}{}\n", console.bold(), INFOMAP_VERSION, console.reset());
    Log().print("{}{} started {}{}\n", console.dim(), console.branch(), startDate, console.reset());
    Log().print("{} Input      {}\n", console.bullet(), config.networkFile);
    Log().print("{} Output     {}\n", console.bullet(), config.noFileOutput ? "disabled" : config.outDirectory);
    if (!config.parsedOptions.empty()) {
      Log().print("{} Options    ", console.bullet());
      for (unsigned int i = 0; i < config.parsedOptions.size(); ++i) {
        if (i > 0)
          Log() << ", ";
        Log() << config.parsedOptions[i];
      }
      Log() << "\n";
    }
    if (numInitialModuleIds > 0) {
      Log().print("{} Initial    {} module ids\n", console.bullet(), numInitialModuleIds);
    }
    Log() << "\n";
  }

  void printPrettyEnd(const Date& endDate, const Stopwatch& elapsedTime)
  {
    Console console;
    Log().print("\n{}Finished {} in {}{}\n", console.dim(), endDate, elapsedTime, console.reset());
  }

  void printPrettyTrialTable(const std::vector<double>& codelengths, const std::vector<unsigned int>& numTopModules)
  {
    Console console;
    // Mark the single global best (the first trial reaching the minimum
    // codelength); other trials sharing that minimum are "tie". The previous
    // running comparison flagged every new minimum — including trial 1 — as
    // "best".
    const double bestCodelength = codelengths.empty() ? 0.0 : *std::min_element(codelengths.begin(), codelengths.end());
    bool markedBest = false;
    Log().print("  {}Trial  Codelength     Top modules  Best{}\n", console.dim(), console.reset());
    for (unsigned int i = 0; i < codelengths.size(); ++i) {
      const bool isMin = std::abs(codelengths[i] - bestCodelength) < 1e-10;
      const bool isBest = isMin && !markedBest;
      const bool isTie = isMin && markedBest;
      if (isBest)
        markedBest = true;
      Log().print("  {:>5}  {:>12}  {:>11}  {}{}{}\n",
                  i + 1,
                  io::toPrecision(codelengths[i]),
                  numTopModules[i],
                  isBest ? console.green() : isTie ? console.yellow()
                                                   : "",
                  isBest ? "best" : isTie ? "tie"
                                          : "",
                  console.reset());
    }
  }

#ifdef _OPENMP
  class ScopedOpenMpMaxActiveLevels {
  public:
    explicit ScopedOpenMpMaxActiveLevels(int value)
        : m_previous(omp_get_max_active_levels())
    {
      omp_set_max_active_levels(value);
    }

    ~ScopedOpenMpMaxActiveLevels()
    {
      omp_set_max_active_levels(m_previous);
    }

  private:
    int m_previous;
  };

  // The run's thread budget goes into the OpenMP nthreads-var ICV so every
  // parallel region inherits it, but that ICV belongs to the process, not to us.
  // An embedding host that limits OpenMP programmatically -- threadpoolctl, a BLAS
  // wrapper, torch -- had its limit replaced by one run and never given back, for
  // the rest of the process. Restore it the same way the neighbouring guard above
  // restores max-active-levels.
  class ScopedOpenMpNumThreads {
  public:
    explicit ScopedOpenMpNumThreads(int value)
        : m_previous(omp_get_max_threads())
    {
      omp_set_num_threads(value);
    }

    ~ScopedOpenMpNumThreads()
    {
      omp_set_num_threads(m_previous);
    }

  private:
    int m_previous;
  };

  // Clamp before the signed cast: the budget is unsigned, and a value above
  // INT_MAX would make omp_set_num_threads see a negative thread count.
  int asOpenMpThreadCount(unsigned int threads)
  {
    return static_cast<int>(std::min(threads, static_cast<unsigned int>(std::numeric_limits<int>::max())));
  }
#endif

  // The serial trial loop moves the config seed to the current trial's seed (see
  // RunSession::seedTrial), but everything downstream of the loop -- the run
  // manifest, the config fingerprint, the shard file's base_seed -- must report
  // the seed the user asked for. Give it back on the way out, including when a
  // trial throws or the run is interrupted.
  class ScopedConfigSeed {
  public:
    explicit ScopedConfigSeed(Config& config)
        : m_config(config), m_previous(config.seedToRandomNumberGenerator) {}

    ~ScopedConfigSeed()
    {
      m_config.seedToRandomNumberGenerator = m_previous;
    }

    ScopedConfigSeed(const ScopedConfigSeed&) = delete;
    ScopedConfigSeed& operator=(const ScopedConfigSeed&) = delete;

  private:
    Config& m_config;
    unsigned long m_previous;
  };

} // namespace

class InfomapBase::RunSession {
public:
  struct Result {
    std::ostringstream bestSolutionStatistics;
    unsigned int bestNumLevels = 0;
    double bestHierarchicalCodelength = std::numeric_limits<double>::max();
    NodePaths bestTree;
    unsigned int bestTrialIndex = 0;
    unsigned int threadsUsed = 1;
    bool bestTreeNeedsRestore = false;
  };

  RunSession(InfomapBase& infomap, Network& network, TimingRegistry& timing) : m_infomap(infomap), m_network(network), m_timing(timing) {}

  Result run()
  {
    {
      // Resolve the thread budget once for the whole run. RunSession is only
      // constructed for the main Infomap (run(Network&) requires isMainInfomap),
      // so this fires exactly once, before any OpenMP region.
      ThreadSources threadSources = readThreadSourcesFromEnv();
      threadSources.explicitThreads = m_infomap.numThreads; // 0 = auto
      m_threadBudget = resolveThreadBudget(threadSources);
      m_cpusetCount = threadSources.cpusetCount;
    }
#ifdef _OPENMP
    // Publishing the budget as the process-global max lets every parallel region
    // (recursive partition, parallel trials, inner parallelization) inherit it
    // through its existing omp_get_max_threads() call. Function scope on purpose:
    // the budget has to hold for the whole run, and the guard has to give the ICV
    // back on the way out -- including when a run throws.
    const ScopedOpenMpNumThreads scopedThreadBudget(asOpenMpThreadCount(m_threadBudget.threads));
#endif
    preflightOutputTargets(m_infomap);
    validateNetwork();
    {
      auto timer = m_timing.scope("pre_run_output_s");
      writeOutputArtifacts(m_infomap, m_network, OutputPhase::BeforeFlow);
    }
    {
      auto timer = m_timing.scope("configure_network_s");
      configureNetworkMode();
    }
    calculateFlowAndInitNetwork();
    {
      auto timer = m_timing.scope("post_flow_output_s");
      writeOutputArtifacts(m_infomap, m_network, OutputPhase::AfterFlow);
    }
    m_reportNetwork.nodes = m_network.numNodes();
    m_reportNetwork.links = m_network.numLinks();
    m_reportNetwork.directed = !m_infomap.isUndirectedFlow();
    m_runParallelTrials = selectParallelTrialMode();
    m_threadsUsed = m_runParallelTrials ? parallelTrialWorkers() : 1;
    releaseInputLinksIfCli();
    logRunPartitionStart();
    Result result;
    {
      auto timer = m_timing.scope("trial_optimize_s");
      result = runTrials();
    }
    result.threadsUsed = m_threadsUsed;
    return result;
  }

  void printSummary(const Result& result)
  {
    if (!m_infomap.isMainInfomap()) {
      return;
    }

    if (m_trialsRun > 1) {
      restoreBestResult(result);
    }

    auto summaryTimer = m_timing.scope("summary_s");

    Console console;
    console.section(fmt::format(FMT_STRING("Summary after {}{}"), m_trialsRun, m_trialsRun > 1 ? " trials" : " trial"));
    if (m_autoStopped) {
      const unsigned int patience = Config::convergePatience;
      Console::detail(0, "auto: stopped after {} trials (no improvement in last {})", m_trialsRun, patience);
    }
    std::string codelengthRange;
    std::string topModulesRange;
    if (m_trialsRun > 1) {
      double averageCodelength = 0.0;
      double minCodelength = m_infomap.m_codelengths[0];
      double maxCodelength = m_infomap.m_codelengths[0];
      double averageNumTopModules = 0.0;
      auto minNumTopModules = m_infomap.m_numTopModules[0];
      auto maxNumTopModules = m_infomap.m_numTopModules[0];
      printPrettyTrialTable(m_infomap.m_codelengths, m_infomap.m_numTopModules);
      for (unsigned int i = 0; i < m_trialsRun; ++i) {
        averageCodelength += m_infomap.m_codelengths[i];
        minCodelength = std::min(minCodelength, m_infomap.m_codelengths[i]);
        maxCodelength = std::max(maxCodelength, m_infomap.m_codelengths[i]);
        averageNumTopModules += m_infomap.m_numTopModules[i];
        minNumTopModules = std::min(minNumTopModules, m_infomap.m_numTopModules[i]);
        maxNumTopModules = std::max(maxNumTopModules, m_infomap.m_numTopModules[i]);
      }
      averageCodelength /= m_trialsRun;
      averageNumTopModules /= m_trialsRun;
      codelengthRange = fmt::format(FMT_STRING("{} / {} / {}"), io::toPrecision(minCodelength), io::toPrecision(averageCodelength), io::toPrecision(maxCodelength));
      topModulesRange = fmt::format(FMT_STRING("{} / {} / {}"), minNumTopModules, io::toPrecision(averageNumTopModules, 1, true), maxNumTopModules);
      console.metric("Codelength min/avg/max", codelengthRange);
      console.metric("Top modules min/avg/max", topModulesRange);
    }
    console.metric("Nodes", fmt::to_string(m_infomap.numLeafNodes()));
    console.metric("Links", fmt::to_string(m_network.numLinks()));
    console.metric("Average degree", io::toPrecision(m_network.numLinks() * 2.0 / m_infomap.numLeafNodes(), 1, true));
    console.metric("Top modules", fmt::format(FMT_STRING("{} ({} non-trivial)"), m_infomap.numTopModules(), m_infomap.numNonTrivialTopModules()));
    console.metric("Levels", fmt::to_string(result.bestNumLevels));
    console.metric("One-level codelength", io::toPrecision(m_infomap.getReferenceOneLevelCodelength()));
    console.metric("Best codelength", io::toPrecision(result.bestHierarchicalCodelength));
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
    if (m_infomap.lossy) {
      console.metric("Objective J", io::toPrecision(result.bestHierarchicalCodelength));
      console.metric("Lossy rate", io::toPrecision(m_infomap.getLossyRate()));
      console.metric("Lossy distortion", io::toPrecision(m_infomap.getLossyDistortion()));
      console.metric("Lambda", io::toPrecision(m_infomap.lossyLambda));
      const auto noise = m_infomap.noiseTopModules();
      std::string noiseSummary = std::to_string(noise.size()) + " of " + std::to_string(m_infomap.numTopModules());
      if (!noise.empty()) {
        noiseSummary += " (index:";
        for (auto i : noise)
          noiseSummary += " " + std::to_string(i);
        noiseSummary += ")";
      }
      console.metric("Noise modules", noiseSummary);
    }
#endif
    console.metric("Relative savings", fmt::format(FMT_STRING("{}%"), io::toPrecision(m_infomap.getRelativeCodelengthSavings() * 100, 2, true)));
    Log() << "\n";
    Log().print("{}\n", result.bestSolutionStatistics.str());
  }

private:
  Result runTrials()
  {
    Result result;
    result.bestTree = NodePaths(m_infomap.numLeafNodes());
    m_infomap.m_codelengths.clear();
    m_infomap.m_numTopModules.clear();
    m_timing.resetTrials(m_numTrials);

    if (m_runParallelTrials) {
      return runTrialsInParallel();
    }

    // The serial paths below reseed the config per trial; give the run's base
    // seed back once the loop is done. Constructed before the first seedTrial,
    // so what it captures is the base seed.
    const ScopedConfigSeed scopedSeed(m_infomap);

    if (m_convergeTrials) {
      runTrialsUntilConverged(result);
    } else {
      for (unsigned int i = 0; i < m_numTrials; ++i) {
        m_infomap.pollInterrupt();
        runTrial(i, result);
      }
    }

    return result;
  }

  // --converge: m_numTrials is a cap. Run trials sequentially and stop once the
  // best codelength has plateaued — no improvement beyond a relative tolerance
  // for `convergePatience` consecutive trials, after a `convergeMinTrials` floor.
  // The decision depends only on the deterministic per-trial codelengths, so the
  // trial count and result are reproducible for a given input + seed.
  void runTrialsUntilConverged(Result& result)
  {
    unsigned int trialsSinceImprovement = 0;
    unsigned int i = 0;
    for (; i < m_numTrials; ++i) {
      m_infomap.pollInterrupt();
      const double bestBefore = result.bestHierarchicalCodelength;
      runTrial(i, result);

      const bool improved = result.bestHierarchicalCodelength < bestBefore * (1.0 - Config::convergeTolerance);
      trialsSinceImprovement = improved ? 0 : trialsSinceImprovement + 1;

      if (i + 1 >= Config::convergeMinTrials && trialsSinceImprovement >= Config::convergePatience) {
        m_autoStopped = true;
        ++i;
        break;
      }
    }
    m_trialsRun = i;
  }

  Result runTrialsInParallel()
  {
    Result result;
    result.bestTree = NodePaths(m_infomap.numLeafNodes());
    result.bestTreeNeedsRestore = true;
    const unsigned int numWorkers = parallelTrialWorkers();

    std::vector<double> codelengths(m_numTrials, std::numeric_limits<double>::max());
    std::vector<unsigned int> numTopModules(m_numTrials, 0);
    bool hasError = false;
    unsigned int firstErrorTrial = 0;
    std::string firstError;
    std::mutex bestResultMutex;
    std::mutex errorMutex;
    std::mutex outputMutex;

#ifdef _OPENMP
    ScopedOpenMpMaxActiveLevels scopedMaxActiveLevels(1);
#pragma omp parallel for schedule(dynamic) num_threads(numWorkers)
#endif
    for (int workerIndex = 0; workerIndex < static_cast<int>(numWorkers); ++workerIndex) {
      auto workerConfig = m_infomap.getConfig();
      workerConfig.numTrials = 1;
      workerConfig.parallelTrials = false;
      workerConfig.innerParallelization = false;
      workerConfig.seedToRandomNumberGenerator = m_baseSeed + static_cast<unsigned int>(workerIndex);

      InfomapBase worker(workerConfig);
      worker.m_initialPartition = m_infomap.m_initialPartition;
      // Share the cancel flag so workers observe a cancel; they only read it and
      // unwind, caught by the per-trial handler below (issue #412).
      worker.inheritRuntimeContext(m_infomap);

      for (unsigned int trialIndex = static_cast<unsigned int>(workerIndex); trialIndex < m_numTrials; trialIndex += numWorkers) {
        try {
          Log::ScopedMute muteWorkerLogs;
          const auto seed = trialSeed(trialIndex);
          worker.seedToRandomNumberGenerator = seed;
          worker.reseed(static_cast<unsigned int>(seed));
          int threadNumber = 0;
#ifdef _OPENMP
          threadNumber = omp_get_thread_num();
#endif
          Stopwatch trialTimer(true);

          NodePaths trialTree;
          unsigned int trialNumLevels = 0;
          std::ostringstream trialStatistics;
          worker.initNetwork(m_network);
          initTrialPartition(worker);
          executeTrial(worker);
          worker.root().sortChildrenOnFlow();

          trialNumLevels = printPerLevelCodelength(worker.root(), trialStatistics);
          trialTree.reserve(worker.numLeafNodes());
          for (auto it(worker.iterLeafNodes()); !it.isEnd(); ++it) {
            trialTree.emplace_back(it->stateId, it.path());
          }

          const auto trialCodelength = worker.m_hierarchicalCodelength;
          const auto trialTopModules = worker.numTopModules();
          const auto trialTime = trialTimer.getElapsedTimeInSec();
          codelengths[trialIndex] = trialCodelength;
          numTopModules[trialIndex] = trialTopModules;
          m_timing.recordTrial(trialIndex, threadNumber, seed, trialTime, trialCodelength, trialTopModules, trialNumLevels);

          if (worker.printAllTrials && m_numTrials > 1) {
            std::lock_guard<std::mutex> lock(outputMutex);
            auto outputTimer = m_timing.scope("output_s");
            // Record this worker's one trial before it writes, so the header's
            // "trials N of M" describes the file. The worker never pushes into
            // m_codelengths itself -- that happens on the serial path -- and an
            // artifact that reported 0 trials would be worse than none (#906).
            if (worker.m_codelengths.empty())
              worker.m_codelengths.push_back(trialCodelength);
            worker.writeResult(static_cast<int>(m_infomap.trialOffset + trialIndex + 1));
          }

          std::lock_guard<std::mutex> lock(bestResultMutex);
          const auto bestIndexMissing = result.bestTrialIndex >= m_numTrials;
          const auto isBetter = trialCodelength < result.bestHierarchicalCodelength - 1e-10;
          const auto isEarlierTie = std::abs(trialCodelength - result.bestHierarchicalCodelength) < 1e-10 && trialIndex < result.bestTrialIndex;
          if (bestIndexMissing || isBetter || isEarlierTie) {
            result.bestSolutionStatistics.clear();
            result.bestSolutionStatistics.str(trialStatistics.str());
            result.bestNumLevels = trialNumLevels;
            result.bestHierarchicalCodelength = trialCodelength;
            result.bestTrialIndex = trialIndex;
            result.bestTree = std::move(trialTree);
          }
        } catch (const std::exception& e) {
          std::lock_guard<std::mutex> lock(errorMutex);
          if (!hasError || trialIndex < firstErrorTrial) {
            hasError = true;
            firstErrorTrial = trialIndex;
            firstError = e.what();
          }
        } catch (...) {
          std::lock_guard<std::mutex> lock(errorMutex);
          if (!hasError || trialIndex < firstErrorTrial) {
            hasError = true;
            firstErrorTrial = trialIndex;
            firstError = "unknown error";
          }
        }
      }
    }

    // Surface a cancel as InterruptionError, not the generic "Parallel trial
    // failed" the std::exception handler above would otherwise report (#412).
    if (m_infomap.interruptRequested()) {
      throw InterruptionError();
    }

    if (hasError) {
      throw std::runtime_error(fmt::format(FMT_STRING("Parallel trial {}/{} failed: {}"), firstErrorTrial + 1, m_numTrials, firstError));
    }

    m_infomap.m_codelengths = std::move(codelengths);
    m_infomap.m_numTopModules = std::move(numTopModules);
    return result;
  }

  unsigned int parallelTrialWorkers() const
  {
#ifdef _OPENMP
    // Worker count is driven directly by the OpenMP thread count (OMP_NUM_THREADS),
    // giving HPC users explicit control with no hidden heuristic. Clamp to numTrials
    // so we never spawn idle workers. Peak memory scales with the worker count (each
    // worker holds its own leaf network) — reduce OMP_NUM_THREADS if memory-constrained.
    const unsigned int maxOpenMpThreads = static_cast<unsigned int>(std::max(1, omp_get_max_threads()));
    return std::max(1u, std::min(m_numTrials, maxOpenMpThreads));
#else
    return 1;
#endif
  }

  void restoreBestResult(const Result& result)
  {
    // Compare against the actual executed trial count (m_trialsRun), not the cap
    // (m_numTrials), so --converge does not force a redundant restore + rewrite
    // when the best trial was the last one executed.
    if (m_trialsRun > 1 && (result.bestTreeNeedsRestore || result.bestTrialIndex < m_trialsRun - 1)) {
      // This restore + rewrite only refreshes the output file's elapsed-time
      // header; the user already saw the "Output" line when the best trial ran.
      // Mute it so the redundant "Initial …" (from initTree) and second "Output"
      // line don't leak into the console. The file still writes (via SafeFile,
      // not Log).
      Log::ScopedMute mute;
      // Restore Infomap tree to best solution.
      {
        auto timer = m_timing.scope("best_restore_s");
        m_infomap.initTree(result.bestTree);
      }
      if (!m_infomap.noFinalOutput) {
        auto outputTimer = m_timing.scope("output_s");
        m_infomap.writeResult(); // Overwrite result to get total elapsed time in output file header.
      }
    }
  }

  bool selectParallelTrialMode()
  {
    if (!m_infomap.parallelTrials) {
      return false;
    }

    if (m_numTrials < 2) {
      Console::warn(0, "--parallel-trials requires --num-trials > 1; running trials serially.");
      return false;
    }

#ifndef _OPENMP
    Console::warn(0, "--parallel-trials requires an OpenMP build; running trials serially.");
    return false;
#else
    if (m_infomap.innerParallelization) {
      Console::warn(0, "--parallel-trials ignores --inner-parallelization inside trial workers.");
    }
    const unsigned int workers = parallelTrialWorkers();
    Console console;
    Log() << "\n"
          << console.dim() << "  Parallel trials: " << workers << " workers from "
          << omp_get_max_threads() << " OpenMP threads (memory scales with workers; inner parallelization off)"
          << console.reset() << "\n";
    return true;
#endif
  }

  void validateNetwork()
  {
    if (m_network.numNodes() == 0) {
      m_network.postProcessInputData();
      if (m_network.numNodes() == 0) {
        throw std::domain_error("Network is empty");
      }
    }
  }

  void configureNetworkMode()
  {
    if (m_network.haveMemoryInput()) {
      Console::detail(1, "found higher-order network input, using the Map Equation for higher-order flows");
      m_infomap.setStateInput();
      m_infomap.setStateOutput();

      if (m_network.isMultilayerNetwork() && !m_infomap.isMultilayerNetwork()) {
        m_infomap.setMultilayerInput();
      }
    } else {
      if (m_infomap.haveMemory() || m_network.higherOrderInputMethodCalled()) {
        Console::warn(0, "Higher-order network specified but no higher-order input found.");
        // Use state output anyway for consistency even in the special case when input is first order
        m_infomap.setStateOutput();
      }
      Console::detail(1, "ordinary network input, using the Map Equation for first-order flows");
    }

    if (m_network.haveDirectedInput() && m_infomap.isUndirectedFlow()) {
      Console::note(1, "Directed input found, changing flow model from '{}' to '{}'.", m_infomap.flowModel, FlowModel(FlowModel::directed));
      m_infomap.flowModel = FlowModel::directed;
    }
    m_network.setConfig(m_infomap);
  }

  void calculateFlowAndInitNetwork()
  {
    // Build the consumed CSR link storage before flow: FlowCalculator reads CSR
    // weights and writes computed flows directly into the CSR arrays, which
    // initNetwork and the output writers then consume.
    m_network.finalizeLinks();

    {
      auto timer = m_timing.scope("flow_calculation_s");
      calculateFlow(m_network, m_infomap);
    }

    if (m_network.isBipartite()) {
      m_infomap.bipartite = true;
    }

    {
      auto timer = m_timing.scope("init_network_s");
      m_infomap.initNetwork(m_network);
    }

    if (m_infomap.numLeafNodes() == 0)
      throw std::domain_error("No nodes to partition");

    checkFlowPostCondition();
    warnOnObjectiveOptionsNotApplied();
  }

  // Options that are accepted, echoed in the run banner and recorded in the run metadata as
  // active, but that the objective this input selects does not act on. Silence made them read as
  // applied: the reported codelength and partition were simply those of a run without them (#904).
  //
  // Warnings rather than errors, and the objective dispatch is left as it is -- combining these
  // objectives is a separate design question. The run path is the right home for this, for the
  // same reason checkFlowPostCondition lives here: initNetwork is public C++ API that embedders
  // call directly, and it runs again per parallel-trial worker.
  void warnOnObjectiveOptionsNotApplied()
  {
    if (!m_infomap.m_optimizer->implementsObjectiveParameters()) {
      // Worded from the objective, not from the input: the condition is which objective
      // initOptimizer selected, and state, multilayer and meta-data input are the common way to
      // select another one but not the only way -- --lossy does it on an ordinary network too, and
      // its own validation does not reject these options. The banner above the warning names the
      // objective in use, so the message does not have to.
      if (m_infomap.entropyBiasCorrection)
        Console::warn(0, "--entropy-corrected has no effect on this run: the entropy bias correction is implemented by the ordinary map equation only, and this run uses a different objective. The run optimizes and reports the uncorrected codelength.");
      if (m_infomap.preferredNumberOfModules != 0)
        Console::warn(0, "--preferred-number-of-modules has no effect on this run: the module-count preference is implemented by the ordinary map equation only, and this run uses a different objective.");
    }

    if (m_infomap.haveMetaData() && (m_infomap.haveMemory() || m_infomap.isMultilayerNetwork())) {
      // initOptimizer tests haveMetaData() before haveMemory(), so meta-data wins and the run is
      // scored by a first-order objective with no physical-node codebook -- the aggregation that
      // makes this a memory network does not enter the objective at all.
      Console::warn(0, "--meta-data takes precedence over higher-order input: the run optimizes the meta-data objective over state nodes, without the physical-node codebook of the higher-order map equation. The higher-order structure affects the network, not the objective.");
    }
  }

  // The map equation is defined on a visit-rate distribution, so a non-finite
  // or empty flow distribution leaves every codelength downstream meaningless.
  // Such a run used to complete successfully and report "codelength 0" -- a
  // value lower than any true codelength, so a script comparing flow models
  // would select the broken run as the best model. Reject it instead.
  //
  // Only the run path checks this: initNetwork() is public C++ API that callers
  // legitimately invoke before any flow has been calculated, so the condition
  // does not belong inside it.
  //
  // A total that merely differs from 1.0 stays a warning. That is what
  // --skip-adjust-bipartite-flow asks for by design ("Keep flow on bipartite
  // nodes instead of distributing it to primary nodes"), and other flow models
  // may have their own reasons; the run is still interpretable. It is warned at
  // the default verbosity so it is no longer invisible.
  void checkFlowPostCondition() const
  {
    const double sumNodeFlow = m_infomap.root().data.flow;

    // A link-free network has no flow to distribute and no structure to find;
    // the trivial partition at zero bits is its correct answer, so an empty
    // total is only degenerate when there were links to carry flow.
    const bool zeroFlowIsExpected = m_network.numLinks() == 0;

    if (!std::isfinite(sumNodeFlow) || (sumNodeFlow <= 0 && !zeroFlowIsExpected))
      throw InfomapError(ExitCode::InputError,
                         fmt::format(FMT_STRING("Degenerate flow: total flow on nodes is {:g} for a network with {} links, so no "
                                                "codelength is defined. This happens when the flow model cannot distribute flow "
                                                "over the network -- for example a directed flow model with recorded teleportation "
                                                "on a bipartite network whose links all point from the primary side to the feature side."),
                                     sumNodeFlow,
                                     m_network.numLinks()));

    if (zeroFlowIsExpected)
      return;

    checkNoNegativeFlow();

    if (std::abs(sumNodeFlow - 1.0) > 1e-10)
      Console::warn(0, "Total flow on nodes is {:g}, not 1. The reported codelength is scaled accordingly.", sumNodeFlow);
  }

  // A visit rate and a boundary-crossing rate are probabilities, so a negative one leaves every
  // codelength downstream meaningless -- and it does not announce itself: plogp of a negative
  // flow is not a number, but where the negative terms cancel the map equation quietly charges
  // too little for a boundary that exists. The bipartite flow adjustment used to leave uncoded
  // feature nodes at flow 0 with negative enter/exit flow, and the two-level index codelength
  // came out at 1.1e-16 -- no cost at all for entering either of two modules (#957). An error
  // rather than a warning: there is no reading of the map equation under which this is a result.
  // Not an exact zero bound: the flow models accumulate, so a value a few epsilons below zero is
  // rounding, not a broken distribution.
  static constexpr double negativeFlowTolerance = -1e-10;

  static const char* negativeFlowQuantity(const FlowData& data) noexcept
  {
    if (data.flow < negativeFlowTolerance)
      return "flow";
    if (data.enterFlow < negativeFlowTolerance)
      return "enter flow";
    if (data.exitFlow < negativeFlowTolerance)
      return "exit flow";
    return nullptr;
  }

  void checkNoNegativeFlow() const
  {
    for (const auto* leafNode : m_infomap.leafNodes()) {
      const auto& data = leafNode->data;
      const char* quantity = negativeFlowQuantity(data);
      if (quantity == nullptr)
        continue;

      throw InfomapError(ExitCode::InternalError,
                         fmt::format(FMT_STRING("Negative {} on node {}: flow {:g}, enter flow {:g}, exit flow {:g}. "
                                                "These are probabilities, so no codelength is defined. This is a bug in the "
                                                "flow calculation for this combination of input and flow model, not something "
                                                "the input can be wrong about."),
                                     quantity,
                                     leafNode->stateId,
                                     data.flow,
                                     data.enterFlow,
                                     data.exitFlow));
    }
  }

  // The same invariant one level up. A module's enter/exit flow is not the sum of its leaves': the
  // optimizer maintains it incrementally and removes the teleportation internal to the module, and
  // that removal can over-subtract even when every leaf is sound. On a bipartite network under
  // --regularized this left modules at enter/exit -0.00121 and the resulting per-module codelength
  // negative, which showed up only as a small disagreement between two evaluations of the same
  // partition rather than as a failure (#958).
  //
  // Checked once per trial on the finished tree rather than after every consolidation: a walk over
  // the modules costs nothing next to the search, and a codelength is only reported from here on.
  //
  // Takes the tree to walk rather than reading m_infomap: under --parallel-trials the trial runs on
  // a worker instance, and checking the main instance's tree would check an unpartitioned one.
  //
  // Running after restoreHardPartition() rather than before is deliberate, and does not check a
  // different tree than the codelength was computed from. A collapsed hard-partition super-node has
  // its child chain swapped out, so firstChild is null and it is a leaf while collapsed -- skipped
  // here either way -- and the restore only splices the original leaves into its place. Measured on
  // a hard partition of ninetriangles.net: same four modules at the same depths with bit-identical
  // flow, enter and exit flow before and after; only childDegree() changes (1 super-node child to 9
  // real ones), and reporting the real count in the message is the more useful of the two.
  void checkNoNegativeModuleFlow(InfoNode& root) const
  {
    for (auto it = root.begin_post_depth_first(); !it.isEnd(); ++it) {
      const auto& node = *it;
      if (node.isLeaf())
        continue;

      const char* quantity = negativeFlowQuantity(node.data);
      if (quantity == nullptr)
        continue;

      throw InfomapError(ExitCode::InternalError,
                         fmt::format(FMT_STRING("Negative {} on a module at depth {} with {} children: flow {:g}, "
                                                "enter flow {:g}, exit flow {:g}. These are probabilities, so no codelength "
                                                "is defined. Every leaf node's flow is sound here, so this is a bug in how "
                                                "module flow is accumulated for this combination of input and flow model."),
                                     quantity,
                                     it.depth(),
                                     node.childDegree(),
                                     node.data.flow,
                                     node.data.enterFlow,
                                     node.data.exitFlow));
    }
  }

  void releaseInputLinksIfCli()
  {
    // If used as a library, we may want to reuse the network instance, else clear to use less memory
    // TODO: May have to use some meta data for output?
    if (m_infomap.isCLI && !m_runParallelTrials) {
      m_network.clearLinks();
    }
  }

  void logRunPartitionStart()
  {
    if (m_infomap.haveMemory())
      Console::detail(2, "run Infomap with memory");
    else
      Console::detail(2, "run Infomap");
#ifdef _OPENMP
    {
      // Prefer the cpuset size (the actual allocation under Linux/SLURM) over
      // omp_get_num_procs(), which ignores affinity; fall back when unknown.
      const unsigned int cpus = m_cpusetCount > 0 ? m_cpusetCount : static_cast<unsigned int>(std::max(1, omp_get_num_procs()));
      const char* source = threadSourceName(m_threadBudget.source);
      const char* innerState = m_infomap.innerParallelization ? "enabled" : "disabled";
      Console console;
      console.section("Threads");
      console.metric("CPUs available", fmt::to_string(cpus));
      console.metric("Thread budget", fmt::format(FMT_STRING("{} ({})"), m_threadBudget.threads, source));
      console.metric("Inner parallelization", innerState);

      // Oversubscription warnings (only meaningful when cpuset is known, i.e. Linux).
      // Warn whenever the *resolved* budget exceeds the cpuset, regardless of source
      // (explicit --num-threads, INFOMAP_NUM_THREADS, SLURM_CPUS_PER_TASK, OMP_NUM_THREADS).
      if (m_cpusetCount > 0 && m_threadBudget.threads > m_cpusetCount) {
        Console::warn(0, "Thread budget {} ({}) exceeds the available cpuset of {} CPUs; this may oversubscribe the node.", m_threadBudget.threads, source, m_cpusetCount);
      }
      if (m_cpusetCount > 0) {
        const unsigned int innerThreads = (m_infomap.innerParallelization && !m_runParallelTrials) ? m_threadBudget.threads : 1;
        const unsigned int demand = m_threadsUsed * innerThreads;
        if (demand > m_cpusetCount) {
          Console::warn(0, "{} trial workers × {} inner threads = {} exceeds the available cpuset of {} CPUs.", m_threadsUsed, innerThreads, demand, m_cpusetCount);
        }
      }
    }
#endif
  }

  // The seed of the trial at global index `trialOffset + trialIndex`. Every mode
  // derives it the same way, so "the trial with seed s" names one run whether it
  // was produced serially, by --parallel-trials, or by a shard at an offset.
  //
  // Kept at the config field's own width. Every route into the engine takes
  // unsigned int (`BasicRandom(unsigned int)` on construction,
  // `Random::seed(unsigned int)` from reseed and from setNonMainConfig), so the
  // RNG only ever sees the low 32 bits of the seed.
  //
  // Whether that discards anything depends on the platform:
  // Config::seedToRandomNumberGenerator is `unsigned long`, which is 64-bit on
  // LP64 (Linux, macOS) and 32-bit on LLP64 (Windows). On Windows the field is
  // exactly the engine's width and no seed can outrun it. On LP64 the narrowing
  // is real but uniform, which is what keeps the reported seed usable: a
  // standalone run given the value recorded here narrows to the same engine seed
  // and reproduces the trial, and the recursion reading the config field cannot
  // land on a different stream than the engine. The one visible consequence
  // there is aliasing -- seeds congruent mod 2^32 name the same run -- which
  // follows from the width mismatch, not from anything about per-trial seeding.
  unsigned long trialSeed(unsigned int trialIndex) const
  {
    return m_baseSeed + m_infomap.trialOffset + trialIndex;
  }

  // Reseed for every trial, in every mode.
  //
  // This used to happen only in "sharding mode" (trialOffset > 0 or a
  // --trial-results path given), which made the presence of an output flag
  // decide whether the trials were independently seeded or one continuous RNG
  // stream. Adding --trial-results therefore changed the published partition,
  // and the per-trial seeds reported in --timing-json named seeds that never
  // ran their trial. Both are fixed by seeding unconditionally.
  //
  // Set the config field as well as the engine: the recursive partition builds
  // its sub-Infomaps with setNonMainConfig, which seeds each one from
  // Config::seedToRandomNumberGenerator, so an engine-only reseed would leave
  // the whole recursion on the base seed and make shard trial i differ from the
  // same trial run any other way. runTrials restores the base seed afterwards.
  void seedTrial(unsigned int trialIndex)
  {
    const auto seed = trialSeed(trialIndex);
    m_infomap.seedToRandomNumberGenerator = seed;
    m_infomap.reseed(static_cast<unsigned int>(seed));
  }

  void runTrial(unsigned int trialIndex, Result& result)
  {
    seedTrial(trialIndex);
    m_infomap.removeModules();
    auto startDate = Date();
    Stopwatch timer(true);

    logTrialStart(trialIndex, startDate);
    initTrialPartition(m_infomap);
    executeTrial(m_infomap);
    finishTrial(trialIndex, timer, result);
  }

  void logTrialStart(unsigned int trialIndex, const Date& startDate)
  {
    if (!m_infomap.isMainInfomap()) {
      return;
    }

    Console console;
    Log() << "\n"
          << console.cyan() << "Trial " << (trialIndex + 1) << "/" << m_numTrials << console.reset()
          << console.dim() << " started " << startDate << console.reset() << "\n";
  }

  void initTrialPartition(InfomapBase& infomap)
  {
    if (!infomap.isMainInfomap()) {
      return;
    }

    if (!infomap.clusterDataFile.empty() && !m_network.initialPartitionPaths().empty())
      Console::note(0, "--cluster-data overrides embedded JSON initial-partition paths.");

    if (!infomap.clusterDataFile.empty())
      infomap.initPartition(infomap.clusterDataFile, infomap.clusterDataIsHard, &m_network);
    else if (!infomap.m_multilayerInitialPartition.empty()) {
      // Resolve the physical (layer_id, node_id) keys to the generated state
      // ids now that the network is built, then init as a normal partition.
      const auto& layerNodeToStateId = m_network.layerNodeToStateId();
      InfomapBase::InitialPartition statePartition;
      for (const auto& row : infomap.m_multilayerInitialPartition) {
        const auto layerId = row[0];
        const auto nodeId = row[1];
        const auto moduleId = row[2];
        auto layerIt = layerNodeToStateId.find(layerId);
        if (layerIt == layerNodeToStateId.end())
          throw std::out_of_range("Initial partition references a multilayer node not in the network: (layer " + std::to_string(layerId) + ", node " + std::to_string(nodeId) + ")");
        auto nodeIt = layerIt->second.find(nodeId);
        if (nodeIt == layerIt->second.end())
          throw std::out_of_range("Initial partition references a multilayer node not in the network: (layer " + std::to_string(layerId) + ", node " + std::to_string(nodeId) + ")");
        statePartition[nodeIt->second] = moduleId;
      }
      infomap.initPartition(statePartition, infomap.clusterDataIsHard);
    } else if (!infomap.m_initialPartition.empty())
      infomap.initPartition(infomap.m_initialPartition, infomap.clusterDataIsHard);
    else if (!m_network.initialPartitionPaths().empty()) {
      // Embedded JSON nodes[].path / states[].path, reusing the .tree channel.
      unsigned int numNodesNotInNetwork = 0;
      const auto normalizedTree = infomap.normalizeTreePaths(m_network.initialPartitionPaths(), numNodesNotInNetwork);
      if (numNodesNotInNetwork > 0)
        Console::note(1, "{} nodes in embedded initial-partition paths not found in network.", numNodesNotInNetwork);
      infomap.initTree(normalizedTree);
    }
  }

  void executeTrial(InfomapBase& infomap)
  {
    if (!infomap.noInfomap)
      infomap.runPartition();
    else
      infomap.m_hierarchicalCodelength = infomap.calcCodelengthOnTree(infomap.root(), true);

    if (infomap.haveHardPartition())
      infomap.restoreHardPartition();

    checkNoNegativeModuleFlow(infomap.root());
  }

  void finishTrial(unsigned int trialIndex, Stopwatch& timer, Result& result)
  {
    if (!m_infomap.isMainInfomap()) {
      return;
    }

    Console console;
    Log() << "\n"
          << console.green() << "Done" << console.reset()
          << " trial " << (trialIndex + 1) << "/" << m_numTrials
          << " in " << timer.getElapsedTimeInSec() << "s"
          << " | codelength " << io::toPrecision(m_infomap.m_hierarchicalCodelength) << "\n";
    m_infomap.m_codelengths.push_back(m_infomap.m_hierarchicalCodelength);
    m_infomap.m_numTopModules.push_back(m_infomap.numTopModules());
    // From trialSeed, not from the live config field: seedTrial has already moved
    // that field to this trial's seed, so recomputing from it would double-count
    // the global index.
    m_timing.recordTrial(trialIndex, 0, trialSeed(trialIndex), timer.getElapsedTimeInSec(), m_infomap.m_hierarchicalCodelength, m_infomap.numTopModules(), m_infomap.numLevels());

    if (m_infomap.printAllTrials && m_numTrials > 1) {
      auto outputTimer = m_timing.scope("output_s");
      m_infomap.writeResult(static_cast<int>(m_infomap.trialOffset + trialIndex + 1));
    }

    if (m_infomap.m_hierarchicalCodelength < result.bestHierarchicalCodelength - 1e-10) {
      updateBestResult(trialIndex, result);
    }
  }

  void updateBestResult(unsigned int trialIndex, Result& result)
  {
    result.bestSolutionStatistics.clear();
    result.bestSolutionStatistics.str("");
    result.bestNumLevels = printPerLevelCodelength(m_infomap.root(), result.bestSolutionStatistics);
    result.bestHierarchicalCodelength = m_infomap.m_hierarchicalCodelength;
    result.bestTrialIndex = trialIndex;
    m_infomap.root().sortChildrenOnFlow();
    auto outputTimer = m_timing.scope("output_s");
    if (!m_infomap.noFinalOutput) {
      m_infomap.writeResult();
    }
    if (m_numTrials > 1) {
      result.bestTree.clear();
      for (auto it(m_infomap.iterLeafNodes()); !it.isEnd(); ++it) {
        result.bestTree.emplace_back(it->stateId, it.path());
      }
    }
  }

public:
  void writeRunReports(const Result& result)
  {
    if (!m_infomap.summaryJsonPath.empty()) {
      RunSummaryReport report;
      report.version = INFOMAP_VERSION;
      report.codelength = result.bestHierarchicalCodelength;
      report.topModules = m_infomap.numTopModules();
      report.levels = result.bestNumLevels;
      report.trials = m_trialsRun;
      report.bestTrial = result.bestTrialIndex + 1;
      report.autoStopped = m_autoStopped;
      report.haveFlowConvergence = m_infomap.m_network.haveFlowConvergence();
      report.flowConverged = m_infomap.m_network.flowConverged();
      report.flowIterations = m_infomap.m_network.flowIterations();
      report.flowError = m_infomap.m_network.flowError();
      report.trialCodelengths = m_infomap.m_codelengths;
      report.trialTopModules = m_infomap.m_numTopModules;
      writeJsonReport(m_infomap.summaryJsonPath, runSummaryReportJson(report), m_infomap.overwriteOutput());
    }

    if (!m_infomap.timingJsonPath.empty()) {
      RunTimingReport report;
      report.version = INFOMAP_VERSION;
#ifdef _OPENMP
      report.openmp = true;
      report.threadsRequested = m_threadBudget.threads;
      report.threadSource = threadSourceName(m_threadBudget.source);
#else
      report.openmp = false;
      report.threadsRequested = 1;
      report.threadSource = "serial";
#endif
      report.threadsUsed = result.threadsUsed;
      report.network = m_reportNetwork;
      report.timing = m_timing.phases();
      report.trials = m_timing.trials();
      report.includeMemory = m_infomap.memoryReport;
      if (report.includeMemory) {
        report.memory = currentMemoryReport(report.network.nodes, report.network.links);
      }
      writeJsonReport(m_infomap.timingJsonPath, runTimingReportJson(report), m_infomap.overwriteOutput());
    }

    if (!m_infomap.runManifestPath.empty()) {
      writeJsonReport(m_infomap.runManifestPath, runManifestJson(m_infomap), m_infomap.overwriteOutput());
    }

    if (!m_infomap.trialResultsPath.empty()) {
      TrialResultsFile trialResultsFile;
      trialResultsFile.networkFingerprint = networkFingerprint(m_infomap.networkFile);
      trialResultsFile.configFingerprint = configFingerprint(m_infomap.getConfig());
      trialResultsFile.infomapVersion = INFOMAP_VERSION;
      trialResultsFile.baseSeed = m_baseSeed;
      trialResultsFile.trialOffset = m_infomap.trialOffset;
      trialResultsFile.numTrials = m_numTrials;

      // Build per-trial entries from timing records.
      const auto timingTrials = m_timing.trials();
      for (std::size_t s = 0; s < timingTrials.size(); ++s) {
        const auto& rec = timingTrials[s];
        if (!rec.valid) continue;
        TrialResultEntry entry;
        entry.trial = m_infomap.trialOffset + static_cast<unsigned int>(s);
        entry.seed = rec.seed;
        entry.codelength = rec.codelength;
        entry.numTopModules = rec.topModules;
        entry.numLevels = rec.numLevels;
        entry.thread = rec.thread;
        entry.timeSec = rec.timeSec;
        trialResultsFile.trials.push_back(entry);
      }

      // Determine best tree file.
      // After printSummary / restoreBestResult, m_infomap holds the best tree in memory.
      // We write the best tree with an explicit per-trial-numbered filename (always, not
      // conditional on printAllTrials) so the merge step has a stable reference.
      // If printAllTrials already wrote it and --no-overwrite is set, skip the write to
      // avoid a collision on the already-existing file.
      const unsigned int globalBestIndex = m_infomap.trialOffset + result.bestTrialIndex;
      if (m_infomap.printTree && !m_infomap.noFileOutput) {
        // Always use the trial-numbered basename for the shard best tree, regardless of
        // printAllTrials, so the merge step always finds a per-shard-named file.
        const auto treeBasename = m_infomap.outDirectory + m_infomap.outName
            + "_trial_" + std::to_string(globalBestIndex + 1);
        const auto treeAbsPath = outputFilenameForResultKey(treeBasename, "tree");

        const bool alreadyWritten = pathExists(treeAbsPath);
        if (!alreadyWritten || m_infomap.overwriteOutput()) {
          auto outputTimer = m_timing.scope("output_s");
          m_infomap.writeTree(treeAbsPath, false);
        }

        // For higher-order input the physical tree is a projection: the same node id
        // appears under several modules, so it cannot express the partition that was
        // found. Write the state-level companion too, since this file is the only
        // artifact infomap.merge reads and --no-final-output (which the HPC recipe
        // prescribes) means nothing else writes it (#906).
        std::string statesTreeAbsPath;
        if (m_infomap.haveMemory()) {
          statesTreeAbsPath = outputFilenameForResultKey(treeBasename + "_states", "tree");
          if (!pathExists(statesTreeAbsPath) || m_infomap.overwriteOutput()) {
            auto outputTimer = m_timing.scope("output_s");
            m_infomap.writeTree(statesTreeAbsPath, true);
          }
        }

        // Store the tree path relative to the results file's directory, so the
        // merge tool can resolve it from wherever the results file lives.
        std::string resultsDir;
        const auto lastSlash = m_infomap.trialResultsPath.find_last_of("/\\");
        if (lastSlash != std::string::npos) {
          resultsDir = m_infomap.trialResultsPath.substr(0, lastSlash + 1);
        }
        const auto relativeToResults = [&resultsDir](const std::string& path) {
          if (!resultsDir.empty() && path.substr(0, resultsDir.size()) == resultsDir) {
            // File lives under the results directory — store the relative remainder.
            return path.substr(resultsDir.size());
          }
          // Results file is in the CWD (resultsDir empty) or an unrelated directory: the
          // path as written (relative to the CWD) is correct. Stripping to a basename
          // here would drop the output-directory prefix.
          return path;
        };
        trialResultsFile.bestTreeFile = relativeToResults(treeAbsPath);
        if (!statesTreeAbsPath.empty())
          trialResultsFile.bestStatesTreeFile = relativeToResults(statesTreeAbsPath);
      }

      writeJsonReport(m_infomap.trialResultsPath, serializeTrialResults(trialResultsFile), m_infomap.overwriteOutput());
    }
  }

private:
  InfomapBase& m_infomap;
  Network& m_network;
  TimingRegistry& m_timing;
  RunReportNetwork m_reportNetwork;
  const unsigned int m_numTrials = m_infomap.numTrials;
  // Captured before the first trial reseeds the config field, so the run keeps a
  // stable notion of "the seed the user asked for".
  const unsigned long m_baseSeed = m_infomap.seedToRandomNumberGenerator;
  const bool m_convergeTrials = m_infomap.convergeTrials;
  unsigned int m_trialsRun = m_infomap.numTrials; // actual count; equals m_numTrials unless --converge stops early
  bool m_autoStopped = false;
  bool m_runParallelTrials = false;
  unsigned int m_threadsUsed = 1;
  ThreadBudget m_threadBudget;
  unsigned int m_cpusetCount = 0;
};

std::map<unsigned int, std::vector<unsigned int>> InfomapBase::getMultilevelModules(bool states)
{
  if (haveMemory() && !states) {
    throw std::runtime_error("Cannot get multilevel modules on higher-order network without states.");
  }
  unsigned int maxDepth = maxTreeDepth();
  unsigned int numModuleLevels = maxDepth - 1;
  std::map<unsigned int, std::vector<unsigned int>> modules;
  if (maxDepth < 2) return modules;
  for (unsigned int level = 1; level <= numModuleLevels; ++level) {
    for (auto it(iterLeafNodes(static_cast<int>(level))); !it.isEnd(); ++it) {
      auto& node = *it;
      auto nodeId = states ? node.stateId : node.physicalId;
      modules[nodeId].push_back(it.moduleId());
    }
  }
  return modules;
}

// ===================================================
// IO
// ===================================================

std::ostream& operator<<(std::ostream& out, const InfomapBase& infomap)
{
  return infomap.toString(out);
}

// ===================================================
// Getters
// ===================================================

unsigned int InfomapBase::numLevels() const
{
  // TODO: Make sure this is not called unless tree is guaranteed to have even depth!
  unsigned int depth = 0;
  InfoNode* n = m_root.firstChild;
  while (n != nullptr) {
    ++depth;
    n = n->firstChild;
  }
  return depth;
}

unsigned int InfomapBase::maxTreeDepth() const
{
  unsigned int maxDepth = 0;
  for (auto it = root().begin_tree(); !it.isEnd(); ++it) {
    const InfoNode& node = *it;
    if (node.isLeaf()) {
      if (it.depth() > maxDepth) {
        maxDepth = it.depth();
      }
    }
  }
  return maxDepth;
}

// ===================================================
// Run
// ===================================================

void InfomapBase::run(const std::string& parameters)
{
  if (!isMainInfomap()) {
    runPartition();
    return;
  }

  m_elapsedTime = Stopwatch(true);
  m_startDate = Date();

  std::string currentParameters = fmt::format(FMT_STRING("{}{}{}"), m_initialParameters, parameters.empty() ? "" : " ", parameters);
  if (currentParameters != m_currentParameters) {
    m_currentParameters = currentParameters;
    setConfig(Config(m_currentParameters, isCLI));
    m_network.setConfig(*this);
  }

  Log::init(verbosity, silent, verboseNumberPrecision);

  printPrettyStart(*this, m_startDate, m_initialPartition.size());
#ifdef _OPENMP
#pragma omp parallel
#pragma omp master
  {
    Console console;
    Log().print("{} Runtime    OpenMP {} with {} threads\n", console.bullet(), _OPENMP, omp_get_num_threads());
  }
#endif

  m_network.postProcessInputData();
  if (m_network.numNodes() == 0) {
    m_network.readInputData(networkFile);
  }

  if (!metaDataFile.empty()) {
    initMetaData(metaDataFile);
  }

  run(m_network);

  printPrettyEnd(m_endDate, m_elapsedTime);
}

void InfomapBase::run(Network& network)
{
  if (!isMainInfomap())
    throw std::logic_error("Can't run a non-main Infomap with an input network");

  // Stamp the owner thread (gates the callback) and clear any leftover cancel so
  // the instance is reusable after an interrupted run (issue #412).
  m_ownerThreadId = std::this_thread::get_id();
  m_cancelRequested.store(false, std::memory_order_relaxed);
  pollInterrupt();

  TimingRegistry timing;
  Stopwatch totalTimer(true);
  RunSession runSession(*this, network, timing);
  auto runResult = runSession.run();
  m_elapsedTime.stop();
  m_endDate = Date();
  runSession.printSummary(runResult);
  totalTimer.stop();
  timing.setPhase("total_s", totalTimer.getElapsedTimeInSec());
  runSession.writeRunReports(runResult);
}

// ===================================================
// Run: Init: *
// ===================================================

InfomapBase& InfomapBase::initMetaData(const std::string& metaDataFile)
{
  m_network.readMetaData(metaDataFile);
  return *this;
}

InfomapBase& InfomapBase::initNetwork(Network& network)
{
  if (network.numNodes() == 0)
    throw std::domain_error("No nodes in network");
  // A fresh network init invalidates any previous hard partition restore buffer.
  m_originalLeafNodes.clear();
  if (m_root.firstChild != nullptr || m_root.collapsedFirstChild != nullptr) {
    m_root.deleteChildren();
    m_leafNodes.clear();
  }
  generateSubNetwork(network);

  initOptimizer();

  // Per-instance network properties for entropy bias correction must be set after the
  // optimizer/objective exists and before init() reads them via calcCodelength().
  if (entropyBiasCorrection)
    m_optimizer->setNetworkProperties(network);

  init();
  return *this;
}

InfomapBase& InfomapBase::initNetwork(InfoNode& parent, bool asSuperNetwork)
{
  m_originalLeafNodes.clear();
  generateSubNetwork(parent);

  if (asSuperNetwork)
    transformNodeFlowToEnterFlow(m_root);

  init();
  return *this;
}

InfomapBase& InfomapBase::initPartition(const std::string& clusterDataFile, bool hard, const Network* network)
{
  ClusterMap clusterMap;
  if (network != nullptr && network->isMultilayerNetwork()) {
    const auto& map = network->layerNodeToStateId();
    clusterMap.readClusterData(clusterDataFile, false, &map);
  } else {
    clusterMap.readClusterData(clusterDataFile);
  }

  Console().status("Initial", fmt::format(FMT_STRING("reading {}"), clusterDataFile));
  Console::detail(1, "init partition from file '{}'", clusterDataFile);

  const auto& ext = clusterMap.extension();

  if (ext == "tree" || ext == "ftree") {
    unsigned int numTreeNodesNotInNetwork = 0;
    const auto normalizedTree = normalizeTreePaths(clusterMap.treePaths(), numTreeNodesNotInNetwork);
    if (numTreeNodesNotInNetwork > 0) {
      Console::note(1, "{} physical nodes in tree not found in network.", numTreeNodesNotInNetwork);
    }
    initTree(normalizedTree);
  } else if (ext == "clu") {
    initPartition(clusterMap.clusterIds(), hard);
  }

  Console().status("Initial", fmt::format(FMT_STRING("generated {} levels, codelength {}"), numLevels(), io::toPrecision(m_hierarchicalCodelength)));
  Console::detail(1, "generated {} levels, codelength {:g} + {:g} = {}", numLevels(), getIndexCodelength(), m_hierarchicalCodelength - getIndexCodelength(), io::toPrecision(m_hierarchicalCodelength));

  return *this;
}

NodePaths InfomapBase::normalizeTreePaths(const TreePaths& tree, unsigned int& numNodesNotInNetwork) const
{
  numNodesNotInNetwork = 0;
  NodePaths normalized;
  normalized.reserve(tree.size());

  // Fast path: every row is a state-level leaf — pass through.
  bool allState = std::all_of(tree.begin(), tree.end(), [](const auto& tp) { return tp.idType == TreeLeafIdType::state; });
  if (allState) {
    for (const auto& tp : tree)
      normalized.emplace_back(tp.nodeId, tp.path);
    return normalized;
  }

  // Build physical -> [state ids] index across the active leaf set.
  std::map<unsigned int, std::vector<unsigned int>> physicalToStateIds;
  for (const auto* leafNode : m_leafNodes) {
    physicalToStateIds[leafNode->physicalId].push_back(leafNode->stateId);
  }

  // Count how many tree rows reference each physical id; a physical id may
  // appear multiple times (one row per state) on higher-order outputs.
  std::map<unsigned int, unsigned int> physicalRowCount;
  for (const auto& tp : tree) {
    if (tp.idType == TreeLeafIdType::physical)
      ++physicalRowCount[tp.nodeId];
  }

  std::map<unsigned int, unsigned int> physicalRowsConsumed;
  std::map<unsigned int, unsigned int> physicalStatesConsumed;
  // Physical nodes whose states are spread over more than one module. A physical .tree
  // has one row per (module, physical node) and no state ids, so which state went to
  // which module simply is not in the file -- the pairing below is ascending state id
  // against file order, which is a guess (#908).
  unsigned int numSplitPhysicalNodes = 0;

  for (const auto& tp : tree) {
    if (tp.idType == TreeLeafIdType::state) {
      normalized.emplace_back(tp.nodeId, tp.path);
      continue;
    }

    const auto it = physicalToStateIds.find(tp.nodeId);
    if (it == physicalToStateIds.end()) {
      ++numNodesNotInNetwork;
      continue;
    }

    auto& rowsConsumed = physicalRowsConsumed[tp.nodeId];
    auto& statesConsumed = physicalStatesConsumed[tp.nodeId];
    const auto totalRows = physicalRowCount[tp.nodeId];
    const auto remainingRows = totalRows - rowsConsumed;
    const auto remainingStates = it->second.size() - statesConsumed;

    if (remainingStates == 0) {
      ++numNodesNotInNetwork;
      ++rowsConsumed;
      continue;
    }

    if (rowsConsumed == 0 && totalRows > 1 && it->second.size() > 1)
      ++numSplitPhysicalNodes;

    // If this is the last row for this physical id, claim every remaining state.
    // Otherwise consume one state per row, preserving 1:1 alignment when the
    // original output wrote one row per state.
    const auto statesToAssign = (remainingRows <= 1) ? remainingStates : 1u;
    for (unsigned int i = 0; i < statesToAssign; ++i) {
      normalized.emplace_back(it->second[statesConsumed + i], tp.path);
    }
    statesConsumed += statesToAssign;
    ++rowsConsumed;
  }

  if (numSplitPhysicalNodes > 0) {
    Console::warn(0,
                  "{} {} states split across modules in this tree. A physical "
                  "tree cannot express which state belongs to which module, so the states were paired "
                  "with the rows in ascending state-id order and the partition read back is likely not "
                  "the one that was written. Use the states tree (_states.tree) for an exact round trip.",
                  numSplitPhysicalNodes,
                  numSplitPhysicalNodes == 1 ? "physical node has its" : "physical nodes have their");
  }

  return normalized;
}

InfomapBase& InfomapBase::initTree(const NodePaths& tree)
{
  Log(4) << "Init tree... ";
  int maxDepth = 2;
  // If only two-level partition, we can directly use initPartition
  for (const auto& nodePath : tree) {
    if (nodePath.second.size() > 2) {
      maxDepth = std::max(maxDepth, static_cast<int>(nodePath.second.size()));
    }
  }
  // A deeper tree under --two-level takes the same route, keeping each leaf's top-level
  // module and dropping the sub-module structure below it. That is the collapse the
  // hierarchical path performs with removeSubModules, and doing it here means the
  // objective is initialised by the branch that works: building the full tree and then
  // running a two-level search over it left the search with the deeper tree's
  // index/module codelengths, which it reported as 0 -- or, with the tuning loop
  // enabled, aborted on "fineTune() called but numLevels != 2" (#898).
  if (maxDepth == 2 || twoLevel) {
    std::map<unsigned int, unsigned int> clusterIds;
    for (const auto& nodePath : tree) {
      const auto nodeId = nodePath.first;
      const auto clusterId = nodePath.second[0]; // First level module
      clusterIds[nodeId] = clusterId;
    }
    return initPartition(clusterIds, false);
  }

  // Refuse a module that would get both a leaf and a sub-module as children. The
  // search never produces that shape -- checked across 205k leaves at depths 2 to 8 --
  // and the code downstream assumes it cannot happen: the per-level aggregation used to
  // dereference the leaf's null firstChild, the Newick writer emits an empty taxon for
  // the leaf, and the objective gives the same partition two different codelengths
  // depending on the order the rows appear in the file. Supporting the shape means
  // deciding what the map equation is for a module that codes a leaf and a sub-module
  // side by side, which is a modelling question rather than a fix, so the input is
  // named and rejected instead of silently reinterpreted (#898). Note this runs after
  // the two-level branch above, where the sub-module structure is dropped anyway.
  std::map<Path, bool> prefixHasLeafChild;
  std::map<Path, bool> prefixHasModuleChild;
  for (const auto& nodePath : tree) {
    const auto& path = nodePath.second;
    if (path.size() < 2)
      continue;
    Path prefix;
    for (size_t i = 0; i + 1 < path.size(); ++i) {
      prefix.push_back(path[i]);
      if (i + 2 == path.size())
        prefixHasLeafChild[prefix] = true;
      else
        prefixHasModuleChild[prefix] = true;
    }
  }
  for (const auto& it : prefixHasLeafChild) {
    if (prefixHasModuleChild.count(it.first) == 0)
      continue;
    std::string modulePath;
    for (const auto level : it.first) {
      if (!modulePath.empty())
        modulePath += ':';
      modulePath += io::stringify(level);
    }
    throw InfomapError(ExitCode::InputError,
                       fmt::format(FMT_STRING("Cluster data mixes depths under module {}: it has both a leaf node and a "
                                              "sub-module as children. Infomap's own output never has that shape, and the "
                                              "codelength of such a module is not defined, so the tree cannot be used as an "
                                              "initial partition. Give every node under module {} the same depth, or pass "
                                              "--two-level to keep only the top-level assignment."),
                                   modulePath,
                                   modulePath));
  }

  // Tear down the existing partition. Detach the leaves first so that the
  // module destruction below does not also destroy them, then collect every
  // released leaf for re-attachment.
  std::map<unsigned int, unsigned int> nodeIdToIndex;
  unsigned int leafIndex = 0;
  for (auto* leafNode : m_leafNodes) {
    if (leafNode->parent != nullptr)
      leafNode->parent->releaseChildren();
    leafNode->parent = nullptr;
    leafNode->next = nullptr;
    leafNode->previous = nullptr;
    nodeIdToIndex[leafNode->stateId] = leafIndex;
    ++leafIndex;
  }
  m_root.deleteChildren();

  // Build the module tree by walking each path top-down and creating only
  // missing prefixes — this guarantees that two leaves sharing a path prefix
  // share the same ancestor modules, which is what determines codelength.
  std::map<Path, InfoNode*> moduleByPathPrefix;

  auto numNodesFound = 0;
  auto numNodesNotInNetwork = 0;
  for (const auto& nodePath : tree) {
    const auto& path = nodePath.second;
    if (path.empty())
      continue;

    ++numNodesFound;
    InfoNode* leafNode = nullptr;
    const auto nodeIdIt = nodeIdToIndex.find(nodePath.first);
    if (nodeIdIt == nodeIdToIndex.end()) {
      ++numNodesNotInNetwork;
      continue;
    }
    leafNode = m_leafNodes[nodeIdIt->second];

    InfoNode* parent = &m_root;
    Path prefix;
    for (size_t i = 0; i + 1 < path.size(); ++i) {
      prefix.push_back(path[i]);
      auto inserted = moduleByPathPrefix.emplace(prefix, nullptr);
      if (inserted.second) {
        auto* module = &allocNode();
        parent->addChild(module);
        inserted.first->second = module;
      }
      parent = inserted.first->second;
    }

    parent->addChild(leafNode);
    if (static_cast<int>(path.size()) > maxDepth)
      maxDepth = static_cast<int>(path.size());
  }

  if (numNodesNotInNetwork > 0) {
    Console::note(1, "{}/{} nodes in tree not found in network.", numNodesNotInNetwork, numNodesFound);
  }

  auto numNodesAddedToNeighbouringModules = 0;
  auto numNodesWithoutClusterInfo = 0;

  // Set orphaned nodes to their own or neighbouring module
  for (auto& leafNode : m_leafNodes) {
    if (leafNode->parent != nullptr) {
      continue;
    }

    ++numNodesWithoutClusterInfo;
    if (assignToNeighbouringModule) {
      // Take first neighbour that has a module assigned
      for (auto& edge : leafNode->outEdges()) {
        if (edge->target->parent != nullptr) {
          edge->target->parent->addChild(leafNode);
          ++numNodesAddedToNeighbouringModules;
          break;
        }
      }
      // Check incoming links if still orphan
      if (leafNode->parent == nullptr) {
        for (auto& edge : leafNode->inEdges()) {
          if (edge->source->parent != nullptr) {
            edge->source->parent->addChild(leafNode);
            ++numNodesAddedToNeighbouringModules;
            break;
          }
        }
      }
      // Set to own module if no neighbour module available
      if (leafNode->parent == nullptr) {
        auto* module = &allocNode();
        root().addChild(module);
        module->addChild(leafNode);
      }
    } else {
      // Set to own module if no neighbour module available
      auto* module = &allocNode();
      root().addChild(module);
      module->addChild(leafNode);
    }
  }
  if (numNodesWithoutClusterInfo > 0) {
    if (assignToNeighbouringModule) {
      Console::note(1, "{} nodes not found in tree, {} are put into neighbouring modules and {} into separate modules.", numNodesWithoutClusterInfo, numNodesAddedToNeighbouringModules, numNodesWithoutClusterInfo - numNodesAddedToNeighbouringModules);
    } else {
      Console::note(1, "{} nodes not found in tree are put into separate modules.", numNodesWithoutClusterInfo);
    }
  }

  aggregateFlowValuesFromLeafToRoot();
  initTree();
  initNetwork();
  m_numNonTrivialTopModules = calculateNumNonTrivialTopModules();
  // Init partition to calculate indexCodelength and moduleCodelength
  if (haveModules())
    setActiveNetworkFromChildrenOfRoot();
  else
    setActiveNetworkFromLeafs();
  initPartition();

  m_hierarchicalCodelength = calcCodelengthOnTree(root(), true);
  Console().status("Initial", fmt::format(FMT_STRING("codelength {:g}, {} levels, {} top modules"), m_hierarchicalCodelength, maxDepth, numTopModules()));
  Log(1) << std::setprecision(6);
  return *this;
}

/**
 * clusterIds is a map from nodeId -> clusterId
 */
InfomapBase& InfomapBase::initPartition(const std::map<unsigned int, unsigned int>& clusterIds, bool hard)
{
  // Generate map from node id to index in leaf node vector
  std::map<unsigned int, unsigned int> nodeIdToIndex;
  unsigned int leafIndex = 0;
  for (auto& nodePtr : m_leafNodes) {
    auto& leafNode = *nodePtr;
    nodeIdToIndex[leafNode.stateId] = leafIndex;
    ++leafIndex;
  }

  // Remap clusterIds to vector from leaf node index to module index instead of nodeId -> clusterId
  std::map<unsigned int, std::set<unsigned int>> clusterIdToNodeIds;
  for (auto it : clusterIds) {
    clusterIdToNodeIds[it.second].insert(it.first);
  }
  unsigned int numNodes = numLeafNodes();
  std::vector<unsigned int> modules(numNodes);
  std::vector<unsigned int> selectedNodes(numNodes, 0);
  unsigned int moduleIndex = 0;
  unsigned int numClusterNodeIdsNotInNetwork = 0;
  for (const auto& it : clusterIdToNodeIds) {
    auto& nodes = it.second;
    bool moduleHasMember = false;
    for (auto& nodeId : nodes) {
      // find, not operator[]: a node id the network does not have would otherwise be
      // inserted with a default index of 0, quietly reassigning the *first* leaf node
      // to this module and marking it as covered -- which also suppressed the note
      // below about the nodes that really are missing an assignment (#898).
      const auto nodeIndexIt = nodeIdToIndex.find(nodeId);
      if (nodeIndexIt == nodeIdToIndex.end()) {
        ++numClusterNodeIdsNotInNetwork;
        continue;
      }
      modules[nodeIndexIt->second] = moduleIndex;
      ++selectedNodes[nodeIndexIt->second];
      moduleHasMember = true;
    }
    // Only claim an index for a module that got a member, so a cluster whose ids are
    // all unknown does not leave an empty module behind.
    if (moduleHasMember)
      ++moduleIndex;
  }

  if (numClusterNodeIdsNotInNetwork > 0) {
    Console::note(1, "{} node ids in cluster file not found in network.", numClusterNodeIdsNotInNetwork);
  }

  unsigned int numNodesWithoutClusterInfo = static_cast<unsigned int>(
      std::count_if(selectedNodes.begin(), selectedNodes.end(), [](unsigned int count) { return count == 0; }));

  if (numNodesWithoutClusterInfo == 0) {
    return initPartition(modules, hard);
  }

  if (assignToNeighbouringModule) {
    Console::note(1, "{} nodes not found in cluster file are put into neighbouring modules if possible.", numNodesWithoutClusterInfo);
    for (unsigned int i = 0; i < numNodes; ++i) {
      if (selectedNodes[i] == 0) {
        // Check out edges greedily for connected modules
        auto& node = *m_leafNodes[i];
        for (auto& e : node.outEdges()) {
          auto& edge = *e;
          auto targetNodeIndex = nodeIdToIndex[edge.target->stateId];
          if (selectedNodes[targetNodeIndex] != 0) {
            selectedNodes[i] = 1;
            modules[i] = modules[targetNodeIndex];
            break;
          }
        }
        if (selectedNodes[i] == 0) {
          // Check in edges greedily for connected modules
          for (auto& e : node.inEdges()) {
            auto& edge = *e;
            auto sourceNodeIndex = nodeIdToIndex[edge.source->stateId];
            if (selectedNodes[sourceNodeIndex] != 0) {
              selectedNodes[i] = 1;
              modules[i] = modules[sourceNodeIndex];
              break;
            }
          }
        }
        if (selectedNodes[i] == 0) {
          // No connected node with a module index, fall back to create new module
          selectedNodes[i] = 1;
          modules[i] = moduleIndex;
          ++moduleIndex;
        }
      }
    }
  } else {
    Console::note(1, "{} nodes not found in cluster file are put into separate modules.", numNodesWithoutClusterInfo);
    // Put non-selected nodes in its own module
    for (unsigned int i = 0; i < numNodes; ++i) {
      if (selectedNodes[i] == 0) {
        modules[i] = moduleIndex;
        ++moduleIndex;
      }
    }
  }

  return initPartition(modules);
}

InfomapBase& InfomapBase::initPartition(std::vector<std::vector<unsigned int>>& clusters, bool hard)
{
  Log(3).print("\n -> Init {} partition... ", hard ? "hard" : "soft") << std::flush;
  unsigned int numNodes = numLeafNodes();
  // Get module indices from the merge structure
  std::vector<unsigned int> modules(numNodes);
  std::vector<unsigned int> selectedNodes(numNodes, 0);
  unsigned int moduleIndex = 0;

  for (auto& cluster : clusters) {
    if (cluster.empty())
      continue;
    for (auto id : cluster) {
      ++selectedNodes[id];
      modules[id] = moduleIndex;
    }
    ++moduleIndex;
  }

  // Put non-selected nodes in its own module
  unsigned int numNodesWithoutClusterInfo = 0;
  for (unsigned int i = 0; i < numNodes; ++i) {
    if (selectedNodes[i] == 0) {
      modules[i] = moduleIndex;
      ++moduleIndex;
      ++numNodesWithoutClusterInfo;
    }
  }
  if (numNodesWithoutClusterInfo > 0)
    Console::note(1, "{} nodes not found in cluster file are put into separate modules.", numNodesWithoutClusterInfo);

  return initPartition(modules, hard);
}

InfomapBase& InfomapBase::initPartition(std::vector<unsigned int>& modules, bool hard)
{
  removeModules();
  setActiveNetworkFromLeafs();
  initPartition(); // TODO: confusing same name, should be able to init default without arguments here too
  moveActiveNodesToPredefinedModules(modules);
  consolidateModules(false);

  if (hard) {
    // Save the original network
    m_originalLeafNodes.swap(m_leafNodes);

    // Use the pre-partitioned network as the new leaf network
    unsigned int numNodesInNewNetwork = numTopModules();
    m_leafNodes.resize(numNodesInNewNetwork);
    unsigned int nodeIndex = 0;
    unsigned int numLinksInNewNetwork = 0;
    for (InfoNode& node : m_root) {
      m_leafNodes[nodeIndex] = &node;
      // Collapse children
      node.collapseChildren();
      numLinksInNewNetwork += node.outDegree();
      ++nodeIndex;
    }

    Console::detail(1, "hard-partitioned the network to {} nodes and {} links with codelength {}", numNodesInNewNetwork, numLinksInNewNetwork, io::stringify(*this));
  } else {
    Console().status("Initial", fmt::format(FMT_STRING("codelength {}, {} top modules"), io::stringify(*this), numTopModules()));
  }
  m_hierarchicalCodelength = getCodelength();

  return *this;
}

void InfomapBase::generateSubNetwork(Network& network)
{
  unsigned int numNodes = network.numNodes();
  auto& metaData = network.metaData();
  numMetaDataDimensions = network.numMetaDataColumns();

  Console::detail(1, "build internal network with {} nodes and {} links", numNodes, network.numLinks());
  if (!metaData.empty()) {
    Console::detail(1, "with {} meta-data records in {} dimensions", metaData.size(), numMetaDataDimensions);
  }

  m_leafNodes.reserve(numNodes);
  m_nodePool.reserve(numNodes);
  m_edgePool.reserve(network.numLinks());
  double sumNodeFlow = 0.0;
  double sumTeleFlow = 0.0;
  std::unordered_map<unsigned int, unsigned int> nodeIndexMap;
  nodeIndexMap.reserve(numNodes);
  for (auto& nodeIt : network.nodes()) {
    auto& networkNode = nodeIt.second;
    auto* node = &allocNode(networkNode.flow, networkNode.id, networkNode.physicalId, networkNode.layerId);
    node->data.teleportWeight = networkNode.weight;
    node->data.teleportFlow = networkNode.teleFlow;
    node->data.exitFlow = networkNode.exitFlow;
    node->data.enterFlow = networkNode.enterFlow;
    // Only regularized-multilayer flow consumes layerTeleFlowData; populating it
    // unconditionally would allocate an InfoNodeExtras for every node and defeat
    // the out-of-line shrink (extras must stay null in the common case).
    if (isRegularizedMultilayerFlow())
      node->ensureExtras().layerTeleFlowData.emplace_back(networkNode.layerId, networkNode.intraLayerTeleFlow, networkNode.intraLayerTeleWeight);
    // Log(1) << "Node " << node->stateId << " (" << node->layerId << "," << node->physicalId << ") data: " << node->data << ", tele-flow: " << networkNode.teleFlow << ", intra-tele: " << networkNode.intraLayerTeleFlow << "\n";
    if (haveMetaData()) {
      auto meta = metaData.find(networkNode.id);
      if (meta != metaData.end()) {
        node->ensureExtras().metaData = meta->second;
      } else {
        node->ensureExtras().metaData = std::vector<int>(numMetaDataDimensions, -1);
      }
    }
    sumNodeFlow += networkNode.flow;
    sumTeleFlow += networkNode.teleFlow;
    m_root.addChild(node);
    nodeIndexMap[networkNode.id] = m_leafNodes.size();
    m_leafNodes.push_back(node);
  }
  m_root.data.flow = sumNodeFlow;
  m_root.data.teleportFlow = sumTeleFlow;
  // m_calculateEnterExitFlow = true; //TODO: Implement always in flow calculation
  if (!this->regularized) {
    m_calculateEnterExitFlow = true;
  }
  if (std::abs(sumNodeFlow - 1.0) > 1e-10)
    Console::warn(1, "Total flow on nodes differs from 1.0 by {:g}.", sumNodeFlow - 1.0);

  bool changeMarkovTime = std::abs(markovTime - 1) > 1e-3;
  if (changeMarkovTime) {
    Console::detail(1, "rescale link flow with global Markov time {:g}", markovTime);
  }

  network.forEachLink([&](unsigned int sourceIndex, unsigned int targetIndex, double weight, double& flow) {
    // Ignore self-links in optimization as it doesn't change enter/exit flow on modular level
    if (sourceIndex != targetIndex) {
      m_leafNodes[sourceIndex]->addOutEdge(*m_leafNodes[targetIndex], weight, flow * markovTime);
    }
  });

  if (variableMarkovTime) {
    Console::detail(1, "rescale link flow with variable Markov time");
    if (std::abs(variableMarkovTimeDamping - 1) > 1e-9) {
      Console::detail(1, "use variable Markov time strength {:g}", variableMarkovTimeDamping);
    }
  }

  double maxEntropy = 0.0;
  double maxFlow = 0.0;
  double entropyRate = 0.0;
  unsigned int maxDegree = 0;
  unsigned int maxOutDegree = 0;
  unsigned int maxInDegree = 0;
  double totDegree = network.sumDegree();
  std::vector<double> entropies(numNodes, 0);

  for (unsigned i = 0; i < numNodes; ++i) {
    InfoNode& node = *m_leafNodes[i];
    maxDegree = std::max(maxDegree, node.degree());
    maxOutDegree = std::max(maxOutDegree, node.outDegree());
    maxInDegree = std::max(maxInDegree, node.inDegree());
    double entropy = 0;
    double sumOut = 0;
    for (InfoEdge* e : node.outEdges()) {
      entropy -= infomath::plogp(e->data.weight);
      sumOut += e->data.weight;
    }
    if (isUndirectedFlow()) {
      for (InfoEdge* e : node.inEdges()) {
        entropy -= infomath::plogp(e->data.weight);
        sumOut += e->data.weight;
      }
    }
    entropy = sumOut > 1e-9 ? (entropy + infomath::plogp(sumOut)) / sumOut : 0;
    maxEntropy = std::max(maxEntropy, entropy);
    entropyRate += m_leafNodes[i]->data.flow * entropy;
    maxFlow = std::max(maxFlow, node.data.flow);
    entropies[i] = entropy; // Store for undirected networks
  }

  m_entropyRate = entropyRate;
  m_maxEntropy = maxEntropy;
  m_maxFlow = maxFlow;

  double minLocalScale = variableMarkovTimeMinLocalScale;
  double damping = variableMarkovTimeDamping;

  double maxScale = infomath::linlog(maxFlow * totDegree, damping);

  if (variableMarkovTime) {
    if (damping < 0) {
      maxScale = infomath::linlog(pow(2.0, maxEntropy), -damping);
    }
    for (unsigned i = 0; i < numNodes; ++i) {
      InfoNode& node = *m_leafNodes[i];
      double localScale = damping < 0 ? infomath::linlog(pow(2.0, entropies[i]), -damping) : infomath::linlog(std::max(minLocalScale, node.data.flow * totDegree), damping);
      for (InfoEdge* e : node.outEdges()) {
        if (isUndirectedFlow()) {
          double oppositeLocalScale = damping < 0 ? infomath::linlog(pow(2.0, entropies[nodeIndexMap[e->target->stateId]]), -damping) : infomath::linlog(std::max(minLocalScale, e->target->data.flow * totDegree), damping);
          localScale = std::max(localScale, oppositeLocalScale);
        }
        double localMarkovTimeScale = maxScale / std::max(minLocalScale, localScale);
        e->data.flow *= localMarkovTimeScale;
        // Note: deliberately do NOT write the scaled flow back into the shared StateNetwork
        // (network.nodeLinkMap()). The optimizer and output use only the per-instance InfoEdge
        // flow (e->data.flow), so the write-back was a non-idempotent side effect that made
        // --parallel-trials unsafe (workers share one Network). Library getLinks(flow=True) /
        // getLinkResults() now return the input link flow rather than the VMT-scaled flow.
      }
    }
  }
}

void InfomapBase::generateSubNetwork(InfoNode& parent)
{
  // Store parent module
  m_root.owner = &parent;
  m_root.data = parent.data;
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
  m_root.lossyEntropy = parent.lossyEntropy;
  m_root.lossyFlowLogFlow = parent.lossyFlowLogFlow;
#endif

  unsigned int numNodes = parent.childDegree();
  m_leafNodes.resize(numNodes);
  // Reserve both pools to the exact clone counts. The per-pool ramp slack is
  // small, but the recursive phase keeps many sub-Infomaps alive at once and
  // the summed slack dominates the pool RSS overhead on higher-order networks;
  // exact first chunks remove it. Module nodes allocated later still ramp from
  // kInitialChunk, which stays tight for the tiny refined modules.
  m_nodePool.reserve(numNodes);
  unsigned int numInternalEdges = 0;
  for (InfoNode& node : parent) {
    for (InfoEdge* e : node.outEdges()) {
      if (e->target->parent == &parent)
        ++numInternalEdges;
    }
  }
  m_edgePool.reserve(numInternalEdges);

  Console::detail(1, "generate sub network with {} nodes", numNodes);

  unsigned int childIndex = 0;
  for (InfoNode& node : parent) {
    auto* clonedNode = &allocNode(node);
    clonedNode->initClean();
    m_root.addChild(clonedNode);
    node.index = childIndex; // Set index to its place in this subnetwork to be able to find edge target below
    m_leafNodes[childIndex] = clonedNode;
    ++childIndex;
  }

  InfoNode* parentPtr = &parent;
  // Clone edges
  for (InfoNode& node : parent) {
    for (InfoEdge* e : node.outEdges()) {
      InfoEdge& edge = *e;
      // If neighbour node is within the same module, add the link to this subnetwork.
      if (edge.target->parent == parentPtr) {
        m_leafNodes[edge.source->index]->addOutEdge(*m_leafNodes[edge.target->index], edge.data.weight, edge.data.flow);
      }
    }
  }
}

void InfomapBase::init()
{
  Log(3) << "InfomapBase::init()...\n";

  if (m_calculateEnterExitFlow && isMainInfomap())
    initEnterExitFlow();

  initNetwork();

  m_oneLevelCodelength = calcCodelength(m_root);
  Console::detail(1, "one-level codelength: {}", io::toPrecision(m_oneLevelCodelength));
}

#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
std::vector<unsigned int> InfomapBase::noiseTopModules() const
{
  std::vector<unsigned int> noise;
  if (!lossy)
    return noise;
  // A module is noise iff its naming overhead exceeds the priced dynamics it covers:
  // l_i = plogp(P_i) - sum_alpha plogp(p_alpha) > lambda * H_i. The per-module leaf
  // sums (P_i, F_i, H_i) are carried on the module node by consolidateModules.
  unsigned int index = 1;
  for (const InfoNode& module : m_root) {
    const double loss = infomath::plogp(module.data.flow) - module.lossyFlowLogFlow;
    if (loss > lossyLambda * module.lossyEntropy + 1e-12)
      noise.push_back(index);
    ++index;
  }
  return noise;
}
#endif

// ===================================================
// Run: *
// ===================================================

void InfomapBase::hierarchicalPartition()
{
  Console::detail(1, "hierarchical partition");

  const auto depth = maxTreeDepth();
  if (depth > 2) {
    Console::detail(1, "continuing from a tree with {} levels", depth);

    if (fastHierarchicalSolution == 0) {
      Console::detail(1, "removing sub modules");
      removeSubModules(true);
      m_hierarchicalCodelength = calcCodelengthOnTree(root(), true);
    }

    else if (fastHierarchicalSolution == 1) {
      const bool shouldRestoreSilent = isMainInfomap() && !Log::isThreadMuted();
      bool isSilent = shouldRestoreSilent && Log::isSilent();

      double codelengthBefore = 0.0;
      double codelengthAfter = 0.0;

      for (InfoNode& module : m_root.infomapTree()) {
        if (!module.isLeaf() && module.firstChild->isLeafModule()) {
          codelengthBefore += module.codelength;
          auto numLeafs = 0;
          for (auto& leafModule : module) {
            numLeafs += leafModule.childDegree();
          }

          std::vector<unsigned int> modules(numLeafs);
          std::vector<InfoNode*> leafs(numLeafs);

          auto i = 0;
          for (auto it = module.begin_tree(); !it.isEnd(); ++it) {
            if (it->isLeaf()) {
              modules[i] = it.moduleIndex();
              leafs[i] = it.current();
              ++i;
            }
          }

          module.replaceChildrenWithGrandChildren();

          auto& subInfomap = getSubInfomap(module).initNetwork(module);

          subInfomap.initPartition(modules);

          // Run two-level partition + find hierarchically super modules (skip recursion).
          // Free the sub-Infomap on interrupt before unwinding (issue #412).
          try {
            subInfomap.setOnlySuperModules(true).run();
          } catch (const InterruptionError&) {
            module.disposeInfomap();
            throw;
          }

          // Collect sub Infomap modules
          i = 0;
          for (auto& subLeafPtr : subInfomap.leafNodes()) {
            modules[i] = subLeafPtr->index;
            ++i;
          }

          // Create new sub modules
          std::vector<InfoNode*> subModules(numLeafs, nullptr);
          module.releaseChildren();

          for (auto j = 0; j < numLeafs; ++j) {
            InfoNode* leaf = leafs[j];
            unsigned int moduleIndex = modules[j];
            if (subModules[moduleIndex] == nullptr) {
              subModules[moduleIndex] = &allocNode(subInfomap.leafNodes()[j]->parent->data);
              subModules[moduleIndex]->index = moduleIndex;
              module.addChild(subModules[moduleIndex]);
            }
            subModules[moduleIndex]->addChild(leaf);
          }
          module.disposeInfomap();
          module.codelength = calcCodelength(module);
          codelengthAfter += module.codelength;
        }
      }

      if (shouldRestoreSilent)
        Log::setSilent(isSilent);

      const double diffCodelength = codelengthBefore - codelengthAfter;
      Console::detail(1, "fine-tune bottom modules: improvement {:g}% {} codelength {}", (diffCodelength / codelengthBefore) * 100, Console::arrow(), io::toPrecision(codelengthAfter));
    }

    recursivePartition();
    return;
  }

  partition();

  if (numTopModules() == 1 || numTopModules() == numLeafNodes()) {
    Console::detail(1, "trivial partition, skipping hierarchical search");
    return;
  }

  if (numTopModules() > preferredNumberOfModules) {
    findHierarchicalSuperModules();
  }

  if (onlySuperModules) {
    removeSubModules(true);
    m_hierarchicalCodelength = calcCodelengthOnTree(root(), true);
    return;
  }

  if (fastHierarchicalSolution >= 2) {
    // FIXME Calculate individual module codelengths and store on the modules?
    return;
  }

  if (fastHierarchicalSolution == 0) {
    removeSubModules(true);
    m_hierarchicalCodelength = calcCodelengthOnTree(root(), true);
  }

  recursivePartition();
}

void InfomapBase::partition()
{
  std::vector<std::string> prettyCompression;
  Console().section("Optimization");
  auto progress = Console().progress("Two-level");
  progress.update("working…");
  double initialCodelength = m_oneLevelCodelength;
  double oldCodelength = initialCodelength;

  m_tuneIterationIndex = 0;
  findTopModulesRepeatedly(levelAggregationLimit);

  double newCodelength = getCodelength();
  double compression = oldCodelength < 1e-16 ? 0.0 : (oldCodelength - newCodelength) / oldCodelength;
  prettyCompression.push_back(Console::percent(compression * 100));
  oldCodelength = newCodelength;
  progress.update(fmt::format(FMT_STRING("working… codelength {} · {} modules"), io::toPrecision(getCodelength()), numTopModules()));

  bool doFineTune = true;
  bool coarseTuned = false;
  while (numTopModules() > 1 && (m_tuneIterationIndex + 1) != tuneIterationLimit) {
    ++m_tuneIterationIndex;
    if (doFineTune) {
      Log(3) << "\n";
      unsigned int numEffectiveLoops = fineTune();
      if (numEffectiveLoops > 0) {
        Console::detail(1, 3, "fine-tune: done in {} effective loops {} codelength {}, {} modules", numEffectiveLoops, Console::arrow(), io::stringify(*this), numTopModules());
        // Continue to merge fine-tuned modules
        findTopModulesRepeatedly(levelAggregationLimit);
      } else {
        Console::detail(1, 3, "fine-tune: no improvement");
      }
    } else {
      coarseTuned = true;
      if (!noCoarseTune) {
        Log(3) << "\n";
        unsigned int numEffectiveLoops = coarseTune();
        if (numEffectiveLoops > 0) {
          Console::detail(1, 3, "coarse-tune: done in {} effective loops {} codelength {}, {} modules", numEffectiveLoops, Console::arrow(), io::stringify(*this), numTopModules());
          // Continue to merge fine-tuned modules
          findTopModulesRepeatedly(levelAggregationLimit);
        } else {
          Console::detail(1, 3, "coarse-tune: no improvement");
        }
      }
    }
    newCodelength = getCodelength();
    compression = oldCodelength < 1e-16 ? 0.0 : (oldCodelength - newCodelength) / oldCodelength;
    bool isImprovement = newCodelength <= oldCodelength - minimumCodelengthImprovement && newCodelength < oldCodelength - initialCodelength * minimumRelativeTuneIterationImprovement;
    if (!isImprovement) {
      // Make sure coarse-tuning have been tried
      if (coarseTuned)
        break;
    } else {
      oldCodelength = newCodelength;
    }
    prettyCompression.push_back(Console::percent(compression * 100));
    progress.update(fmt::format(FMT_STRING("working… codelength {} · {} modules"), io::toPrecision(getCodelength()), numTopModules()));
    doFineTune = !doFineTune;
  }

  std::string nonTrivialSummary;
  if (m_numNonTrivialTopModules != numTopModules())
    nonTrivialSummary = fmt::format(FMT_STRING(" ({} non-trivial)"), m_numNonTrivialTopModules);
  progress.finish(fmt::format(FMT_STRING("{} {} codelength {}, {} modules{}"), compressionSummary(prettyCompression), Console::arrow(), io::toPrecision(getCodelength()), numTopModules(), nonTrivialSummary));
  const bool regularizedPriorOnly = regularized && network().numLinks() == 0;
  if (!preferModularSolution && preferredNumberOfModules == 0 && (haveNonTrivialModules() || regularizedPriorOnly) && getCodelength() > getOneLevelCodelength()) {
    Console::detail(1, "worse codelength than one-level, putting all nodes in one module");

    // Create new single module between modules and root
    auto& module = root().replaceChildrenWithOneNode();
    // TODO: Extract copying from root and resetting index to a method, also copy metadata?
    module.data = m_root.data;
    module.physicalNodes = m_root.physicalNodes;
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
    module.lossyEntropy = m_root.lossyEntropy;
    module.lossyFlowLogFlow = m_root.lossyFlowLogFlow;
#endif
    module.index = 0;
    for (auto& node : module) {
      node.index = 0;
    }
    module.codelength = getOneLevelCodelength();
    m_hierarchicalCodelength = getOneLevelCodelength();
    // calcCodelengthOnTree(root(), true);
    m_root.codelength = 0.0;
    m_numNonTrivialTopModules = calculateNumNonTrivialTopModules();

  } else {
    // Set consolidated cluster index on nodes and modules
    for (auto clusterIt = m_root.begin_tree(); !clusterIt.isEnd(); ++clusterIt) {
      clusterIt->index = clusterIt.moduleIndex();
    }

    // Calculate individual module codelengths and store on the modules
    calcCodelengthOnTree(root(), false);
    m_root.codelength = getIndexCodelength();
    m_hierarchicalCodelength = getCodelength();
  }
}

void InfomapBase::restoreHardPartition()
{
  // Collect all leaf nodes in a separate sequence to be able to iterate independent of tree structure
  std::vector<InfoNode*> leafNodes(numLeafNodes());
  unsigned int leafIndex = 0;
  for (InfoNode& node : m_root.infomapTree()) {
    if (node.isLeaf()) {
      leafNodes[leafIndex] = &node;
      ++leafIndex;
    }
  }
  // Expand all children
  unsigned int numExpandedChildren = 0;
  unsigned int numExpandedNodes = 0;
  for (InfoNode* node : leafNodes) {
    ++numExpandedNodes;
    numExpandedChildren += node->expandChildren();
    node->replaceWithChildren();
  }

  // Swap back original leaf nodes
  m_originalLeafNodes.swap(m_leafNodes);

  Console::detail(1, "expanded {} hard modules to {} original nodes", numExpandedNodes, numExpandedChildren);
}

// ===================================================
// Run: Init: *
// ===================================================

void InfomapBase::initEnterExitFlow()
{
  // Not done in Bayesian
  // TODO: Skip this, always add enter/exit/tele flow from flow calculator
  // Calculate enter/exit
  double alpha = teleportationProbability;
  double beta = 1.0 - alpha;

  for (auto* n : m_leafNodes) {
    n->data.enterFlow = n->data.exitFlow = 0.0;
  }
  if (!isUndirectedClustering()) {
    for (auto* n : m_leafNodes) {
      auto& node = *n;
      node.data.teleportFlow = alpha * node.data.flow;
      node.data.teleportSourceFlow = node.data.flow;
      if (node.isDangling()) {
        node.data.teleportFlow = node.data.flow;
        node.data.danglingFlow = node.data.flow;
        m_sumDanglingFlow += node.data.flow;
      }
    }
    for (auto* n : m_leafNodes) {
      auto& node = *n;
      for (InfoEdge* e : node.outEdges()) {
        InfoEdge& edge = *e;
        // Self-links not included here, should not add to enter and exit flow in its enclosing module
        edge.source->data.exitFlow += edge.data.flow;
        edge.target->data.enterFlow += edge.data.flow;
      }
      if (recordedTeleportation) {
        // Don't let self-teleportation add to the enter/exit flow (i.e. multiply with (1.0 - node.data.teleportWeight))
        node.data.exitFlow += (alpha * node.data.flow + beta * node.data.danglingFlow) * (1.0 - node.data.teleportWeight);
        node.data.enterFlow += (alpha * (1.0 - node.data.flow) + beta * (m_sumDanglingFlow - node.data.danglingFlow)) * node.data.teleportWeight;
      }
    }
  } else {
    for (auto* n : m_leafNodes) {
      auto& node = *n;
      node.data.teleportFlow = alpha * node.data.flow;
      node.data.teleportSourceFlow = node.data.flow;
      if (node.isDangling()) {
        node.data.teleportFlow = node.data.flow;
        node.data.danglingFlow = node.data.flow;
      }
      for (InfoEdge* e : node.outEdges()) {
        InfoEdge& edge = *e;
        // Self-links not included here, should not add to enter and exit flow in its enclosing module
        double halfFlow = edge.data.flow / 2;
        edge.source->data.exitFlow += halfFlow;
        edge.target->data.exitFlow += halfFlow;
        edge.source->data.enterFlow += halfFlow;
        edge.target->data.enterFlow += halfFlow;
      }
    }
  }
}

// Aggregate node and enter/exit flow to all tree nodes
void InfomapBase::aggregateFlowValuesFromLeafToRoot()
{
  for (auto& node : root().infomapTree()) {
    if (!node.isLeaf()) {
      node.data = {};
    }
  }

  // Aggregate flow from leaf nodes to root node
  unsigned int numLevels = 0;
  for (auto it = root().begin_post_depth_first(); !it.isEnd(); ++it) {
    auto& node = *it;
    if (!node.isRoot()) {
      node.parent->data.flow += node.data.flow;
      node.parent->data.teleportFlow += node.data.teleportFlow;
      node.parent->data.teleportSourceFlow += node.data.teleportSourceFlow;
      node.parent->data.teleportWeight += node.data.teleportWeight;
      node.parent->data.danglingFlow += node.data.danglingFlow;
    }
    // Don't aggregate enter and exit flow
    if (!node.isLeaf()) {
      node.index = it.depth(); // Use index to store the depth on modules
      node.data.enterFlow = 0.0;
      node.data.exitFlow = 0.0;
    } else
      numLevels = std::max(numLevels, it.depth());
  }

  if (std::abs(root().data.flow - 1.0) > 1e-10) {
    Console::warn(1, "Aggregated flow is not exactly 1.0, but {:.10g}.", root().data.flow);
  }

  // Aggregate enter and exit flow between modules
  for (auto& leafNode : m_leafNodes) {
    auto& leafNodeSource = *leafNode;
    for (InfoEdge* e : leafNodeSource.outEdges()) {
      auto& edge = *e;
      auto& leafNodeTarget = edge.target;
      double linkFlow = edge.data.flow;
      double halfFlow = linkFlow / 2;

      InfoNode* node1 = leafNodeSource.parent;
      InfoNode* node2 = leafNodeTarget->parent;

      if (node1 == node2)
        continue;

      // First aggregate link flow until equal depth
      while (node1->index > node2->index) {
        if (isUndirectedClustering()) {
          node1->data.exitFlow += halfFlow;
          node1->data.enterFlow += halfFlow;
        } else {
          node1->data.exitFlow += linkFlow;
        }
        node1 = node1->parent;
      }
      while (node2->index > node1->index) {
        if (isUndirectedClustering()) {
          node2->data.enterFlow += halfFlow;
          node2->data.exitFlow += halfFlow;
        } else {
          node2->data.enterFlow += linkFlow;
        }
        node2 = node2->parent;
      }

      // Then aggregate link flow until equal parent
      while (node1 != node2) {
        if (isUndirectedClustering()) {
          node1->data.exitFlow += halfFlow;
          node1->data.enterFlow += halfFlow;
          node2->data.enterFlow += halfFlow;
          node2->data.exitFlow += halfFlow;
        } else {
          node1->data.exitFlow += linkFlow;
          node2->data.enterFlow += linkFlow;
        }
        node1 = node1->parent;
        node2 = node2->parent;
      }
    }
  }
  if (recordedTeleportation) {
    Console::detail(1, "aggregating enter/exit flow for recorded teleportation, sum teleFlow {:g}", m_root.data.teleportFlow);

    for (auto& node : m_root.infomapTree()) {
      if (!node.isLeaf()) {
        // Don't code self-teleportation

        node.data.enterFlow += (m_root.data.teleportFlow - node.data.teleportFlow) * node.data.teleportWeight;
        node.data.exitFlow += node.data.teleportFlow * (1.0 - node.data.teleportWeight);
      }
    }
  }
}

double InfomapBase::calcCodelengthOnTree(InfoNode& root, bool includeRoot) const
{
  double totalCodelength = 0.0;
  // Calculate enter/exit flow (assuming 0 by default)
  for (auto& node : root.infomapTree()) {
    if (node.isLeaf() || (!includeRoot && node.isRoot()))
      continue;
    node.codelength = calcCodelength(node);
    totalCodelength += node.codelength;
  }
  return totalCodelength;
}

// ===================================================
// Run: Partition: *
// ===================================================

void InfomapBase::setActiveNetworkFromChildrenOfRoot()
{
  m_moduleNodes.resize(m_root.childDegree());
  unsigned int i = 0;
  for (auto& node : m_root)
    m_moduleNodes[i++] = &node;
  m_activeNetwork = &m_moduleNodes;
}

void InfomapBase::findTopModulesRepeatedly(unsigned int maxLevels)
{
  std::vector<std::string> passList; // "<nodes>*<loops>" per aggregation pass, for the -v trace
  Log(3).print("\nIteration {}:\n", m_tuneIterationIndex + 1);
  m_aggregationLevel = 0;
  unsigned int numLevelsConsolidated = numLevels() - 1;
  if (maxLevels == 0)
    maxLevels = std::numeric_limits<unsigned int>::max();

  std::string initialCodelength;

  // Reapply core algorithm on modular network, replacing modules with super modules
  while (numTopModules() > 1 && numLevelsConsolidated != maxLevels) {
    if (haveModules())
      setActiveNetworkFromChildrenOfRoot();
    else
      setActiveNetworkFromLeafs();
    initPartition();
    if (m_aggregationLevel == 0) {
      initialCodelength = io::stringify(*this);
    }

    const auto numActiveBeforeMove = activeNetwork().size();
    Log(3).print("Level {} (codelength: {}): Moving {} nodes... ", numLevelsConsolidated + 1, io::stringify(*this), activeNetwork().size()) << std::flush;
    // Core loop, merging modules
    unsigned int numOptimizationLoops = optimizeActiveNetwork();

    passList.push_back(fmt::format(FMT_STRING("{}*{}"), numActiveBeforeMove, numOptimizationLoops));
    Log(3).print("done! -> codelength {} in {} modules.\n", io::stringify(*this), numActiveModules());

    // If no improvement, revert codelength terms to the last consolidated state
    if (haveModules() && restoreConsolidatedOptimizationPointIfNoImprovement()) {
      Log(3).print("-> Restored to codelength {} in {} modules.\n", io::stringify(*this), numTopModules());
      break;
    }

    // Consolidate modules
    bool replaceExistingModules = haveModules();
    consolidateModules(replaceExistingModules);
    ++numLevelsConsolidated;
    ++m_aggregationLevel;
  }
  m_aggregationLevel = 0;

  m_numNonTrivialTopModules = calculateNumNonTrivialTopModules();

  Console::detail(1, 2, "iter {}: codelength {} {} {}, {} modules ({} non-trivial); passes ({}*loops): {}", m_tuneIterationIndex + 1, initialCodelength, Console::arrow(), io::stringify(*this), numTopModules(), m_numNonTrivialTopModules, m_isCoarseTune ? "modules" : "nodes", io::stringify(passList, ", "));
}

unsigned int InfomapBase::fineTune()
{
  if (numLevels() != 2)
    throw std::logic_error("InfomapBase::fineTune() called but numLevels != 2");

  setActiveNetworkFromLeafs();
  initPartition();

  // Make sure module nodes have correct index //TODO: Assume always correct here?
  unsigned int moduleIndex = 0;
  for (auto& module : m_root) {
    module.index = moduleIndex++;
  }

  // Put leaf modules in existing modules
  std::vector<unsigned int> modules(numLeafNodes());
  for (unsigned int i = 0; i < numLeafNodes(); ++i) {
    modules[i] = m_leafNodes[i]->parent->index;
  }
  moveActiveNodesToPredefinedModules(modules);

  Log(3).print(" -> moved to codelength {} in {} existing modules. Try tuning...\n", io::stringify(*this), numActiveModules());

  // Continue to optimize from there to tune leaf nodes
  unsigned int numEffectiveLoops = optimizeActiveNetwork();
  if (numEffectiveLoops == 0) {
    restoreConsolidatedOptimizationPointIfNoImprovement();
    Log(4) << "Fine-tune didn't improve solution, restoring last.\n";
  } else {
    // Delete existing modules and consolidate fine-tuned modules
    root().replaceChildrenWithGrandChildren();
    consolidateModules(false);
    Log(4).print("Fine-tune done in {} effective loops to codelength {} in {} modules.\n", numEffectiveLoops, io::stringify(*this), numTopModules());
  }

  return numEffectiveLoops;
}

unsigned int InfomapBase::coarseTune()
{
  if (numLevels() != 2)
    throw std::logic_error("InfomapBase::coarseTune() called but numLevels != 2");

  Log(4) << "Coarse-tune...\nPartition each module in sub-modules for coarse tune...\n";

  Log::ScopedMute muteNestedMainRun(isMainInfomap());

  unsigned int moduleIndexOffset = 0;
  for (auto& node : m_root) {
    // Don't search for sub-modules in too small modules
    if (node.childDegree() < 2) {
      for (auto& child : node) {
        child.index = moduleIndexOffset;
      }
      moduleIndexOffset += 1;
      continue;
    } else {
      InfomapBase& subInfomap = getSubInfomap(node)
                                    .setTwoLevel(true)
                                    .setTuneIterationLimit(1);
      // Free the sub-Infomap on interrupt before unwinding (issue #412).
      try {
        subInfomap.initNetwork(node).run();
      } catch (const InterruptionError&) {
        node.disposeInfomap();
        throw;
      }

      auto originalLeafIt = node.begin_child();
      for (auto& subLeafPtr : subInfomap.leafNodes()) {
        InfoNode& subLeaf = *subLeafPtr;
        originalLeafIt->index = subLeaf.index + moduleIndexOffset;
        ++originalLeafIt;
      }
      moduleIndexOffset += subInfomap.numTopModules();

      node.disposeInfomap();
    }
  }

  Log(4).print("Move leaf nodes to {} sub-modules... \n", moduleIndexOffset);
  // Put leaf modules in the calculated sub-modules
  std::vector<unsigned int> subModules(numLeafNodes());
  for (unsigned int i = 0; i < numLeafNodes(); ++i) {
    subModules[i] = m_leafNodes[i]->index;
  }

  setActiveNetworkFromLeafs();
  initPartition();
  moveActiveNodesToPredefinedModules(subModules);

  Log(4).print("Moved to {} modules...\n", numActiveModules());

  // Replace existing modules but store that structure on the sub-modules
  consolidateModules(true);

  Log(4).print("Consolidated {} sub-modules to codelength {}\n", numTopModules(), io::stringify(*this));

  Log(4) << "Move sub-modules to former modular structure... \n";
  // Put sub-modules in the former modular structure
  std::vector<unsigned int> modules(numTopModules());
  unsigned int iModule = 0;
  for (auto& subModule : m_root) {
    modules[iModule++] = subModule.index;
  }

  setActiveNetworkFromChildrenOfRoot();
  initPartition();
  // Move sub-modules to former modular structure
  moveActiveNodesToPredefinedModules(modules);

  Log(4).print("Tune sub-modules from codelength {} in {} modules... \n", io::stringify(*this), numActiveModules());
  // Continue to optimize from there to tune sub-modules
  unsigned int numEffectiveLoops = optimizeActiveNetwork();
  Log(4).print("Tuned sub-modules in {} effective loops to codelength {} in {} modules.\n", numEffectiveLoops, io::stringify(*this), numActiveModules());
  consolidateModules(true);
  Log(4).print("Consolidated tuned sub-modules to codelength {} in {} modules.\n", io::stringify(*this), numTopModules());
  return numEffectiveLoops;
}

unsigned int InfomapBase::calculateNumNonTrivialTopModules() const
{
  if (root().childDegree() == 1)
    return 0;
  unsigned int numNonTrivialTopModules = 0;
  for (auto& module : m_root) {
    if (module.childDegree() > 1)
      ++numNonTrivialTopModules;
  }
  return numNonTrivialTopModules;
}

unsigned int InfomapBase::calculateMaxDepth() const
{
  unsigned int maxDepth = 0;
  for (auto it(m_root.begin_tree()); !it.isEnd(); ++it) {
    if (it->isLeaf()) {
      maxDepth = std::max(maxDepth, it.depth());
    }
  }
  return maxDepth;
}

unsigned int InfomapBase::findHierarchicalSuperModules(unsigned int superLevelLimit)
{
  if (superLevelLimit == 0)
    return 0;

  unsigned int numLevelsCreated = 0;
  double oldIndexLength = getIndexCodelength();
  double hierarchicalCodelength = getCodelength();
  double workingHierarchicalCodelength = hierarchicalCodelength;

  const unsigned int prefBaseDepth = numLevels(); // #308: depth before super levels are added

  if (!haveModules())
    throw std::logic_error("Trying to find hierarchical super modules without any modules");

  std::vector<std::string> prettyCompression;
  auto progress = Console().progress("Super-level");
  progress.update("working…");

  // Add index codebooks as long as the code gets shorter
  do {
    Console::detail(1, "super iter {}: searching {} modules", numLevelsCreated + 1, numTopModules());
    const bool shouldMuteNestedMainRun = isMainInfomap() && !Log::isThreadMuted();
    bool isSilent = false;
    if (shouldMuteNestedMainRun) {
      isSilent = Log::isSilent();
      Log::setSilent(true);
    }
    InfoNode tmp(m_root.data); // Temporary owner for the super Infomap instance
    auto& superInfomap = getSuperInfomap(tmp)
                             .setTwoLevel(true);
    superInfomap.initNetwork(m_root, true); //.initSuperNetwork();
    // No interrupt catch needed: local `tmp` owns the super-Infomap, so
    // ~InfoNode frees it if run() throws InterruptionError (#412).
    superInfomap.run();
    if (shouldMuteNestedMainRun) {
      Log::setSilent(isSilent);
    }
    double superCodelength = superInfomap.getCodelength();
    double superIndexCodelength = superInfomap.getIndexCodelength();

    unsigned int numSuperModules = superInfomap.numTopModules();
    bool trivialSolution = numSuperModules == 1 || numSuperModules == numTopModules();

    if (trivialSolution) {
      Console::detail(1, "super: no non-trivial super modules found");
      break;
    }

    // #308: bias the super-level accept test by the depth it would create.
    const double superPrefBias = twoLevel
        ? 0.0
        : prefLevelsBias(preferredNumberOfLevels, preferredNumberOfLevelsStrength, prefBaseDepth + numLevelsCreated + 1);
    bool improvedCodelength = superCodelength < oldIndexLength - minimumCodelengthImprovement + superPrefBias;
    if (!improvedCodelength) {
      Console::detail(1, "super: index codebook not improved over one-level");
      break;
    }

    workingHierarchicalCodelength += superCodelength - oldIndexLength;

    const double superCompression = ((hierarchicalCodelength - workingHierarchicalCodelength)
                                     / hierarchicalCodelength * 100);
    prettyCompression.push_back(Console::percent(superCompression));
    progress.update(fmt::format(FMT_STRING("working… codelength {} · {} top modules"), io::toPrecision(workingHierarchicalCodelength), numSuperModules));
    Console::detail(1, "super iter {}: found {} super modules, codelength {:g} + {:g} + {:g} = {:g}", numLevelsCreated + 1, numSuperModules, superIndexCodelength, superCodelength - superIndexCodelength, workingHierarchicalCodelength - superCodelength, workingHierarchicalCodelength);

    // Merge current top modules to two-level solution found in the super Infomap instance.
    std::vector<unsigned int> modules(numTopModules());
    unsigned int iModule = 0;
    for (InfoNode* n : superInfomap.leafNodes()) {
      modules[iModule++] = n->index;
    }

    setActiveNetworkFromChildrenOfRoot();
    initPartition();
    // Move top modules to super modules
    moveActiveNodesToPredefinedModules(modules);

    // Consolidate the dynamic modules without replacing any existing ones.
    consolidateModules(false);

    Console::detail(1, "super: consolidated {} modules, index codelength {:g} {} {}", numTopModules(), oldIndexLength, Console::arrow(), io::stringify(*this));

    hierarchicalCodelength = workingHierarchicalCodelength;
    oldIndexLength = superIndexCodelength;

    ++numLevelsCreated;

    m_numNonTrivialTopModules = calculateNumNonTrivialTopModules();

  } while (m_numNonTrivialTopModules > 1 && numLevelsCreated != superLevelLimit);

  // Restore the temporary transformation of flow to enter flow on super modules
  resetFlowOnModules();

  progress.finish(fmt::format(FMT_STRING("{} {} codelength {}, {} top modules"), compressionSummary(prettyCompression), Console::arrow(), io::toPrecision(hierarchicalCodelength), numTopModules()));
  Console::detail(1, "super-level: added {} levels, codelength {}", numLevelsCreated, io::toPrecision(hierarchicalCodelength));

  // Calculate and store hierarchical codelengths
  m_hierarchicalCodelength = calcCodelengthOnTree(root(), true);

  return numLevelsCreated;
}

unsigned int InfomapBase::findHierarchicalSuperModulesFast(unsigned int superLevelLimit)
{
  if (superLevelLimit == 0)
    return 0;

  unsigned int numLevelsCreated = 0;
  double oldIndexLength = getIndexCodelength();
  double hierarchicalCodelength = getCodelength();
  double workingHierarchicalCodelength = hierarchicalCodelength;

  Log(1) << "\nFind hierarchical super modules iteratively...\n";
  Log(1) << "Fast super-level compression: " << std::setprecision(2) << std::flush;

  // Add index codebooks as long as the code gets shorter
  do {
    Log(1).print("Iteration {}, finding super modules fast to {}{}... \n", numLevelsCreated + 1, numTopModules(), haveModules() ? " modules" : " nodes");

    if (haveModules()) {
      setActiveNetworkFromChildrenOfRoot();
      transformNodeFlowToEnterFlow(m_root);
      initSuperNetwork();
    } else {
      setActiveNetworkFromLeafs();
    }

    initPartition();

    unsigned int numEffectiveLoops = optimizeActiveNetwork();

    double codelength = getCodelength();
    double indexCodelength = getIndexCodelength();
    unsigned int numSuperModules = numActiveModules();
    bool trivialSolution = numSuperModules == 1 || numSuperModules == activeNetwork().size();

    bool acceptSolution = !trivialSolution && codelength < oldIndexLength - minimumCodelengthImprovement;
    // Force at least one modular level!
    bool acceptByForce = !acceptSolution && !haveModules();
    if (acceptByForce)
      acceptSolution = true;

    workingHierarchicalCodelength += codelength - oldIndexLength;

    Log(1) << hideIf(!acceptSolution) << ((hierarchicalCodelength - workingHierarchicalCodelength) / hierarchicalCodelength * 100) << "% " << std::flush;
    Console::note(1, "Found {} super modules in {} effective loops with hierarchical codelength {:g} + {:g} = {:g}{}", numActiveModules(), numEffectiveLoops, indexCodelength, workingHierarchicalCodelength - indexCodelength, workingHierarchicalCodelength, acceptSolution ? "" : ", discarding the solution.");
    Log(1) << std::flush;
    Log(1).print("{:g} -> {}\n", oldIndexLength, io::stringify(*this));

    if (!acceptSolution) {
      restoreConsolidatedOptimizationPointIfNoImprovement(true);
      break;
    }

    // Consolidate the dynamic modules without replacing any existing ones.
    consolidateModules(false);

    hierarchicalCodelength = workingHierarchicalCodelength;
    oldIndexLength = indexCodelength;

    ++numLevelsCreated;

    m_numNonTrivialTopModules = calculateNumNonTrivialTopModules();

  } while (m_numNonTrivialTopModules > 1 && numLevelsCreated != superLevelLimit);

  // Restore the temporary transformation of flow to enter flow on super modules
  resetFlowOnModules();

  Log(1).print("to codelength {} in {} top modules.\n", io::toPrecision(hierarchicalCodelength), numTopModules());
  Log(1).print("Finding super modules done! Added {} levels with hierarchical codelength {} in {} top modules.\n", numLevelsCreated, io::toPrecision(hierarchicalCodelength), numTopModules());

  // Calculate and store hierarchical codelengths
  m_hierarchicalCodelength = calcCodelengthOnTree(root(), true);

  return numLevelsCreated;
}

void InfomapBase::transformNodeFlowToEnterFlow(InfoNode& parent)
{
  double sumFlow = 0.0;
  for (auto& module : m_root) {
    module.data.flow = module.data.enterFlow;
    sumFlow += module.data.flow;
  }
  parent.data.flow = sumFlow;
}

void InfomapBase::resetFlowOnModules()
{
  // Reset flow on all modules
  for (auto& module : m_root.infomapTree()) {
    if (!module.isLeaf())
      module.data.flow = 0.0;
  }
  // Aggregate flow from leaf nodes up in the tree
  for (InfoNode* n : m_leafNodes) {
    double leafNodeFlow = n->data.flow;
    InfoNode* module = n->parent;
    do {
      module->data.flow += leafNodeFlow;
      module = module->parent;
    } while (module != nullptr);
  }
}

unsigned int InfomapBase::removeModules()
{
  // Flatten until every child of root is a leaf. numLevels() only follows the
  // firstChild chain (it assumes a uniform-depth tree, see its TODO), so it
  // under-counts a ragged multilevel tree and would stop early, leaving deeper
  // module subtrees behind. Loop on "any non-leaf child" instead so ragged trees
  // (e.g. from materialize-and-free) flatten completely. For the uniform trees
  // produced elsewhere this does the same number of passes as before.
  auto rootHasModuleChild = [this]() {
    for (auto& child : m_root)
      if (!child.isLeaf())
        return true;
    return false;
  };

  unsigned int numLevelsDeleted = 0;
  while (rootHasModuleChild()) {
    m_root.replaceChildrenWithGrandChildren();
    ++numLevelsDeleted;
  }
  if (numLevelsDeleted > 0) {
    Console::detail(2, "removed {} levels of modules", numLevelsDeleted);
  }
  return numLevelsDeleted;
}

unsigned int InfomapBase::removeSubModules(bool recalculateCodelengthOnTree)
{
  // Loop on "any top module still has a module child" for the reason removeModules
  // gives above: numLevels() follows the firstChild chain only, so on a ragged tree it
  // reports the leftmost branch's depth and this loop stopped while deeper branches
  // were still nested. What came out was neither the two-level tree the caller asked
  // for nor a codelength of it -- calcCodelengthOnTree below then scored a tree that
  // still had sub-modules in it (#898).
  auto anyModuleHasModuleChild = [this]() {
    for (auto& module : m_root)
      for (auto& child : module)
        if (!child.isLeaf())
          return true;
    return false;
  };

  unsigned int numLevelsDeleted = 0;
  while (anyModuleHasModuleChild()) {
    for (InfoNode& module : m_root)
      module.replaceChildrenWithGrandChildren();
    ++numLevelsDeleted;
  }
  if (numLevelsDeleted > 0) {
    Console::detail(2, "removed {} levels of sub-modules", numLevelsDeleted);
    if (recalculateCodelengthOnTree)
      calcCodelengthOnTree(root(), false);
  }
  return numLevelsDeleted;
}

// Sub-modules with fewer children than this are partitioned inline by the
// task that found them instead of as their own task, keeping the task count
// proportional to the work that can actually run in parallel.
constexpr unsigned int minChildDegreeForPartitionTask = 16;

namespace detail {
  // Per-module result of the recursive partitioning task graph. Children are
  // sized once before their tasks are spawned, so each task writes only its
  // own record through a stable address.
  struct PartitionTaskRecord {
    double indexCodelength = 0.0;
    double moduleCodelength = 0.0;
    double leafCodelength = 0.0;
    double queueFlow = 0.0; // flow of the sub-queue spawned from this module
    std::vector<PartitionTaskRecord> children;
  };
} // namespace detail

unsigned int InfomapBase::recursivePartition()
{
  double indexCodelength = getIndexCodelength();
  double hierarchicalCodelength = m_hierarchicalCodelength;

  std::vector<std::string> prettyCompression;
  auto progress = Console().progress("Recursive");
  progress.update("working…");

  PartitionQueue partitionQueue;
  std::string queueSource;
  if (fastHierarchicalSolution > 0) {
    queueLeafModules(partitionQueue);
    queueSource = fmt::format(FMT_STRING("{} sub modules on level {}"), partitionQueue.size(), numLevels() - 2);
  } else {
    queueTopModules(partitionQueue);
    queueSource = fmt::format(FMT_STRING("{} top modules"), partitionQueue.size());
    double moduleCodelength = 0.0;
    for (auto& module : m_root)
      moduleCodelength += module.codelength;
    hierarchicalCodelength = indexCodelength + moduleCodelength;
  }
  Console::detail(1, "recursive: from {}, codelength {:g} + {:g} = {}", queueSource, indexCodelength, hierarchicalCodelength - indexCodelength, io::toPrecision(hierarchicalCodelength));

  double sumConsolidatedCodelength = hierarchicalCodelength - partitionQueue.moduleCodelength;

  const bool shouldMuteNestedMainRun = isMainInfomap() && !Log::isThreadMuted();
  bool isSilent = false;
  if (shouldMuteNestedMainRun) {
    isSilent = Log::isSilent();
  }

  // Run the whole recursion as a task graph: every module is one task that
  // spawns tasks for its accepted sub-modules, so threads flow from finished
  // branches into deeper levels instead of waiting at per-level barriers.
  // Each sub-Infomap re-seeds from the config seed, so per-module results are
  // independent of execution order. Per-module statistics are recorded and
  // aggregated afterwards in the same order as the former level-synchronous
  // queue, reproducing its codelengths and log output exactly.
  const unsigned int startLevel = partitionQueue.level;
  std::vector<detail::PartitionTaskRecord> rootRecords(partitionQueue.size());

  if (shouldMuteNestedMainRun)
    Log::setSilent(true);

  {
    // Spawn the largest modules first so the dominant subtrees start immediately
    std::vector<PartitionQueue::size_t> spawnOrder(partitionQueue.size());
    std::iota(spawnOrder.begin(), spawnOrder.end(), PartitionQueue::size_t { 0 });
    std::sort(spawnOrder.begin(), spawnOrder.end(), [&partitionQueue](PartitionQueue::size_t a, PartitionQueue::size_t b) {
      return partitionQueue[a]->data.flow > partitionQueue[b]->data.flow;
    });

#ifdef _OPENMP
#pragma omp parallel
#pragma omp single
#endif
    {
      for (auto moduleIndex : spawnOrder) {
        // Stop spawning once cancelled; never throw inside the task region (#412).
        if (interruptRequested())
          break;
        InfoNode* module = partitionQueue[moduleIndex];
        detail::PartitionTaskRecord* record = &rootRecords[moduleIndex];
#ifdef _OPENMP
#pragma omp task firstprivate(module, record)
#endif
        partitionModuleRecursively(*module, startLevel, *record);
      }
#ifdef _OPENMP
#pragma omp taskwait
#endif
    }
  }

  if (shouldMuteNestedMainRun)
    Log::setSilent(isSilent);

  // Tasks have joined: throwing here is safe (outside any OpenMP block) and
  // aborts the recursion if a cancel was requested while tasks ran (#412).
  checkCancelled();

  // Aggregate the recorded statistics level by level, visiting the records in
  // the same order as the former level-synchronous queue concatenated them, so
  // all floating-point sums and printed values match it exactly.
  std::vector<const detail::PartitionTaskRecord*> levelRecords;
  levelRecords.reserve(rootRecords.size());
  for (auto& record : rootRecords)
    levelRecords.push_back(&record);

  unsigned int level = startLevel;
  double queueFlow = partitionQueue.flow;
  while (!levelRecords.empty()) {
    const auto levelNum = level;
    const auto levelFlowPct = queueFlow * 100;
    const auto levelModules = levelRecords.size();

    double sumIndexCodelength = 0.0;
    double sumModuleCodelengths = 0.0;
    double sumLeafCodelength = 0.0;
    std::vector<const detail::PartitionTaskRecord*> nextLevelRecords;
    double nextQueueFlow = 0.0;
    for (const auto* record : levelRecords) {
      sumIndexCodelength += record->indexCodelength;
      sumModuleCodelengths += record->moduleCodelength;
      sumLeafCodelength += record->leafCodelength;
      if (!record->children.empty()) {
        nextQueueFlow += record->queueFlow;
        for (const auto& child : record->children)
          nextLevelRecords.push_back(&child);
      }
    }

    double leftToImprove = sumModuleCodelengths;
    sumConsolidatedCodelength += sumIndexCodelength + sumLeafCodelength;
    double limitCodelength = sumConsolidatedCodelength + leftToImprove;

    const double recursiveCompression = ((hierarchicalCodelength - limitCodelength) / hierarchicalCodelength) * 100;
    prettyCompression.push_back(Console::percent(recursiveCompression));
    progress.update(fmt::format(FMT_STRING("working… {} levels · codelength {}"), levelNum, io::toPrecision(limitCodelength)));
    Console::detail(1, "level {}: {} modules ({:g}% flow) {} codelength {:.6g} + {:.6g} (+{:.6g} to improve) = {} bits", levelNum, levelModules, levelFlowPct, Console::arrow(), sumIndexCodelength, sumLeafCodelength, leftToImprove, io::toPrecision(limitCodelength));

    hierarchicalCodelength = limitCodelength;

    levelRecords = std::move(nextLevelRecords);
    queueFlow = nextQueueFlow;
    ++level;
  }

  // Store resulting hierarchical codelength
  m_hierarchicalCodelength = hierarchicalCodelength;

  progress.finish(fmt::format(FMT_STRING("{} {} {} levels, codelength {}"), compressionSummary(prettyCompression), Console::arrow(), level, io::toPrecision(hierarchicalCodelength)));

  return level;
}

void InfomapBase::queueTopModules(PartitionQueue& partitionQueue)
{
  // Add modules to partition queue
  unsigned int numNonTrivialModules = 0;
  partitionQueue.resize(numTopModules());
  double sumFlow = 0.0;
  double sumNonTrivialFlow = 0.0;
  double sumModuleCodelength = 0.0;
  unsigned int moduleIndex = 0;
  for (auto& module : m_root) {
    partitionQueue[moduleIndex] = &module;
    sumFlow += module.data.flow;
    sumModuleCodelength += module.codelength;
    if (module.childDegree() > 1) {
      ++numNonTrivialModules;
      sumNonTrivialFlow += module.data.flow;
    }
    ++moduleIndex;
  }
  partitionQueue.flow = sumFlow;
  partitionQueue.numNonTrivialModules = numNonTrivialModules;
  partitionQueue.nonTrivialFlow = sumNonTrivialFlow;
  partitionQueue.indexCodelength = getIndexCodelength();
  partitionQueue.moduleCodelength = sumModuleCodelength;
  partitionQueue.level = 1;
}

void InfomapBase::queueLeafModules(PartitionQueue& partitionQueue)
{
  unsigned int numLeafModules = 0;
  for (auto& node : m_root.infomapTree()) {
    if (node.isLeafModule())
      ++numLeafModules;
  }

  // Add modules to partition queue
  partitionQueue.resize(numLeafModules);
  unsigned int numNonTrivialModules = 0;
  double sumFlow = 0.0;
  double sumNonTrivialFlow = 0.0;
  double sumModuleCodelength = 0.0;
  unsigned int moduleIndex = 0;
  unsigned int maxDepth = 0;
  for (auto it = m_root.begin_tree(); !it.isEnd(); ++it) {
    if (it->isLeafModule()) {
      auto& module = *it;
      partitionQueue[moduleIndex] = it.current();
      sumFlow += module.data.flow;
      sumModuleCodelength += module.codelength;
      if (module.childDegree() > 1) {
        ++numNonTrivialModules;
        sumNonTrivialFlow += module.data.flow;
      }
      maxDepth = std::max(maxDepth, it.depth());
      ++moduleIndex;
    }
  }
  partitionQueue.flow = sumFlow;
  partitionQueue.numNonTrivialModules = numNonTrivialModules;
  partitionQueue.nonTrivialFlow = sumNonTrivialFlow;
  partitionQueue.indexCodelength = getIndexCodelength();
  partitionQueue.moduleCodelength = sumModuleCodelength;
  partitionQueue.level = maxDepth;
}

void InfomapBase::partitionModuleRecursively(InfoNode& module, unsigned int level, detail::PartitionTaskRecord& record) const
{
  module.codelength = calcCodelength(module);
  // Delete former sub-structure if exists
  if (module.disposeInfomap())
    module.codelength = calcCodelength(module);

  // Inside an OpenMP task: bail without throwing; recursivePartition throws once
  // after the region (#412).
  if (interruptRequested()) {
    record.leafCodelength = module.codelength;
    return;
  }

  // If only trivial substructure is to be found, no need to create infomap instance to find sub-module structures.
  if (module.childDegree() <= 2) {
    module.codelength = calcCodelength(module);
    record.leafCodelength = module.codelength;
    return;
  }

  double oldModuleCodelength = module.codelength;

  auto& subInfomap = getSubInfomap(module)
                         .initNetwork(module);
  // Run two-level partition + find hierarchically super modules (skip recursion).
  // Contain an interrupt from the sub-run so it can't escape the enclosing
  // OpenMP task; free the sub-Infomap and bail (the flag stays set) (#412).
  try {
    subInfomap.setOnlySuperModules(true).run();
  } catch (const InterruptionError&) {
    module.disposeInfomap();
    module.codelength = oldModuleCodelength;
    record.leafCodelength = module.codelength;
    return;
  }

  double subCodelength = subInfomap.getHierarchicalCodelength();
  double subIndexCodelength = subInfomap.root().codelength;
  double subModuleCodelength = subCodelength - subIndexCodelength;
  InfoNode& subRoot = *module.getInfomapRoot();
  unsigned int numSubModules = subRoot.childDegree();
  bool trivialSubPartition = numSubModules == 1 || numSubModules == module.childDegree();
  // #308: bias the refinement accept test by the local recursion depth (level + 1).
  const double subPrefBias = twoLevel
      ? 0.0
      : prefLevelsBias(preferredNumberOfLevels, preferredNumberOfLevelsStrength, level + 1);
  bool improvedCodelength = subCodelength < oldModuleCodelength - minimumCodelengthImprovement + subPrefBias;

  if (trivialSubPartition || !improvedCodelength) {
    Console::detail(1, "disposing unaccepted sub Infomap instance");
    module.disposeInfomap();
    module.codelength = oldModuleCodelength;
    record.leafCodelength = module.codelength;
    return;
  }

  // Improvement
  PartitionQueue subQueue;
  subInfomap.queueTopModules(subQueue);
  record.indexCodelength = subIndexCodelength;
  record.moduleCodelength = subModuleCodelength;
  record.queueFlow = subQueue.flow;

  const auto numSubModulesQueued = subQueue.size();
  // Sized once before spawning so child records have stable addresses
  record.children.resize(numSubModulesQueued);

  // Materialize this sub-partition into real children of `module`, then dispose
  // the sub-Infomap so its working set (cloned active network + pools + maps) is
  // freed immediately instead of being held live until the whole recursion ends.
  // Top-down: `module`'s children are original main-tree leaves, so the new
  // sub-module nodes reparent originals (never clones-of-clones). NOTE: this does
  // NOT guarantee byte-identical output — recursing on original leaves sees a
  // different edge order than the lazy clone path, so the order-sensitive greedy
  // optimizer can settle on a different (quality-equivalent) optimum on deep
  // networks; output is byte-identical only where recursion stays shallow.
  // New module nodes use raw `new` (NOT the instance pool): the
  // task graph runs this function concurrently and ObjectPool is per-instance,
  // not thread-safe; malloc is. They are few, and InfoNode teardown deletes
  // non-pool nodes via `delete`.
  std::vector<InfoNode*> materializedSubModules(numSubModulesQueued, nullptr);
  {
    std::unordered_map<const InfoNode*, InfoNode*> cloneToMaterialized;
    cloneToMaterialized.reserve(numSubModulesQueued);
    for (unsigned int k = 0; k < numSubModulesQueued; ++k) {
      const InfoNode* cloneSubModule = subQueue[k];
      auto* newSubModule = new InfoNode(cloneSubModule->data);
      newSubModule->codelength = cloneSubModule->codelength;
      // Copy the objective-specific module aggregates the lazy sub-Infomap module
      // node carried, so the materialized node is equivalent: calcCodelength and
      // the deeper recursion read these on module nodes (e.g. MemMapEquation reads
      // parent.physicalNodes). Without them, feature/higher-order builds compute
      // wrong codelengths.
      newSubModule->physicalNodes = cloneSubModule->physicalNodes;
      newSubModule->stateNodes = cloneSubModule->stateNodes;
#ifndef SWIG
      // #659 moved these feature-only fields out-of-line into InfoNodeExtras;
      // route the materialize-copy through the extras API and only allocate
      // extras when the clone actually carries them (preserve the shrink).
      if (!cloneSubModule->layerTeleFlowData().empty())
        newSubModule->ensureExtras().layerTeleFlowData = cloneSubModule->layerTeleFlowData();
      if (cloneSubModule->hasMetaCollection())
        newSubModule->ensureMetaCollection() = cloneSubModule->metaCollection();
#endif
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
      newSubModule->lossyEntropy = cloneSubModule->lossyEntropy;
      newSubModule->lossyFlowLogFlow = cloneSubModule->lossyFlowLogFlow;
#endif
      cloneToMaterialized[cloneSubModule] = newSubModule;
      materializedSubModules[k] = newSubModule;
    }

    std::vector<InfoNode*> originalLeaves;
    originalLeaves.reserve(module.childDegree());
    for (auto& child : module)
      originalLeaves.push_back(&child);

    module.releaseChildren();

    // The sub-Infomap clones `module`'s leaves 1:1 in the same order, so clone
    // leaf i corresponds to original leaf i. Assert that invariant explicitly:
    // an indexed loop with a size check fails fast with a clear error if future
    // subnetwork-generation changes ever break it, instead of walking the clone
    // iterator past end() (UB / SIGSEGV).
    const auto& cloneLeaves = subInfomap.leafNodes();
    if (cloneLeaves.size() != originalLeaves.size())
      throw std::logic_error("InfomapBase::partitionModuleRecursively(): sub-Infomap leaf count does not match the module's original leaves");
    for (std::size_t i = 0; i < originalLeaves.size(); ++i) {
      const InfoNode* cloneSubModule = cloneLeaves[i]->parent;
      cloneToMaterialized.at(cloneSubModule)->addChild(originalLeaves[i]);
    }

    for (InfoNode* newSubModule : materializedSubModules)
      module.addChild(newSubModule);

    // module now plays the role the sub-Infomap root did under getInfomapRoot():
    // its index codelength is the sub-partition's index codelength.
    module.codelength = subInfomap.root().codelength;
    module.disposeInfomap();
  }

  // Spawn the largest sub-modules first so the dominant subtrees start immediately.
  // Children recurse into the materialized real nodes (disjoint subtrees, no shared
  // sub-Infomap), so the task graph is unchanged and no taskwait is needed.
  std::vector<PartitionQueue::size_t> spawnOrder(numSubModulesQueued);
  std::iota(spawnOrder.begin(), spawnOrder.end(), PartitionQueue::size_t { 0 });
  std::sort(spawnOrder.begin(), spawnOrder.end(), [&materializedSubModules](PartitionQueue::size_t a, PartitionQueue::size_t b) {
    return materializedSubModules[a]->data.flow > materializedSubModules[b]->data.flow;
  });

  for (auto subModuleIndex : spawnOrder) {
    // Stop spawning once cancelled (in-flight children bail at their own check).
    if (interruptRequested())
      break;
    InfoNode* subModule = materializedSubModules[subModuleIndex];
    detail::PartitionTaskRecord* childRecord = &record.children[subModuleIndex];
    // Small subtrees run inline: spawning them as tasks costs more in task
    // scheduling and idle-thread stealing than their work is worth.
    // NOLINTNEXTLINE(bugprone-branch-clone) -- the branches differ by the task pragma
    if (subModule->childDegree() < minChildDegreeForPartitionTask) {
      partitionModuleRecursively(*subModule, level + 1, *childRecord);
    } else {
#ifdef _OPENMP
#pragma omp task firstprivate(subModule, childRecord)
#endif
      partitionModuleRecursively(*subModule, level + 1, *childRecord);
    }
  }
  // No taskwait needed: child tasks only write their own records, which outlive
  // them, and all tasks complete before the spawning parallel region ends.
}

// ===================================================
// Write output
// ===================================================

void InfomapBase::writeResult(int trial)
{
  writeOutputArtifacts(*this, m_network, OutputPhase::AfterPartition, trial);
}

unsigned int printPerLevelCodelength(const InfoNode& parent, std::ostream& out)
{
  std::vector<detail::PerLevelStat> perLevelStats;
  aggregatePerLevelCodelength(parent, perLevelStats);

  unsigned int numLevels = perLevelStats.size();

  std::vector<std::vector<std::string>> rows;
  rows.reserve(numLevels + 1);
  unsigned int sumNumModules = 0;
  unsigned int sumNumLeafNodes = 0;
  double sumAverageChildDegree = 0.0;
  double sumIndexLengths = 0.0;
  double sumLeafLengths = 0.0;
  double sumCodelengths = 0.0;
  for (unsigned int i = 0; i < numLevels; ++i) {
    double averageChildDegree = perLevelStats[i].numNodes();
    if (i > 0 && perLevelStats[i - 1].numModules > 0) {
      averageChildDegree = perLevelStats[i].numNodes() * 1.0 / perLevelStats[i - 1].numModules;
    }
    sumNumModules += perLevelStats[i].numModules;
    sumNumLeafNodes += perLevelStats[i].numLeafNodes;
    sumAverageChildDegree += averageChildDegree * perLevelStats[i].numNodes();
    sumIndexLengths += perLevelStats[i].indexLength;
    sumLeafLengths += perLevelStats[i].leafLength;
    sumCodelengths += perLevelStats[i].codelength();
    rows.push_back({
        fmt::to_string(i + 1),
        fmt::to_string(perLevelStats[i].numModules),
        fmt::to_string(perLevelStats[i].numLeafNodes),
        io::toPrecision(averageChildDegree, 2, true),
        Console::fixed(perLevelStats[i].indexLength),
        Console::fixed(perLevelStats[i].leafLength),
        Console::fixed(perLevelStats[i].codelength()),
    });
  }
  const unsigned int sumNumNodes = sumNumModules + sumNumLeafNodes;
  rows.push_back({
      "Total",
      fmt::to_string(sumNumModules),
      fmt::to_string(sumNumLeafNodes),
      sumNumNodes > 0 ? io::toPrecision(sumAverageChildDegree / sumNumNodes, 2, true) : "0",
      Console::fixed(sumIndexLengths),
      Console::fixed(sumLeafLengths),
      Console::fixed(sumCodelengths),
  });
  out << "\nLevels\n";
  out << "  Level  Modules  Leaves  Avg children  Module bits  Leaf bits  Total bits\n";
  for (const auto& row : rows) {
    out << "  " << std::setw(5) << row[0]
        << "  " << std::setw(7) << row[1]
        << "  " << std::setw(6) << row[2]
        << "  " << std::setw(12) << row[3]
        << "  " << std::setw(11) << row[4]
        << "  " << std::setw(9) << row[5]
        << "  " << std::setw(10) << row[6] << "\n";
  }
  return numLevels;
}

void aggregatePerLevelCodelength(const InfoNode& parent, std::vector<detail::PerLevelStat>& perLevelStat, unsigned int level)
{
  if (perLevelStat.size() < level + 1)
    perLevelStat.resize(level + 1);

  // Classify each child on its own rather than reading the first one and assuming the
  // rest match. A module with both kinds of child made this crash or lie depending on
  // which came first (#898): a leaf after a sub-module was recursed into, dereferencing
  // its null firstChild, while a leaf *before* one returned early and dropped every
  // remaining sub-module subtree from the report.
  unsigned int numLeafChildren = 0;
  unsigned int numModuleChildren = 0;
  for (const auto& child : parent) {
    if (child.isLeaf())
      ++numLeafChildren;
    else
      ++numModuleChildren;
  }

  if (numLeafChildren > 0) {
    perLevelStat[level].numLeafNodes += numLeafChildren;
    // The codelength of coding these leaves belongs to this level. Charged once even
    // when sub-modules sit alongside them, since it is the parent's own codelength.
    perLevelStat[level].leafLength += parent.codelength;
  }

  if (numModuleChildren == 0)
    return;

  perLevelStat[level].numModules += numModuleChildren;
  if (numLeafChildren == 0)
    perLevelStat[level].indexLength += parent.codelength;

  for (auto& module : parent) {
    if (module.isLeaf())
      continue;
    if (module.getInfomapRoot() != nullptr)
      aggregatePerLevelCodelength(*module.getInfomapRoot(), perLevelStat, level + 1);
    else
      aggregatePerLevelCodelength(module, perLevelStat, level + 1);
  }
}

void InfomapBase::initOptimizer(bool forceNoMemory)
{
  // The value-member root is never pool-freed (isRoot() guards every free path),
  // but it can allocate (e.g. the middle node in replaceChildrenWithOneNode), so
  // it needs this instance's pool. Idempotent: initOptimizer may run twice.
  m_root.m_pool = &m_nodePool;
  m_root.m_edgePool = &m_edgePool;

#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
  // Config::adaptDefaults validates --lossy against the parsed flags, but input
  // parsing can flip state/multilayer input and auto-switch the flow model to
  // directed afterwards. Re-validate here so lossy never silently falls through
  // to another objective.
  if (lossy) {
    if (haveMetaData() || haveMemory() || isMultilayerNetwork())
      throw std::runtime_error("--lossy does not support memory, multilayer or meta-data input");
    if (!isUndirectedFlow())
      throw std::runtime_error("--lossy requires undirected flow");
  }
#endif
  if (haveMetaData()) {
    m_optimizer = std::make_unique<InfomapOptimizer<MetaMapEquation>>();
  } else if (haveMemory() && !forceNoMemory) {
#if INFOMAP_FEATURE_REGULARIZED_MULTILAYER
    if (isRegularizedMultilayerFlow()) {
      m_optimizer = std::make_unique<InfomapOptimizer<RegularizedMultilayerMapEquation>>();
    } else {
      m_optimizer = std::make_unique<InfomapOptimizer<MemMapEquation>>();
    }
#else
    m_optimizer = std::make_unique<InfomapOptimizer<MemMapEquation>>();
#endif
  } else {
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
    if (lossy) {
      m_optimizer = std::make_unique<InfomapOptimizer<LossyMapEquation>>();
    } else {
      m_optimizer = std::make_unique<InfomapOptimizer<BiasedMapEquation>>();
    }
#else
    m_optimizer = std::make_unique<InfomapOptimizer<BiasedMapEquation>>();
#endif
  }
  m_optimizer->init(this);
}

} // namespace infomap
